import numpy as np
import pandas as pd
import sqlite3
import sqlalchemy as sqlac
import sqlalchemy.types as sqlacts
import typing


def load_table_df(conn, table_name):
    return pd.read_sql("SELECT * from '{}'".format(table_name), conn)


class DBSubsetBuilder:
    def __init__(self, src_db: str, dst_db, genome_lst: typing.List[str]) -> None:
        self.src_engine = sqlac.create_engine("sqlite:///" + src_db)
        self.src_conn = self.src_engine.connect()
        self.dst_engine = sqlac.create_engine("sqlite:///" + dst_db)
        self.dst_conn = self.dst_engine.connect()
        self.genomes = genome_lst
        self.genome_ids = {gx: ix for ix, gx in enumerate(genome_lst)}
        self.tetras_dtypes = {"tetramer": sqlacts.INTEGER, "genomes": sqlacts.BLOB}
        self.genomes_dtypes = {"genome_id": sqlacts.INTEGER, "tetramers": sqlacts.BLOB}
        self.meta_dypes = {
            "genome_name": sqlacts.TEXT,
            "genome_id": sqlacts.INTEGER,
            "genome_length": sqlacts.INTEGER,
            "genome_class": sqlacts.INTEGER,
            "SCP_count": sqlacts.INTEGER,
        }
        self.scp_dtypes = {
            "genome_id": sqlacts.INTEGER,
            "SCP_acc": sqlacts.TEXT,
            "SCP_score": sqlacts.REAL,
            "tetra_count": sqlacts.INTEGER,
        }
        self.get_metadata()

    # generate ids for the subset genes starting with 0
    # Step 0:  Acquire data
    #   1. Get protein names, genome names
    #   2. Current ids  for the genomes
    #   3. r_genome_ids : Genome ids of the genomes need to be removed
    def get_metadata(self):
        self.tab_names = pd.read_sql(
            "SELECT name FROM sqlite_schema WHERE type='table'", self.src_conn
        )
        self.protein_df = pd.read_sql("SELECT scp_acc from 'scp_data'",
                                      self.src_conn)
        self.protien_names = set(self.protein_df["SCP_acc"].tolist())
        print("No. of Proteins : ", len(self.protien_names))
        #
        self.genome_rdf = pd.read_sql("SELECT * from genome_metadata",
                                      self.src_conn)
        self.scp_df = pd.read_sql("SELECT * from scp_data", self.src_conn)
        #
        self.gselect_df = self.genome_rdf[
            self.genome_rdf["genome_name"].isin(set(self.genomes))  # type: ignore
        ]
        self.gremove_df = self.genome_rdf[
            self.genome_rdf["genome_name"].isin(set(self.genomes)) == False  # type: ignore
        ]
        #
        self.select_ids = self.gselect_df["genome_id"].tolist()
        grecords = self.gselect_df[["genome_name", "genome_id"]].to_dict(  # type:ignore
            "records"
        )
        self.c2s_map = {}
        for rcd in grecords:
            gname = rcd["genome_name"]
            gid = rcd["genome_id"]
            self.c2s_map[gid] = self.genome_ids[gname]
        self.rgenomes_str = ",".join(str(rx) for rx in self.gremove_df["genome_id"])

    # Step 1: Update _tetras
    # for each ${protien_id}_tetras table
    #    - Remove the rows with r_genome_ids
    #    - Update the existing with the
    #
    def subset_tetras_df(self, tetras_df, ids_select, ids_map):
        def filter_ids(row_nx):
            return np.extract(np.isin(row_nx, ids_select), row_nx)

        def map_to_ids(row_nx):
            return np.array(
                sorted([ids_map[y] for y in row_nx if y in ids_map]),
                dtype=np.int32
            )

        #
        tetras_df["genomes"] = (
            tetras_df["genomes"]
            .apply(np.frombuffer, dtype=np.int32)
            .apply(filter_ids)
            .apply(map_to_ids)
            .apply(lambda x: x.tobytes())
        )
        tetras_df = tetras_df[tetras_df["genomes"] != b""]
        tetras_df.reset_index(drop=True)
        return tetras_df

    def create_tetras_table(self, table_name, tetras_df):
        # drop_stmt = "DROP TABLE '{}'".format(table_name)
        create_stmt = "CREATE TABLE IF NOT EXISTS '{}' (tetramer INTEGER PRIMARY KEY, genomes BLOB)".format(
            table_name
        )
        idx_stmt = "CREATE INDEX `{}_index` ON `{}` (tetramer)".format(
            table_name, table_name
        )
        self.dst_conn.execute(sqlac.text(create_stmt))
        tetras_df.to_sql(
            table_name,
            self.dst_conn,
            index=False,
            # index_label="tetramer",
            if_exists="append",
            dtype=self.tetras_dtypes,  # type:ignore
        )
        self.dst_conn.execute(sqlac.text(idx_stmt))
        self.dst_conn.commit()

    def subset_tetras_table(self, table_name):
        tetras_df = pd.read_sql("SELECT * from '{}'".format(table_name),
                                self.src_conn)
        tetras_df = self.subset_tetras_df(tetras_df, self.select_ids,
                                          self.c2s_map)
        self.create_tetras_table(table_name, tetras_df)

    # Step 2: Update genomes table
    #   1.  Delete the genomes with ids to be removed
    #   2.  Update the table with new genome ids
    def create_genomes_table(self, table_name, genomes_df):
        # delete_stmt = "DELETE FROM '{}' where genome_id in ({})".format(
        #     table_name, self.rgenomes_str
        # )
        # self.src_conn.execute(sqlac.text(drop_stmt))
        # self.src_conn.commit()
        # update_stmt = "UPDATE '{}' SET genome_id = {} WHERE genome_id == {}"
        # drop_stmt = "DROP TABLE '{}'".format(table_name)
        create_stmt = "CREATE TABLE IF NOT EXISTS '{}' (genome_id INTEGER PRIMARY KEY, tetramers BLOB)".format(
            table_name
        )
        self.dst_conn.execute(sqlac.text(create_stmt))
        genomes_df.to_sql(
            table_name,
            self.dst_conn,
            index=False,
            # index_label="genome_id",
            if_exists="append",
            dtype=self.genomes_dtypes,  # type:ignore
        )
        self.dst_conn.commit()
        # self.conn.execute(sqlac.text(delete_stmt))
        # self.conn.commit()
        # for rx, ux in sorted(self.c2s_map.items()):
        #     print("UPD: ", rx, ux)
        #     self.conn.execute(
        #         sqlac.text(update_stmt.format(table_name, str(rx), str(ux)))
        #     )
        #     self.conn.commit()
        # self.conn.commit()

    def subset_genomes_table(self, table_name):
        select_stmt = "SELECT * FROM '{}'".format(table_name)
        genomes_df = pd.read_sql(select_stmt, self.src_conn)
        genomes_df = genomes_df[genomes_df["genome_id"].isin(self.select_ids)]
        genomes_df["genome_id"] = genomes_df["genome_id"].map(  # type:ignore
            self.c2s_map  # type:ignore
        )
        genomes_df = genomes_df.reset_index(drop=True)
        self.create_genomes_table(table_name, genomes_df)

    #
    # Step 3: Update scp_data table
    def create_scp_data(self):
        uscp_df = self.scp_df[self.scp_df["genome_id"].isin(self.select_ids)].copy()
        uscp_df["genome_id"] = uscp_df["genome_id"].map(  # type:ignore
            self.c2s_map  # type:ignore
        )
        uscp_df.reset_index(drop=True)
        # print(self.scp_df.shape, uscp_df.shape)
        # drop_stmt = "DROP TABLE '{}'".format("scp_data")
        # self.src_conn.execute(sqlac.text(drop_stmt))
        # self.src_conn.commit()
        create_stmt = "CREATE TABLE IF NOT EXISTS 'scp_data' (genome_id INTEGER, SCP_acc TEXT, SCP_score REAL, tetra_count INTEGER)"
        self.dst_conn.execute(sqlac.text(create_stmt))
        uscp_df.to_sql(
            "scp_data",
            self.dst_conn,
            index=False,
            # index_label="genome_id",
            if_exists="append",
            dtype=self.scp_dtypes,  # type:ignore
        )
        self.dst_conn.commit()

    # Step 4: Update genome_metadata table
    def create_genome_meta_data(self):
        ugenomes_df = self.genome_rdf[
            self.genome_rdf["genome_id"].isin(self.select_ids)
        ].copy()
        ugenomes_df["genome_id"] = ugenomes_df["genome_id"].map(  # type:ignore
            self.c2s_map  # type:ignore
        )
        ugenomes_df = ugenomes_df.reset_index(drop=True)
        print(ugenomes_df)
        # drop_stmt = "DROP TABLE '{}'".format("genome_metadata")
        # self.src_conn.execute(sqlac.text(drop_stmt))
        # self.src_conn.commit()

        create_stmt = "CREATE TABLE IF NOT EXISTS 'genome_metadata' (genome_name TEXT, genome_id INTEGER PRIMARY KEY, genome_length INTEGER, genome_class INTEGER, SCP_count INTEGER);"
        self.dst_conn.execute(sqlac.text(create_stmt))
        self.dst_conn.commit()
        ugenomes_df.to_sql(
            "genome_metadata",
            self.dst_conn,
            index=False,
            # index_label="genome_id",
            if_exists="append",
            dtype=self.meta_dypes,  # type:ignore
        )
        self.dst_conn.commit()

    def create_prtein_indices(self):
        idxp_cstmt = """CREATE TABLE index_protein (
                                protein_number INTEGER PRIMARY KEY,
                                protein_string VARCHAR(255) NOT NULL
                        )"""
        protix_cstmt = """CREATE TABLE protein_index (
                                protein_string VARCHAR(255) NOT NULL PRIMARY KEY,
                                protein_number INTEGER
                        )"""
        idxp_df = pd.read_sql("SELECT * FROM 'index_protein'", self.src_conn)
        protix_df = pd.read_sql("SELECT * FROM 'protein_index'", self.src_conn)
        self.dst_conn.execute(sqlac.text(idxp_cstmt))
        idxp_df.to_sql(
            "index_protein",
            self.dst_conn,
            index=False,
            # index_label="genome_id",
            if_exists="append",
            dtype={
                "protein_number": sqlacts.INTEGER,
                "protein_string": sqlacts.VARCHAR,
            },  # type:ignore
        )
        self.dst_conn.commit()
        #
        self.dst_conn.execute(sqlac.text(protix_cstmt))
        protix_df.to_sql(
            "protein_index",
            self.dst_conn,
            index=False,
            # index_label="genome_id",
            if_exists="append",
            dtype={
                "protein_string": sqlacts.VARCHAR,
                "protein_number": sqlacts.INTEGER,
            },
        )
        self.dst_conn.commit()

    def run(self):
        self.create_genome_meta_data()
        self.create_scp_data()
        self.create_prtein_indices()
        for ix, px in enumerate(self.protien_names):
            print("Processing : ", ix, px)
            self.subset_tetras_table(px + "_tetras")
            self.subset_genomes_table(px + "_genomes")

    def close(self):
        self.src_conn.close()
        self.dst_conn.close()


def gen_subset_db(src_db, subset_db, subset):
    dsub = DBSubsetBuilder(src_db, subset_db, subset)
    dsub.run()
    dsub.close()


source_db = "./modified_xantho_fastaai2.db"
rsubset1_loc = "./xdb_subset1.db"
rsubset1 = [
    "Xanthomonas_albilineans_GCA_000962915_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000962945_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963065_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963195_1.fna.gz",
]
rsubset2_loc = "./xdb_subset2.db"
rsubset2 = [
    "Xanthomonas_albilineans_GCA_000962925_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963025_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963155_1.fna.gz",
    "_Pseudomonas__cissicola_GCA_008801575_1.fna.gz",
]
rsubset12_loc = "./xdb_subset_combo12.db"
rsubset12 = [
    "Xanthomonas_albilineans_GCA_000962915_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000962945_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963065_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963195_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000962925_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963025_1.fna.gz",
    "Xanthomonas_albilineans_GCA_000963155_1.fna.gz",
    "_Pseudomonas__cissicola_GCA_008801575_1.fna.gz",
]
if __name__ == "__main__":
    gen_subset_db(source_db, rsubset1_loc, rsubset1)
    gen_subset_db(source_db, rsubset2_loc, rsubset2)
    gen_subset_db(source_db, rsubset12_loc, rsubset12)
