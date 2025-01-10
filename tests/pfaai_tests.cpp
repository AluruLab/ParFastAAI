///
// @file parfaai_test.cpp
// @brief Main Entry function for tests for Parallel Fast AAI
// @author Sriram P C <srirampc@gatech.edu>
//
// Copyright 2024 Georgia Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include "pfaai/algorithm_impl.hpp"
#include "pfaai/ds_helper.hpp"
#include "pfaai/ds_impl.hpp"
#include "pfaai/interface.hpp"
#include "pfaai/scp_db.hpp"
#include "pfaai_tests.hpp"  // NOLINT

using IdType = int;
using ValueType = double;
using IdPairType = DPair<IdType, IdType>;
using IdMatType = DMatrix<IdType>;
using SQLiteDB = SQLiteSCPDataBase<IdType, DatabaseNames>;
using QTSQLiteDB = QTSQLiteSCPDataBase<IdType, DatabaseNames>;
using PFImpl = ParFAAIImpl<IdType, ValueType>;
using PFDataT = ParFAAIData<IdType>;
using PFQSubDataT = ParFAAIQSubData<IdType>;
using PFQTData = ParFAAIQryTgtData<IdType>;

static constexpr char G_DB_PATH[] = "data/modified_xantho_fastaai2.db";
static constexpr char G_DB_SUBSET1_PATH[] = "data/xdb_subset1.db";
static constexpr char G_DB_SUBSET2_PATH[] = "data/xdb_subset2.db";
static constexpr char G_DB_SUBSET_COMBO_12_PATH[] =
    "data/xdb_subset_combo12.db";
static constexpr IdType G_NTETRAMERS = (20 * 20 * 20 * 20);
static constexpr IdType G_GENOMESET_SIZE = 20;
static constexpr IdType G_PROTEINSET_SIZE = 80;
static constexpr IdType G_QT_GENOMESET_SIZE = 8;
static const std::vector<std::string> G_PROTEINSET = TESTDB_PROTEIN_SET;
static const std::vector<std::string> G_GENOMESET = TESTDB_GENOME_SET;
static const std::vector<std::string> G_QT_TARGET_GENOMESET =
    TESTDB_SUBSET1_GENOME_SET;
static const std::vector<std::string> G_QT_QUERY_GENOMESET =
    TESTDB_SUBSET2_GENOME_SET;
static const std::vector<std::string> G_QT_PROTEINSET = QT_TESTDB_PROTEIN_SET;
// Testing of genomes query for protein
static const std::vector<std::vector<IdType>> G_QRY_PST_TETRA_CTS = {
    {260, 260, 260, 260, 260, 260, 260, 260, 260, 260,
     260, 260, 260, 260, 260, 260, 260, 263, 263, 263},
    {244, 244, 244, 244, 244, 244, 244, 244, 244, 244,
     244, 244, 244, 244, 244, 244, 244, 246, 246, 246},
    {386, 386, 386, 386, 386, 386, 386, 386, 386, 386,
     386, 386, 386, 386, 386, 386, 386, 383, 384, 384},
    {121, 121, 121, 121, 121, 121, 121, 121, 121, 121,
     121, 121, 121, 121, 121, 121, 121, 121, 121, 121}};

// TETRAMER, # genomes, source table
static constexpr IdType G_QRT_PST_GP_PAIRS[388][3] = TESTDB_PSET_GP_PAIRS;
static std::vector<IdPairType> G_QT_LC_VALUES = TEST_QT_LC_VALUES;

// Reference data files to compare
static constexpr char REF_LC_ARRAY[] = "data/xanthodb_lc_array.bin";
static constexpr char REF_LP_ARRAY[] = "data/xanthodb_lp_array.bin";
static constexpr char REF_F_ARRAY[] = "data/xanthodb_f_array.bin";
static constexpr char REF_T_MATRIX[] = "data/xanthodb_t_matrix.bin";
static constexpr char REF_SRTD_E_ARRAY[] = "data/xanthodb_sorted_e_array.bin";
static constexpr char REF_JAC_DATA[] = "data/xanthodb_jac.bin";
static constexpr char REF_AJI_DATA[] = "data/xanthodb_aji.bin";

static constexpr char REF_SUBSET1_LC_ARRAY[] = "data/xdb_subset1_lc_array.bin";
static constexpr char REF_SUBSET1_LP_ARRAY[] = "data/xdb_subset1_lp_array.bin";
static constexpr char REF_SUBSET1_F_ARRAY[] = "data/xdb_subset1_f_array.bin";
static constexpr char REF_SUBSET1_T_MATRIX[] = "data/xdb_subset1_t_matrix.bin";
static constexpr char REF_SUBSET1_SRTD_E_ARRAY[] =
    "data/xdb_subset1_sorted_e_array.bin";
static constexpr char REF_SUBSET1_JAC_DATA[] = "data/xdb_subset1_jac.bin";
static constexpr char REF_SUBSET1_AJI_DATA[] = "data/xdb_subset1_aji.bin";

static constexpr char REF_SUBSET2_LC_ARRAY[] = "data/xdb_subset2_lc_array.bin";
static constexpr char REF_SUBSET2_LP_ARRAY[] = "data/xdb_subset2_lp_array.bin";
static constexpr char REF_SUBSET2_F_ARRAY[] = "data/xdb_subset2_f_array.bin";
static constexpr char REF_SUBSET2_T_MATRIX[] = "data/xdb_subset2_t_matrix.bin";
static constexpr char REF_SUBSET2_SRTD_E_ARRAY[] =
    "data/xdb_subset2_sorted_e_array.bin";
static constexpr char REF_SUBSET2_JAC_DATA[] = "data/xdb_subset2_jac.bin";
static constexpr char REF_SUBSET2_AJI_DATA[] = "data/xdb_subset2_aji.bin";

static constexpr char DATA_QSUB_TEST_INPUT[] = "data/qsub_test_input.txt";
static constexpr char REF_QSUB_JAC_DATA[] = "data/xdb_qry_subset_jac.bin";
static constexpr char REF_QSUB_AJI_DATA[] = "data/xdb_qry_subset_aji.bin";

static constexpr char REF_QT_LC_ARRAY[] = "data/xdb_qt_lc_array.bin";
static constexpr char REF_QT_LP_ARRAY[] = "data/xdb_qt_lp_array.bin";
static constexpr char REF_QT_F_ARRAY[] = "data/xdb_qt_f_array.bin";
static constexpr char REF_QT_T_MATRIX[] = "data/xdb_qt_t_matrix.bin";
static constexpr char REF_QT_COMBO_T_MATRIX[] =
    "data/xdb_qt_combo_t_matrix.bin";

TEST_CASE("Query Genome Metadata", "[db_meta_data]") {
    SQLiteDB sqlt_if(G_DB_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    REQUIRE(G_PROTEINSET == dbMeta.proteinSet);
    REQUIRE(G_GENOMESET == dbMeta.genomeSet);
}

TEST_CASE("Query No. of Tetramers", "[ntetramers_query]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0);
    SQLiteDB sqlt_if(G_DB_PATH);
    sqlt_if.tetramerOccCounts("PF00119.20", IdPairType(2060, 2144), Lc);
    REQUIRE(Lc[2060] == 20);
    REQUIRE(Lc[2100] == 20);
    REQUIRE(Lc[2140] == 3);
    REQUIRE(Lc[2144] == 17);
}

TEST_CASE("Query Protein Tetramer Counts", "[prot_tetramer_counts]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    IdMatType T(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    SQLiteDB sqlt_if(G_DB_PATH);
    sqlt_if.proteinTetramerCounts(IdPairType(0, 3), T);
    REQUIRE(G_QRY_PST_TETRA_CTS[0] ==
            std::vector<IdType>(T.row_begin(0), T.row_end(0)));
    REQUIRE(G_QRY_PST_TETRA_CTS[1] ==
            std::vector<IdType>(T.row_begin(1), T.row_end(1)));
    REQUIRE(G_QRY_PST_TETRA_CTS[2] ==
            std::vector<IdType>(T.row_begin(2), T.row_end(2)));
    REQUIRE(G_QRY_PST_TETRA_CTS[3] ==
            std::vector<IdType>(T.row_begin(3), T.row_end(3)));
}

TEST_CASE("Query Protein Set Tetramers", "[prot_set_tetramers]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    for (auto& rx : G_QRT_PST_GP_PAIRS) {
        Lc[rx[0]] += rx[1] / sizeof(IdType);
    }
    // Parallel prefix sum on Lc to construct Lp
    DefaultDSHelper<IdType>::parallelPrefixSum(Lc, Lp);
    REQUIRE(Lp[2000] == 0);

    std::vector<IdPairType> F(Lp.back() + Lc.back(), IdPairType(-1, -1));
    IdType fct = 0;
    SQLiteDB sqlt_if(G_DB_PATH);
    int rc = sqlt_if.proteinSetGPPairs(IdPairType(2000, 3000), F.begin(), &fct);
    REQUIRE(rc == SQLITE_OK);
    REQUIRE(F[Lp[2000]] == IdPairType(5, 0));
    REQUIRE(F[Lp[2415] - 1] == IdPairType(35, 0x13));
    REQUIRE(F[Lp[2415]] == IdPairType(74, 0x11));
}

TEST_CASE("Test Data Structures Construction", "[construct_LFTE]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);

    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Test Data Subset 1 Structures Construction",
          "[construct_subset1_LFTE]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET1_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_SUBSET1_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_SUBSET1_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_SUBSET1_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_SUBSET1_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SUBSET1_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Test Data Subset 2 Structures Construction",
          "[construct_subset2_LFTE]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET2_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_SUBSET2_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_SUBSET2_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_SUBSET2_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_SUBSET2_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    //{
    //    std::ofstream lcx1(REF_SUBSET1_SRTD_E_ARRAY, std::ios::binary);
    //    cereal::BinaryOutputArchive orachive(lcx1);
    //    orachive(E);
    //}
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SUBSET2_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Implementation Test", "[compute_JAC_AJI]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //

    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;
    {
        std::ifstream tcx(REF_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;
    {
        std::ifstream tcx(REF_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("Subset 1 Implementation Test", "[compute_subset1_JAC_AJI]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET1_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //

    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;

    {
        std::ifstream tcx(REF_SUBSET1_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;
    {
        std::ifstream tcx(REF_SUBSET1_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("Subset 2 Implementation Test", "[compute_subset2_JAC_AJI]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET2_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;

    {
        std::ifstream tcx(REF_SUBSET2_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;

    {
        std::ifstream tcx(REF_SUBSET2_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("PF Query Subset Test", "[query_subset]") {
    SQLiteSCPDataBase<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    std::vector<std::string> qryGenomeSet;
    std::ifstream in_stream(DATA_QSUB_TEST_INPUT);
    std::string str;
    while (in_stream >> str) {
        qryGenomeSet.push_back(str);
    }
    //
    PFQSubDataT pfaaiData(sqlt_if, dbMeta, qryGenomeSet);
    int errorCode = pfaaiData.construct();

    sqlt_if.closeDB();
    //
    assert(errorCode == PFAAI_OK);
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;
    {
        std::ifstream tcx(REF_QSUB_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;

    {
        std::ifstream tcx(REF_QSUB_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
}

TEST_CASE("QT Genome Metadata", "[qt_db_meta_data]") {
    QTSQLiteDB sqlt_if(G_DB_SUBSET1_PATH, G_DB_SUBSET2_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    REQUIRE(G_QT_PROTEINSET == dbMeta.proteinSet);
    REQUIRE(G_QT_TARGET_GENOMESET == dbMeta.genomeSet);
    REQUIRE(G_QT_QUERY_GENOMESET == dbMeta.qyGenomeSet);
}

TEST_CASE("QT No. of Tetramers", "[qt_ntetramers_query]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0);
    QTSQLiteDB sqlt_if(G_DB_SUBSET1_PATH, G_DB_SUBSET2_PATH);
    sqlt_if.tetramerOccCounts("PF00119.20", IdPairType(1000, 2500), Lc);
    REQUIRE(Lc[1151] == 8);
    REQUIRE(Lc[1255] == 8);
    REQUIRE(Lc[1704] == 7);
    REQUIRE(Lc[2144] == 7);
}

TEST_CASE("QT Protein Tetramer Counts", "[qt_prot_tetramer_counts]") {
    QTSQLiteDB qtsqlt_if(G_DB_SUBSET1_PATH, G_DB_SUBSET2_PATH);
    const DBMetaData& dbMeta = qtsqlt_if.getMeta();
    IdMatType T(G_PROTEINSET_SIZE, G_QT_GENOMESET_SIZE);
    qtsqlt_if.proteinTetramerCounts(IdPairType(0, dbMeta.proteinSet.size() - 1),
                                    T);

    SQLiteDB sqlt_if(G_DB_SUBSET_COMBO_12_PATH);
    const DBMetaData& dbMeta2 = sqlt_if.getMeta();
    IdMatType T2(dbMeta2.proteinSet.size(), dbMeta2.genomeSet.size());
    sqlt_if.proteinTetramerCounts(IdPairType(0, dbMeta.proteinSet.size()), T2);

    //
    REQUIRE(dbMeta.proteinSet ==
            std::vector<std::string>(dbMeta2.proteinSet.begin(),
                                     dbMeta2.proteinSet.begin() +
                                         dbMeta.proteinSet.size()));
    for (std::size_t ix = 0; ix < dbMeta.proteinSet.size(); ix++) {
        REQUIRE(std::vector<IdType>(T2.row_begin(ix), T2.row_end(ix)) ==
                std::vector<IdType>(T.row_begin(ix), T.row_end(ix)));
    }
}

void constructL(const std::vector<std::string>& protSet,
                const DefaultDBInterface<IdType>& sqlt_if,
                std::vector<IdType>& Lc, std::vector<IdType>& Lp) {  // NOLINT
#pragma omp parallel default(none) shared(protSet, sqlt_if, Lc)
    {
        int nThreads = omp_get_num_threads(), threadID = omp_get_thread_num();
        IdPairType tetraRange(BLOCK_LOW(threadID, nThreads, G_NTETRAMERS),
                              BLOCK_HIGH(threadID, nThreads, G_NTETRAMERS));
        for (const std::string& protein : protSet) {
            sqlt_if.tetramerOccCounts(protein, tetraRange, Lc);
        }
    }
    // Parallel prefix sum on Lc to construct
    DefaultDSHelper<IdType>::parallelPrefixSum(Lc, Lp);
}

TEST_CASE("QT Protein Set Tetramers", "[qt_prot_set_tetramers]") {
    QTSQLiteDB qtsqlt_if(G_DB_SUBSET1_PATH, G_DB_SUBSET2_PATH);
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    // construct L
    constructL(qtsqlt_if.getMeta().proteinSet, qtsqlt_if, Lc, Lp);
    std::vector<IdType> tLc(3001, 0);
    for (auto& lx : G_QT_LC_VALUES) {
        tLc[lx.first] = lx.second;
    }

    // test first 3000 value
    REQUIRE(std::vector(Lc.begin(), Lc.begin() + 3001) == tLc);

    // construct F
    int fct = 0;
    std::vector<IdPairType> F(Lp.back() + Lc.back(), IdPairType(-1, -1));
    int rc = qtsqlt_if.proteinSetGPPairs(IdPairType(0, 3000), F.begin(), &fct);
    REQUIRE(rc == SQLITE_OK);
    REQUIRE(fct == Lp[3001]);
    REQUIRE(F[Lp[12]] == IdPairType(17, 0));
    REQUIRE(F[Lp[16]] == IdPairType(2, 0));
}

TEST_CASE("QT Test Data Structures Construction", "[qt_construct_LFTE]") {
    QTSQLiteDB sqlt_if(G_DB_SUBSET1_PATH, G_DB_SUBSET2_PATH);
    const DBMetaData& dbMeta = sqlt_if.getMeta();
    //
    PFQTData pfaaiQTData(sqlt_if, dbMeta, dbMeta.proteinSet);
    pfaaiQTData.constructL();
    std::vector<IdType> Lc = pfaaiQTData.refLc();
    std::vector<IdType> Lp = pfaaiQTData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_QT_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_QT_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    std::vector<IdType> tLc(3001, 0);
    for (auto& lx : G_QT_LC_VALUES) {
        tLc[lx.first] = lx.second;
    }

    // test first 3000 value
    REQUIRE(std::vector(Lc.begin(), Lc.begin() + 3001) == tLc);
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiQTData.constructF();
    std::vector<IdPairType> F = pfaaiQTData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_QT_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);

    pfaaiQTData.constructT();
    IdMatType T = pfaaiQTData.refT();
    IdMatType refT(dbMeta.proteinSet.size(),
                   dbMeta.genomeSet.size() + dbMeta.qyGenomeSet.size()),
        comboT(dbMeta.proteinSet.size(),
               dbMeta.genomeSet.size() + dbMeta.qyGenomeSet.size());
    {
        std::ifstream tcx(REF_QT_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
    {
        std::ifstream tcx(REF_QT_COMBO_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(comboT);
    }
    REQUIRE(IdMatType(comboT.rows() - 1, comboT.cols(), comboT.data().data()) ==
            T);
    // TODO(x)

    // //
    // pfaaiQTData.constructE();
    // const std::vector<ETriple<IdType>>& E = pfaaiQTData.refE().E;
    // std::vector<ETriple<IdType>> refE;
    // {
    //     std::ifstream tcx(REF_SRTD_E_ARRAY, std::ios::binary);
    //     cereal::BinaryInputArchive iarchive(tcx);
    //     iarchive(refE);
    // }
    // REQUIRE(E == refE);
    // IdType nTotal = 0;
    // for (const auto& ex : E)
    //     if (ex.genomeA == 0 || ex.genomeA == 2)
    //         nTotal++;
    // std::cout << "N total : " << nTotal << std::endl;
}
