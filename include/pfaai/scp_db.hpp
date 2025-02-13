///
// @file scp_db.hpp
// @brief The classes to query SCP database (sqlite) and populate the
//        data structure arrays.
// @author Sriram P C <srirampc@gatech.edu>, Hoang Le <hanh9@gatech.edu>
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

#ifndef DATABASE_INTERFACE_HPP
#define DATABASE_INTERFACE_HPP

#include <cstddef>
#include <fmt/format.h>
#include <iostream>
#include <sqlite3.h>
#include <sstream>
#include <string>
#include <vector>

#include "pfaai/db_helper.hpp"
#include "pfaai/interface.hpp"

// Names definition of the database
struct DatabaseNames {
    const std::string TMTAB_SUFFIX = "_tetras";
    const std::string TMTAB_COLUMN_TT = "tetramer";
    const std::string TMTAB_COLUMN_GM = "genomes";
    //
    const std::string GNMTAB_SUFFIX = "_genomes";
    const std::string GNMTAB_COLUMN_GID = "genome_id";
    const std::string GNMTAB_COLUMN_TMS = "tetramers";
    //
    const std::string GMTTAB = "genome_metadata";
    const std::string GMTTAB_COLUMN_GNAME = "genome_name";
    const std::string GMTTAB_COLUMN_GID = "genome_id";
    const std::string GMTTAB_COLUMN_GCLASS = "genome_class";
    const std::string GMTTAB_COLUMN_GLEN = "genome_length";
    const std::string GMTTAB_COLUMN_SCPC = "scp_count";
    //
    const std::string SCPDTAB = "scp_data";
    const std::string SCPDTAB_COLUMN_ACC = "scp_acc";
};

//
// Interface to a SQLite database, which acts as both query and target database
template <typename IdType, typename DBNameType>
class SQLiteSCPDataBase : public DefaultDBInterface<IdType> {
    std::string m_pathToDb;
    sqlite3* m_sqltDbPtr;
    int m_dbErrorCode;
    DBNameType m_dbNames;
    DBMetaData m_dbMeta;

    // Function to load meta data : genome and protein sets
    DBMetaData dbMetaData() {
        DBMetaData dbMeta;

        m_dbErrorCode = SQLiteHelper::dbProteinSet(m_dbNames, m_sqltDbPtr,
                                                   dbMeta.proteinSet);
        if (m_dbErrorCode == SQLITE_OK) {
            m_dbErrorCode = SQLiteHelper::dbGenomeSet(
                m_dbNames, m_sqltDbPtr, m_dbNames.GMTTAB, dbMeta.genomeSet);
        }
        return dbMeta;
    }

  public:
    using ParentT = DefaultDBInterface<IdType>;
    using IdPairType = typename ParentT::IdPairType;
    using IdMatType = typename ParentT::IdMatType;
    using PairIterT = typename ParentT::PairIterT;

    explicit SQLiteSCPDataBase(const std::string dbPath,
                               DBNameType dbn = DBNameType())
        : m_pathToDb(dbPath),
          m_sqltDbPtr(SQLiteHelper::initDB(dbPath, &m_dbErrorCode)),
          m_dbNames(dbn), m_dbMeta(dbMetaData()) {}

    virtual ~SQLiteSCPDataBase() {
        if (m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr) {
            closeDB();
        }
    }

    // getter functions
    virtual inline int getDBErrorCode() const { return m_dbErrorCode; }
    virtual inline std::string getDBPath() const { return m_pathToDb; }
    const DBMetaData& getMeta() const { return m_dbMeta; }

    // Functions to close database connections and validate
    virtual inline int closeDB() {
        m_dbErrorCode = SQLiteHelper::closeDB(m_pathToDb, m_sqltDbPtr);
        m_sqltDbPtr = nullptr;
        m_dbErrorCode = SQLITE_OK;
        return m_dbErrorCode;
    }

    virtual inline PFAAI_ERROR_CODE validate() const {
        if (m_dbErrorCode != SQLITE_OK || m_sqltDbPtr == nullptr) {
            std::cerr << "Error in opening " << getDBPath() << std::endl;
            SQLiteHelper::displayDBError(m_dbErrorCode);
            return PFAAI_ERR_SQLITE_DB;
        }
        return PFAAI_OK;
    }

    // Function to query and load counts in the tetramer counts array, Lc
    virtual int tetramerOccCounts(const std::string protein,
                                  IdPairType tetraRange,
                                  std::vector<IdType>& Lc) const {  // NOLINT
        IdType tetramerStart = tetraRange.first;
        IdType tetramerEnd = tetraRange.second;

        const std::string genomeTetramerQueryFmt =
            "SELECT {}, {} FROM `{}{}` WHERE {} BETWEEN ? AND ?";

        std::string sqlQuery =
            fmt::format(genomeTetramerQueryFmt, m_dbNames.TMTAB_COLUMN_TT,
                        m_dbNames.TMTAB_COLUMN_GM, protein,
                        m_dbNames.TMTAB_SUFFIX, m_dbNames.TMTAB_COLUMN_TT);

        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        sqlite3_bind_int(statement, 1, tetramerStart);
        sqlite3_bind_int(statement, 2, tetramerEnd);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement" << sqlQuery
                      << std::endl;
            std::cerr << "SQL error : " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        }

        while (sqlite3_step(statement) == SQLITE_ROW) {
            const int tetraID = sqlite3_column_int(statement, 0);
            size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
            IdType countGenomes = sizeOfBlobInBytes / sizeof(IdType);

            Lc[tetraID] += countGenomes;
        }
        sqlite3_finalize(statement);
        return errorCode;
    }

    // Function to query and load (protein, genome) pairs in the pairs array F
    virtual int proteinSetGPPairs(IdPairType tetraRange, PairIterT iterF,
                                  IdType* fCount) const {
        const std::vector<std::string>& proteinSet = m_dbMeta.proteinSet;
        IdType tetramerStart = tetraRange.first;
        IdType tetramerEnd = tetraRange.second;

        std::ostringstream oss;
        const std::string proteinSetTetramersQueryFmt =
            "SELECT {}, {}, {} as source_table FROM `{}{}` WHERE {} BETWEEN "
            "{} and {} ";
        for (std::size_t pIndex = 0; pIndex < proteinSet.size(); pIndex++) {
            const std::string protein = proteinSet[pIndex];
            oss << fmt::format(
                proteinSetTetramersQueryFmt, m_dbNames.TMTAB_COLUMN_TT,
                m_dbNames.TMTAB_COLUMN_GM, pIndex, protein,
                m_dbNames.TMTAB_SUFFIX, m_dbNames.TMTAB_COLUMN_TT,
                tetramerStart, tetramerEnd);
            if (pIndex < proteinSet.size() - 1) {
                oss << " UNION ALL ";
            }
        }
        oss << " ORDER BY " << m_dbNames.TMTAB_COLUMN_TT << ", source_table";

        std::string sqlQuery = oss.str();
        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement" << sqlQuery
                      << std::endl;
            std::cerr << "SQL Error: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        } else {
            while (sqlite3_step(statement) == SQLITE_ROW) {
                // const int tetraID =
                sqlite3_column_int(statement, 0);
                const void* genomeBlob = sqlite3_column_blob(statement, 1);
                const IdType proteinIndex = sqlite3_column_int(statement, 2);

                const int* genomeArray = static_cast<const int*>(genomeBlob);
                size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
                int countGenomes = sizeOfBlobInBytes / sizeof(int);

                for (int i = 0; i < countGenomes; i++) {
                    int genomeID = genomeArray[i];
                    *iterF = IdPairType(proteinIndex, genomeID);
                    (*fCount)++;
                    iterF++;
                }
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    // Function to query and load tetramer counts in the pairs matrix T
    virtual int proteinTetramerCounts(IdPairType proteinRange,
                                      IdMatType& T) const {  // NOLINT
        const std::vector<std::string>& proteinSet = m_dbMeta.proteinSet;
        IdType proteinStart = proteinRange.first;
        IdType proteinEnd = proteinRange.second;

        std::ostringstream oss;
        const std::string proteinSetTetramerCtQueryFmt =
            "SELECT {}, length({}), {} as source_table from `{}{}` ";
        for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd;
             proteinIndex++) {
            const std::string& protein = proteinSet[proteinIndex];
            oss << fmt::format(proteinSetTetramerCtQueryFmt,
                               m_dbNames.GNMTAB_COLUMN_GID,
                               m_dbNames.GNMTAB_COLUMN_TMS, proteinIndex,
                               protein, m_dbNames.GNMTAB_SUFFIX);
            if (proteinIndex < proteinEnd) {
                oss << " UNION ALL ";
            }
        }
        std::string sqlQuery = oss.str();

        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "SQLite error: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        } else {
            while (sqlite3_step(statement) == SQLITE_ROW) {
                const int genomeID = sqlite3_column_int(statement, 0);
                int sizeOfBlobInBytes = sqlite3_column_int(statement, 1);
                int countTetras = sizeOfBlobInBytes / sizeof(int);
                int proteinIndex = sqlite3_column_int(statement, 2);
                T(proteinIndex, genomeID) = countTetras;
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }
};

//
// Interface to two different SQLite databases -- one query and one target
template <typename IdType, typename DBNameType>
class QTSQLiteSCPDataBase : public DefaultDBInterface<IdType> {
  public:
    using ParentT = DefaultDBInterface<IdType>;
    using IdPairType = typename ParentT::IdPairType;
    using IdMatType = typename ParentT::IdMatType;
    enum DBIndicator {
        target_db = 0,
        query_db = 1,
    };

  private:
    std::string m_pathToTgtDb, m_pathToQryDb;
    sqlite3* m_sqltDbPtr;
    DBNameType m_dbNames;
    int m_dbErrorCode;
    DBMetaData m_dbMeta;
    IdType m_qryGenomeIdOffset;

    DBMetaData dbMetaData() {
        DBMetaData dbMeta;
        m_dbErrorCode = SQLiteHelper::qtDBProteinSet(
            m_dbNames, m_sqltDbPtr,
            fmt::format("`{}`.{}", c_dbAlias[DBIndicator::target_db],
                        m_dbNames.SCPDTAB),
            fmt::format("`{}`.{}", c_dbAlias[DBIndicator::query_db],
                        m_dbNames.SCPDTAB),
            dbMeta.proteinSet);
        if (m_dbErrorCode == SQLITE_OK) {
            m_dbErrorCode = SQLiteHelper::dbGenomeSet(
                m_dbNames, m_sqltDbPtr,
                fmt::format("`{}`.{}", c_dbAlias[DBIndicator::target_db],
                            m_dbNames.GMTTAB),
                dbMeta.genomeSet);
        }
        if (m_dbErrorCode == SQLITE_OK) {
            m_dbErrorCode = SQLiteHelper::dbGenomeSet(
                m_dbNames, m_sqltDbPtr,
                fmt::format("`{}`.{}", c_dbAlias[DBIndicator::query_db],
                            m_dbNames.GMTTAB),
                dbMeta.qyGenomeSet);
        }
        return dbMeta;
    }
    inline std::string tetraTableName(std::string protein,
                                      DBIndicator dbind) const {
        static const std::string tabNameFmt = "{}.`{}{}`";
        return fmt::format(tabNameFmt, c_dbAlias[dbind], protein,
                           m_dbNames.TMTAB_SUFFIX);
    }

    inline std::string genomeTableName(std::string protein,
                                       DBIndicator dbind) const {
        static const std::string tabNameFmt = "{}.`{}{}`";
        return fmt::format(tabNameFmt, c_dbAlias[dbind], protein,
                           m_dbNames.GNMTAB_SUFFIX);
    }

    const char* c_dbAlias[2] = {"main", "QueryDB"};

  public:
    explicit QTSQLiteSCPDataBase(const std::string& tgtDBPath,
                                 const std::string& qryDBPath,
                                 DBNameType dbn = DBNameType())
        : m_pathToTgtDb(tgtDBPath), m_pathToQryDb(qryDBPath),
          m_sqltDbPtr(SQLiteHelper::initDB(tgtDBPath, &m_dbErrorCode)),
          m_dbNames(dbn) {
        if (m_dbErrorCode) {
            SQLiteHelper::displayDBError(m_dbErrorCode);
            return;
        }

        // Attach database
        const std::string attachSQL =
            fmt::format("ATTACH DATABASE '{}' as {} ;", m_pathToQryDb,
                        c_dbAlias[DBIndicator::query_db]);
        m_dbErrorCode = sqlite3_exec(m_sqltDbPtr, attachSQL.c_str(), nullptr,
                                     nullptr, nullptr);
        if (m_dbErrorCode) {
            std::cerr << "Error in attaching query database : " << m_pathToQryDb
                      << std::endl;
            std::cerr << "SQL Error: " << sqlite3_errmsg(m_sqltDbPtr);
            return;
        }

        m_dbMeta = dbMetaData();
        m_qryGenomeIdOffset = m_dbMeta.genomeSet.size();
    }

    virtual ~QTSQLiteSCPDataBase() {
        if (m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr) {
            closeDB();
        }
    }

    // getter functions
    virtual inline int getDBErrorCode() const { return m_dbErrorCode; }
    virtual inline std::string getDBPath() const { return m_pathToTgtDb; }
    const DBMetaData& getMeta() const { return m_dbMeta; }

    // Functions to close database connections
    virtual inline int closeDB() {
        if (m_sqltDbPtr == nullptr) {
            std::cerr << "Database already closed : " << m_pathToQryDb
                      << std::endl;
            return m_dbErrorCode;
        }
        // Detach database
        const std::string detachSQL = fmt::format(
            "DETACH DATABASE '{}' ;", c_dbAlias[DBIndicator::query_db]);
        m_dbErrorCode = sqlite3_exec(m_sqltDbPtr, detachSQL.c_str(), nullptr,
                                     nullptr, nullptr);
        if (m_dbErrorCode) {
            std::cerr << "Error in deatach database : " << m_pathToQryDb
                      << std::endl;
            std::cerr << "SQL Error: " << sqlite3_errmsg(m_sqltDbPtr);
            return m_dbErrorCode;
        }

        m_dbErrorCode = sqlite3_close(m_sqltDbPtr);
        m_sqltDbPtr = nullptr;
        m_dbErrorCode = SQLITE_OK;
        return m_dbErrorCode;
    }

    PFAAI_ERROR_CODE validate() const {
        if (m_dbErrorCode != SQLITE_OK || m_sqltDbPtr == nullptr) {
            std::cerr << "Error in opening " << getDBPath() << std::endl;
            SQLiteHelper::displayDBError(m_dbErrorCode);
            return PFAAI_ERR_SQLITE_DB;
        }
        return PFAAI_OK;
    }

    // Function to query and load counts in the tetramer counts array, Lc
    int tetramerOccCounts(const std::string protein, IdPairType tetraRange,
                          std::vector<IdType>& Lc) const {  // NOLINT

        IdType tetramerStart = tetraRange.first;
        IdType tetramerEnd = tetraRange.second;

        const char* queryFormat = "SELECT target_table.{} as tetramer, \n"
                                  "    target_table.{} as target_genomes, \n"
                                  "    query_table.{} as query_genomes  \n"
                                  "FROM {} as target_table,  \n"
                                  "     {} as query_table  \n"
                                  "WHERE target_table.{} = query_table.{} \n"
                                  "  AND target_table.{} BETWEEN ? AND ?";

        std::string sqlQuery =
            fmt::format(queryFormat, m_dbNames.TMTAB_COLUMN_TT,
                        m_dbNames.TMTAB_COLUMN_GM, m_dbNames.TMTAB_COLUMN_GM,
                        tetraTableName(protein, DBIndicator::target_db),
                        tetraTableName(protein, DBIndicator::query_db),
                        m_dbNames.TMTAB_COLUMN_TT, m_dbNames.TMTAB_COLUMN_TT,
                        m_dbNames.TMTAB_COLUMN_TT);

        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        sqlite3_bind_int(statement, 1, tetramerStart);
        sqlite3_bind_int(statement, 2, tetramerEnd);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement" << sqlQuery
                      << std::endl;
            std::cerr << "SQL error : " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        }

        while (sqlite3_step(statement) == SQLITE_ROW) {
            const int tetraID = sqlite3_column_int(statement, 0);
            size_t sizeOfTgtBlob = sqlite3_column_bytes(statement, 1);
            size_t sizeOfQryBlob = sqlite3_column_bytes(statement, 2);
            IdType countGenomes = (sizeOfTgtBlob / sizeof(IdType)) +
                                  (sizeOfQryBlob / sizeof(IdType));
            Lc[tetraID] += countGenomes;
        }
        sqlite3_finalize(statement);
        return errorCode;
    }

    virtual int proteinSetGPPairs(IdPairType tetraRange,
                                  typename ParentT::PairIterT iterF,
                                  IdType* fCount) const {
        const std::vector<std::string>& proteinSet = m_dbMeta.proteinSet;
        IdType tetramerStart = tetraRange.first;
        IdType tetramerEnd = tetraRange.second;

        std::ostringstream oss;
        // Query to jointly from both databases
        const char* queryFormat = "SELECT target_table.{0} as tetramer, \n"
                                  "    target_table.{1} as target_genomes,\n"
                                  "    query_table.{1} as query_genomes, \n"
                                  "    {2} as source_table \n"
                                  "FROM {3} as target_table, \n"
                                  "     {4} as query_table \n"
                                  "WHERE target_table.{0} = query_table.{0} \n"
                                  "  AND target_table.{0} BETWEEN {5} AND {6} ";

        for (std::size_t pIndex = 0; pIndex < proteinSet.size(); pIndex++) {
            const std::string protein = proteinSet[pIndex];

            // Qualified table names
            oss << fmt::format(queryFormat, m_dbNames.TMTAB_COLUMN_TT,
                               m_dbNames.TMTAB_COLUMN_GM, pIndex,
                               tetraTableName(protein, DBIndicator::target_db),
                               tetraTableName(protein, DBIndicator::query_db),
                               tetramerStart, tetramerEnd);
            if (pIndex < proteinSet.size() - 1) {
                oss << "\n UNION ALL \n";
            }
        }
        oss << "\n ORDER BY " << m_dbNames.TMTAB_COLUMN_TT << ", source_table";

        std::string sqlQuery = oss.str();
        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement" << sqlQuery
                      << std::endl;
            std::cerr << "SQL Error: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        } else {
            while (sqlite3_step(statement) == SQLITE_ROW) {
                // const int tetraID =
                sqlite3_column_int(statement, 0);
                const void* tgtBlob = sqlite3_column_blob(statement, 1);
                const void* qryBlob = sqlite3_column_blob(statement, 2);
                const IdType proteinIndex = sqlite3_column_int(statement, 3);

                size_t sizeOfTargetBlob = sqlite3_column_bytes(statement, 1);
                size_t sizeOfQryBlob = sqlite3_column_bytes(statement, 2);
                int tgtGenomes = sizeOfTargetBlob / sizeof(int);
                int qryGenomes = sizeOfQryBlob / sizeof(int);

                // genome array from target table
                const int* genomeArray = static_cast<const int*>(tgtBlob);
                for (int i = 0; i < tgtGenomes; i++) {
                    int genomeID = genomeArray[i];
                    *iterF = IdPairType(proteinIndex, genomeID);
                    (*fCount)++;
                    iterF++;
                }
                // genome array from query table
                genomeArray = static_cast<const int*>(qryBlob);
                for (int i = 0; i < qryGenomes; i++) {
                    int genomeID = genomeArray[i];
                    *iterF = IdPairType(proteinIndex,
                                        m_qryGenomeIdOffset + genomeID);
                    (*fCount)++;
                    iterF++;
                }
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    // Function to query and load tetramer counts in the pairs matrix T
    template <typename MapFn>
    int dbProteinTetramerCounts(IdPairType protRange, DBIndicator dbind,
                                MapFn&& gidMapFn,
                                IdMatType& T) const {  // NOLINT
        //
        const std::vector<std::string>& proteinSet = m_dbMeta.proteinSet;
        IdType proteinStart = protRange.first;
        IdType proteinEnd = protRange.second;
        //
        std::ostringstream oss;
        const std::string proteinSetTetramerCtQueryFmt =
            "SELECT {}, length({}), {} as source_table from {} ";
        for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd;
             proteinIndex++) {
            const std::string& protein = proteinSet[proteinIndex];
            oss << fmt::format(proteinSetTetramerCtQueryFmt,
                               m_dbNames.GNMTAB_COLUMN_GID,
                               m_dbNames.GNMTAB_COLUMN_TMS, proteinIndex,
                               genomeTableName(protein, dbind));
            if (proteinIndex < proteinEnd) {
                oss << " UNION ALL ";
            }
        }
        std::string sqlQuery = oss.str();

        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "SQLite error: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        } else {
            while (sqlite3_step(statement) == SQLITE_ROW) {
                const int genomeID = sqlite3_column_int(statement, 0);
                int sizeOfBlobInBytes = sqlite3_column_int(statement, 1);
                int countTetras = sizeOfBlobInBytes / sizeof(int);
                int proteinIndex = sqlite3_column_int(statement, 2);
                T(proteinIndex, gidMapFn(genomeID)) += countTetras;
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    virtual int proteinTetramerCounts(IdPairType protRange,
                                      IdMatType& T) const {  // NOLINT
        //  Get counts from both databases
        dbProteinTetramerCounts(
            protRange, DBIndicator::target_db, [](IdType gid) { return gid; },
            T);
        dbProteinTetramerCounts(
            protRange, DBIndicator::query_db,
            [&](IdType gid) { return m_dbMeta.genomeSet.size() + gid; }, T);
        return SQLITE_OK;
    }
};

#endif  // !DATABASE_INTERFACE_HPP
