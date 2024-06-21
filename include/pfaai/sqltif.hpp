#ifndef SQLITE_INTERFACE_HPP
#define SQLITE_INTERFACE_HPP

#include <fmt/format.h>
#include <iostream>
#include <sqlite3.h>
#include <sstream>
#include <string>
#include <vector>

struct DatabaseNames {
    const std::string TMTAB_SUFFIX = "_tetras";
    const std::string TMTAB_COLUMN_TT = "tetra";
    const std::string TMTAB_COLUMN_GM = "genomes";
    const std::string TMTAB_TABLE_COLUMNS[2] = {TMTAB_COLUMN_TT,
                                                TMTAB_COLUMN_GM};
    //
    const std::string GNMTAB_SUFFIX = "_genomes";
    const std::string GNMTAB_COLUMN_GID = "genome_id";
    const std::string GNMTAB_COLUMN_TMS = "tetramers";
    //
    const std::string MTTAB = "genome_metadata";
    //
    const std::string SCPDTAB = "scp_data";
    const std::string SCPDTAB_COLUMN_ACC = "scp_acc";
};

template <typename DBNameType> class SQLiteInterface {
    std::string m_pathToDb;
    sqlite3* m_sqltDbPtr;
    int m_dbErrorCode;
    DBNameType m_dbNames;

  public:
    inline bool isDBOpenError() const { return m_dbErrorCode != SQLITE_OK; }

    inline bool isDBNull() const { return m_sqltDbPtr == nullptr; }

    const char* getDBError() const {
        if (m_dbErrorCode == SQLITE_OK)
            return "NO ERROR AVAILABLE";
        return sqlite3_errstr(m_dbErrorCode);
    }

    static sqlite3* initDB(const std::string dbPath, int* errorCode) {
        sqlite3* db;
        *errorCode = sqlite3_open(dbPath.c_str(), &db);
        return db;
    }

    explicit SQLiteInterface(const std::string dbPath, DBNameType dbn)
        : m_pathToDb(dbPath), m_sqltDbPtr(initDB(dbPath, &m_dbErrorCode)),
          m_dbNames(dbn) {}

    inline bool isDBOpen() const {
        return m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr;
    }

    std::string getDBPath() const { return m_pathToDb; }

    inline int closeDB() {
        int errorCode = sqlite3_close(m_sqltDbPtr);
        m_sqltDbPtr = nullptr;
        m_dbErrorCode = SQLITE_OK;
        return errorCode;
    }

    ~SQLiteInterface() {
        if (m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr) {
            closeDB();
        }
    }

    template <typename IdType>
    int queryGenomeTetramers(const std::string protein, IdType tetramerStart,
                             IdType tetramerEnd, std::vector<IdType>& Lc) {
        // std::string sqlQuery = "SELECT tetra, genomes FROM `" + protein +
        //                        "_tetras` WHERE tetra BETWEEN ? AND ?";
        //
        // std::string sqlQuery0 = "SELECT " + _dbNames.TMTAB_COLUMN_TT + "," +
        //                       _dbNames.TMTAB_COLUMN_GM + " FROM `" + protein
        //                       + _dbNames.TMTAB_SUFFIX +
        //                       "` WHERE tetra BETWEEN ? AND ?";
        //
        const std::string genomeTetramerQueryFmt =
            "SELECT {}, {} FROM `{}_{}` WHERE {} BETWEEN ? AND ?";

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

    template <typename IdType, typename IdPairType>
    int queryProtienSetTetramers(const std::vector<std::string>& proteinSet,
                                 IdType tetramerStart, IdType tetramerEnd,
                                 std::vector<IdType>& Lp,
                                 std::vector<IdPairType>& F) {
        std::ostringstream oss;
        const std::string proteinSetTetramersQueryFmt =
            "SELECT {}, {}, {} as source_table FROM `{}_{}` WHERE {} BETWEEN "
            "{} and {} ";
        for (int pIndex = 0; pIndex < proteinSet.size(); pIndex++) {
            const std::string protein = proteinSet[pIndex];
            oss << fmt::format(
                proteinSetTetramersQueryFmt, m_dbNames.TMTAB_COLUMN_TT,
                m_dbNames.TMTAB_COLUMN_GM, pIndex, protein,
                m_dbNames.TMTAB_SUFFIX, tetramerStart, tetramerEnd);
            // oss << "SELECT " << _dbNames.TMTAB_COLUMN_TT << ","
            //     << _dbNames.TMTAB_COLUMN_GM << ", " << pIndex
            //     << " as source_table FROM `" + protein +
            //     _dbNames.TMTAB_SUFFIX +
            //            "` WHERE "
            //     << _dbNames.TMTAB_COLUMN_TT << " BETWEEN " << tetramerStart
            //     << " AND " << tetramerEnd << " ";
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
            // Starting location in F
            IdType indexInF = Lp[tetramerStart];
            while (sqlite3_step(statement) == SQLITE_ROW) {
                const int tetraID = sqlite3_column_int(statement, 0);
                const void* genomeBlob = sqlite3_column_blob(statement, 1);
                const IdType proteinIndex = sqlite3_column_int(statement, 2);

                const int* genomeArray = static_cast<const int*>(genomeBlob);
                size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
                int countGenomes = sizeOfBlobInBytes / sizeof(int);

                for (int i = 0; i < countGenomes; i++) {
                    int genomeID = genomeArray[i];
                    F[indexInF] = IdPairType(proteinIndex, genomeID);
                    indexInF++;
                }
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    template <typename IdType>
    int queryProtienSetTetramerCount(const std::vector<std::string>& proteinSet,
                                     IdType proteinStart, IdType proteinEnd,
                                     std::vector<std::vector<IdType>>& T) {
        std::ostringstream oss;
        const std::string proteinSetTetramerCtQueryFmt =
            "SELECT {}, length({}), {} as source_table from `{}_{}` ";
        for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd;
             proteinIndex++) {
            std::string protein = proteinSet[proteinIndex];
            oss << fmt::format(proteinSetTetramerCtQueryFmt,
                               m_dbNames.GNMTAB_COLUMN_GID,
                               m_dbNames.GNMTAB_COLUMN_TMS, proteinIndex,
                               protein, m_dbNames.GNMTAB_SUFFIX);
            // oss << "SELECT " << _dbNames.GNMTAB_COLUMN_GID << ", length("
            //     << _dbNames.GNMTAB_COLUMN_TMS << "), " << proteinIndex
            //     << " as source_table from `" << protein
            //     << _dbNames.GNMTAB_SUFFIX + "` ";
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
                T[proteinIndex][genomeID] = countTetras;
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    int queryMetaData(std::vector<std::string>& proteinSet, int& nGenomes) {

        const std::string countQueryFmt =
            "SELECT count(*) as count_genome FROM {}";
        // std::string sqlQuery =
        //     "SELECT count(*) as count_genome FROM " + _dbNames.MTTAB;
        std::string sqlQuery = fmt::format(countQueryFmt, m_dbNames.MTTAB);
        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        }

        if (sqlite3_step(statement) == SQLITE_ROW) {
            nGenomes = sqlite3_column_int(statement, 0);
        }

        const std::string protCountQueryFmt = "SELECT count(*) FROM {}";
        // sqlQuery = "SELECT count(*) from " + _dbNames.SCPDTAB;
        sqlQuery = fmt::format(protCountQueryFmt, m_dbNames.SCPDTAB);
        errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                       &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        }

        if (sqlite3_step(statement) == SQLITE_ROW) {
            int scpCount = sqlite3_column_int(statement, 0);
            proteinSet.resize(scpCount);
        }

        const std::string protQueryFmt = "SELECT {} FROM {}";
        // sqlQuery = "SELECT scp_acc from " + _dbNames.SCPDTAB;
        sqlQuery = fmt::format(protQueryFmt, m_dbNames.SCPDTAB_COLUMN_ACC,
                               m_dbNames.SCPDTAB);
        errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
                                       &statement, nullptr);
        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(m_sqltDbPtr);
            return errorCode;
        }

        int index = 0;
        while (sqlite3_step(statement) == SQLITE_ROW) {
            const unsigned char* scp_acc = sqlite3_column_text(statement, 0);
            std::string proteinName(reinterpret_cast<const char*>(scp_acc));
            proteinSet[index] = proteinName;
            index++;
        }

        return SQLITE_OK;
    }
};

#endif  // !SQLITE_INTERFACE_HPP
