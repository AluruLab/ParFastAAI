#ifndef SQLITE_INTERFACE_HPP
#define SQLITE_INTERFACE_HPP

#include "pfaai/interface.hpp"
#include <cstddef>
#include <fmt/format.h>
#include <iostream>
#include <iterator>
#include <sqlite3.h>
#include <sstream>
#include <string>
#include <vector>

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

template <typename IdType, typename IdPairType, typename IdMatType,
          typename DBNameType>
class SQLiteInterface
    : public DataBaseInterface<IdType, IdPairType, IdMatType> {
    std::string m_pathToDb;
    sqlite3* m_sqltDbPtr;
    int m_dbErrorCode;
    DBNameType m_dbNames;

  public:
    using ParentT = DataBaseInterface<IdType, IdPairType, IdMatType>;
    explicit SQLiteInterface(const std::string dbPath,
                             DBNameType dbn = DBNameType())
        : m_pathToDb(dbPath), m_sqltDbPtr(initDB(dbPath, &m_dbErrorCode)),
          m_dbNames(dbn) {}

    virtual inline bool isDBOpen() const {
        return m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr;
    }

    virtual inline std::string getDBPath() const { return m_pathToDb; }

    virtual inline const char* getDBError() const {
        if (m_dbErrorCode == SQLITE_OK) {
            if (m_sqltDbPtr == NULL) {
                return "SQLite is unable to allocate memory for the database ";
            }
            return "NO Error Message Available";
        }
        return sqlite3_errstr(m_dbErrorCode);
    }

    virtual inline int closeDB() {
        int errorCode = sqlite3_close(m_sqltDbPtr);
        m_sqltDbPtr = nullptr;
        m_dbErrorCode = SQLITE_OK;
        return errorCode;
    }

    virtual int getDBErrorCode() const { return m_dbErrorCode; }

    PFAAI_ERROR_CODE validate() const {
        if (!isDBOpen()) {
            std::cerr << "Error in opening " << getDBPath() << std::endl;
            std::cerr << "DB ERROR: " << getDBError() << std::endl;

            if (getDBError())
                return PFAAI_ERR_SQLITE_MEM_ALLOC;
            else
                return PFAAI_ERR_SQLITE_DB;
        }
        return PFAAI_OK;
    }

    virtual ~SQLiteInterface() {
        if (m_dbErrorCode == SQLITE_OK && m_sqltDbPtr != nullptr) {
            closeDB();
        }
    }

    static sqlite3* initDB(const std::string dbPath, int* errorCode) {
        sqlite3* db;
        *errorCode = sqlite3_open(dbPath.c_str(), &db);
        return db;
    }

    int queryGenomeTetramers(const std::string protein, IdType tetramerStart,
                             IdType tetramerEnd,
                             std::vector<IdType>& Lc) const {  // NOLINT
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

    int queryProteinSetGPPairs(const std::vector<std::string>& proteinSet,
                               IdType tetramerStart, IdType tetramerEnd,
                               typename ParentT::PIterT iterF) const {
        std::ostringstream oss;
        const std::string proteinSetTetramersQueryFmt =
            "SELECT {}, {}, {} as source_table FROM `{}{}` WHERE {} BETWEEN "
            "{} and {} ";
        for (int pIndex = 0; pIndex < proteinSet.size(); pIndex++) {
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
                const int tetraID = sqlite3_column_int(statement, 0);
                const void* genomeBlob = sqlite3_column_blob(statement, 1);
                const IdType proteinIndex = sqlite3_column_int(statement, 2);

                const int* genomeArray = static_cast<const int*>(genomeBlob);
                size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
                int countGenomes = sizeOfBlobInBytes / sizeof(int);

                for (int i = 0; i < countGenomes; i++) {
                    int genomeID = genomeArray[i];
                    *iterF =  IdPairType(proteinIndex, genomeID);
                    iterF++;
                }
            }
        }

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    int queryProteinTetramerCounts(const std::vector<std::string>& proteinSet,
                                   IdType proteinStart, IdType proteinEnd,
                                   IdMatType& T) const {  // NOLINT
        std::ostringstream oss;
        const std::string proteinSetTetramerCtQueryFmt =
            "SELECT {}, length({}), {} as source_table from `{}{}` ";
        for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd;
             proteinIndex++) {
            std::string protein = proteinSet[proteinIndex];
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

    int queryGenomeSet(std::vector<std::string>& genomeSet) const {  // NOLINT
        const std::string countQueryFmt =
            "SELECT count(*) as count_genome FROM {}";
        // std::string sqlQuery =
        //     "SELECT count(*) as count_genome FROM " + _dbNames.GMTTAB;
        std::string sqlQuery = fmt::format(countQueryFmt, m_dbNames.GMTTAB);
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
            int nGenomes = sqlite3_column_int(statement, 0);
            genomeSet.resize(nGenomes);
        }
        sqlite3_finalize(statement);
        statement = NULL;

        const std::string genomeQueryFmt = "SELECT {} FROM {}";
        sqlQuery = fmt::format(genomeQueryFmt, m_dbNames.GMTTAB_COLUMN_GNAME,
                               m_dbNames.GMTTAB);
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
            const unsigned char* genome_name =
                sqlite3_column_text(statement, 0);
            std::string genomeName(reinterpret_cast<const char*>(genome_name));
            genomeSet[index] = genomeName;
            index++;
        }

        return SQLITE_OK;
    }

    int queryProteinSet(std::vector<std::string>& proteinSet) const {  // NOLINT
        sqlite3_stmt* statement;
        const std::string protCountQueryFmt =
            "SELECT count(DISTINCT {}) FROM {}";
        // sqlQuery = "SELECT count(*) from " + _dbNames.SCPDTAB;
        std::string sqlQuery = fmt::format(
            protCountQueryFmt, m_dbNames.SCPDTAB_COLUMN_ACC, m_dbNames.SCPDTAB);
        int errorCode = sqlite3_prepare_v2(m_sqltDbPtr, sqlQuery.c_str(), -1,
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
        sqlite3_finalize(statement);
        statement = NULL;

        const std::string protQueryFmt = "SELECT DISTINCT {} FROM {}";
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

        sqlite3_finalize(statement);
        return SQLITE_OK;
    }

    int queryMetaData(std::vector<std::string>& proteinSet,         // NOLINT
                      std::vector<std::string>& genomeSet) const {  // NOLINT
        int errCode = queryProteinSet(proteinSet);
        if (errCode == SQLITE_OK) {
            return queryGenomeSet(genomeSet);
        } else {
            return errCode;
        }
    }
};

#endif  // !SQLITE_INTERFACE_HPP
