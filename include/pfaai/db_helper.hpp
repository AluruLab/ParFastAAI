///
// @file db_helper.hpp
// @brief Helper classes to query data base (sqlite) and populate the
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

#ifndef DATABASE_HELPER_HPP
#define DATABASE_HELPER_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <sqlite3.h>

struct SQLiteHelper {
    static inline sqlite3* initDB(const std::string& dbPath, int* errorCode) {
        sqlite3* db;
        *errorCode = sqlite3_open(dbPath.c_str(), &db);
        return db;
    }

    static inline const char* displayDBError(int dbErrorCode) {
        if (dbErrorCode == SQLITE_OK) {
            return "NO Error Message Available";
        }
        const char* errMsg = sqlite3_errstr(dbErrorCode);
        std::cerr << "DB ERROR: " << errMsg << std::endl;
        return errMsg;
    }

    static inline int closeDB(const std::string& dbPath, sqlite3* sqltDbPtr) {
        if (sqltDbPtr == nullptr) {
            std::cerr << "Database already closed : " << dbPath
                      << std::endl;
            return -1;
        }
        return sqlite3_close(sqltDbPtr);
    }

    // Function to query set of genomes in sqlite db
    template <typename DBNameType>
    static int dbGenomeSet(const DBNameType& dbNames, sqlite3* sqltDbPtr,
                           const std::string& genomeTableName,
                           std::vector<std::string>& genomeSet) {  // NOLINT
        const std::string countQueryFmt =
            "SELECT count(*) as count_genome FROM {}";
        // std::string sqlQuery =
        //     "SELECT count(*) as count_genome FROM " + _dbNames.GMTTAB;
        std::string sqlQuery = fmt::format(countQueryFmt, genomeTableName);
        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
            return errorCode;
        }

        if (sqlite3_step(statement) == SQLITE_ROW) {
            int nGenomes = sqlite3_column_int(statement, 0);
            genomeSet.resize(nGenomes);
        }
        sqlite3_finalize(statement);
        statement = NULL;

        const std::string genomeQueryFmt = "SELECT {} FROM {}";
        sqlQuery = fmt::format(genomeQueryFmt, dbNames.GMTTAB_COLUMN_GNAME,
                               genomeTableName);
        errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                       &statement, nullptr);
        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
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

    template <typename DBNameType>
    static int qtDBProteinSet(const DBNameType& dbNames, sqlite3* sqltDbPtr,
                              const std::string& mainTableName,
                              const std::string& qryTableName,
                              std::vector<std::string>& proteinSet) {  // NOLINT
        const char* ctQueryFormat =
            "SELECT COUNT(DISTINCT target_table.{0}) \n"
            "  FROM {1} as target_table, {2} as query_table \n"
            "  WHERE target_table.{0} = query_table.{0};";
        std::string sqlQuery =
            fmt::format(ctQueryFormat, dbNames.SCPDTAB_COLUMN_ACC,
                        mainTableName, qryTableName);
        sqlite3_stmt* statement;
        int errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
            return errorCode;
        }

        if (sqlite3_step(statement) == SQLITE_ROW) {
            int scpCount = sqlite3_column_int(statement, 0);
            proteinSet.resize(scpCount);
        }
        sqlite3_finalize(statement);
        statement = NULL;

        const char* queryFormat =
            "SELECT DISTINCT target_table.{0} \n"
            "  FROM {1} as target_table, {2} as query_table \n"
            "  WHERE target_table.{0} = query_table.{0};";

        // sqlQuery = "SELECT scp_acc from " + _dbNames.SCPDTAB;
        sqlQuery = fmt::format(queryFormat, dbNames.SCPDTAB_COLUMN_ACC,
                               mainTableName, qryTableName);
        errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                       &statement, nullptr);
        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
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
        return errorCode;
    }

    // Function to query set of proteins
    template <typename DBNameType>
    static int dbProteinSet(const DBNameType& dbNames, sqlite3* sqltDbPtr,
                            std::vector<std::string>& proteinSet) {  // NOLINT
        sqlite3_stmt* statement;
        const std::string protCountQueryFmt =
            "SELECT count(DISTINCT {}) FROM {}";
        // sqlQuery = "SELECT count(*) from " + ;
        std::string sqlQuery = fmt::format(
            protCountQueryFmt, dbNames.SCPDTAB_COLUMN_ACC, dbNames.SCPDTAB);
        int errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                           &statement, nullptr);

        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
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
        sqlQuery = fmt::format(protQueryFmt, dbNames.SCPDTAB_COLUMN_ACC,
                               dbNames.SCPDTAB);
        errorCode = sqlite3_prepare_v2(sqltDbPtr, sqlQuery.c_str(), -1,
                                       &statement, nullptr);
        if (errorCode != SQLITE_OK) {
            std::cerr << "Error in preparing sql statement " << sqlQuery
                      << std::endl;
            std::cerr << "The error was: " << sqlite3_errmsg(sqltDbPtr);
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
};

#endif  // DATABASE_HELPER_HPP
