///
// @file main.cpp
// @brief Main Entry function for Parallel Fast AAI
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

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <omp.h>
#include <sqlite3.h>

#include "pfaai/algorithm_impl.hpp"
#include "pfaai/ds_impl.hpp"
#include "pfaai/interface.hpp"
#include "pfaai/scp_db.hpp"

using IdType = int;
using ValueType = double;
using IdPairType = DPair<IdType, IdType>;
using IdMatrixType = DMatrix<IdType>;

// database classes
using SQLiteDB = SQLiteSCPDataBase<IdType, DatabaseNames>;
using QTSQLiteDB = QTSQLiteSCPDataBase<IdType, DatabaseNames>;

// data structure classes
using PFDSInterface = DefaultDataStructInterface<IdType>;
using PFData = ParFAAIData<IdType>;
using PFQSubData = ParFAAIQSubData<IdType>;
using PFQTData = ParFAAIQryTgtData<IdType>;

// algorithm implementation
using PFImpl = ParFAAIImpl<IdType, ValueType>;

struct AppParams {
    CLI::App app;
    std::string pathToDatabase;
    std::string pathToQryDatabase;
    std::string pathToQrySubsetFile;
    std::string pathToOutputFile;
    std::string outFieldSeparator;
    std::vector<std::string> qryGenomeSet;

    AppParams()
        : app(), pathToDatabase(""), pathToQryDatabase(""),
          pathToQrySubsetFile(""), pathToOutputFile(""),
          outFieldSeparator(",") {
        app.add_option("path_to_input_db", pathToDatabase,
                       "Path to the Input Database")
            ->check(CLI::ExistingFile)
            ->required(true);
        app.add_option("path_to_output_file", pathToOutputFile,
                       "Path to output csv file.")
            ->required(true);
        app.add_option("-r,--query_db", pathToQryDatabase,
                       "Path to the Query Database [Optional (default: Same "
                       "as the Input DB)]")
            ->check(CLI::ExistingFile);
        app.add_option(
               "-s,--separator", outFieldSeparator,
               "Field Separator in the output file [Optional (default: ,)].")
            ->capture_default_str();
        app.add_option("-q,--query_subset", pathToQrySubsetFile,
                       "Path to Query List (Should be subset of genomoes in "
                       "the input DB.)")
            ->check(CLI::ExistingFile);
    }

    void print() const {
        // Display arguments
        std::string arg1 =
            fmt::format(" Input Database  : {} ", pathToDatabase);
        std::string arg2 =
            fmt::format(" Query Database  : {} ", pathToQryDatabase);
        std::string arg3 =
            fmt::format(" Query Subset    : {} ", pathToQrySubsetFile);
        std::string arg4 =
            fmt::format(" Output File     : {} ", pathToOutputFile);
        std::string arg5 =
            fmt::format(" Field Separator : {} ", outFieldSeparator);
        std::size_t argsz = std::max(
            {arg1.size(), arg2.size(), arg3.size(), arg4.size(), arg5.size()});
        fmt::print(" ┌{0:─^{1}}┐\n"
                   " │{2: <{1}}│\n"
                   " │{3: <{1}}│\n"
                   " │{4: <{1}}│\n"
                   " │{5: <{1}}│\n"
                   " │{6: <{1}}│\n"
                   " └{0:─^{1}}┘\n",
                   "", argsz, arg1, arg2, arg3, arg4, arg5);
    }

    void load_query_genomes() {
        std::ifstream in_stream(pathToQrySubsetFile);
        std::string str;
        while (in_stream >> str) {
            qryGenomeSet.push_back(str);
        }

        // Printing each new line separately
        // std::copy(v.begin(), v.end(),
        //          std::ostream_iterator<std::string>(std::cout, ","));
    }

    void print_query_genomes() {
        std::cout << "Query Set: [" << qryGenomeSet.size() << ", "
                  << fmt::format("({})", fmt::join(qryGenomeSet, ", ")) << "]"
                  << std::endl;
    }
};

void printOutput(const PFDSInterface& pfdata, const PFImpl& impl,
                 const AppParams& appArgs, bool isSubset = true) {
    IdType nQryGenomes = pfdata.qrySetSize();
    IdType nTgtGeneomes = pfdata.tgtSetSize();
    //
    const std::vector<PFImpl::JACType>& cJAC = impl.getJAC();
    const std::vector<ValueType>& cAJI = impl.getAJI();

    fmt::print("Writing output with {} query genomes and {} target genomes. \n",
               nQryGenomes, nTgtGeneomes);
    DMatrix<ValueType> ajiMatrix(nQryGenomes, nTgtGeneomes, 0.0);
    for (std::size_t i = 0; i < cJAC.size(); i++) {
        IdType fga = cJAC[i].genomeA;
        IdType fgb = cJAC[i].genomeB;
        ajiMatrix(pfdata.mapQueryId(fga), pfdata.mapTargetId(fgb)) = cAJI[i];
        //
        // Other side of the diagonal if genomeB is also a query genome
        if (isSubset && pfdata.isQryGenome(cJAC[i].genomeB)) {
            ajiMatrix(pfdata.mapQueryId(fgb), pfdata.mapTargetId(fga)) =
                cAJI[i];
        }
    }

    FILE* fpx = fopen(appArgs.pathToOutputFile.c_str(), "w");

    const auto& tgtRef = pfdata.refTargetSet();
    // Header
    fmt::print(
        fpx, "{}{}\n", appArgs.outFieldSeparator,
        fmt::join(tgtRef.begin(), tgtRef.end(), appArgs.outFieldSeparator));

    // Data
    const std::vector<ValueType>& ajiMatData = ajiMatrix.data();
    for (std::size_t i = 0; i < ajiMatrix.rows(); i++) {
        fmt::print(fpx, "{}{}", pfdata.refQuerySet()[i],
                   appArgs.outFieldSeparator);
        fmt::print(fpx, "{}\n",
                   fmt::join(ajiMatData.begin() + i * nTgtGeneomes,
                             ajiMatData.begin() + (i + 1) * nTgtGeneomes,
                             appArgs.outFieldSeparator));
    }
    fclose(fpx);
}

int parallel_fastaai(const AppParams& appArgs) {
    // Initialize databases
    SQLiteDB sqltIf(appArgs.pathToDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    //
    PFData pfaaiData(sqltIf, sqltIf.getMeta());
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    sqltIf.closeDB();
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiData);
    pfaaiImpl.run();
#ifndef NDEBUG
    pfaaiImpl.print_aji();
#endif
    if (appArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, appArgs);
    }
    return PFAAI_OK;
}

int validate_subset(const AppParams& appArgs, const DBMetaData& dbMeta) {
    std::unordered_set<std::string> dbGenomeSet(dbMeta.genomeSet.begin(),
                                                dbMeta.genomeSet.end());

    std::vector<std::string> missingGenomes;
    for (auto& genome : appArgs.qryGenomeSet) {
        if (dbGenomeSet.find(genome) == dbGenomeSet.end()) {
            missingGenomes.push_back(genome);
        }
    }
    if (missingGenomes.size() > 0) {
        std::cout << "--------------------ERROR-----------------------------"
                  << std::endl
                  << " In the query subset file, the following genomes are "
                  << std::endl
                  << " missing from the database : " << std::endl
                  << fmt::format("    {}",
                                 fmt::join(missingGenomes.begin(),
                                           missingGenomes.end(), "\n    "))
                  << std::endl
                  << std::endl
                  << " Please remove them and before running Fast AAI. "
                  << std::endl
                  << "------------------------------------------------------"
                  << std::endl;
        return PFAAI_ERR_CONSTRUCT;
    }
    return PFAAI_OK;
}

int parallel_subset_fastaai(const AppParams& appArgs) {
    // Initialize databases
    SQLiteDB sqltIf(appArgs.pathToDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    // Validate Input
    const DBMetaData& dbMeta = sqltIf.getMeta();
    if (validate_subset(appArgs, dbMeta) != PFAAI_OK) {
        return PFAAI_ERR_CONSTRUCT;
    }

    // PHASE 1: Construction of the data structures
    PFQSubData pfaaiData(sqltIf, dbMeta, appArgs.qryGenomeSet);
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    sqltIf.closeDB();
    //
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiData);
    pfaaiImpl.run();
#ifndef NDEBUG
    // pfaaiImpl.print_e();
    pfaaiImpl.print_aji();
#endif
    if (appArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, appArgs);
    }
    return PFAAI_OK;
}

int validate_qry2tgt(const DBMetaData& dbMeta) {
    std::unordered_set<std::string> dbGenomeSet(dbMeta.genomeSet.begin(),
                                                dbMeta.genomeSet.end());

    std::vector<std::string> commonGenomes;
    for (auto& genome : dbMeta.qyGenomeSet) {
        if (dbGenomeSet.find(genome) != dbGenomeSet.end()) {
            commonGenomes.push_back(genome);
        }
    }
    if (commonGenomes.size() > 0) {
        std::cout << "---------------------ERROR------------------------------"
                  << std::endl
                  << "  Query database should have no intersecting genes with "
                  << std::endl
                  << "  the main data base." << std::endl
                  << std::endl
                  << "  In the query data base, the following genomes are "
                  << std::endl
                  << "  overlapping with the target database:" << std::endl
                  << fmt::format("    {}",
                                 fmt::join(commonGenomes.begin(),
                                           commonGenomes.end(), "\n    "))
                  << std::endl
                  << std::endl
                  << "  Please remove them before running Fast AAI."
                  << std::endl
                  << "--------------------------------------------------------"
                  << std::endl;
        return PFAAI_ERR_CONSTRUCT;
    }
    return PFAAI_OK;
}

int parallel_qry2tgt_fastaai(const AppParams& appArgs) {
    //
    // Initialize databases and Meta data query
    QTSQLiteDB qtDBIf(appArgs.pathToDatabase, appArgs.pathToQryDatabase);
    PFAAI_ERROR_CODE qErrCode = qtDBIf.validate();
    if (qErrCode != PFAAI_OK) {
        return qErrCode;
    }
    //
    const DBMetaData& dbMeta = qtDBIf.getMeta();
    if (validate_qry2tgt(dbMeta) != PFAAI_OK) {
        return PFAAI_ERR_CONSTRUCT;
    }
    // PHASE 1: Construction of the data structures
    PFQTData pfaaiQTData(qtDBIf, dbMeta, dbMeta.proteinSet);
    PFAAI_ERROR_CODE pfErrorCode = pfaaiQTData.construct();
    qtDBIf.closeDB();
    if (pfErrorCode != PFAAI_OK) {
        return pfErrorCode;
    }

    // PHASE 2: Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiQTData);
    pfaaiImpl.run();
#ifndef NDEBUG
    pfaaiImpl.print_aji();
#endif
    if (appArgs.pathToOutputFile.size() > 0) {
        // TODO(x)::
        printOutput(pfaaiQTData, pfaaiImpl, appArgs, false);
    }

    return PFAAI_OK;
}

int main(int argc, char* argv[]) {
    AppParams pfaaiAppArgs;
    CLI11_PARSE(pfaaiAppArgs.app, argc, argv);
    pfaaiAppArgs.print();
    //
    if (pfaaiAppArgs.pathToQryDatabase.size() == 0 ||
        pfaaiAppArgs.pathToQryDatabase == pfaaiAppArgs.pathToDatabase) {
        if (pfaaiAppArgs.pathToQrySubsetFile.size() == 0) {
            return parallel_fastaai(pfaaiAppArgs);
        } else {
            pfaaiAppArgs.load_query_genomes();
#ifndef NDEBUG
            pfaaiAppArgs.print_query_genomes();
#endif
            return parallel_subset_fastaai(pfaaiAppArgs);
        }
    } else {
        return parallel_qry2tgt_fastaai(pfaaiAppArgs);
    }
}
