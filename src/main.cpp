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

#include "fmt/core.h"
#include "pfaai/algorithm_impl.hpp"
#include "pfaai/data_impl.hpp"
#include "pfaai/database.hpp"
#include "pfaai/interface.hpp"

using IdType = int;
using ValueType = double;

using IdPairType = DPair<IdType, IdType>;
using IdMatrixType = DMatrix<IdType>;
using SQLiteIfT = SQLiteInterface<IdType, DatabaseNames>;
using QTSQLiteIfT = QTSQLiteInterface<IdType, DatabaseNames>;
using PFImpl = ParFAAIImpl<IdType, ValueType>;

using PFDSInterface = DefaultDataStructInterface<IdType>;
using PFData = ParFAAIData<IdType>;
using PFQSubData = ParFAAIQSubData<IdType>;
using PFQTData = ParFAAIQryTgtData<IdType>;

struct AppParams {
    CLI::App app;
    std::string pathToQryDatabase;
    std::string pathToTgtDatabase;
    std::string pathToQrySubsetFile;
    std::string pathToOutputFile;
    std::string outFieldSeparator;
    std::vector<std::string> qryGenomeSet;

    AppParams()
        : app(), pathToQryDatabase(""), pathToTgtDatabase(""),
          pathToQrySubsetFile(""), pathToOutputFile(""),
          outFieldSeparator(",") {
        app.add_option("path_to_query_db", pathToQryDatabase,
                       "Path to the Query Database")
            ->check(CLI::ExistingFile);
        app.add_option("path_to_output_file", pathToOutputFile,
                       "Path to output csv file.");
        app.add_option("-t,--target", pathToTgtDatabase,
                       "Path to the Target Database [Optional (default: Same "
                       "as the Query DB)] *NOT IMPLEMENTED YET*")
            ->check(CLI::ExistingFile);
        app.add_option(
               "-s,--separator", outFieldSeparator,
               "Field Separator in the output file [Optional (default: ,)].")
            ->capture_default_str();
        app.add_option(
               "-q,--query_subset", pathToQrySubsetFile,
               "Path to Query List (Should be subset of genomoes in input DB.)")
            ->check(CLI::ExistingFile);
    }

    void print() const {
        // Display arguments
        std::string arg1 =
            fmt::format(" Query Database  : {} ", pathToQryDatabase);
        std::string arg2 =
            fmt::format(" Target Database : {} ", pathToTgtDatabase);
        std::string arg3 =
            fmt::format(" Output File     : {} ", pathToOutputFile);
        std::string arg4 =
            fmt::format(" Field Separator : {} ", outFieldSeparator);
        std::string arg5 =
            fmt::format(" Query Subset    : {} ", pathToQrySubsetFile);
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

void print_aji(const PFDSInterface& dsif,
               const std::vector<PFImpl::JACType>& cJAC,
               const std::vector<ValueType>& cAJI) {
    std::cout << "AJI Ouput : " << std::endl;
    std::cout << " [(GP1, GP2,   SUM, NCP) ->  AJI]" << std::endl;
    for (std::size_t i = 0; i < cAJI.size(); i++) {
        fmt::print(" [{}, {} -> {:03.2f}] \n",
                   dsif.genomePairToIndex(cJAC[i].genomeA, cJAC[i].genomeB),
                   cJAC[i], cAJI[i]);
    }
}

void printOutput(const PFDSInterface& pfdata, const PFImpl& impl,
                 const AppParams& pfaaiAppArgs) {
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
        if (pfdata.isQryGenome(cJAC[i].genomeB)) {
            ajiMatrix(pfdata.mapQueryId(fgb), pfdata.mapTargetId(fga)) =
                cAJI[i];
        }
    }

    // std::ofstream ofx(pathToOutputFile);
    FILE* fpx = fopen(pfaaiAppArgs.pathToOutputFile.c_str(), "w");

    const auto& tgtRef = pfdata.refTargetSet();
    // Header
    fmt::print(fpx, "{}{}\n", pfaaiAppArgs.outFieldSeparator,
               fmt::join(tgtRef.begin(), tgtRef.end(),
                         pfaaiAppArgs.outFieldSeparator));

    // Data
    const std::vector<ValueType>& ajiMatData = ajiMatrix.data();
    for (std::size_t i = 0; i < ajiMatrix.rows(); i++) {
        fmt::print(fpx, "{}{}", pfdata.refQuerySet()[i],
                   pfaaiAppArgs.outFieldSeparator);
        fmt::print(fpx, "{}\n",
                   fmt::join(ajiMatData.begin() + i * nTgtGeneomes,
                             ajiMatData.begin() + (i + 1) * nTgtGeneomes,
                             pfaaiAppArgs.outFieldSeparator));
    }
    fclose(fpx);
}

int parallel_fastaai(const AppParams& pfaaiAppArgs) {
    // Initialize databases
    SQLiteIfT sqltIf(pfaaiAppArgs.pathToQryDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    const DBMetaData& dbMeta = sqltIf.getMeta();
    //
    PFData pfaaiData(sqltIf, dbMeta);
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
    print_aji(pfaaiData, pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, pfaaiAppArgs);
    }
    return PFAAI_OK;
}

int parallel_subset_fastaai(const AppParams& pfaaiAppArgs) {
    // Initialize databases
    SQLiteIfT sqltIf(pfaaiAppArgs.pathToQryDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    const DBMetaData& dbMeta = sqltIf.getMeta();
    //
    PFQSubData pfaaiData(sqltIf, dbMeta, pfaaiAppArgs.qryGenomeSet);
    // PHASE 1: Construction of the data structures
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
    print_aji(pfaaiData, pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, pfaaiAppArgs);
    }
    return PFAAI_OK;
}

int parallel_qry2tgt_fastaai(const AppParams& pfaaiAppArgs) {
    //
    std::cout << "Implementation NOT Complete " << std::endl;
    return 0;
    // Initialize databases
    QTSQLiteIfT qtDBIf(pfaaiAppArgs.pathToQryDatabase,
                       pfaaiAppArgs.pathToTgtDatabase);
    PFAAI_ERROR_CODE qErrCode = qtDBIf.validate();
    if (qErrCode != PFAAI_OK) {
        return qErrCode;
    }
    // Meta data query
    const DBMetaData& dbMeta = qtDBIf.getMeta();
    //
    std::unordered_set<std::string> unionProtiens;
    std::vector<std::string> sharedProtiens;

    // TODO(x):: get shared proteins
    //
    // for (const auto& sx : qryDbMeta.proteinSet) {
    //     unionProtiens.insert(sx);
    // }
    // for (const auto& sx : tgtDbMeta.proteinSet) {
    //     if (unionProtiens.find(sx) != unionProtiens.end()) {
    //         sharedProtiens.emplace_back(sx);
    //     } else {
    //         unionProtiens.insert(sx);
    //     }
    // }

    PFQTData pfaaiQTData(qtDBIf, dbMeta, sharedProtiens);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiQTData.construct();
    qtDBIf.closeDB();
    if (pfErrorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiQTData);
    pfaaiImpl.run();
#ifndef NDEBUG
    print_aji(pfaaiQTData, pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        // TODO(x)::
        printOutput(pfaaiQTData, pfaaiImpl, pfaaiAppArgs);
    }

    return PFAAI_OK;
}

int main(int argc, char* argv[]) {
    AppParams pfaaiAppArgs;
    CLI11_PARSE(pfaaiAppArgs.app, argc, argv);
    pfaaiAppArgs.print();
    //
    if (pfaaiAppArgs.pathToTgtDatabase.size() == 0 ||
        pfaaiAppArgs.pathToTgtDatabase == pfaaiAppArgs.pathToQryDatabase) {
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
