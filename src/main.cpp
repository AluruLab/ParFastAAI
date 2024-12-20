/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <unordered_set>
#include <vector>

#include "fmt/format.h"
#include "pfaai/algorithm_impl.hpp"
#include "pfaai/data_impl.hpp"
#include "pfaai/database.hpp"
#include <CLI/CLI.hpp>

using IdType = int;
using ValueType = double;

using IdPairType = DPair<IdType, IdType>;
using IdMatrixType = DMatrix<IdType>;
using SQLiteIfT = SQLiteInterface<IdType, DatabaseNames>;
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
        app.add_option("-t,--target", pathToTgtDatabase,
                       "Path to the Target Database [Optional, Same as Query, "
                       "if not given]")
            ->check(CLI::ExistingFile);
        app.add_option("-s,--separator", outFieldSeparator,
                       "Field Separator in the output file")
            ->capture_default_str();
        app.add_option("-q,--query_subset", pathToQrySubsetFile,
                       "Path to output csv file")
            ->check(CLI::ExistingFile);
        app.add_option("-o,--output_file", pathToOutputFile,
                       "Path to output csv file");
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
};

void print_aji(const std::vector<PFImpl::JACType>& cJAC,
               const std::vector<ValueType>& cAJI) {
    std::cout << "AJI Ouput : " << std::endl;
    std::cout << " [(GP1, GP2,   SUM, NCP) ->  AJI]" << std::endl;
    for (std::size_t i = 0; i < cAJI.size(); i++) {
        fmt::print(" [{} -> {:03.2f}] \n", cJAC[i], cAJI[i]);
    }
}

void printOutput(const PFDSInterface& pfdata, const PFImpl& impl,
                 const AppParams& pfaaiAppArgs) {
    IdType nQryGenomes = pfdata.qrySetSize();
    IdType nTgtGeneomes = pfdata.tgtSetSize();
    //
    const std::vector<PFImpl::JACType>& cJAC = impl.getJAC();
    const std::vector<ValueType>& cAJI = impl.getAJI();

    DMatrix<ValueType> ajiMatrix(nQryGenomes, nTgtGeneomes, 0.0);
    for (std::size_t i = 0; i < cJAC.size(); i++) {
        ajiMatrix(cJAC[i].genomeA, cJAC[i].genomeB) = cAJI[i];
        ajiMatrix(cJAC[i].genomeB, cJAC[i].genomeA) = cAJI[i];
    }

    // std::ofstream ofx(pathToOutputFile);
    FILE* fpx = fopen(pfaaiAppArgs.pathToOutputFile.c_str(), "w");
    const std::vector<ValueType>& ajiMatData = ajiMatrix.data();
    for (std::size_t i = 0; i < ajiMatrix.rows(); i++) {
        fmt::print(fpx, "{}\n",
                   fmt::join(ajiMatData.begin() + i * nTgtGeneomes,
                             ajiMatData.begin() + (i + 1) * nTgtGeneomes,
                             pfaaiAppArgs.outFieldSeparator));
    }
    fclose(fpx);
}

void printQTOutput(const PFDSInterface& pfdata, const PFImpl& impl,
                   const AppParams& pfaaiAppArgs) {
    IdType nQryGenomes = pfdata.qrySetSize();
    IdType nTgtGeneomes = pfdata.tgtSetSize();
    //
    const std::vector<PFImpl::JACType>& cJAC = impl.getJAC();
    const std::vector<ValueType>& cAJI = impl.getAJI();

    DMatrix<ValueType> ajiMatrix(nQryGenomes, nTgtGeneomes, 0.0);
    for (std::size_t i = 0; i < cJAC.size(); i++) {
        ajiMatrix(cJAC[i].genomeA, cJAC[i].genomeB - nTgtGeneomes) = cAJI[i];
    }

    // std::ofstream ofx(pathToOutputFile);
    FILE* fpx = fopen(pfaaiAppArgs.pathToOutputFile.c_str(), "w");
    const std::vector<ValueType>& ajiMatData = ajiMatrix.data();
    for (std::size_t i = 0; i < ajiMatrix.rows(); i++) {
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
    DBMetaData dbMeta;
    sqltIf.queryMetaData(dbMeta);
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
    print_aji(pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
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
    DBMetaData dbMeta;
    sqltIf.queryMetaData(dbMeta);
    std::cout << "{" << dbMeta.genomeSet.size() << ", "
              << fmt::format("[{}]", fmt::join(dbMeta.genomeSet, ", ")) << "}"
              << std::endl;
    //
    PFQSubData pfaaiData(sqltIf, dbMeta, pfaaiAppArgs.qryGenomeSet,
                         dbMeta.proteinSet);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    sqltIf.closeDB();
    std::cout << "{" << pfaaiData.refF().size() << ", "
              << pfaaiData.refE().E.size() << "}" << std::endl;
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiData);
    pfaaiImpl.run();
    return 0;
#ifndef NDEBUG
    print_aji(pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, pfaaiAppArgs);
    }
    return PFAAI_OK;
}

int parallel_qry2tgt_fastaai(const AppParams& pfaaiAppArgs) {
    // Initialize databases
    SQLiteIfT qryDBIf(pfaaiAppArgs.pathToQryDatabase);
    PFAAI_ERROR_CODE qErrCode = qryDBIf.validate();
    if (qErrCode != PFAAI_OK) {
        return qErrCode;
    }
    SQLiteIfT tgtDBIf(pfaaiAppArgs.pathToQryDatabase);
    PFAAI_ERROR_CODE tgtErrCode = tgtDBIf.validate();
    if (tgtErrCode != PFAAI_OK) {
        return tgtErrCode;
    }
    // Meta data query
    DBMetaData qryDbMeta, tgtDbMeta;
    qryDBIf.queryMetaData(qryDbMeta);
    tgtDBIf.queryMetaData(tgtDbMeta);
    //
    std::unordered_set<std::string> unionProtiens;
    std::vector<std::string> sharedProtiens;
    for (const auto& sx : qryDbMeta.proteinSet) {
        unionProtiens.insert(sx);
    }
    for (const auto& sx : tgtDbMeta.proteinSet) {
        if (unionProtiens.find(sx) != unionProtiens.end()) {
            sharedProtiens.emplace_back(sx);
        } else {
            unionProtiens.insert(sx);
        }
    }

    PFQTData pfaaiQTData(qryDBIf, tgtDBIf, qryDbMeta, tgtDbMeta,
                         sharedProtiens);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiQTData.construct();
    qryDBIf.closeDB();
    tgtDBIf.closeDB();
    if (pfErrorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiQTData);
    pfaaiImpl.run();
#ifndef NDEBUG
    print_aji(pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        // TODO(x)::
        printQTOutput(pfaaiQTData, pfaaiImpl, pfaaiAppArgs);
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
            std::cout << "{" << pfaaiAppArgs.qryGenomeSet.size() << ", "
                      << fmt::format("[{}]",
                                     fmt::join(pfaaiAppArgs.qryGenomeSet, ", "))
                      << "}" << std::endl;
            return parallel_subset_fastaai(pfaaiAppArgs);
        }
    } else {
        return parallel_qry2tgt_fastaai(pfaaiAppArgs);
    }
}
