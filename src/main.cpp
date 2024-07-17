/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <vector>

#include "fmt/format.h"
#include "pfaai/data.hpp"
#include "pfaai/impl.hpp"
#include "pfaai/sqltif.hpp"
#include <CLI/CLI.hpp>


using IdType = int;
using ValueType = double;
using IdPairType = DPair<IdType, IdType>;
using IdMatrixType = DMatrix<IdType>;
using SQLiteIfT =
    SQLiteInterface<IdType, IdPairType, IdMatrixType, DatabaseNames>;
using PFImplT = ParFAAIImpl<IdType, IdPairType, IdMatrixType, ValueType>;
using PFDataT = ParFAAIData<IdType, IdPairType, IdMatrixType, PFImplT::JACType>;

struct AppParams {
    CLI::App app;
    std::string pathToDatabase;
    std::string pathToOutputFile;
    std::string outFieldSeparator;

    AppParams()
        : app(), pathToDatabase(""), pathToOutputFile(""),
          outFieldSeparator(",") {
        app.add_option("path_to_db", pathToDatabase,
                       "Path to the Input Database")
            ->check(CLI::ExistingFile);
        app.add_option("-s,--separator", outFieldSeparator,
                       "Field Separator in the output file")
            ->capture_default_str();
        app.add_option("-o,--output_file", pathToOutputFile,
                       "Path to output csv file");
    }

    void print() const {
        // Display arguments
        std::string arg1 = fmt::format(" Input Database  : {} ", pathToDatabase);
        std::string arg2 = fmt::format(" Output File     : {} ", pathToOutputFile);
        std::string arg3 = fmt::format(" Field Separator : {} ", outFieldSeparator);
        std::size_t argsz = std::max({arg1.size(), arg2.size(), arg3.size()});
        fmt::print(" ┌{0:─^{1}}┐\n"
                   " │{2: <{1}}│\n"
                   " │{3: <{1}}│\n"
                   " │{4: <{1}}│\n"
                   " └{0:─^{1}}┘\n",
                   "", argsz, arg1, arg2, arg3);
    }
};



void print_aji(const std::vector<PFImplT::JACType>& cJAC,
               const std::vector<ValueType>& cAJI) {
    std::cout << "AJI Ouput : " << std::endl;
    std::cout << " [(GP1, GP2,   SUM, NCP) ->  AJI]" << std::endl;
    for (std::size_t i = 0; i < cAJI.size(); i++) {
        fmt::print(" [{} -> {:03.2f}] \n", cJAC[i], cAJI[i]);
    }
}

void printOutput(const PFDataT& pfdata, const PFImplT& impl,
                 const AppParams& pfaaiAppArgs) {
    IdType nQryGenomes = pfdata.getQryGenomeCount();
    IdType nTgtGeneomes = pfdata.getTgtGenomeCount();
    //
    const std::vector<PFImplT::JACType>& cJAC = impl.getJAC();
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

int parallel_fastaai(const AppParams& pfaaiAppArgs) {
    // Initialize databases
    SQLiteIfT sqltIf(pfaaiAppArgs.pathToDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    DBMetaData dbMeta;
    sqltIf.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqltIf, dbMeta);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    sqltIf.closeDB();
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImplT pfaaiImpl(pfaaiData);
    pfaaiImpl.run();
#ifndef NDEBUG
    print_aji(pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
#endif
    if (pfaaiAppArgs.pathToOutputFile.size() > 0) {
        printOutput(pfaaiData, pfaaiImpl, pfaaiAppArgs);
    }
    return PFAAI_OK;
}


int main(int argc, char* argv[]) {

    AppParams pfaaiAppArgs;
    CLI11_PARSE(pfaaiAppArgs.app, argc, argv);
    pfaaiAppArgs.print();
    //
    return parallel_fastaai(pfaaiAppArgs);
}
