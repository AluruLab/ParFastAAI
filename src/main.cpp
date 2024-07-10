/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <vector>

#include "pfaai/data.hpp"
#include "pfaai/impl.hpp"
#include "pfaai/sqltif.hpp"
#include <CLI/CLI.hpp>
#include "fmt/format.h"

using IdType = int;
using ValueType = double;
using IdPairType = IDPair<IdType, IdType>;
using IdMatrixType = DMatrix<IdType>;
using SQLiteIfT =
    SQLiteInterface<IdType, IdPairType, IdMatrixType, DatabaseNames>;
using PFImplT = ParFAAIImpl<IdType, IdPairType, IdMatrixType, ValueType>;
using PFDataT = ParFAAIData<IdType, IdPairType, IdMatrixType, PFImplT::JACType>;

void print_aji(const std::vector<PFImplT::JACType>& cJAC,
               const std::vector<ValueType>& cAJI) {
    std::cout << "AJI Ouput : " << std::endl;
    std::cout << " [(GP1, GP2,   SUM, NCP) ->  AJI]" << std::endl;
    for (int i = 0; i < cAJI.size(); i++) {
        fmt::print(" [{} -> {:03.2f}] \n", cJAC[i], cAJI[i]);
    }
}

void printOutput(const PFDataT& pfdata, const PFImplT& impl,
                 const std::string pathToOutputFile) {
    IdType nQryGenomes = pfdata.getQryGenomeCount();
    IdType nTgtGeneomes = pfdata.getTgtGenomeCount();
    //
    const std::vector<PFImplT::JACType>& cJAC = impl.getJAC();
    const std::vector<ValueType>& cAJI = impl.getAJI();

    DMatrix<ValueType> ajiMatrix(nQryGenomes, nTgtGeneomes, 0.0);
    for (IdType i = 0; i < cJAC.size(); i++) {
        ajiMatrix(cJAC[i].genomeA, cJAC[i].genomeB) = cAJI[i];
        ajiMatrix(cJAC[i].genomeB, cJAC[i].genomeA) = cAJI[i];
    }

    //std::ofstream ofx(pathToOutputFile);
    FILE *fpx = fopen(pathToOutputFile.c_str(), "w");
    const std::vector<ValueType>& ajiMatData = ajiMatrix.data();
    for (IdType i = 0; i < ajiMatrix.rows(); i++) {
        fmt::print(fpx, "{}\n",
                   fmt::join(ajiMatData.begin() + i * nTgtGeneomes,
                             ajiMatData.begin() + (i + 1) * nTgtGeneomes, ","));
    }
    fclose(fpx);
}

int parallel_fastaai(const std::string pathToDatabase,
                     const std::string pathToOutputFile) {
    // Initialize databases
    SQLiteIfT sqltIf(pathToDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    std::vector<std::string> proteinSet;
    std::vector<std::string> genomeSet;
    sqltIf.queryMetaData(proteinSet, genomeSet);
    //
    PFDataT pfaaiData(sqltIf, proteinSet, genomeSet);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    sqltIf.closeDB();
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    PFImplT pfaaiImpl(pfaaiData);
    pfaaiImpl.run();
    print_aji(pfaaiImpl.getJAC(), pfaaiImpl.getAJI());
    if(pathToOutputFile.size() > 0){
        printOutput(pfaaiData, pfaaiImpl, pathToOutputFile);
    }
    return PFAAI_OK;
}

int main(int argc, char* argv[]) {
    CLI::App app;
    std::string pathToDatabase, pathToOutputFile("");

    app.add_option("path_to_db", pathToDatabase, "Path to the Input Database")
        ->check(CLI::ExistingFile);
    app.add_option("-o,--output_file", pathToOutputFile, "Path to output csv file");
    CLI11_PARSE(app, argc, argv);
    //
    // Display arguments
    std::string arg1 = fmt::format(" Input Database : {} ", pathToDatabase);
    std::string arg2 = fmt::format(" Output File    : {} ", pathToOutputFile);
    std::size_t argsz = std::max(arg1.size(), arg2.size());
    fmt::print(" ┌{0:─^{3}}┐\n"
               " │{1: <{3}}│\n"
               " │{2: <{3}}│\n"
               " └{0:─^{3}}┘\n",
               "", arg1, arg2, argsz);
    return parallel_fastaai(pathToDatabase, pathToOutputFile);
}
