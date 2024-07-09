/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <vector>

#include "pfaai/data.hpp"
#include "pfaai/sqltif.hpp"
#include "pfaai/runner.hpp"

using IdType = int;
using IdPairType = IDPair<int, int>;
using IdMatrixType = IDMatrix<int>;

int parallel_fastaai(const std::string pathToDatabase) {
    SQLiteInterface<IdType, IdPairType, IdMatrixType, DatabaseNames> sqltIf(pathToDatabase);
    PFAAI_ERROR_CODE errorCode = sqltIf.validate();
    if (errorCode != PFAAI_OK) {
        return errorCode;
    }
    std::vector<std::string> proteinSet;
    std::vector<std::string> genomeSet;
    sqltIf.queryMetaData(proteinSet, genomeSet);
    //
    ParFAAIData<IdType, IdPairType, IdMatrixType> pfaaiData(sqltIf, proteinSet, genomeSet);
    // PHASE 1: Construction of the data structures
    PFAAI_ERROR_CODE pfErrorCode = pfaaiData.construct();
    if (errorCode != PFAAI_OK) {
        return pfErrorCode;
    }
    // Run parallel Fast AAI algorithm
    ParFAAIRunner<IdType, IdPairType, IdMatrixType> pfaaiRunner(pfaaiData);
    pfaaiRunner.run();
    pfaaiRunner.print_aji();
    return PFAAI_OK;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Missing Database location " << std::endl;
        std::cout << "Usage: " << argv[0] << " Database Location " << std::endl;
        return 1;
    }
    // const std::string pathToDatabse = "modified_xantho_fastaai2.db";
    const std::string pathToDatabase = argv[1];
    return parallel_fastaai(pathToDatabase);
}
