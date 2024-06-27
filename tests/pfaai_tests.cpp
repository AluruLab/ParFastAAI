
#include "pfaai/sqltif.hpp"
#include "pfaai_tests.hpp"
#include <catch2/catch_all.hpp>

TEST_CASE( "META DATA TEST", "[meta data test]" ) {
    const std::string db_path = "data/modified_xantho_fastaai2.db";
    std::vector<std::string> testProtienSet = TESTDB_PROTEIN_SET;
    SQLiteInterface<DatabaseNames> sqlt_if(db_path);
    std::vector<std::string> protienSet;
    int nGenomes;
    sqlt_if.queryMetaData(protienSet, nGenomes);
    REQUIRE(nGenomes == 20);
    REQUIRE(testProtienSet == protienSet);
}
