
#include "pfaai_tests.hpp"
#include "catch2/catch_test_macros.hpp"
#include "pfaai/data.hpp"
#include "pfaai/sqltif.hpp"
#include "pfaai/runner.hpp"
#include <catch2/catch_all.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>

using IdType = int;
using IdPairType = IDPair<IdType, IdType>;
// using IdMatType = std::vector<std::vector<IdType>>;
using IdMatType = IDMatrix<IdType>;

static constexpr char G_DB_PATH[] = "data/modified_xantho_fastaai2.db";
static constexpr IdType G_NTETRAMERS = (20 * 20 * 20 * 20);
static constexpr IdType G_GENOMESET_SIZE = 20;
static constexpr IdType G_PROTIENSET_SIZE = 80;
static const std::vector<std::string> G_PROTIENSET = TESTDB_PROTEIN_SET;
static const std::vector<std::string> G_GENOMESET = TESTDB_GENOME_SET;
// Testing of genomes query for protien
static const std::vector<std::vector<IdType>> G_QRY_PST_TETRA_CTS = {
    {260, 260, 260, 260, 260, 260, 260, 260, 260, 260,
     260, 260, 260, 260, 260, 260, 260, 263, 263, 263},
    {244, 244, 244, 244, 244, 244, 244, 244, 244, 244,
     244, 244, 244, 244, 244, 244, 244, 246, 246, 246},
    {386, 386, 386, 386, 386, 386, 386, 386, 386, 386,
     386, 386, 386, 386, 386, 386, 386, 383, 384, 384},
    {121, 121, 121, 121, 121, 121, 121, 121, 121, 121,
     121, 121, 121, 121, 121, 121, 121, 121, 121, 121}};
//

// TETRAMER, # genomes, source table
static constexpr IdType G_QRT_PST_GP_PAIRS[388][3] = TESTDB_PSET_GP_PAIRS;

TEST_CASE("Query Genome Metadata", "[meta data test]") {
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    std::vector<std::string> protienSet, genomeSet;
    IdType nGenomes;
    sqlt_if.queryMetaData(protienSet, genomeSet);
    REQUIRE(G_PROTIENSET == protienSet);
    REQUIRE(G_GENOMESET == genomeSet);
}

TEST_CASE("Query No. of Tetramers", "[ntetramers query test]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0);
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    sqlt_if.queryGenomeTetramers("PF00119.20", 2060, 2144, Lc);
    REQUIRE(Lc[2060] == 20);
    REQUIRE(Lc[2100] == 20);
    REQUIRE(Lc[2140] == 3);
    REQUIRE(Lc[2144] == 17);
}

TEST_CASE("Query Protien Tetramer Counts", "[prot teteramer counts]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    IdMatType T(G_PROTIENSET_SIZE, G_GENOMESET_SIZE);
    // IdMatType T(G_PROTIENSET.size());
    // for (auto& tx : T) {
    //     tx.resize(G_NGENOMES);
    // }
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    sqlt_if.queryProtienTetramerCounts(G_PROTIENSET, 0, 3, T);
    REQUIRE(G_QRY_PST_TETRA_CTS[0] == T.row(0));
    REQUIRE(G_QRY_PST_TETRA_CTS[1] == T.row(1));
    REQUIRE(G_QRY_PST_TETRA_CTS[2] == T.row(2));
    REQUIRE(G_QRY_PST_TETRA_CTS[3] == T.row(3));
}

TEST_CASE("Query Protien Set Tetramers", "[prot set teteramers]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    for (auto& rx : G_QRT_PST_GP_PAIRS) {
        Lc[rx[0]] += rx[1] / sizeof(IdType);
    }
    //
    // Parallel prefix sum on Lc to construct Lp
    int cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan, + : cumulativeSum)
    for (int i = 0; i < Lc.size(); i++) {
        Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
        cumulativeSum += Lc[i];
    }

    std::vector<IdPairType> F(Lp[G_NTETRAMERS - 1] + Lc[G_NTETRAMERS - 1],
                              IdPairType(-1, -1));
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    int rc = sqlt_if.queryProtienSetGPPairs(G_PROTIENSET, 2000, 3000, Lp, F);
    REQUIRE(rc == SQLITE_OK);
    REQUIRE(Lp[2000] == 0);
    REQUIRE(F[Lp[2000]] == IdPairType(5, 0));
    REQUIRE(F[Lp[2415] - 1] == IdPairType(35, 0x13));
    REQUIRE(F[Lp[2415]] == IdPairType(74, 0x11));
}

TEST_CASE("Test Data Structures Construction", "[construct Lc Lp F T]") {
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    std::vector<std::string> proteinSet;
    std::vector<std::string> genomeSet;
    sqlt_if.queryMetaData(proteinSet, genomeSet);
    //
    ParFAAIData<IdType, IdPairType, IdMatType> pfaaiData(sqlt_if, proteinSet,
                                                         genomeSet);
    pfaaiData.constructLcandLp();
    std::vector<IdType> Lc = pfaaiData.getLc();
    std::vector<IdType> Lp = pfaaiData.getLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx("data/xanthodb_lc_array.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx); 
        iarchive(refLc);
    }
    {
        std::ifstream lcx("data/xanthodb_lp_array.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx); 
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.getF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx("data/xanthodb_f_array.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.getT();
    IdMatType refT(G_PROTIENSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx("data/xanthodb_t_matrix.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
}

TEST_CASE("Test Runner ouputs", "[construct E JAC and AJI]") {
    SQLiteInterface<IdType, IdPairType, IdMatType, DatabaseNames> sqlt_if(
        G_DB_PATH);
    std::vector<std::string> proteinSet;
    std::vector<std::string> genomeSet;
    sqlt_if.queryMetaData(proteinSet, genomeSet);
    //
    ParFAAIData<IdType, IdPairType, IdMatType> pfaaiData(sqlt_if, proteinSet,
                                                         genomeSet);
    pfaaiData.construct();
    //
    ParFAAIRunner<IdType, IdPairType, IdMatType, double> pfaaiRunner(pfaaiData);

    pfaaiRunner.generateTetramerTuples();
    std::vector<ETriple<IdType>> E = pfaaiRunner.getE();
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx("data/xanthodb_sorted_e_array.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    //
    pfaaiRunner.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiRunner.getJAC(), refJAC;
    {
        std::ifstream tcx("data/xanthodb_jac.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC.size() == JAC.size());
    for(auto i = 0; i < JAC.size(); i++)
        REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiRunner.computeAJI();
    std::vector<double> AJI = pfaaiRunner.getAJI(), refAJI;
    {
        std::ifstream tcx("data/xanthodb_aji.bin", std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    for(auto i = 0; i < refAJI.size(); i++)
        REQUIRE(refAJI[i] == AJI[i]);
}
