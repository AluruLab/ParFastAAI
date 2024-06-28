
#include "pfaai_tests.hpp"
#include "catch2/catch_test_macros.hpp"
#include "pfaai/sqltif.hpp"
#include <catch2/catch_all.hpp>

template <typename IT> struct GPPair {
    IT genome_id;
    IT prot_id;

    explicit GPPair(int g, int p) {
        genome_id = g;
        prot_id = p;
    }
    bool operator==(const GPPair& other) const {
        return genome_id == other.genome_id && prot_id == other.prot_id;
    }
};
template<typename IT>
std::ostream& operator<<(std::ostream& ox, GPPair<IT> const& cx) {
    ox << "{" << cx.genome_id << ", " << cx.prot_id << "}";
    return ox;
}

using IdType = int;
//using IdPairType = std::pair<int, int>;
using IdPairType = GPPair<IdType>;
using IdMatType = std::vector<std::vector<IdType>>;

constexpr char G_DB_PATH[] = "data/modified_xantho_fastaai2.db";
constexpr static IdType G_NTETRAMERS = (20 * 20 * 20 * 20);
constexpr static IdType G_NGENOMES = 20;
constexpr IdType G_PROTIENSET_SIZE = 80;
const std::vector<std::string> G_PROTIENSET = TESTDB_PROTEIN_SET;
// Testing of genomes query for protien
const IdMatType G_QRY_PST_TETRA_CTS = {
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
constexpr IdType G_QRT_PST_GP_PAIRS[388][3] = TESTDB_PSET_GP_PAIRS;

TEST_CASE("Query Genome Metadata", "[meta data test]") {
    SQLiteInterface<DatabaseNames> sqlt_if(G_DB_PATH);
    std::vector<std::string> protienSet;
    IdType nGenomes;
    sqlt_if.queryMetaData(protienSet, nGenomes);
    REQUIRE(nGenomes == 20);
    REQUIRE(G_PROTIENSET == protienSet);
}

TEST_CASE("Query No. of Tetramers", "[ntetramers query test]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0);
    SQLiteInterface<DatabaseNames> sqlt_if(G_DB_PATH);
    sqlt_if.queryGenomeTetramers("PF00119.20", 2060, 2144, Lc);
    REQUIRE(Lc[2060] == 20);
    REQUIRE(Lc[2100] == 20);
    REQUIRE(Lc[2140] == 3);
    REQUIRE(Lc[2144] == 17);
}

TEST_CASE("Query Protien Tetramer Counts", "[prot teteramer counts]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    IdMatType T(G_PROTIENSET.size());
    for (auto& tx : T) {
        tx.resize(G_NGENOMES);
    }
    SQLiteInterface<DatabaseNames> sqlt_if(G_DB_PATH);
    sqlt_if.queryProtienTetramerCounts(G_PROTIENSET, 0, 3, T);
    REQUIRE(G_QRY_PST_TETRA_CTS[0] == T[0]);
    REQUIRE(G_QRY_PST_TETRA_CTS[1] == T[1]);
    REQUIRE(G_QRY_PST_TETRA_CTS[2] == T[2]);
    REQUIRE(G_QRY_PST_TETRA_CTS[3] == T[3]);
}

TEST_CASE("Query Protien Set Tetramers", "[prot set teteramers]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    for (auto& rx : G_QRT_PST_GP_PAIRS) {
        Lc[rx[0]] += rx[1]/sizeof(IdType);
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
    SQLiteInterface<DatabaseNames> sqlt_if(G_DB_PATH);
    int rc = sqlt_if.queryProtienSetGPPairs(G_PROTIENSET, 2000, 3000, Lp, F);
    REQUIRE(rc == SQLITE_OK);
    REQUIRE(Lp[2000] == 0);
    REQUIRE(F[Lp[2000]] == IdPairType(5, 0));
    REQUIRE(F[Lp[2415]-1] == IdPairType(35, 0x13));
    REQUIRE(F[Lp[2415]] == IdPairType(74, 0x11));
}
