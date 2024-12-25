
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include <catch2/catch_all.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <pfaai/algorithm_impl.hpp>
#include <pfaai/data_impl.hpp>
#include <pfaai/database.hpp>

#include "pfaai_tests.hpp"  // NOLINT

using IdType = int;
using ValueType = double;
using IdPairType = DPair<IdType, IdType>;
using IdMatType = DMatrix<IdType>;
using SQLiteIfT = SQLiteInterface<IdType, DatabaseNames>;
using PFImpl = ParFAAIImpl<IdType, ValueType>;
using PFDataT = ParFAAIData<IdType>;
using PFQSubDataT = ParFAAIQSubData<IdType>;
using PFQTData = ParFAAIQryTgtData<IdType>;

static constexpr char G_DB_PATH[] = "data/modified_xantho_fastaai2.db";
static constexpr char G_DB_SUBSET1_PATH[] = "data/xdb_subset1.db";
static constexpr char G_DB_SUBSET2_PATH[] = "data/xdb_subset2.db";
static constexpr IdType G_NTETRAMERS = (20 * 20 * 20 * 20);
static constexpr IdType G_GENOMESET_SIZE = 20;
static constexpr IdType G_PROTEINSET_SIZE = 80;
static const std::vector<std::string> G_PROTEINSET = TESTDB_PROTEIN_SET;
static const std::vector<std::string> G_GENOMESET = TESTDB_GENOME_SET;
// Testing of genomes query for protein
static const std::vector<std::vector<IdType>> G_QRY_PST_TETRA_CTS = {
    {260, 260, 260, 260, 260, 260, 260, 260, 260, 260,
     260, 260, 260, 260, 260, 260, 260, 263, 263, 263},
    {244, 244, 244, 244, 244, 244, 244, 244, 244, 244,
     244, 244, 244, 244, 244, 244, 244, 246, 246, 246},
    {386, 386, 386, 386, 386, 386, 386, 386, 386, 386,
     386, 386, 386, 386, 386, 386, 386, 383, 384, 384},
    {121, 121, 121, 121, 121, 121, 121, 121, 121, 121,
     121, 121, 121, 121, 121, 121, 121, 121, 121, 121}};

// TETRAMER, # genomes, source table
static constexpr IdType G_QRT_PST_GP_PAIRS[388][3] = TESTDB_PSET_GP_PAIRS;

// Reference data files to compare
static constexpr char REF_LC_ARRAY[] = "data/xanthodb_lc_array.bin";
static constexpr char REF_LP_ARRAY[] = "data/xanthodb_lp_array.bin";
static constexpr char REF_F_ARRAY[] = "data/xanthodb_f_array.bin";
static constexpr char REF_T_MATRIX[] = "data/xanthodb_t_matrix.bin";
static constexpr char REF_SRTD_E_ARRAY[] = "data/xanthodb_sorted_e_array.bin";
static constexpr char REF_JAC_DATA[] = "data/xanthodb_jac.bin";
static constexpr char REF_AJI_DATA[] = "data/xanthodb_aji.bin";

static constexpr char REF_SUBSET1_LC_ARRAY[] = "data/xdb_subset1_lc_array.bin";
static constexpr char REF_SUBSET1_LP_ARRAY[] = "data/xdb_subset1_lp_array.bin";
static constexpr char REF_SUBSET1_F_ARRAY[] = "data/xdb_subset1_f_array.bin";
static constexpr char REF_SUBSET1_T_MATRIX[] = "data/xdb_subset1_t_matrix.bin";
static constexpr char REF_SUBSET1_SRTD_E_ARRAY[] =
    "data/xdb_subset1_sorted_e_array.bin";
static constexpr char REF_SUBSET1_JAC_DATA[] = "data/xdb_subset1_jac.bin";
static constexpr char REF_SUBSET1_AJI_DATA[] = "data/xdb_subset1_aji.bin";

static constexpr char REF_SUBSET2_LC_ARRAY[] = "data/xdb_subset2_lc_array.bin";
static constexpr char REF_SUBSET2_LP_ARRAY[] = "data/xdb_subset2_lp_array.bin";
static constexpr char REF_SUBSET2_F_ARRAY[] = "data/xdb_subset2_f_array.bin";
static constexpr char REF_SUBSET2_T_MATRIX[] = "data/xdb_subset2_t_matrix.bin";
static constexpr char REF_SUBSET2_SRTD_E_ARRAY[] =
    "data/xdb_subset2_sorted_e_array.bin";
static constexpr char REF_SUBSET2_JAC_DATA[] = "data/xdb_subset2_jac.bin";
static constexpr char REF_SUBSET2_AJI_DATA[] = "data/xdb_subset2_aji.bin";

static constexpr char DATA_QSUB_TEST_INPUT[] = "data/qsub_test_input.txt";
static constexpr char REF_QSUB_JAC_DATA[] = "data/xdb_qry_subset_jac.bin";
static constexpr char REF_QSUB_AJI_DATA[] = "data/xdb_qry_subset_aji.bin";

TEST_CASE("Query Genome Metadata", "[db_meta_data]") {
    SQLiteIfT sqlt_if(G_DB_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    REQUIRE(G_PROTEINSET == dbMeta.proteinSet);
    REQUIRE(G_GENOMESET == dbMeta.genomeSet);
}

TEST_CASE("Query No. of Tetramers", "[ntetramers_query]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0);
    SQLiteIfT sqlt_if(G_DB_PATH);
    sqlt_if.queryTetramerOccCounts("PF00119.20", 2060, 2144, Lc);
    REQUIRE(Lc[2060] == 20);
    REQUIRE(Lc[2100] == 20);
    REQUIRE(Lc[2140] == 3);
    REQUIRE(Lc[2144] == 17);
}

TEST_CASE("Query Protein Tetramer Counts", "[prot_teteramer_counts]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    IdMatType T(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    SQLiteIfT sqlt_if(G_DB_PATH);
    sqlt_if.queryProteinTetramerCounts(G_PROTEINSET, 0, 3, T);
    REQUIRE(G_QRY_PST_TETRA_CTS[0] ==
            std::vector<IdType>(T.row_begin(0), T.row_end(0)));
    REQUIRE(G_QRY_PST_TETRA_CTS[1] ==
            std::vector<IdType>(T.row_begin(1), T.row_end(1)));
    REQUIRE(G_QRY_PST_TETRA_CTS[2] ==
            std::vector<IdType>(T.row_begin(2), T.row_end(2)));
    REQUIRE(G_QRY_PST_TETRA_CTS[3] ==
            std::vector<IdType>(T.row_begin(3), T.row_end(3)));
}

TEST_CASE("Query Protein Set Tetramers", "[prot_set_teteramers]") {
    std::vector<IdType> Lc(G_NTETRAMERS, 0), Lp(G_NTETRAMERS, 0);
    for (auto& rx : G_QRT_PST_GP_PAIRS) {
        Lc[rx[0]] += rx[1] / sizeof(IdType);
    }
    //
    // Parallel prefix sum on Lc to construct Lp
    int cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan, + : cumulativeSum)
    for (std::size_t i = 0; i < Lc.size(); i++) {
        Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
        cumulativeSum += Lc[i];
    }

    REQUIRE(Lp[2000] == 0);
    std::vector<IdPairType> F(Lp.back() + Lc.back(), IdPairType(-1, -1));
    SQLiteIfT sqlt_if(G_DB_PATH);
    int fct = 0;
    int rc = sqlt_if.queryProteinSetGPPairs(G_PROTEINSET, 2000, 3000, F.begin(),
                                            &fct);
    REQUIRE(rc == SQLITE_OK);
    REQUIRE(F[Lp[2000]] == IdPairType(5, 0));
    REQUIRE(F[Lp[2415] - 1] == IdPairType(35, 0x13));
    REQUIRE(F[Lp[2415]] == IdPairType(74, 0x11));
}

TEST_CASE("Test Data Structures Construction", "[construct_LFTE]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);

    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Test Data Subset 1 Structures Construction",
          "[construct_subset1_LFTE]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET1_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_SUBSET1_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_SUBSET1_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_SUBSET1_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_SUBSET1_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SUBSET1_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Test Data Subset 2 Structures Construction",
          "[construct_subset2_LFTE]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET2_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.constructL();
    std::vector<IdType> Lc = pfaaiData.refLc();
    std::vector<IdType> Lp = pfaaiData.refLp();
    std::vector<IdType> refLc, refLp;
    {
        std::ifstream lcx(REF_SUBSET2_LC_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLc);
    }
    {
        std::ifstream lcx(REF_SUBSET2_LP_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(lcx);
        iarchive(refLp);
    }
    REQUIRE(refLc == Lc);
    REQUIRE(refLp == Lp);
    //
    pfaaiData.constructF();
    std::vector<IdPairType> F = pfaaiData.refF();
    std::vector<IdPairType> refF;
    {
        std::ifstream fcx(REF_SUBSET2_F_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(fcx);
        iarchive(refF);
    }
    REQUIRE(refF == F);
    //
    pfaaiData.constructT();
    IdMatType T = pfaaiData.refT();
    IdMatType refT(G_PROTEINSET_SIZE, G_GENOMESET_SIZE);
    {
        std::ifstream tcx(REF_SUBSET2_T_MATRIX, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refT);
    }
    REQUIRE(refT == T);
    //
    pfaaiData.constructE();
    const std::vector<ETriple<IdType>>& E = pfaaiData.refE().E;
    //{
    //    std::ofstream lcx1(REF_SUBSET1_SRTD_E_ARRAY, std::ios::binary);
    //    cereal::BinaryOutputArchive orachive(lcx1);
    //    orachive(E);
    //}
    std::vector<ETriple<IdType>> refE;
    {
        std::ifstream tcx(REF_SUBSET2_SRTD_E_ARRAY, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refE);
    }
    REQUIRE(E == refE);
    IdType nTotal = 0;
    for (const auto& ex : E)
        if (ex.genomeA == 0 || ex.genomeA == 2)
            nTotal++;
    std::cout << "N total : " << nTotal << std::endl;
}

TEST_CASE("Implementation Test", "[compute_JAC_AJI]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //

    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;
    {
        std::ifstream tcx(REF_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;
    {
        std::ifstream tcx(REF_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("Subset 1 Implementation Test", "[compute_subset1_JAC_AJI]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET1_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //

    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;

    {
        std::ifstream tcx(REF_SUBSET1_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;
    {
        std::ifstream tcx(REF_SUBSET1_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("Subset 2 Implementation Test", "[compute_subset2_JAC_AJI]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_SUBSET2_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    PFDataT pfaaiData(sqlt_if, dbMeta);
    pfaaiData.construct();
    //
    ParFAAIImpl<IdType, double> pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;

    {
        std::ifstream tcx(REF_SUBSET2_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;

    {
        std::ifstream tcx(REF_SUBSET2_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
    // for(auto i = 0; i < refAJI.size(); i++)
    //    REQUIRE(refAJI[i] == AJI[i]);
}

TEST_CASE("PF Query Subset Test", "[query_subset]") {
    SQLiteInterface<IdType, DatabaseNames> sqlt_if(G_DB_PATH);
    DBMetaData dbMeta;
    sqlt_if.queryMetaData(dbMeta);
    //
    std::vector<std::string> qryGenomeSet;
    std::ifstream in_stream(DATA_QSUB_TEST_INPUT);
    std::string str;
    while (in_stream >> str) {
        qryGenomeSet.push_back(str);
    }
    //
    PFQSubDataT pfaaiData(sqlt_if, dbMeta, qryGenomeSet);
    int errorCode = pfaaiData.construct();

    sqlt_if.closeDB();
    //
    assert(errorCode == PFAAI_OK);
    // Run parallel Fast AAI algorithm
    PFImpl pfaaiImpl(pfaaiData);
    pfaaiImpl.computeJAC();
    std::vector<JACTuple<IdType, double>> JAC = pfaaiImpl.getJAC(), refJAC;
    {
        std::ifstream tcx(REF_QSUB_JAC_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refJAC);
    }
    REQUIRE(refJAC == JAC);
    // for(auto i = 0; i < JAC.size(); i++)
    //    REQUIRE(refJAC[i] == JAC[i]);
    //
    pfaaiImpl.computeAJI();
    std::vector<double> AJI = pfaaiImpl.getAJI(), refAJI;

    {
        std::ifstream tcx(REF_QSUB_AJI_DATA, std::ios::binary);
        cereal::BinaryInputArchive iarchive(tcx);
        iarchive(refAJI);
    }
    REQUIRE(refAJI == AJI);
}

// TEST_CASE("PF QueryxTarget Test", "[qry_tgt]") {
//     SQLiteInterface<IdType, DatabaseNames> qryDBIf(G_DB_SUBSET1_PATH);
//     DBMetaData qryDbMeta;
//     qryDBIf.queryMetaData(qryDbMeta);
//     //
//     SQLiteInterface<IdType, DatabaseNames> tgtDBIf(G_DB_SUBSET2_PATH);
//     DBMetaData tgtDbMeta;
//     tgtDBIf.queryMetaData(tgtDbMeta);
// 
//     //
//     std::unordered_set<std::string> unionProtiens;
//     std::vector<std::string> sharedProtiens;
//     for (const auto& sx : qryDbMeta.proteinSet) {
//         unionProtiens.insert(sx);
//     }
//     for (const auto& sx : tgtDbMeta.proteinSet) {
//         if (unionProtiens.find(sx) != unionProtiens.end()) {
//             sharedProtiens.emplace_back(sx);
//         } else {
//             unionProtiens.insert(sx);
//         }
//     }
// 
//     PFQTData pfaaiQTData(qryDBIf, tgtDBIf, qryDbMeta, tgtDbMeta,
//                          sharedProtiens);
//     // PHASE 1: Construction of the data structures
//     PFAAI_ERROR_CODE pfErrorCode = pfaaiQTData.construct();
//     qryDBIf.closeDB();
//     tgtDBIf.closeDB();
// }
