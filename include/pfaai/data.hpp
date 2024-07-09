//
//
//

#ifndef PAR_FAST_AAI_H
#define PAR_FAST_AAI_H
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <sqlite3.h>

#include "fmt/core.h"
#include "fmt/format.h"
#include "pfaai/interface.hpp"
#include "pfaai/utils.hpp"

using time_point_t =
    std::chrono::time_point<std::chrono::high_resolution_clock>;

template <typename DT1, typename DT2> struct IDPair {
    DT1 first;
    DT2 second;
    explicit IDPair(DT1 a, DT2 b) : first(a), second(b) {}
    explicit IDPair() : first(DT1(-1)), second(DT2(-1)) {}
    IDPair(const IDPair& other) : first(other.first), second(other.second) {}
    const IDPair& operator=(const IDPair& other) {
        first = other.first;
        second = other.second;
        return *this;
    }

    inline bool operator==(const IDPair<DT1, DT2>& other) const {
        return (first == other.first) && (second == other.second);
    }

    template <class Archive> void serialize(Archive& archive) {
        archive(first, second);
    }
};

template <typename DT1, typename DT2>
std::ostream& operator<<(std::ostream& ox, IDPair<DT1, DT2> const& cx) {
    ox << "(" << cx.first << ", " << cx.second << ")";
    return ox;
}

template <typename DT> class IDMatrix {
    std::size_t m_nrows, m_ncols;
    std::vector<DT> m_data;

  public:
    explicit IDMatrix(std::size_t nrows, std::size_t ncols, DT ivx = DT(0))
        : m_nrows(nrows), m_ncols(ncols), m_data(nrows * ncols, ivx) {}

    explicit IDMatrix(std::size_t nrows, std::size_t ncols, const DT* arr)
        : m_nrows(nrows), m_ncols(ncols), m_data(arr, arr + (nrows * ncols)) {}

    inline DT& operator()(std::size_t i, std::size_t j) {
        assert(i < m_nrows);
        assert(j < m_ncols);
        return m_data[i * m_ncols + j];
    }

    inline DT at(std::size_t i, std::size_t j) const {
        assert(i < m_nrows);
        assert(j < m_ncols);
        return m_data[i * m_ncols + j];
    }

    void resize(std::size_t nrows, std::size_t ncols) {
        m_nrows = nrows;
        m_ncols = ncols;
        m_data.resize(m_nrows * m_ncols);
    }

    inline bool operator==(const IDMatrix<DT>& other) const {
        return (m_ncols == other.m_ncols && m_nrows == other.m_nrows &&
                m_data == other.m_data);
    }
    std::vector<DT> row(std::size_t i) {
        return std::vector(m_data.begin() + i * m_ncols,
                           m_data.begin() + (i + 1) * m_ncols);
    }
    std::size_t rows() const { return m_nrows; }
    std::size_t cols() const { return m_ncols; }
    const std::vector<DT>& data() const { return m_data; }
    template <class Archive> void serialize(Archive& archive) {
        archive(m_nrows, m_ncols, m_data);
    }
};

template <typename IT>
std::ostream& operator<<(std::ostream& ox, IDMatrix<IT> const& cx) {
    ox << "{" << cx.rows() << ", " << cx.cols() << ", "
       << fmt::format("{}", fmt::join(cx.data(), ",")) << "}";
    return ox;
}

template <typename IdType, typename PairType, typename MatrixType>
class ParFAAIData : public DataStructInterface<IdType, PairType, MatrixType> {
    // Inputs
    DataBaseInterface<IdType, PairType, MatrixType>& m_DBIf;
    const std::vector<std::string>& m_proteinSet;
    const std::vector<std::string>& m_genomeSet;
    IdType m_nProtiens, m_nGenomes;
    float m_slack;
    // Error Codes
    PFAAI_ERROR_CODE m_errorCode;
    // Data structure
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    std::vector<PairType> m_F;
    MatrixType m_T;
    bool m_flagInitL;

  public:
    constexpr static IdType NTETRAMERS = (20 * 20 * 20 * 20);
    constexpr static float DEFAULT_SLACK_PCT = 0.0;

    explicit ParFAAIData(DataBaseInterface<IdType, PairType, MatrixType>& dbif,
                         const std::vector<std::string>& protSet,
                         const std::vector<std::string>& gnmSet,
                         float slack = DEFAULT_SLACK_PCT)
        : m_DBIf(dbif), m_proteinSet(protSet), m_genomeSet(gnmSet),
          m_nProtiens(m_proteinSet.size()), m_nGenomes(m_genomeSet.size()),
          m_slack(slack), m_errorCode(PFAAI_OK), m_Lc(NTETRAMERS, 0),
          m_Lp(NTETRAMERS, 0), m_T(m_nProtiens, m_nGenomes),
          m_flagInitL(false) {}

    inline const std::vector<IdType>& getLc() const { return m_Lc; }
    inline const std::vector<IdType>& getLp() const { return m_Lp; }
    inline const MatrixType& getT() const { return m_T; }
    inline const std::vector<PairType>& getF() const { return m_F; }
    inline float getSlackPercentage() const { return m_slack; }
    inline IdType getTetramerCount() const { return NTETRAMERS; }
    inline IdType getGenomeCount() const { return m_nGenomes; }
    inline IdType getGPCount() const {
        return (m_nGenomes * (m_nGenomes - 1) / 2);
    }

    ~ParFAAIData() { m_DBIf.closeDB(); }

  private:
    void printThreadTimes(const std::string& prt_prefix,
                          const std::vector<timer>& threadTimers) const {
        for (int threadID = 0; threadID < threadTimers.size(); threadID++) {
            std::cout << "Thread " << threadID;
            threadTimers[threadID].print_elapsed(prt_prefix, std::cout);
        }
    }

    inline IdType genomePairToJACIndex(IdType genomeCount, IdType genomeA,
                                       IdType genomeB) const {
        return (genomeCount * genomeA) + genomeB -
               static_cast<int>((genomeA + 2) * (genomeA + 1) / 2);
    }

  public:
    inline IdType genomePairToJACIndex(IdType genomeA, IdType genomeB) const {
        // (37/2) * a - a^2 / 2 + b - 1
        // return (37 * genomeA - genomeA * genomeA) / 2 + genomeB - 1;
        return genomePairToJACIndex(m_nGenomes, genomeA, genomeB);
    }

    PFAAI_ERROR_CODE constructLcandLp() {
        std::vector<int> errorCodes;
        // Lc construction
#pragma omp parallel default(none) shared(errorCodes)
        {
            int totalNumThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            int tetramerStart =
                BLOCK_LOW(threadID, totalNumThreads, NTETRAMERS);
            int tetramerEnd = BLOCK_HIGH(threadID, totalNumThreads, NTETRAMERS);
#pragma omp single
            { errorCodes.resize(totalNumThreads, 0); }

            for (const std::string protein : m_proteinSet) {
                int qryErrCode = m_DBIf.queryGenomeTetramers(
                    protein, tetramerStart, tetramerEnd, m_Lc);
                if (qryErrCode != SQLITE_OK) {
                    errorCodes[threadID] = qryErrCode;
                }
            }
        }
        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (m_errorCode = PFAAI_ERR_CONSTRUCT);
        // Parallel prefix sum on Lc to construct Lp
        IdType cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan, + : cumulativeSum)
        for (int i = 0; i < m_Lc.size(); i++) {
            m_Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
            cumulativeSum += m_Lc[i];
        }
        m_flagInitL = true;
        return PFAAI_OK;
    }

  private:  // construction functions
    inline void distributeFTasks(int nThreads, std::vector<int>& tetramerStart,
                                 std::vector<int>& tetramerEnd) {
        int nTasks = m_F.size();
        std::vector<int> nTasksDist(nThreads, 0);
        int nTasksPerThread =
            (static_cast<float>(nTasks) / nThreads) * (1 + m_slack);

        for (int tetraId = 0, thid = 0; tetraId < m_Lc.size(); tetraId++) {
            IdType nTetraTasks = m_Lc[tetraId];
            if (nTasksDist[thid] + nTetraTasks <= nTasksPerThread ||
                thid == nThreads - 1) {
                nTasksDist[thid] += nTetraTasks;
                if (tetramerStart[thid] == -1) {
                    tetramerStart[thid] = tetraId;
                }
                tetramerEnd[thid] = tetraId;
            } else {
                // Move to the next processor and assign the task there
                thid++;
                nTasksDist[thid] += nTetraTasks;
                tetramerStart[thid] = tetraId;
                tetramerEnd[thid] = tetraId;
            }
        }
    }

  public:
    PFAAI_ERROR_CODE constructF() {
        std::vector<int> tetramerStart, tetramerEnd, errorCodes;
        std::vector<timer> threadTimers;
        assert(m_Lp[NTETRAMERS - 1] > 0);
        if(m_flagInitL == false){
            std::cout << "Lc and Lp are not Initialized" << std::endl;
            return PFAAI_ERR_CONSTRUCT;
        }
        //
        m_F.resize(m_Lp[NTETRAMERS - 1] + m_Lc[NTETRAMERS - 1],
                   PairType(-1, -1));
#pragma omp parallel default(none)                                             \
    shared(tetramerStart, tetramerEnd, errorCodes, threadTimers)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
#pragma omp single
            {
                tetramerStart.resize(nThreads, -1);
                tetramerEnd.resize(nThreads, -1);
                errorCodes.resize(nThreads, 0);
                threadTimers.resize(nThreads);
                distributeFTasks(nThreads, tetramerStart, tetramerEnd);
            }
            threadTimers[threadID].reset();
            errorCodes[threadID] = m_DBIf.queryProtienSetGPPairs(
                m_proteinSet, tetramerStart[threadID], tetramerEnd[threadID],
                m_Lp, m_F);
            threadTimers[threadID].elapsed();
        }
        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (m_errorCode = PFAAI_ERR_CONSTRUCT);
        //
        printThreadTimes(" F construction : ", threadTimers);
        return PFAAI_OK;
    }

    PFAAI_ERROR_CODE constructT() {
        std::vector<timer> threadTimers;
        int nProteins = m_proteinSet.size();
        std::vector<int> errorCodes;
        //
        // m_T.resize(nProteins);
        // for (auto& row : m_T) {
        //     row.resize(m_nGenomes, 0);
        // }
        //
#pragma omp parallel default(none) shared(nProteins, errorCodes, threadTimers)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            //
#pragma omp single
            {
                threadTimers.resize(nThreads);
                errorCodes.resize(nThreads);
            }
            //
            threadTimers[threadID].reset();
            int proteinStart = BLOCK_LOW(threadID, nThreads, nProteins);
            int proteinEnd = BLOCK_HIGH(threadID, nThreads, nProteins);

            errorCodes[threadID] = m_DBIf.queryProtienTetramerCounts(
                m_proteinSet, proteinStart, proteinEnd, m_T);
            threadTimers[threadID].elapsed();
        }

        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (m_errorCode = PFAAI_ERR_CONSTRUCT);
        //
        printThreadTimes("T construction : ", threadTimers);
        return PFAAI_OK;
    }

    PFAAI_ERROR_CODE construct() {
        timer run_timer;
        PFAAI_ERROR_CODE errorCode = constructLcandLp();
        run_timer.elapsed();
        run_timer.print_elapsed("Lc, Lp contruction time: ", std::cout);
        if (errorCode != PFAAI_OK) {
            std::cerr << "Error in constructing Lc and Lp, error code = "
                      << errorCode << std::endl;
            return errorCode;
        }
        //
        run_timer.reset();
        errorCode = constructF();
        run_timer.elapsed();
        run_timer.print_elapsed("F contruction time: ", std::cout);
        if (errorCode != PFAAI_OK) {
            std::cerr << "Error in constructing F, error code = " << errorCode
                      << std::endl;
            return errorCode;
        }
        //
        run_timer.reset();
        errorCode = constructT();
        run_timer.elapsed();
        run_timer.print_elapsed("T construction time: ", std::cout);
        if (errorCode != PFAAI_OK) {
            std::cerr << "Error in constructing T, error code = " << errorCode
                      << std::endl;
            return errorCode;
        }
        return PFAAI_OK;
    }
};

#endif  // !PAR_FAST_AAI_H
#define PAR_FAST_AAI_H
