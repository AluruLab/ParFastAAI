//
//
//

#ifndef PAR_FAST_AAI_H
#define PAR_FAST_AAI_H
#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <utility>
#include <vector>

#include "pfaai/sqltif.hpp"
#include "pfaai/utils.hpp"

enum PFAAI_ERROR_CODE {
    PFAAI_OK = 0,
    PFAAI_ERR_SQLITE_DB = 1,
    PFAAI_ERR_SQLITE_MEM_ALLOC = 2,
    PFAAI_ERR_CONSTRUCT = 3
};

using time_point_t =
    std::chrono::time_point<std::chrono::high_resolution_clock>;

template <typename IdType, typename IdPairType = std::pair<IdType, IdType>>
class ParFAAIData {
    // Inputs
    SQLiteInterface<DatabaseNames> m_sqltIf;
    std::vector<std::string> m_proteinSet;
    IdType m_nGenomes;
    float m_slack;
    // Error Codes
    PFAAI_ERROR_CODE m_errorCode;
    // Data structure
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    std::vector<IdPairType> m_F;
    std::vector<std::vector<IdType>> m_T;

  public:
    // constexpr static IdType _tetramerCount = 160000;
    constexpr static IdType NTETRAMERS = (20 * 20 * 20 * 20);
    constexpr static float DEFAULT_SLACK_PCT = 0.0;

    explicit ParFAAIData(const std::string dbPath,
                         float slack = DEFAULT_SLACK_PCT)
        : m_sqltIf(dbPath, DatabaseNames()), m_slack(slack),
          m_errorCode(PFAAI_OK) {
        if (m_sqltIf.isDBOpen() == PFAAI_OK) {
            m_sqltIf.queryMetaData(m_proteinSet, m_nGenomes);
        }
    }

    int validateDBOpen() {
        if (m_sqltIf.isDBOpenError()) {
            std::cerr << "Error in opening " << m_sqltIf.getDBPath()
                      << std::endl;
            std::cerr << "The error was: " << m_sqltIf.getDBError()
                      << std::endl;
            return PFAAI_ERR_SQLITE_DB;
        }
        if (m_sqltIf.isDBNull()) {
            std::cerr << "SQLite is unable to allocate memory for the database "
                      << m_sqltIf.getDBPath() << std::endl;
            return PFAAI_ERR_SQLITE_MEM_ALLOC;
        }
        return PFAAI_OK;
    }

    inline const std::vector<IdType>& getLc() const { return m_Lc; }
    inline const std::vector<IdType>& getLp() const { return m_Lp; }
    inline const std::vector<std::vector<IdType>>& getT() const { return m_T; }
    inline const std::vector<IdPairType>& getF() const { return m_F; }
    inline float getSlackPercentage() const { return m_slack; }
    inline IdType getTetramerCount() const { return NTETRAMERS; }
    inline IdType getGenomeCount() const { return m_nGenomes; }

    ~ParFAAIData() { m_sqltIf.closeDB(); }

    void printThreadTimes(const std::string& prt_prefix,
                          const std::vector<timer>& threadTimers) const {
        for (int threadID = 0; threadID < threadTimers.size(); threadID++) {
            std::cout << "Thread " << threadID;
            threadTimers[threadID].print_elapsed(prt_prefix, std::cout);
        }
    }

  public:  // construction functions
    PFAAI_ERROR_CODE constructLcandLp() {
        m_Lc.resize(NTETRAMERS, 0);
        m_Lp.resize(NTETRAMERS, 0);
        std::vector<int> errorCodes;
        // Lc construction
#pragma omp parallel default(none) shared(errorCodes)
        {
            int totalNumThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            int tetramerStart =
                BLOCK_LOW(threadID, totalNumThreads, NTETRAMERS);
            int tetramerEnd = BLOCK_HIGH(threadID, totalNumThreads, NTETRAMERS);

            for (const std::string protein : m_proteinSet) {
                int qryErrCode = m_sqltIf.queryGenomeTetramers(
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
        return PFAAI_OK;
    }

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

    PFAAI_ERROR_CODE constructF() {
        std::vector<int> tetramerStart, tetramerEnd, errorCodes;
        std::vector<timer> threadTimers;
        //
        m_F.resize(m_Lp[NTETRAMERS - 1] + m_Lc[NTETRAMERS - 1],
                   IdPairType(-1, -1));
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
            errorCodes[threadID] = m_sqltIf.queryProtienSetGPPairs(
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
        m_T.resize(nProteins);
        for (auto& row : m_T) {
            row.resize(m_nGenomes, 0);
        }
        //
#pragma omp parallel default(none) shared(nProteins, errorCodes, threadTimers)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            //
#pragma omp single
            { threadTimers.resize(nThreads); }
            //
            threadTimers[threadID].reset();
            int proteinStart = BLOCK_LOW(threadID, nThreads, nProteins);
            int proteinEnd = BLOCK_HIGH(threadID, nThreads, nProteins);

            errorCodes[threadID] = m_sqltIf.queryProtienTetramerCounts(
                m_proteinSet, proteinStart, proteinEnd, m_T);
            threadTimers[threadID].elapsed();
        }

        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (m_errorCode = PFAAI_ERR_CONSTRUCT);
        //
        printThreadTimes(" T construction : ", threadTimers);
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
