//
//
//

#ifndef PAR_FAST_AAI_DATA_H
#define PAR_FAST_AAI_DATA_H
#include <algorithm>
#include <cassert>
#include <chrono>  // NOLINT
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <vector>

#include "pfaai/interface.hpp"
#include "pfaai/utils.hpp"

using time_point_t =
    std::chrono::time_point<std::chrono::high_resolution_clock>;

template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename JACType>
class ParFAAIData
    : public DataStructInterface<IdType, IdPairType, IdMatrixType, JACType> {

  public:
    using ParentT =
        DataStructInterface<IdType, IdPairType, IdMatrixType, JACType>;
    using DBIfxT = DataBaseInterface<IdType, IdPairType, IdMatrixType>;

  private:
    // Inputs
    const DBIfxT& m_DBIf;
    const DBMetaData& m_dbMeta;
    const std::vector<std::string>& m_proteinSet;
    const std::vector<std::string>& m_genomeSet;
    IdType m_nProteins, m_nGenomes;
    float m_slack;
    // Data structure
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    std::vector<IdPairType> m_F;
    IdMatrixType m_T;
    bool m_flagInitL;

  public:
    explicit ParFAAIData(const DBIfxT& dbif, const DBMetaData& dbMeta,
                         float slack = ParentT::DEFAULT_SLACK_PCT)
        : m_DBIf(dbif), m_dbMeta(dbMeta), m_proteinSet(dbMeta.proteinSet),
          m_genomeSet(dbMeta.genomeSet), m_nProteins(m_proteinSet.size()),
          m_nGenomes(m_genomeSet.size()), m_slack(slack), 
          m_Lc(ParentT::NTETRAMERS, 0), m_Lp(ParentT::NTETRAMERS, 0),
          m_T(m_nProteins, m_nGenomes), m_flagInitL(false) {}

    inline const std::vector<IdType>& getLc() const { return m_Lc; }
    inline const std::vector<IdType>& getLp() const { return m_Lp; }
    inline const IdMatrixType& getT() const { return m_T; }
    inline const std::vector<IdPairType>& getF() const { return m_F; }
    inline float getSlackPercentage() const { return m_slack; }
    inline IdType getQryGenomeCount() const { return m_nGenomes; }
    inline IdType getTgtGenomeCount() const { return m_nGenomes; }
    inline IdType getGPCount() const {
        return (m_nGenomes * (m_nGenomes - 1) / 2);
    }

    inline IdType genomePairToJACIndex(IdType genomeA, IdType genomeB) const {
        return ParentT::genomePairToJACIndex(m_nGenomes, genomeA, genomeB);
    }

    void initJAC(std::vector<JACType>& jac_tuples) const {  //  NOLINT
        ParentT::initJAC(m_nGenomes, getGPCount(), jac_tuples); 
    }
    
    ~ParFAAIData() {}

    PFAAI_ERROR_CODE constructLcandLp() {
        this->m_errorCode = ParentT::constructLc(m_proteinSet, m_DBIf, m_Lc);
        // Parallel prefix sum on Lc to construct Lp
        ParentT::parallelPrefixSum(m_Lc, m_Lp);
        m_flagInitL = true;
        return this->m_errorCode;
    }

    PFAAI_ERROR_CODE constructF() {
        assert(m_Lp.back() > 0);
        if (m_flagInitL == false) {
            std::cout << "Lc and Lp are not Initialized" << std::endl;
            return PFAAI_ERR_CONSTRUCT;
        }
        //
        m_F.resize(m_Lp.back() + m_Lc.back(), IdPairType(-1, -1));
        std::vector<IdType> tetramerStart, tetramerEnd;
        std::vector<int> errorCodes;
        std::vector<timer> threadTimers;
#pragma omp parallel default(none)                                             \
    shared(tetramerStart, tetramerEnd, errorCodes, threadTimers)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            IdType fCount = 0;
#pragma omp single
            {
                errorCodes.resize(nThreads, 0);
                threadTimers.resize(nThreads);
                distribute_bags_of_tasks(nThreads, IdType(m_F.size()), m_Lc,
                                         m_slack, tetramerStart, tetramerEnd);
            }
            threadTimers[threadID].reset();
            errorCodes[threadID] = m_DBIf.queryProteinSetGPPairs(
                m_proteinSet, tetramerStart[threadID], tetramerEnd[threadID],
                m_F.begin() + m_Lp[tetramerStart[threadID]], &fCount);
            threadTimers[threadID].elapsed();
        }
        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (this->m_errorCode = PFAAI_ERR_CONSTRUCT);
        //
        // printThreadTimes(" F construction : ", threadTimers);
        return PFAAI_OK;
    }

    PFAAI_ERROR_CODE constructT() {
        return (this->m_errorCode = ParentT::constructT(m_proteinSet, m_DBIf, m_T));
    }


};

#endif  // !PAR_FAST_AAI_DATA_H
#define PAR_FAST_AAI_DATA_H
