
#ifndef PAR_FAST_AAI_QT_DATA_H
#define PAR_FAST_AAI_QT_DATA_H

#include "pfaai/interface.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <typename IdType> struct QTGenomeSetData {
    //
    const std::vector<std::string>&c_qryGenomes, c_tgtGenomes;
    //
    // Flags is 1 if the genome is common
    std::vector<bool> m_qryFlag, m_tgtFlag;
    //
    // Output matrix is |Q| x |T|
    // Following
    std::unordered_set<std::string> m_intersectSet, m_unionSet;

    explicit QTGenomeSetData(const std::vector<std::string>& qryGenomes,
                             const std::vector<std::string>& tgtGenomes)
        : c_qryGenomes(qryGenomes), c_tgtGenomes(tgtGenomes),
          m_qryFlag(qryGenomes.size(), 0), m_tgtFlag(tgtGenomes.size(), 0),
          m_intersectSet(std::min(qryGenomes.size(), tgtGenomes.size())),
          m_unionSet(std::max(qryGenomes.size(), tgtGenomes.size())) {
        // Map for target genes
        std::unordered_map<std::string, IdType> tgtIdMap(tgtGenomes.size());
        ;
        // First, add all the target genes to union map
        // Assuming all these genome ids are unique
        for (IdType jx = 0; jx < tgtGenomes.size(); jx++) {
            const std::string& gx = tgtGenomes[jx];
            tgtIdMap[gx] = jx;
            m_unionSet.insert(gx);
        }

        // Query set
        for (IdType jx = 0; jx < qryGenomes.size(); jx++) {
            const std::string& gx = qryGenomes[jx];
            auto gitr = tgtIdMap.find(gx);
            if (gitr != tgtIdMap.end()) {
                m_intersectSet.insert(gx);
            } else {
                m_unionSet.insert(gx);
                // For a row in the matrix, if the row belongs to a gene in
                //    intersection set, then
                //     - Set intersection flags
                m_tgtFlag[gitr->second] = true;
                m_qryFlag[jx] = true;
            }
        }
    }

    inline IdType genomePairToJACIndex(IdType qryGenome,
                                       IdType tgtGenome) const {
        return qryGenome * c_qryGenomes.size() + tgtGenome;
    }
};

template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename JACType>
class ParFAAIQryTgtData
    : public DataStructInterface<IdType, IdPairType, IdMatrixType, JACType> {

  public:
    using ParentT =
        DataStructInterface<IdType, IdPairType, IdMatrixType, JACType>;
    using DBIfxT = DataBaseInterface<IdType, IdPairType, IdMatrixType>;
    using CHelperT = ConstructHelper<IdType, IdPairType, IdMatrixType, JACType>;

  private:
    // Inputs
    const DBIfxT &m_qryDBIf, &m_tgtDBIf;
    const DBMetaData &m_qryDBMeta, &m_tgtDBMeta;
    const std::vector<std::string>& m_proteinSet;
    float m_slack;
    //
    QTGenomeSetData<IdType> m_qtGenomeSetData;
    //
    IdType m_nSharedProteins;
    // Data structure
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    //
    IdMatrixType m_qryT, m_tgtT, m_T;
    std::vector<IdPairType> m_F;
    //
    bool m_flagInitL;

  public:
    explicit ParFAAIQryTgtData(const DBIfxT& qryDbif, const DBIfxT& tgtDbif,
                               const DBMetaData& qryDbMeta,
                               const DBMetaData& tgtDbMeta,
                               const std::vector<std::string>& protSet,
                               float slack = ParentT::DEFAULT_SLACK_PCT)
        : m_qryDBIf(qryDbif), m_tgtDBIf(tgtDbif), m_qryDBMeta(qryDbMeta),
          m_tgtDBMeta(tgtDbMeta), m_proteinSet(protSet), m_slack(slack),
          m_qtGenomeSetData(qryDbMeta.genomeSet, tgtDbMeta.genomeSet),
          m_nSharedProteins(m_proteinSet.size()), m_Lc(ParentT::NTETRAMERS, 0),
          m_Lp(ParentT::NTETRAMERS, 0),
          m_qryT(m_nSharedProteins, m_qryDBMeta.genomeSet.size()),
          m_tgtT(m_nSharedProteins, m_tgtDBMeta.genomeSet.size()),
          m_T(m_nSharedProteins, m_qtGenomeSetData.m_unionSet.size()), m_F(),
          m_flagInitL(false) {}

    inline const std::vector<IdType>& getLc() const { return m_Lc; }
    inline const std::vector<IdType>& getLp() const { return m_Lp; }
    inline const IdMatrixType& getT() const { return m_T; }
    inline const std::vector<IdPairType>& getF() const { return m_F; }
    inline float getSlackPercentage() const { return m_slack; }
    inline IdType getQryGenomeCount() const {
        return m_qryDBMeta.genomeSet.size();
    }
    inline IdType getTgtGenomeCount() const {
        return m_tgtDBMeta.genomeSet.size();
    }
    inline IdType getGPCount() const {
        return getQryGenomeCount() * getTgtGenomeCount();
    }
    inline IdType getIntersectGenomesCount() const {
        return m_qtGenomeSetData.m_intersectSet.size();
    }
    inline IdType getUnionGenomesCount() const {
        return m_qtGenomeSetData.m_intersectSet.size();
    }

    inline IdType genomePairToJACIndex(IdType genomeA, IdType genomeB) const {
        // TODO(srirampc)
        IdType tgtSize = m_tgtDBMeta.genomeSet.size();
        return genomeA * tgtSize + (genomeB - tgtSize);
    }

    void initJAC(std::vector<JACType>& jac_tuples) const {  //  NOLINT
        jac_tuples.resize(getGPCount());
        // Matrix is of size |Q| x |T|
        //   add tgtSize for each
#pragma omp parallel for
        for (std::size_t i = 0; i < jac_tuples.size(); i++) {
            //
            IdType tgtSize = m_tgtDBMeta.genomeSet.size();
            jac_tuples[i].genomeA = i / tgtSize;
            jac_tuples[i].genomeB = tgtSize + (i % tgtSize);
        }
    }

    ~ParFAAIQryTgtData() {}

  private:
  public:
    //
    virtual PFAAI_ERROR_CODE constructT() {
        this->m_errorCode =
            CHelperT::constructT(m_proteinSet, m_qryDBIf, m_qryT);
        if (this->m_errorCode != PFAAI_OK) {
            return this->m_errorCode;
        }
        this->m_errorCode =
            CHelperT::constructT(m_proteinSet, m_tgtDBIf, m_tgtT);
        if (this->m_errorCode != PFAAI_OK) {
            return this->m_errorCode;
        }
        //
        // merge qry T and tgt T
        //  First target columns, followed by query data size
#pragma omp parallel for
        for (std::size_t i = 0; i < m_tgtT.rows(); i++) {
#pragma omp parallel for
            for (std::size_t j = 0; j < m_tgtT.cols(); j++) {
                m_T(i, j) = m_tgtT(i, j);
            }
        }
        //
#pragma omp parallel for
        for (std::size_t i = 0; i < m_qryT.rows(); i++) {
#pragma omp parallel for
            for (std::size_t j = 0; j < m_qryT.cols(); j++) {
                m_T(i, m_tgtDBMeta.genomeSet.size() + j) = m_qryT(i, j);
            }
        }
        return this->m_errorCode;
    }

    //
    virtual PFAAI_ERROR_CODE constructLcandLp() {
        //
        // Tetramer counts  summed up from both target and query databases
        this->m_errorCode = CHelperT::constructLc(m_proteinSet, m_tgtDBIf, m_Lc);
        if (this->m_errorCode != PFAAI_OK)
            return this->m_errorCode;
        //
        this->m_errorCode = CHelperT::constructLc(m_proteinSet, m_qryDBIf, m_Lc);

        // Parallel prefix sum on Lc to construct Lp
        CHelperT::parallelPrefixSum(m_Lc, m_Lp);
        m_flagInitL = true;

        return this->m_errorCode;
    }

    virtual PFAAI_ERROR_CODE constructF() {
        // TODO(x): verify construction
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
#pragma omp single
            {
                errorCodes.resize(nThreads, 0);
                threadTimers.resize(nThreads);
                distribute_bags_of_tasks(nThreads, IdType(m_F.size()), m_Lc,
                                         m_slack, tetramerStart, tetramerEnd);
            }
            //
            threadTimers[threadID].reset();
            IdType tgtFCount = 0, qryFCount = 0;
            auto threadIter = m_F.begin() + m_Lp[tetramerStart[threadID]];
            errorCodes[threadID] = m_tgtDBIf.queryProteinSetGPPairs(
                m_proteinSet, tetramerStart[threadID], tetramerEnd[threadID],
                threadIter, &tgtFCount);
            //
            threadIter += tgtFCount;
            errorCodes[threadID] = m_qryDBIf.queryProteinSetGPPairs(
                m_proteinSet, tetramerStart[threadID], tetramerEnd[threadID],
                threadIter, &qryFCount);
            //
            // Increment the genome ids of query db with target Genome Size
            for (auto qIter = threadIter; qIter != threadIter + qryFCount;
                 qIter++) {
                (*qIter).second += m_tgtDBMeta.genomeSet.size();
            }
            threadTimers[threadID].elapsed();
        }
        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return (this->m_errorCode = PFAAI_ERR_CONSTRUCT);
        //
        // printThreadTimes(" F construction : ", threadTimers);

        return this->m_errorCode;
    }
};

#endif  // !PAR_FAST_AAI_QT_DATA_H
#define PAR_FAST_AAI_QT_DATA_H
