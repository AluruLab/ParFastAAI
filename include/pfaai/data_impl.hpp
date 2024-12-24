//
//
//
#ifndef PAR_FAST_AAI_DATA_H
#define PAR_FAST_AAI_DATA_H
#include <cassert>
#include <iostream>
#include <omp.h>
#include <sqlite3.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "pfaai/helpers.hpp"
#include "pfaai/interface.hpp"

template <typename IdType>
class ParFAAIData : public DefaultDataStructInterface<IdType> {
  public:
    // types
    using Parent = DefaultDataStructInterface<IdType>;
    using DBIfx = DefaultDBInterface<IdType>;
    using DSHelper = DefaultDSHelper<IdType>;

    using IdPairType = typename Parent::IdPairType;
    using IdMatrixType = typename Parent::IdMatrixType;
    using JACType = typename Parent::JACType;

  private:
    // Inputs
    const DBIfx& m_DBIf;
    const DBMetaData& m_dbMeta;
    const std::vector<std::string>& m_genomeSet;
    // Dimensions
    IdType m_nProteins, m_nGenomes;

  public:
    explicit ParFAAIData(const DBIfx& dbif, const DBMetaData& dbMeta,
                         float slack = Parent::DEFAULT_SLACK_PCT)
        : Parent(dbMeta.proteinSet, dbMeta.proteinSet.size(),
                 dbMeta.genomeSet.size(), slack),
          m_DBIf(dbif), m_dbMeta(dbMeta), m_genomeSet(dbMeta.genomeSet),
          m_nProteins(dbMeta.proteinSet.size()),
          m_nGenomes(m_genomeSet.size()) {}

    ~ParFAAIData() {}

    inline virtual const std::vector<std::string>& refQuerySet() const {
        return m_genomeSet;
    } 
    inline virtual const std::vector<std::string>& refTargetSet() const {
        return m_genomeSet;
    }
    // Dimensions
    inline IdType qrySetSize() const { return m_nGenomes; }
    inline IdType tgtSetSize() const { return m_nGenomes; }
    inline IdType nGenomePairs() const {
        return (m_nGenomes * (m_nGenomes - 1) / 2);
    }

    // Mapping functions
    inline IdType genomePairToIndex(IdType genomeA, IdType genomeB) const {
        return (m_nGenomes * genomeA) + genomeB -
               static_cast<IdType>((genomeA + 2) * (genomeA + 1) / 2);
    }
    inline IdType mapQueryId(IdType qry) const { return qry; }
    inline IdType mapTargetId(IdType tgt) const { return tgt; }
    inline bool isQryGenome(IdType genome) const { return true; }
    inline bool isValidPair(IdType qryId, IdType tgtId) const {
        return qryId < tgtId;
    }
    inline IdType countGenomePairs(IdType nQry, IdType nTgt) const {
        assert(nTgt == 0);
        return (nQry * (nQry - 1) / 2);
    }

    // Construction functions
    std::vector<JACType> initJAC() const {
        std::vector<JACType> jac_tuples(nGenomePairs());
        IdType gA = 0, gB = 1;
        for (std::size_t i = 0; i < jac_tuples.size(); i++) {
            jac_tuples[i].genomeA = gA;
            jac_tuples[i].genomeB = gB;

            if (gB == m_nGenomes - 1) {
                gA += 1;
                gB = gA + 1;
            } else {
                gB += 1;
            }
        }
        return jac_tuples;
    }

    PFAAI_ERROR_CODE constructL() {
        this->m_errorCode = DSHelper::constructLc(*this, m_DBIf, this->m_Lc);
        if (this->m_errorCode != PFAAI_OK) {
            std::cout << "Lc and Lp fail to Initialize" << std::endl;
            return (this->m_errorCode = PFAAI_ERR_CONSTRUCT);
        }
        // Parallel prefix sum on Lc to construct Lp
        DSHelper::parallelPrefixSum(this->m_Lc, this->m_Lp);
        assert(this->m_Lp.back() > 0);
        this->m_initFlags["L"] = true;
        return this->m_errorCode;
    }

    PFAAI_ERROR_CODE constructF() {
        assert(this->refLp().back() > 0);
        if (this->m_errorCode != PFAAI_OK) {
            std::cout << "Lc and Lp are not Initialized" << std::endl;
            return (this->m_errorCode = PFAAI_ERR_CONSTRUCT);
        }
        return (this->m_errorCode =
                    DSHelper::constructF(*this, m_DBIf, this->m_F));
    }

    PFAAI_ERROR_CODE constructT() {
        return (this->m_errorCode =
                    DSHelper::constructT(*this, m_DBIf, this->m_T));
    }

    PFAAI_ERROR_CODE constructE() {
        return this->m_errorCode = DSHelper::constructE(*this, this->m_pE);
    }
};

//
//
template <typename IdType>
class ParFAAIQSubData : public DefaultDataStructInterface<IdType> {
  public:
    using Parent = DefaultDataStructInterface<IdType>;
    using DBIfx = DefaultDBInterface<IdType>;
    using DSHelper = DefaultDSHelper<IdType>;

    using IdPairType = typename Parent::IdPairType;
    using IdMatrixType = typename Parent::IdMatrixType;
    using JACType = typename Parent::JACType;

  private:
    // Input References
    const DBIfx& m_DBIf;
    const DBMetaData& m_dbMeta;
    const std::vector<std::string>& m_qryGenomeSet;
    const std::vector<std::string>& m_genomeSet;

    // Data Sizes
    IdType m_nProteins, m_nGenomes, m_nQryGenomes, m_nTgtGenomes;

    // Data structures
    std::vector<IdType> m_genomeIndexMap;
    std::vector<bool> m_qryIndicator;
    std::vector<IdType> m_qryLookup, m_tgtLookup;
    //
    std::vector<std::string> m_tgtGenomeSet;
    IdType m_qtMatSize;

  public:
    explicit ParFAAIQSubData(const DBIfx& srcDbif, const DBMetaData& dbMeta,
                             const std::vector<std::string>& qryGenomeSet,
                             const std::vector<std::string>& protSet,
                             float slack = Parent::DEFAULT_SLACK_PCT)
        : Parent(protSet, protSet.size(), dbMeta.genomeSet.size(), slack),
          m_DBIf(srcDbif), m_dbMeta(dbMeta), m_qryGenomeSet(qryGenomeSet),
          m_genomeSet(dbMeta.genomeSet), m_nProteins(protSet.size()),
          m_nGenomes(m_genomeSet.size()), m_nQryGenomes(m_qryGenomeSet.size()),
          m_nTgtGenomes(m_nGenomes - m_nQryGenomes),
          m_genomeIndexMap(m_nGenomes), m_qryIndicator(m_nGenomes, false),
          m_qryLookup(m_nQryGenomes), m_tgtLookup(m_nTgtGenomes),
          m_tgtGenomeSet(m_nTgtGenomes),
          m_qtMatSize(m_nQryGenomes * m_nTgtGenomes) {
        //
        // Query Set
        std::unordered_map<std::string, IdType> qSet;
        for (std::size_t ix = 0; ix < m_qryGenomeSet.size(); ix++) {
            qSet[m_qryGenomeSet[ix]] = IdType(ix);
        }
        //
        // Query and Target lookups and index maps
        for (IdType ix = 0, jx = 0; ix < m_nGenomes; ix++) {
            const std::string& gxName = m_genomeSet[ix];
            auto gxItr = qSet.find(m_genomeSet[ix]);
            if (gxItr != qSet.end()) {
                m_qryIndicator[ix] = true;
                m_qryLookup[gxItr->second] = ix;
                m_genomeIndexMap[ix] = gxItr->second;
            } else {
                m_tgtGenomeSet[jx] = gxName;
                m_genomeIndexMap[ix] = jx;
                m_tgtLookup[jx] = ix;
                jx++;
            }
        }

        //
        std::cout << (qSet.find(m_genomeSet[0]) != qSet.end())
                  << fmt::format("[{}]", fmt::join(m_qryGenomeSet, ", "))
                  << fmt::format("[{}]", fmt::join(m_qryIndicator, ", "))
                  << std::endl;
    }

    virtual ~ParFAAIQSubData() {}

    inline virtual const std::vector<std::string>& refQuerySet() const {
        return m_qryGenomeSet;
    } 
    inline virtual const std::vector<std::string>& refTargetSet() const {
        return m_genomeSet;
    }

    //
    inline IdType qrySetSize() const { return m_nQryGenomes; }
    inline IdType tgtSetSize() const { return m_nGenomes; }
    inline IdType nGenomePairs() const {
        // Matrix is of two parts:
        //    1. Full |Q| x |T|
        //    2. Upper triangular matrix of |T| x |T|
        return m_qtMatSize + (m_nQryGenomes * (m_nQryGenomes - 1) / 2);
    }

    inline IdType genomePairToIndex(IdType genomeA, IdType genomeB) const {
        //
        // Matrix is of size |Q| x |NQ| +
        //     each row = number of target genes
        int qIndicator = m_qryIndicator[genomeA] && !m_qryIndicator[genomeB];
        int gia = m_genomeIndexMap[genomeA], gib = m_genomeIndexMap[genomeB];
        return qIndicator * ((gia * m_nTgtGenomes) + gib) +
               // Case 2: If genomeB is a query genome
               (1 - qIndicator) *
                   (m_qtMatSize +
                    ((m_nQryGenomes * gia) + gib -
                     static_cast<IdType>((gia + 2) * (gia + 1) / 2)));
    }
    inline IdType mapQueryId(IdType qry) const { return m_genomeIndexMap[qry]; }
    inline IdType mapTargetId(IdType tgt) const { return tgt; }

    inline bool isQryGenome(IdType genome) const {
        return m_qryIndicator[genome];
    }
    inline bool isValidPair(IdType qry, IdType tgt) const {
        return (m_qryIndicator[qry] && m_qryIndicator[tgt] && qry < tgt) ||
               (m_qryIndicator[qry] && !m_qryIndicator[tgt] && qry != tgt);
    }
    inline IdType countGenomePairs(IdType nQry, IdType nTgt) const {
        return (nQry * nTgt) + (nQry * (nQry - 1) / 2);
    }

    std::vector<JACType> initJAC() const {
        std::vector<JACType> jac_tuples(nGenomePairs());
        //
        // Matrix is of size |Q| x |T|
        //     each row = number of target genes
#pragma omp parallel for
        for (IdType i = 0; i < m_qtMatSize; i++) {
            IdType qid = i / m_nTgtGenomes;
            IdType tid = i % m_nTgtGenomes;
            // TODO(x):: check if correct
            jac_tuples[i].genomeA = m_qryLookup[qid];
            jac_tuples[i].genomeB = m_tgtLookup[tid];
        }

        IdType gA = 0, gB = 1;
        for (std::size_t i = m_qtMatSize; i < jac_tuples.size(); i++) {
            jac_tuples[i].genomeA = m_qryLookup[gA];
            jac_tuples[i].genomeB = m_qryLookup[gB];

            if (gB == m_nQryGenomes - 1) {
                gA += 1;
                gB = gA + 1;
            } else {
                gB += 1;
            }
        }
        return jac_tuples;
    }

    PFAAI_ERROR_CODE constructL() {
        this->m_errorCode = DSHelper::constructLc(*this, m_DBIf, this->m_Lc);
        // Parallel prefix sum on Lc to construct Lp
        DSHelper::parallelPrefixSum(this->m_Lc, this->m_Lp);
        if (this->m_errorCode == PFAAI_OK)
            this->m_initFlags["L"] = true;
        return this->m_errorCode;
    }

    PFAAI_ERROR_CODE constructT() {
        this->m_errorCode = DSHelper::constructT(*this, m_DBIf, this->m_T);
        if (this->m_errorCode == PFAAI_OK)
            this->m_initFlags["T"] = true;
        return this->m_errorCode;
    }

    PFAAI_ERROR_CODE constructF() {
        assert(this->m_Lp.back() > 0);
        this->m_errorCode = DSHelper::constructF(*this, m_DBIf, this->m_F);
        if (this->m_errorCode == PFAAI_OK)
            this->m_initFlags["F"] = true;
        return this->m_errorCode;
    }

    PFAAI_ERROR_CODE constructE() {
        this->m_errorCode = DSHelper::constructE(*this, this->m_pE);
        if (this->m_errorCode == PFAAI_OK)
            this->m_initFlags["E"] = true;
        return this->m_errorCode;
    }
};

// template <typename IdType> struct QTGenomeSetData {
//     //
//     const std::vector<std::string>&c_qryGenomes, c_tgtGenomes;
//     //
//     // Flags is 1 if the genome is common
//     std::vector<bool> m_qryFlag, m_tgtFlag;
//     //
//     // Output matrix is |Q| x |T|
//     // Following
//     std::unordered_set<std::string> m_intersectSet, m_unionSet;
//
//     explicit QTGenomeSetData(const std::vector<std::string>& qryGenomes,
//                              const std::vector<std::string>& tgtGenomes)
//         : c_qryGenomes(qryGenomes), c_tgtGenomes(tgtGenomes),
//           m_qryFlag(qryGenomes.size(), 0), m_tgtFlag(tgtGenomes.size(), 0),
//           m_intersectSet(std::min(qryGenomes.size(), tgtGenomes.size())),
//           m_unionSet(std::max(qryGenomes.size(), tgtGenomes.size())) {
//         // Map for target genes
//         std::unordered_map<std::string, IdType> tgtIdMap(tgtGenomes.size());
//         ;
//         // First, add all the target genes to union map
//         // Assuming all these genome ids are unique
//         for (IdType jx = 0; jx < IdType(tgtGenomes.size()); jx++) {
//             const std::string& gx = tgtGenomes[jx];
//             tgtIdMap[gx] = jx;
//             m_unionSet.insert(gx);
//         }
//
//         // Query set
//         for (IdType jx = 0; jx < IdType(qryGenomes.size()); jx++) {
//             const std::string& gx = qryGenomes[jx];
//             auto gitr = tgtIdMap.find(gx);
//             if (gitr != tgtIdMap.end()) {
//                 m_intersectSet.insert(gx);
//             } else {
//                 m_unionSet.insert(gx);
//                 // For a row in the matrix, if the row belongs to a gene in
//                 //    intersection set, then
//                 //     - Set intersection flags
//                 m_tgtFlag[gitr->second] = true;
//                 m_qryFlag[jx] = true;
//             }
//         }
//     }
//
//     inline IdType genomePairToJACIndex(IdType qryGenome,
//                                        IdType tgtGenome) const {
//         return qryGenome * c_qryGenomes.size() + tgtGenome;
//     }
// };

template <typename IdType>
class ParFAAIQryTgtData : public DefaultDataStructInterface<IdType> {

  public:
    using Parent = DefaultDataStructInterface<IdType>;
    using DBIfx = DefaultDBInterface<IdType>;
    using DSHelper = DefaultDSHelper<IdType>;

    using IdPairType = typename Parent::IdPairType;
    using IdMatrixType = typename Parent::IdMatrixType;
    using JACType = typename Parent::JACType;

  private:
    // Inputs
    const DBIfx &m_qryDBIf, &m_tgtDBIf;
    const DBMetaData &m_qryDBMeta, &m_tgtDBMeta;
    //
    // QTGenomeSetData<IdType> m_qtGenomeSetData;
    //
    IdType m_nSharedProteins;
    IdType m_nUnionSize;
    std::vector<bool> m_qryIndicator;
    // Data structure
    //
    IdMatrixType m_qryT, m_tgtT;
    std::vector<IdPairType> m_F;
    EParData<IdType> m_pE;
    std::vector<IdType> m_queryIndexMap, m_targetIndexMap;

  public:
    explicit ParFAAIQryTgtData(const DBIfx& qryDbif, const DBIfx& tgtDbif,
                               const DBMetaData& qryDbMeta,
                               const DBMetaData& tgtDbMeta,
                               const std::vector<std::string>& protSet,
                               float slack = Parent::DEFAULT_SLACK_PCT)
        : Parent(protSet, protSet.size(),
                 qryDbMeta.genomeSet.size() + tgtDbMeta.genomeSet.size(),
                 slack),
          m_qryDBIf(qryDbif), m_tgtDBIf(tgtDbif), m_qryDBMeta(qryDbMeta),
          m_tgtDBMeta(tgtDbMeta), m_nSharedProteins(protSet.size()),
          m_nUnionSize(m_qryDBMeta.genomeSet.size() +
                       m_tgtDBMeta.genomeSet.size()),
          m_qryT(m_nSharedProteins, m_qryDBMeta.genomeSet.size()),
          m_tgtT(m_nSharedProteins, m_tgtDBMeta.genomeSet.size()) {

        // TODO(x): Initialize qry inicators, lookup, index maps
    }

    ~ParFAAIQryTgtData() {}

    inline virtual const std::vector<std::string>& refQuerySet() const {
        return m_qryDBMeta.genomeSet;
    } 
    inline virtual const std::vector<std::string>& refTargetSet() const {
        return m_tgtDBMeta.genomeSet;
    }

    //
    inline IdType qrySetSize() const { return m_qryDBMeta.genomeSet.size(); }
    inline IdType tgtSetSize() const { return m_tgtDBMeta.genomeSet.size(); }
    inline IdType nGenomePairs() const { return qrySetSize() * tgtSetSize(); }
    inline IdType nUnionGenomes() const { return m_nUnionSize; }

    // Mapper/Validator functions
    inline IdType genomePairToIndex(IdType genomeA, IdType genomeB) const {
        // TODO(x)
        IdType tgtSize = m_tgtDBMeta.genomeSet.size();
        return genomeA * tgtSize + genomeB;
    }
    inline IdType mapQueryId(IdType qry) const { return m_queryIndexMap[qry]; }
    inline IdType mapTargetId(IdType tgt) const { return m_targetIndexMap[tgt]; }
    inline bool isQryGenome(IdType genome) const {
        return m_qryIndicator[genome];
    }
    virtual bool isValidPair(IdType qry, IdType tgt) const {
        return m_qryIndicator[qry] && !m_qryIndicator[tgt];
    }
    virtual IdType countGenomePairs(IdType nQry, IdType nTgt) const {
        return nQry * nTgt;
    }

    std::vector<JACType> initJAC() const {
        std::vector<JACType> jac_tuples(nGenomePairs());
        // Matrix is of size |Q| x |T|
        //   add tgtSize for each
#pragma omp parallel for
        for (std::size_t i = 0; i < jac_tuples.size(); i++) {
            //
            IdType tgtSize = m_tgtDBMeta.genomeSet.size();
            jac_tuples[i].genomeA = i / tgtSize;
            jac_tuples[i].genomeB = i % tgtSize;
        }
        return jac_tuples;
    }

  public:
    //
    virtual PFAAI_ERROR_CODE constructT() {
        this->m_errorCode = DSHelper::constructT(*this, m_qryDBIf, m_qryT);
        if (this->m_errorCode != PFAAI_OK) {
            return this->m_errorCode;
        }
        this->m_errorCode = DSHelper::constructT(*this, m_tgtDBIf, m_tgtT);
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
                this->m_T(i, j) = m_tgtT(i, j);
            }
        }
        //
#pragma omp parallel for
        for (std::size_t i = 0; i < m_qryT.rows(); i++) {
#pragma omp parallel for
            for (std::size_t j = 0; j < m_qryT.cols(); j++) {
                this->m_T(i, m_tgtDBMeta.genomeSet.size() + j) = m_qryT(i, j);
            }
        }
        return this->m_errorCode;
    }

    //
    virtual PFAAI_ERROR_CODE constructL() {
        //
        // Tetramer counts  summed up from both target and query databases
        this->m_errorCode = DSHelper::constructLc(*this, m_tgtDBIf, this->m_Lc);
        if (this->m_errorCode != PFAAI_OK)
            return this->m_errorCode;
        //
        this->m_errorCode = DSHelper::constructLc(*this, m_qryDBIf, this->m_Lc);

        // Parallel prefix sum on Lc to construct Lp
        DSHelper::parallelPrefixSum(this->m_Lc, this->m_Lp);
        this->m_initFlags["L"] = true;

        return this->m_errorCode;
    }

    virtual PFAAI_ERROR_CODE constructF() {
        // TODO(x): construction
        assert(this->m_Lp.back() > 0);
        if (this->m_initFlags["L"] == false) {
            std::cout << "Lc and Lp are not Initialized" << std::endl;
            return PFAAI_ERR_CONSTRUCT;
        }
        //
        m_F.resize(this->m_Lp.back() + this->m_Lc.back(), IdPairType(-1, -1));
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
                distribute_bags_of_tasks(nThreads, IdType(m_F.size()),
                                         this->m_Lc, this->m_slack,
                                         tetramerStart, tetramerEnd);
            }
            //
            threadTimers[threadID].reset();
            IdType tgtFCount = 0, qryFCount = 0;
            auto threadIter = m_F.begin() + this->m_Lp[tetramerStart[threadID]];
            errorCodes[threadID] = m_tgtDBIf.queryProteinSetGPPairs(
                this->refProteinSet(), tetramerStart[threadID],
                tetramerEnd[threadID], threadIter, &tgtFCount);
            //
            threadIter += tgtFCount;
            errorCodes[threadID] = m_qryDBIf.queryProteinSetGPPairs(
                this->refProteinSet(), tetramerStart[threadID],
                tetramerEnd[threadID], threadIter, &qryFCount);
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
        //

        return PFAAI_OK;
    }

    PFAAI_ERROR_CODE constructE() {
        // TODO(x)::
        return PFAAI_OK;
    }
};

#endif  // !PAR_FAST_AAI_DATA_H
#define PAR_FAST_AAI_DATA_H
