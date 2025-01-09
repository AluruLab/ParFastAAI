///
// @file data_impl.hpp
// @brief The helper classes and functions to construct the
//        data structures for different usage cases.
// @author Sriram P C <srirampc@gatech.edu>, Hoang Le <hanh9@gatech.edu>
//
// Copyright 2024 Georgia Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///

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

//
// Class for data structure construction in the simple case where
//  AJI for all pairs of genomes are computed from data of a given database.
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
            std::cerr << "Lc and Lp fail to Initialize" << std::endl;
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
            std::cerr << "Lc and Lp are not Initialized" << std::endl;
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
// Class for data structure construction in the case where
//  AJI is computed for  of genomes are computed for a subset of query genomes
//  in the database.
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
                             float slack = Parent::DEFAULT_SLACK_PCT)
        : Parent(dbMeta.proteinSet, dbMeta.proteinSet.size(),
                 dbMeta.genomeSet.size(), slack),
          m_DBIf(srcDbif), m_dbMeta(dbMeta), m_qryGenomeSet(qryGenomeSet),
          m_genomeSet(dbMeta.genomeSet), m_nProteins(dbMeta.proteinSet.size()),
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
        // std::cout << (qSet.find(m_genomeSet[0]) != qSet.end())
        //           << fmt::format("[{}]", fmt::join(m_qryGenomeSet, ", "))
        //           << fmt::format("[{}]", fmt::join(m_qryIndicator, ", "))
        //           << std::endl;
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

// Class for data structure construction in the case where
//  AJI is computed for  of genomes are computed with query and target are in 
//  two different databases.
//  TODO(X): Implementation is incomplete
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
    const DBIfx& m_qtDBIf;
    const DBMetaData& m_dbMeta;
    //
    IdType m_nSharedProteins;
    IdType m_nUnionGenomes;
    std::vector<bool> m_qryIndicator;
    std::vector<IdType> m_genomeIndexMap;
    //
    // Data structures
    // IdMatrixType m_qryT, m_tgtT;
    // std::vector<IdType> m_qryLc, m_tgtLc;
    // std::vector<IdPairType> m_F;
    // EParData<IdType> m_pE;

  public:
    explicit ParFAAIQryTgtData(const DBIfx& qtDBIf, const DBMetaData& dbMeta,
                               const std::vector<std::string>& protSet,
                               float slack = Parent::DEFAULT_SLACK_PCT)
        : Parent(protSet, dbMeta.qyGenomeSet.size(), dbMeta.genomeSet.size(),
                 slack),
          m_qtDBIf(qtDBIf), m_dbMeta(dbMeta), m_nSharedProteins(protSet.size()),
          m_nUnionGenomes(dbMeta.genomeSet.size() + dbMeta.qyGenomeSet.size()),
          m_qryIndicator(m_nUnionGenomes, false),
          m_genomeIndexMap(m_nUnionGenomes, -1) {
        // Initialize query indicators,  index maps
        IdType nQryGenomes = m_dbMeta.qyGenomeSet.size();
        IdType nTgtGenomes = m_dbMeta.genomeSet.size();
        for (IdType ix = 0; ix < nQryGenomes; ix++) {
            m_qryIndicator[ix] = true;
            m_genomeIndexMap[ix] = ix;
        }
        for (IdType ix = 0; ix < nTgtGenomes; ix++) {
            m_qryIndicator[nQryGenomes + ix] = false;
            m_genomeIndexMap[nQryGenomes + ix] = ix;
        }
    }

    ~ParFAAIQryTgtData() {}

    inline virtual const std::vector<std::string>& refQuerySet() const {
        return m_dbMeta.qyGenomeSet;
    }
    inline virtual const std::vector<std::string>& refTargetSet() const {
        return m_dbMeta.genomeSet;
    }

    //
    inline IdType qrySetSize() const { return m_dbMeta.qyGenomeSet.size(); }
    inline IdType tgtSetSize() const { return m_dbMeta.genomeSet.size(); }
    inline IdType nGenomePairs() const { return qrySetSize() * tgtSetSize(); }
    inline IdType nUnionGenomes() const { return m_nUnionGenomes; }

    // Mapper/Validator functions
    // TODO(x): Verify if these are okay
    inline IdType genomePairToIndex(IdType genomeA, IdType genomeB) const {
        return mapQueryId(genomeA) * tgtSetSize() + mapTargetId(genomeB);
    }
    inline IdType mapQueryId(IdType qry) const { return m_genomeIndexMap[qry]; }
    inline IdType mapTargetId(IdType tgt) const {
        return m_genomeIndexMap[tgt];
    }
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
#pragma omp parallel for
        for (std::size_t i = 0; i < jac_tuples.size(); i++) {
            //
            jac_tuples[i].genomeA = i / tgtSetSize();
            // Shift target genome ids
            jac_tuples[i].genomeB = qrySetSize() + (i % tgtSetSize());
        }
        return jac_tuples;
    }

  public:
    //
    virtual PFAAI_ERROR_CODE constructT() {
        // Older logic : T is constructed seperately for qry and target; merged
        // TODO(x): modify this with  a  joint database query form 2 databases
        //
        //// T is constructed for query and target seperately and merged
        // this->m_errorCode = DSHelper::constructT(*this, m_qryDBIf,
        // m_qryT); if (this->m_errorCode != PFAAI_OK) {
        //     return this->m_errorCode;
        // }
        // this->m_errorCode = DSHelper::constructT(*this, m_tgtDBIf,
        // m_tgtT); if (this->m_errorCode != PFAAI_OK) {
        //     return this->m_errorCode;
        // }
        //// Merge query T and target T
        ////  First query columns, followed by target columns
        // #pragma omp parallel for
        //         for (std::size_t i = 0; i < m_qryT.rows(); i++) {
        // #pragma omp parallel for
        //             for (std::size_t j = 0; j < m_qryT.cols(); j++) {
        //                 this->m_T(i, j) = m_qryT(i, j);
        //             }
        //         }
        // #pragma omp parallel for
        //         for (std::size_t i = 0; i < m_tgtT.rows(); i++) {
        // #pragma omp parallel for
        //             for (std::size_t j = 0; j < m_tgtT.cols(); j++) {
        //                 this->m_T(i, j) = m_tgtT(i, j);
        //             }
        //         }
        //         //
        // #pragma omp parallel for
        //         for (std::size_t i = 0; i < m_qryT.rows(); i++) {
        // #pragma omp parallel for
        //             for (std::size_t j = 0; j < m_qryT.cols(); j++) {
        //                 this->m_T(i, refTargetSet().size() + j) = m_qryT(i,
        //                 j);
        //             }
        //         }
        //         this->m_initFlags["T"] = true;
        return this->m_errorCode;
    }

    //
    virtual PFAAI_ERROR_CODE constructL() {

        // Older logic : L is constructed seperately for qry and target; merged
        // TODO(x): modify this with  a  joint database query form 2 databases
        //
        // // Tetramer counts is summed up from both query and target databases
        // this->m_errorCode = DSHelper::constructLc(*this, m_tgtDBIf, m_qryLc);
        // if (this->m_errorCode != PFAAI_OK)
        //     return this->m_errorCode;
        // //
        // this->m_errorCode = DSHelper::constructLc(*this, m_qryDBIf, m_tgtLc);
        // if (this->m_errorCode != PFAAI_OK)
        //     return this->m_errorCode;

        // #pragma omp parallel for
        //         for (IdType ix = 0; ix < this->nTetramers(); ix++) {
        //             this->m_Lc[ix] = this->m_qryLc[ix] + this->m_tgtLc[ix];
        //         }
        //         // Parallel prefix sum on Lc to build Lp
        //         DSHelper::parallelPrefixSum(this->m_Lc, this->m_Lp);
        //         this->m_initFlags["L"] = true;

        return this->m_errorCode;
    }

    virtual PFAAI_ERROR_CODE constructF() {
        assert(this->m_Lp.back() > 0);
        if (this->m_initFlags["L"] == false) {
            std::cerr << "Lc and Lp are not Initialized" << std::endl;
            return PFAAI_ERR_CONSTRUCT;
        }
        //
        //  Had a logic with F constructed from each database seperately.
        //  Realized doesn't work. Need to have a joint query from both the
        //  databases.
        //
        //
        //        this->m_F.resize(this->m_Lp.back() + this->m_Lc.back(),
        //                         IdPairType(-1, -1));
        //        std::vector<IdType> tetramerStart, tetramerEnd, ntSizes;
        //        std::vector<int> errorCodes;
        //        std::vector<timer> threadTimers;
        // #pragma omp parallel default(none)
        //    shared(tetramerStart, tetramerEnd, errorCodes, threadTimers,
        //    ntSizes)
        //        {
        //            int nThreads = omp_get_num_threads();
        //            int threadID = omp_get_thread_num();
        // #pragma omp single
        //            {
        //                errorCodes.resize(nThreads, 0);
        //                threadTimers.resize(nThreads);
        //                ntSizes = distribute_bags_of_tasks(
        //                    nThreads, IdType(this->m_F.size()), this->m_Lc,
        //                    this->m_slack, tetramerStart, tetramerEnd);
        //            }
        //            //
        //            threadTimers[threadID].reset();
        //            IdType tgtFCount = 0, qryFCount = 0;
        //            auto threadIter =
        //                this->m_F.begin() +
        //                this->m_Lp[tetramerStart[threadID]];
        //            errorCodes[threadID] = m_tgtDBIf.queryProteinSetGPPairs(
        //                this->refProteinSet(), tetramerStart[threadID],
        //                tetramerEnd[threadID], threadIter, &tgtFCount);
        //            //
        //            threadIter += tgtFCount;
        //            errorCodes[threadID] = m_qryDBIf.queryProteinSetGPPairs(
        //                this->refProteinSet(), tetramerStart[threadID],
        //                tetramerEnd[threadID], threadIter, &qryFCount);
        //            assert(qryFCount + tgtFCount == ntSizes[threadID]);
        //            //
        //            // Increment the genome ids of query db with target Genome
        //            Size for (auto qIter = threadIter; qIter != threadIter +
        //            qryFCount;
        //                 qIter++) {
        //                (*qIter).second += m_tgtDBMeta.genomeSet.size();
        //            }
        //            threadTimers[threadID].elapsed();
        //        }
        //        if (std::any_of(errorCodes.begin(), errorCodes.end(),
        //                        [](int rc) { return rc != SQLITE_OK; }))
        //            return (this->m_errorCode = PFAAI_ERR_CONSTRUCT);
        //        //
        //        // printThreadTimes(" F construction : ", threadTimers);
        //        //
        //        this->m_initFlags["F"] = true;
        return PFAAI_OK;
    }

    PFAAI_ERROR_CODE constructE() {
        // this->m_errorCode = DSHelper::constructE(*this, this->m_pE);
        // if (this->m_errorCode == PFAAI_OK)
        //     this->m_initFlags["E"] = true;
        return this->m_errorCode;
    }

    // virtual PFAAI_ERROR_CODE construct() { constructL(); }
};

#endif  // !PAR_FAST_AAI_DATA_H
#define PAR_FAST_AAI_DATA_H
