#ifndef PFAAI_INTERFACE_HPP
#define PFAAI_INTERFACE_HPP

#include <algorithm>
#include <cassert>
#include <fmt/format.h>
#include <omp.h>
#include <ostream>
#include <sqlite3.h>
#include <string>
#include <vector>

#include "pfaai/utils.hpp"

enum PFAAI_ERROR_CODE {
    PFAAI_OK = 0,
    PFAAI_ERR_SQLITE_DB = 1,
    PFAAI_ERR_SQLITE_MEM_ALLOC = 2,
    PFAAI_ERR_CONSTRUCT = 3
};

template <typename DT1, typename DT2> struct DPair {
    DT1 first;
    DT2 second;
    explicit DPair(DT1 a, DT2 b) : first(a), second(b) {}
    DPair() : first(DT1(-1)), second(DT2(-1)) {}
    DPair(const DPair& other) : first(other.first), second(other.second) {}
    const DPair& operator=(const DPair& other) {
        first = other.first;
        second = other.second;
        return *this;
    }

    inline bool operator==(const DPair<DT1, DT2>& other) const {
        return (first == other.first) && (second == other.second);
    }

    template <class Archive> void serialize(Archive& archive) {  // NOLINT
        archive(first, second);
    }
};

template <typename DT1, typename DT2>
std::ostream& operator<<(std::ostream& ox, DPair<DT1, DT2> const& cx) {
    ox << "(" << cx.first << ", " << cx.second << ")";
    return ox;
}

template <typename DT> class DMatrix {
    std::size_t m_nrows, m_ncols;
    std::vector<DT> m_data;

  public:
    explicit DMatrix(std::size_t nrows, std::size_t ncols, DT ivx = DT(0))
        : m_nrows(nrows), m_ncols(ncols), m_data(nrows * ncols, ivx) {}

    explicit DMatrix(std::size_t nrows, std::size_t ncols, const DT* arr)
        : m_nrows(nrows), m_ncols(ncols), m_data(arr, arr + (nrows * ncols)) {}

    inline DT& operator()(std::size_t i, std::size_t j) {
        assert(i < m_nrows);
        assert(j < m_ncols);
        return m_data[i * m_ncols + j];
    }

    inline DT operator()(std::size_t i, std::size_t j) const {
        assert(i < m_nrows);
        assert(j < m_ncols);
        return m_data[i * m_ncols + j];
    }

    void resize(std::size_t nrows, std::size_t ncols) {
        m_nrows = nrows;
        m_ncols = ncols;
        m_data.resize(m_nrows * m_ncols);
    }

    inline bool operator==(const DMatrix<DT>& other) const {
        return (m_ncols == other.m_ncols && m_nrows == other.m_nrows &&
                m_data == other.m_data);
    }
    //
    typename std::vector<DT>::iterator row_begin(std::size_t i){
        return m_data.begin() + i * m_ncols;
    }
    typename std::vector<DT>::iterator row_end(std::size_t i){
            return m_data.begin() + (i + 1) * m_ncols;
    }
    std::size_t rows() const { return m_nrows; }
    std::size_t cols() const { return m_ncols; }
    //
    const std::vector<DT>& data() const { return m_data; }
    //
    template <class Archive> void serialize(Archive& archive) {  // NOLINT
        archive(m_nrows, m_ncols, m_data);
    }
};

template <typename IT>
std::ostream& operator<<(std::ostream& ox, DMatrix<IT> const& cx) {
    ox << "{" << cx.rows() << ", " << cx.cols() << ", "
       << fmt::format("[{}]", fmt::join(cx.data(), ", ")) << "}";
    return ox;
}

template <typename IdType, typename IdPairType, typename IdMatType,
          typename ErrCodeType = int>
class DataBaseInterface {

  public:
    using PIterT = typename std::vector<IdPairType>::iterator;
    //
    DataBaseInterface() {}
    //
    virtual inline bool isDBOpen() const = 0;
    virtual ErrCodeType getDBErrorCode() const = 0;
    virtual const char* getDBError() const = 0;
    virtual std::string getDBPath() const = 0;
    virtual PFAAI_ERROR_CODE validate() const = 0;
    virtual inline int closeDB() = 0;
    //
    virtual int
    queryGenomeTetramers(const std::string protein, IdType tetramerStart,
                         IdType tetramerEnd,
                         std::vector<IdType>& Lc) const = 0;  // NOLINT
    virtual int
    queryProteinSetGPPairs(const std::vector<std::string>& proteinSet,
                           IdType tetramerStart, IdType tetramerEnd,

                           PIterT     iterF) const = 0;
                           //const std::vector<IdType>& Lp,
                           //std::vector<IdPairType>& F) const = 0;  // NOLINT
    virtual int
    queryMetaData(std::vector<std::string>& proteinSet,            // NOLINT
                  std::vector<std::string>& genomeSet) const = 0;  // NOLINT

    virtual int
    queryGenomeSet(std::vector<std::string>& genomeSet) const = 0;  // NOLINT

    virtual int
    queryProteinSet(std::vector<std::string>& proteinSet) const = 0;  // NOLINT

    virtual int
    queryProteinTetramerCounts(const std::vector<std::string>& proteinSet,
                               IdType proteinStart, IdType proteinEnd,
                               IdMatType& T) const = 0;  // NOLINT
    //
    virtual ~DataBaseInterface() {}
};

template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename JACType>
class DataStructInterface {
  public:
    // constants
    constexpr static IdType NTETRAMERS = (20 * 20 * 20 * 20);
    constexpr static float DEFAULT_SLACK_PCT = 0.0;
    //
    DataStructInterface() {}
    //
    virtual const std::vector<IdType>& getLc() const = 0;
    virtual const std::vector<IdType>& getLp() const = 0;
    virtual const IdMatrixType& getT() const = 0;
    virtual const std::vector<IdPairType>& getF() const = 0;
    virtual float getSlackPercentage() const = 0;
    virtual inline IdType getTetramerCount() const { return NTETRAMERS; }
    virtual IdType getQryGenomeCount() const = 0;
    virtual IdType getTgtGenomeCount() const = 0;
    virtual IdType getGPCount() const = 0;
    virtual IdType genomePairToJACIndex(IdType genomeA,
                                        IdType genomeB) const = 0;
    virtual void initJAC(std::vector<JACType>& jac_tuples) const = 0;  // NOLINT
    //
    virtual ~DataStructInterface() {}

    static PFAAI_ERROR_CODE constructT(
        const std::vector<std::string>& proteinSet,
        const DataBaseInterface<IdType, IdPairType, IdMatrixType>& inDBIf,
        IdMatrixType& inT) {  // NOLINT
        // Assuming inT is a initalized matrix of size nProteins x nGenomes
        std::vector<timer> threadTimers;
        std::vector<int> errorCodes;
#pragma omp parallel default(none)                                             \
    shared(proteinSet, errorCodes, threadTimers, inDBIf, inT)
        {
            IdType nProteins = proteinSet.size();
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

            errorCodes[threadID] = inDBIf.queryProteinTetramerCounts(
                proteinSet, proteinStart, proteinEnd, inT);
            threadTimers[threadID].elapsed();
        }

        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return PFAAI_ERR_CONSTRUCT;
        //
        // printThreadTimes("T construction : ", threadTimers);
        return PFAAI_OK;
    }

    static PFAAI_ERROR_CODE constructLc(
        const std::vector<std::string>& proteinSet,
        const DataBaseInterface<IdType, IdPairType, IdMatrixType>& inDBIf,
        std::vector<IdType>& Lc) {  // NOLINT
        std::vector<int> errorCodes;
        // Lc construction
#pragma omp parallel default(none) shared(errorCodes, proteinSet, inDBIf, Lc)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            IdType tetraStart = BLOCK_LOW(threadID, nThreads, NTETRAMERS);
            IdType tetraEnd = BLOCK_HIGH(threadID, nThreads, NTETRAMERS);
#pragma omp single
            { errorCodes.resize(nThreads, 0); }

            for (const std::string protein : proteinSet) {
                int qryErrCode = inDBIf.queryGenomeTetramers(
                    protein, tetraStart, tetraEnd, Lc);
                if (qryErrCode != SQLITE_OK) {
                    errorCodes[threadID] = qryErrCode;
                }
            }
        }
        if (std::any_of(errorCodes.begin(), errorCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return PFAAI_ERR_CONSTRUCT;
        return PFAAI_OK;
    }

    static void parallelPrefixSum(const std::vector<IdType>& Lc,
                                  std::vector<IdType>& Lp) {  // NOLINT
        // Parallel prefix sum on Lc to construct Lp
        IdType cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan, + : cumulativeSum)
        for (int i = 0; i < Lc.size(); i++) {
            Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
            cumulativeSum += Lc[i];
        }
    }
};

#endif  // !PFAAI_INTERFACE_HPP
