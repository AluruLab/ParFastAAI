#ifndef PFAAI_INTERFACE_HPP
#define PFAAI_INTERFACE_HPP

#include <cassert>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <omp.h>
#include <sqlite3.h>

#include "pfaai/utils.hpp"

enum PFAAI_ERROR_CODE {
    PFAAI_OK = 0,
    PFAAI_ERR_SQLITE_DB = 1,
    PFAAI_ERR_SQLITE_MEM_ALLOC = 2,
    PFAAI_ERR_CONSTRUCT = 3
};

#define PFAAI_HANDLE_ERROR(err_code, sprefix, ostream)                         \
    {                                                                          \
        if (err_code != PFAAI_OK) {                                            \
            ostream << "Error in " << sprefix                                  \
                    << " ; Error code : " << err_code << std::endl;            \
            return err_code;                                                   \
        }                                                                      \
    }

//
// Utility Class for Pair with first and second element
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

template <typename DT1, typename DT2>
std::string format_as(DPair<DT1, DT2> const& cx) {
    return fmt::format("({}, {})", cx.first, cx.second);
}

//
// Utility Class for Matrix
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
    typename std::vector<DT>::iterator row_begin(std::size_t i) {
        return m_data.begin() + i * m_ncols;
    }
    typename std::vector<DT>::iterator row_end(std::size_t i) {
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

//
// Utility class for triple
template <typename IdType> struct ETriple {
    IdType proteinIndex;
    IdType genomeA;
    IdType genomeB;

    ETriple() : proteinIndex(-1), genomeA(-1), genomeB(-1) {}
    ETriple(const IdType proteinIndexVal, const IdType genomeAVal,
            const IdType genomeBVal)
        : proteinIndex(proteinIndexVal), genomeA(genomeAVal),
          genomeB(genomeBVal) {}

    inline bool operator<(const ETriple<IdType>& other) const {
        if (genomeA != other.genomeA) {
            return genomeA < other.genomeA;
        }
        if (genomeB != other.genomeB) {
            return genomeB < other.genomeB;
        }
        return proteinIndex < other.proteinIndex;
    }

    inline bool operator==(const ETriple<IdType>& other) const {
        return (genomeA == other.genomeA) && (genomeB == other.genomeB) &&
               (proteinIndex == other.proteinIndex);
    }

    template <class Archive> void serialize(Archive& archive) {  // NOLINT
        archive(proteinIndex, genomeA, genomeB);
    }
};

template <typename IT = int> std::string format_as(const ETriple<IT>& ijx) {
    return fmt::format("({:>3d}, {:>3d}, {:>3d})", ijx.genomeA, ijx.genomeB,
                       ijx.proteinIndex);
}

template <typename IT>
std::ostream& operator<<(std::ostream& ox, ETriple<IT> const& cx) {
    ox << "(" << cx.genomeB << ", " << cx.genomeA << ", " << cx.proteinIndex
       << ")";
    return ox;
}

//
// Tuple class with
//    - genomeA, genomeB
//    - N (No. entries corresponding to the genome pair)
//    - S (Sum of entries corresponding to the genome pair)
//
template <typename IdType, typename ValueType = double> struct JACTuple {
    IdType genomeA;
    IdType genomeB;
    ValueType S;
    IdType N;

    inline bool operator==(const JACTuple<IdType, ValueType>& other) const {
        return (genomeA == other.genomeA) && (genomeB == other.genomeB) &&
               (N == other.N) && (std::abs(S - other.S) < 1e-7);
    }

    template <class Archive> void serialize(Archive& archive) {  // NOLINT
        archive(genomeA, genomeB, S, N);
    }
};

template <typename IT = int, typename VT = double>
std::string format_as(const JACTuple<IT, VT>& ijx) {
    return fmt::format("({:>3d}, {:>3d}, {:>03.2f}, {:>3d})", ijx.genomeA,
                       ijx.genomeB, ijx.S, ijx.N);
}

template <typename IT, typename VT>
std::ostream& operator<<(std::ostream& ox, JACTuple<IT, VT> const& cx) {
    ox << "(" << cx.genomeA << ", " << cx.genomeB << ", " << cx.S << ", "
       << cx.N << ")";
    return ox;
}

//
// Class encapsulating the array of triples
//
template <typename IdType> struct EParData {
    using IdTripleT = ETriple<IdType>;
    // tetramer tuples : Array E construction : Phase 2
    std::vector<IdTripleT> E;
    // thread distributions of E array ~(|E|/p)
    std::vector<IdType> threadEStarts, threadESize;  // size nThreads
};

struct DBMetaData {
    std::vector<std::string> proteinSet;
    std::vector<std::string> genomeSet;
    std::vector<std::string> attGenomeSet;
};

//
// Database Interface
template <typename IdType, typename IdPairT, typename IdMatT,
          typename ErrCodeT = int>
class DataBaseInterface {
  public:
    using IdPairType = IdPairT;
    using IdMatType = IdMatT;
    using ErrCodeType = ErrCodeT;
    using PairIterT = typename std::vector<IdPairType>::iterator;
    //
    DataBaseInterface() {}
    //
    virtual ~DataBaseInterface() {}
    //
    virtual inline bool isDBOpen() const = 0;
    virtual ErrCodeType getDBErrorCode() const = 0;
    virtual const char* getDBError() const = 0;
    virtual std::string getDBPath() const = 0;
    virtual PFAAI_ERROR_CODE validate() const = 0;
    virtual inline int closeDB() = 0;
    //
    virtual int
    queryTetramerOccCounts(const std::string protein, IdType tetramerStart,
                           IdType tetramerEnd,
                           std::vector<IdType>& Lc) const = 0;  // NOLINT
    virtual int
    queryProteinSetGPPairs(const std::vector<std::string>& proteinSet,
                           IdType tetramerStart, IdType tetramerEnd,
                           PairIterT iterF, IdType* fCount) const = 0;
    virtual int queryMetaData(DBMetaData& metaData) const = 0;  // NOLINT

    virtual int
    queryGenomeSet(std::vector<std::string>& genomeSet) const = 0;  // NOLINT

    virtual int
    queryProteinSet(std::vector<std::string>& proteinSet) const = 0;  // NOLINT

    virtual int
    queryProteinTetramerCounts(const std::vector<std::string>& proteinSet,
                               IdType proteinStart, IdType proteinEnd,
                               IdMatType& T) const = 0;  // NOLINT
};

template <typename IdType>
using DefaultDBInterface =
    DataBaseInterface<IdType, DPair<IdType, IdType>, DMatrix<IdType>, int>;

//
// Data Structure interface
template <typename IdType, typename IdPairT, typename IdMatrixT,
          typename JaccardT>
class DataStructInterface {
  public:
    using IdPairType = IdPairT;
    using IdMatrixType = IdMatrixT;
    using JACType = JaccardT;

  protected:
    // Error Codes
    PFAAI_ERROR_CODE m_errorCode;
    std::unordered_map<std::string, bool> m_initFlags{
        {"L", false},
        {"E", false},
        {"F", false},
        {"T", false},
    };
    // Data
    const std::vector<std::string>& m_proteinSet;
    // Data structures
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    std::vector<IdPairType> m_F;
    IdMatrixType m_T;
    EParData<IdType> m_pE;

    // value data
    float m_slack;

  public:
    // constants
    constexpr static IdType NTETRAMERS = (20 * 20 * 20 * 20);
    constexpr static float DEFAULT_SLACK_PCT = 0.0;
    //
    DataStructInterface(const std::vector<std::string>& proteinSet,
                        std::size_t nTRows, std::size_t nTColumns, float slack)
        : m_errorCode(PFAAI_OK), m_proteinSet(proteinSet), m_Lc(NTETRAMERS, 0),
          m_Lp(NTETRAMERS, 0), m_T(nTRows, nTColumns), m_slack(slack) {}
    //
    virtual ~DataStructInterface() {}
    //
    // Reference functions to data structures
    inline virtual const std::vector<IdType>& refLc() const { return m_Lc; }
    inline virtual const std::vector<IdType>& refLp() const { return m_Lp; }
    inline virtual const IdMatrixType& refT() const { return m_T; }
    inline virtual const std::vector<IdPairType>& refF() const { return m_F; }
    inline virtual const EParData<IdType>& refE() const { return m_pE; }
    inline virtual const std::vector<std::string>& refProteinSet() const {
        return m_proteinSet;
    }
    inline virtual const std::vector<std::string>& refQuerySet() const = 0;
    inline virtual const std::vector<std::string>& refTargetSet() const = 0;
    //
    // Getter functions to dimensions
    virtual inline IdType nTetramers() const { return NTETRAMERS; }
    virtual float slack() const { return m_slack; }
    virtual IdType qrySetSize() const = 0;
    virtual IdType tgtSetSize() const = 0;
    virtual IdType nGenomePairs() const = 0;
    //
    // Mapper/Validator functions
    virtual IdType genomePairToIndex(IdType genomeA, IdType genomeB) const = 0;
    virtual inline IdType mapQueryId(IdType qry) const = 0;
    virtual inline IdType mapTargetId(IdType tgt) const = 0;
    virtual bool isQryGenome(IdType genome) const = 0;
    virtual bool isValidPair(IdType qry, IdType tgt) const = 0;
    virtual IdType countGenomePairs(IdType nQry, IdType nTgt) const = 0;
    //
    // Construction functions
    virtual std::vector<JACType> initJAC() const = 0;
    virtual PFAAI_ERROR_CODE constructL() = 0;
    virtual PFAAI_ERROR_CODE constructF() = 0;
    virtual PFAAI_ERROR_CODE constructT() = 0;
    virtual PFAAI_ERROR_CODE constructE() = 0;
    //
    // Main construction function
    virtual PFAAI_ERROR_CODE construct() {
        timer run_timer;
        m_errorCode = constructL();
        PRINT_RUNTIME_MEMUSED(run_timer, "Lc & Lp contruction : ", std::cout);
        PFAAI_HANDLE_ERROR(m_errorCode, "Lc & Lp contruction", std::cerr);
        //
        run_timer.reset();
        m_errorCode = constructF();
        PRINT_RUNTIME_MEMUSED(run_timer, "F contruction       : ", std::cout);
        PFAAI_HANDLE_ERROR(m_errorCode, "F contruction", std::cerr);
        //
        run_timer.reset();
        m_errorCode = constructT();
        PRINT_RUNTIME_MEMUSED(run_timer, "T construction      : ", std::cout);
        PFAAI_HANDLE_ERROR(m_errorCode, "T contruction", std::cerr);
        //
        run_timer.reset();
        m_errorCode = constructE();
        PRINT_RUNTIME_MEMUSED(run_timer, "E constr.   (fin)   : ", std::cout);
        PFAAI_HANDLE_ERROR(m_errorCode, "E constr.   (fin)", std::cerr);
        return PFAAI_OK;
    }
};

template <typename IdType>
using DefaultDataStructInterface =
    DataStructInterface<IdType, DPair<IdType, IdType>, DMatrix<IdType>,
                        JACTuple<IdType>>;

#endif  // !PFAAI_INTERFACE_HPP
