///
// @file interface.hpp
// @brief The interface definition for the database and datastructires,
//        and other utilities and error codes.
// @author Sriram P C <srirampc@gatech.edu>
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

#ifndef PFAAI_INTERFACE_HPP
#define PFAAI_INTERFACE_HPP

#include <cassert>
#include <initializer_list>
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
// Class encapsulating the array E, the array of triples
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
    std::vector<std::string> qyGenomeSet;
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

    DataBaseInterface() {}
    virtual ~DataBaseInterface() {}
    
    //
    // Getter functions
    virtual const DBMetaData& getMeta() const = 0;
    virtual std::string getDBPath() const = 0;
    virtual ErrCodeType getDBErrorCode() const = 0;

    //
    // Functions that open/close validates database
    virtual int closeDB() = 0;
    virtual PFAAI_ERROR_CODE validate() const = 0;

    //
    // Function that Queries the database for the number of occurrences
    // of range tetramers for the given protein.
    // Assumes that the Lc is a vector of length NTETRAMERS.
    virtual int tetramerOccCounts(const std::string protein,
                                  IdPairType tetraRange,
                                  std::vector<IdType>& Lc) const = 0;  // NOLINT
    //
    // For a given set of proteins and a range of tetramers,
    // populates the (protein, genome) pair array.
    // Assumes that the pair iterator has been allocated with enough memory.
    // Also returns the number of populated.
    virtual int proteinSetGPPairs(IdPairType tetraRange, PairIterT iterF,
                                  IdType* fCount) const = 0;
    //
    // For a range of proteins among the set of proteins, populates the
    // 2-D array such that T(i, j)  contains the number of tetramers for the
    // protein i and genome j.
    virtual int proteinTetramerCounts(IdPairType proteinRange,
                                      IdMatType& T) const = 0;  // NOLINT
};

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

    // Input Data
    const std::vector<std::string>& m_proteinSet;

    // Data structures
    std::vector<IdType> m_Lc;
    std::vector<IdType> m_Lp;
    IdMatrixType m_T;
    std::vector<IdPairType> m_F;
    EParData<IdType> m_pE;

    // Parameters
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
    virtual float slack() const { return m_slack; }

    // Reference functions to key data structures
    inline virtual const std::vector<IdType>& refLc() const { return m_Lc; }
    inline virtual const std::vector<IdType>& refLp() const { return m_Lp; }
    inline virtual const IdMatrixType& refT() const { return m_T; }
    inline virtual const std::vector<IdPairType>& refF() const { return m_F; }
    inline virtual const EParData<IdType>& refE() const { return m_pE; }

    // Input data reference
    inline virtual const std::vector<std::string>& refProteinSet() const {
        return m_proteinSet;
    }
    virtual const std::vector<std::string>& refQuerySet() const = 0;
    virtual const std::vector<std::string>& refTargetSet() const = 0;
    virtual IdType qrySetSize() const = 0;
    virtual IdType tgtSetSize() const = 0;

    //
    // Getter functions to the dimensions
    virtual inline IdType nTetramers() const { return NTETRAMERS; }
    virtual IdType nGenomePairs() const = 0;

    //
    // Given that a tetramer occurs nQry times in query genome and nTgt times in
    //   target genome, the function returns the number of tuples that should
    //   in the E array : array of (protein, genome_A, genome_B) tuples
    virtual IdType countGenomePairs(IdType nQry, IdType nTgt) const = 0;

    //
    // Indicator functions for genomes :
    // Accepts genome ids that are used in the 'F' and 'E' array
    //   - isQryGenome : return true if a genome id is a valid query genome
    virtual bool isQryGenome(IdType genome) const = 0;
    //   - isValidPair : returns true if a pair is valid (query, target pair)
    virtual bool isValidPair(IdType qry, IdType tgt) const = 0;

    //
    // Jaccard values are computed in a flat JAC Array, JAC is the
    //     total of the number of genome pairs.
    // genomePairToIndex: maps a genome id pair to the index in JAC array.
    virtual IdType genomePairToIndex(IdType genomeA, IdType genomeB) const = 0;

    //
    // Mapper functions : Used to generate the output Matrix
    //  -  mapQueryId : maps from the 'genome id used in the data structures' to
    //  the index in output matrix, specifically, 'column index of the matrix'
    virtual IdType mapQueryId(IdType qry) const = 0;
    //  -  mapTargetId : maps from the 'genome id used in the data structures'
    //  the index in output matrix, specifically, 'row index of the matrix'
    virtual IdType mapTargetId(IdType tgt) const = 0;

    //
    // Construction functions for Data structures
    virtual std::vector<JACType> initJAC() const = 0;
    virtual PFAAI_ERROR_CODE constructL() = 0;
    virtual PFAAI_ERROR_CODE constructF() = 0;
    virtual PFAAI_ERROR_CODE constructT() = 0;
    virtual PFAAI_ERROR_CODE constructE() = 0;

    //
    // Main construction function, constructs in the following order:
    //    Lc -> Lp -> F -> T -> E
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

// Default interfaces
template <typename IdType>
using DefaultDBInterface =
    DataBaseInterface<IdType, DPair<IdType, IdType>, DMatrix<IdType>, int>;

template <typename IdType>
using DefaultDataStructInterface =
    DataStructInterface<IdType, DPair<IdType, IdType>, DMatrix<IdType>,
                        JACTuple<IdType>>;

#endif  // !PFAAI_INTERFACE_HPP
