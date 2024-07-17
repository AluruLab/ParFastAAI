#ifndef PFAAI_RUNNER_HPP
#define PFAAI_RUNNER_HPP

#include <cmath>
#include <cstdlib>
#include <fmt/format.h>
#include <omp.h>
#include <string>
#include <vector>

#include "pfaai/interface.hpp"
#include "pfaai/psort.hpp"
#include "pfaai/utils.hpp"

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

template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename ValueType>
class ParFAAIImpl {
  public:
    using JACType = JACTuple<IdType, ValueType>;
    using DSType =
        DataStructInterface<IdType, IdPairType, IdMatrixType, JACType>;
    using IdTripleT = ETriple<IdType>;

  private:
    // Data references
    const DSType& c_faaiDataRef;
    const std::vector<int>& c_Lc;
    const std::vector<int>& c_Lp;
    const std::vector<IdPairType>& c_F;
    const IdMatrixType& c_T;
    float m_slack;
    IdType m_nTetramers;
    //
    IdType m_nGenomePairs;

    // tetramer tuples : Array E construction : Phase 2
    std::vector<IdType> m_tetramerStart, m_tetramerEnd;  // size nThreads
    std::vector<IdTripleT> m_E;
    // thread distributions of E array ~(|E|/p)
    std::vector<IdType> m_threadEStarts, m_threadESize;  // size nThreads

    // Starts and end indices for Phase 3
    // Thread distributions for genome pair ~(m_nGenomePairs/p)
    std::vector<IdType> m_threadGPStarts, m_threadGPEnds;  // of size nThreads
    // Start and end corresponding to genome pairs
    std::vector<IdType> m_genomePairEStartIndex,
        m_genomePairEEndIndex;  // of size m_nGenomePairs

    // Jaccard Step
    std::vector<JACType> m_JAC;  // of size m_nGenomePairs

    // AJI step
    std::vector<ValueType> m_AJI;  // of size m_nGenomePairs

  public:
    explicit ParFAAIImpl(const DSType& fDataRef)
        : c_faaiDataRef(fDataRef), c_Lc(c_faaiDataRef.getLc()),
          c_Lp(c_faaiDataRef.getLp()), c_F(c_faaiDataRef.getF()),
          c_T(c_faaiDataRef.getT()),
          m_slack(c_faaiDataRef.getSlackPercentage()),
          m_nTetramers(c_faaiDataRef.getTetramerCount()),
          m_nGenomePairs(c_faaiDataRef.getGPCount()) {}

    const std::vector<IdTripleT>& getE() const { return m_E; }
    const std::vector<JACTuple<IdType>>& getJAC() const { return m_JAC; }
    const std::vector<ValueType>& getAJI() const { return m_AJI; }

    inline IdType genomePairToJACIndex(IdType genomeA, IdType genomeB) const {
        return c_faaiDataRef.genomePairToJACIndex(genomeA, genomeB);
    }

  private:
    void distributeTetramers(const int& nThreads) {
        // Ask 1 thread to carry out the distribution of tasks, i.e.
        // tetramer tuples for all threads
        m_tetramerStart.resize(nThreads, -1);
        m_tetramerEnd.resize(nThreads, -1);
        //
        float nTasks = static_cast<float>(c_F.size());
        std::vector<int> threadTaskSum(nThreads, 0);
        int nTasksPerThread =
            static_cast<int>((nTasks / nThreads) * (1 + m_slack));

        int tid = 0;
        for (IdType tetramer = 0; tetramer < IdType(c_Lc.size()); tetramer++) {
            IdType nTetraTasks = c_Lc[tetramer];
            if (threadTaskSum[tid] + nTetraTasks <= nTasksPerThread ||
                tid == nThreads - 1) {
                threadTaskSum[tid] += nTetraTasks;
                if (m_tetramerStart[tid] == -1) {
                    m_tetramerStart[tid] = tetramer;
                }
                m_tetramerEnd[tid] = tetramer;
            } else {
                // Move to the next processor and assign the task there
                tid++;
                threadTaskSum[tid] += nTetraTasks;
                m_tetramerStart[tid] = tetramer;
                m_tetramerEnd[tid] = tetramer;
            }
        }
        m_threadESize.resize(nThreads);
        m_threadEStarts.resize(nThreads);
    }

    IdType countTetramerTuples(const IdType& tStart, const IdType& tEnd) {
        int nTetraTuples = 0;
        // Ask a thread to compute its own number of elements in its local
        // part of E
        for (IdType tetraID = tStart; tetraID <= tEnd; tetraID++) {
            // INCLUSIVE start and end index of the tetramer block in F
            IdType startIndexInF = c_Lp[tetraID];
            IdType endIndexInF = 0;
            if (tetraID < tEnd) {
                endIndexInF = c_Lp[tetraID + 1] - 1;
            } else {
                // tetramerID == tetramerEnd
                if (tetraID == m_nTetramers - 1) {
                    endIndexInF = c_F.size() - 1;
                } else {
                    endIndexInF = c_Lp[tetraID + 1] - 1;
                }
            }

            int currProteinID = c_F[startIndexInF].first;
            int leftBoundary = startIndexInF;
            int rightBoundary = startIndexInF;

            while (rightBoundary <= endIndexInF) {
                if (c_F[rightBoundary].first == currProteinID) {
                    rightBoundary++;
                } else {
                    int n = rightBoundary - leftBoundary;
                    nTetraTuples += n * (n - 1) / 2;

                    currProteinID = c_F[rightBoundary].first;
                    leftBoundary = rightBoundary;
                }
            }

            // Complete the last calculation for when the while loop
            // finishes because rightBoundary > endIndexInF, but the
            // calculation for currProteinID is not done yet
            int n = rightBoundary - leftBoundary;
            nTetraTuples += n * (n - 1) / 2;
        }
        return nTetraTuples;
    }

    void constructTetramerTuples(const int& threadID, const IdType& tStart,
                                 const IdType& tEnd) {
        IdType currentLocalEIndex = m_threadEStarts[threadID];
        for (IdType tetramerID = tStart; tetramerID <= tEnd; tetramerID++) {
            // INCLUSIVE start and end index of the tetramer block in F
            IdType startIndexInF = c_Lp[tetramerID];
            IdType endIndexInF = 0;
            if (tetramerID < tEnd) {
                endIndexInF = c_Lp[tetramerID + 1] - 1;
            } else {
                // tetramerID == tetramerEnd
                if (tetramerID == m_nTetramers - 1) {
                    endIndexInF = c_F.size() - 1;
                } else {
                    endIndexInF = c_Lp[tetramerID + 1] - 1;
                }
            }

            IdType currProteinID = c_F[startIndexInF].first;

            // leftBoundary and rightBoundary defines the block B_l
            // leftBoundary is inclusive, rightBoundary is exclusive
            IdType leftBoundary = startIndexInF;
            IdType rightBoundary = startIndexInF;

            while (rightBoundary <= endIndexInF) {
                if (c_F[rightBoundary].first == currProteinID) {
                    rightBoundary++;
                } else {
                    for (int i = leftBoundary; i < rightBoundary; i++) {
                        int genomeA_ID = c_F[i].second;
                        for (int j = i + 1; j < rightBoundary; j++) {
                            IdType genomeB_ID = c_F[j].second;
                            IdTripleT newElement(currProteinID, genomeA_ID,
                                                 genomeB_ID);
                            m_E[currentLocalEIndex] = newElement;
                            currentLocalEIndex++;
                        }
                    }

                    currProteinID = c_F[rightBoundary].first;
                    leftBoundary = rightBoundary;
                }
            }

            // Complete the last calculation for when the while loop
            // finishes because rightBoundary > endIndexInF, but the
            // calculation for currProteinID is not done yet
            for (IdType i = leftBoundary; i < rightBoundary; i++) {
                IdType genomeA_ID = c_F[i].second;
                for (IdType j = i + 1; j < rightBoundary; j++) {
                    IdType genomeB_ID = c_F[j].second;
                    IdTripleT newElement(currProteinID, genomeA_ID, genomeB_ID);
                    m_E[currentLocalEIndex] = newElement;
                    currentLocalEIndex++;
                }
            }
        }
    }

  public:
    PFAAI_ERROR_CODE generateTetramerTuples() {
        IdType totalESize(0);
        timer run_timer;
        // PHASE 2: Generate tetramer tuples
#pragma omp parallel default(none) shared(totalESize)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
#pragma omp single
            { distributeTetramers(nThreads); }
            // Tetramers specific to this thread
            int tStart = m_tetramerStart[threadID],
                tEnd = m_tetramerEnd[threadID];
            // Size of tetramer tuples corresponding to this thread
            m_threadESize[threadID] = countTetramerTuples(tStart, tEnd);
#pragma omp barrier
            // Get total length of E
#pragma omp for reduction(+ : totalESize)
            for (std::size_t i = 0; i < m_threadESize.size(); i++) {
                totalESize += m_threadESize[i];
            }
#pragma omp single
            { m_E.resize(totalESize); }
            // Get the start indices of each thread-partions of E by parallel
            // prefix sum (exclusive) on _ESize and store it in _EStartIndex
            int cumulativeSum = 0;
#pragma omp simd reduction(inscan, + : cumulativeSum)
            for (std::size_t i = 0; i < m_threadESize.size(); i++) {
                m_threadEStarts[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
                cumulativeSum += m_threadESize[i];
            }
            // Each thread construct its own part of E
            constructTetramerTuples(threadID, tStart, tEnd);
#pragma omp barrier
        }
        run_timer.elapsed();
        run_timer.print_elapsed("E construction      : ", std::cout);
        // Parallel Sort E TODO(): Sorting speed is inconsistent, why ?
        timer srt_timer;
#pragma omp single
        { parallelMergeSort(m_E, 0, m_E.size() - 1, 5); }
        srt_timer.elapsed();
        srt_timer.print_elapsed("E parallel sorting  : ", std::cout);
        return PFAAI_OK;
    }

  private:
    void prepJAC(const int& nThreads) {
        m_genomePairEStartIndex.resize(m_nGenomePairs);
        m_genomePairEEndIndex.resize(m_nGenomePairs);
        c_faaiDataRef.initJAC(m_JAC);
        //
        m_threadGPStarts.resize(nThreads);
        m_threadGPEnds.resize(nThreads);
    }

  protected:
    void distributeGenomePairs(const int& threadID, const int& nThreads) {
        if (threadID < m_nGenomePairs % nThreads) {
            // Each will be responsible for totalGenomePairs /
            // totalNumThreads + 1 pairs
            m_threadGPStarts[threadID] =
                threadID * (m_nGenomePairs / nThreads + 1);
            m_threadGPEnds[threadID] =
                (threadID + 1) * (m_nGenomePairs / nThreads + 1) - 1;
        } else {
            int buffer =
                (m_nGenomePairs % nThreads) * (m_nGenomePairs / nThreads + 1);
            m_threadGPStarts[threadID] =
                buffer + (threadID - m_nGenomePairs % nThreads) *
                             (m_nGenomePairs / nThreads);
            m_threadGPEnds[threadID] =
                buffer +
                (threadID - m_nGenomePairs % nThreads + 1) *
                    (m_nGenomePairs / nThreads) -
                1;
        }
    }

    void findEBlockExtents(const int& threadID,
                           IdType& currentLocalEIndex) {  // NOLINT
        // We need to know where in the sorted E does each genome pair start
        // and end Each thread can look though its local chunk of E and help
        // fill in the 2 arrays: genomePairEStartIndex, genomePairEEndIndex

        // Edge case - First genome pair seen in the local E chunk. It might
        // not start here
        if (currentLocalEIndex > 0) {
            IdTripleT firstElement = m_E[currentLocalEIndex];
            IdTripleT lastElementPrevThread = m_E[currentLocalEIndex - 1];
            if (firstElement.genomeA != lastElementPrevThread.genomeA ||
                firstElement.genomeB != lastElementPrevThread.genomeB) {
                // Then this means we hold the beginning of this genome pair
                // Safely write this information into the
                // genomePairEStartIndex array
                IdType genomePairIndexInJAC = genomePairToJACIndex(
                    firstElement.genomeA, firstElement.genomeB);
                m_genomePairEStartIndex[genomePairIndexInJAC] =
                    currentLocalEIndex;
            }  // Otherwise, let the previous processor fill in the starting
               // index
        } else {
            // currentLocalEIndex == 0
            // Safely write the beginning of the current genome pair index
            // into genomePairEStartIndex array
            IdTripleT firstElement = m_E[currentLocalEIndex];
            IdType genomePairIndexInJAC = genomePairToJACIndex(
                firstElement.genomeA, firstElement.genomeB);
            m_genomePairEStartIndex[genomePairIndexInJAC] = 0;
        }

        while (currentLocalEIndex <
               m_threadEStarts[threadID] + m_threadESize[threadID]) {
            IdType currentGenomeA = m_E[currentLocalEIndex].genomeA;
            IdType currentGenomeB = m_E[currentLocalEIndex].genomeB;

            while (currentLocalEIndex <
                       m_threadEStarts[threadID] + m_threadESize[threadID] &&
                   m_E[currentLocalEIndex].genomeA == currentGenomeA &&
                   m_E[currentLocalEIndex].genomeB == currentGenomeB) {
                currentLocalEIndex++;
            }

            if (currentLocalEIndex <
                m_threadEStarts[threadID] + m_threadESize[threadID]) {
                // Found the end index of the current genome pair
                // Safely write this information into the
                // genomePairEEndIndex array
                IdType genomePairIndexInJAC =
                    genomePairToJACIndex(currentGenomeA, currentGenomeB);
                m_genomePairEEndIndex[genomePairIndexInJAC] =
                    currentLocalEIndex - 1;

                // Move on to the next pair
                // This means we are seeing the start index of a new genome
                // pair Safely write this information into the
                // genomePairEStartIndex array
                genomePairIndexInJAC =
                    genomePairToJACIndex(m_E[currentLocalEIndex].genomeA,
                                         m_E[currentLocalEIndex].genomeB);
                m_genomePairEStartIndex[genomePairIndexInJAC] =
                    currentLocalEIndex;
            } else {
                // currentLocalEIndex == EChunkStartIndex[threadID] +
                // EChunkSize[threadID] This means we've parsed our local E
                // Chunk fully. However, we still need to check if we have
                // seen the end of the current genome pair. If yes, we must
                // write this information If our local E chunk doesn't store
                // the end of the current genome pair, leave it to the next
                // processor to fill this information

                if (currentLocalEIndex < static_cast<IdType>(m_E.size())) {
                    IdTripleT firstElementOfNextThread = m_E[currentLocalEIndex];

                    if (firstElementOfNextThread.genomeA != currentGenomeA ||
                        firstElementOfNextThread.genomeB != currentGenomeB) {
                        // We have seen the end of the current genome pair
                        // Safely write this information into
                        // genomePairEEndIndex array
                        IdType genomePairIndexInJAC = genomePairToJACIndex(
                            currentGenomeA, currentGenomeB);
                        m_genomePairEEndIndex[genomePairIndexInJAC] =
                            currentLocalEIndex - 1;
                    }  // Otherwise, this information will be filled by the
                       // next processor
                } else {
                    // currentLocalEIndex == E.size()
                    IdType genomePairIndexInJAC =
                        genomePairToJACIndex(currentGenomeA, currentGenomeB);
                    m_genomePairEEndIndex[genomePairIndexInJAC] =
                        m_E.size() - 1;
                }
            }
        }
    }

    void print_e() {
        std::cout << "E array : " << std::endl;
        for (int i = 0; i < m_E.size(); i++) {
            fmt::print("{}", m_E[i]);
        }
        std::cout << std::endl << std::endl;

        std::cout << "Genome Pair Extents:" << std::endl;
        for (int i = 0; i < m_genomePairEStartIndex.size(); i++) {
            fmt::print("GP {} : [{}, {}]", m_E[i], m_genomePairEStartIndex[i],
                       m_genomePairEEndIndex[i]);
        }
        std::cout << std::endl << std::endl;
    }

    void computeEBlockJAC(const int& threadID) {
        for (int genomePair = m_threadGPStarts[threadID];
             genomePair <= m_threadGPEnds[threadID]; genomePair++) {
            int currGenomeA = m_JAC[genomePair].genomeA;
            int currGenomeB = m_JAC[genomePair].genomeB;

            // Block Bl of the same genome pair (Ga, Gb)
            // Note: Again, these are inclusive
            int blockBlStart = m_genomePairEStartIndex[genomePair];
            int blockBlEnd = m_genomePairEEndIndex[genomePair];

            double S = 0.0;
            int N = 0;

            // Subblocks Bk of the same protein Pi inside Bl
            int blockBkStart = blockBlStart;
            int blockBkEnd = blockBkStart;

            while (blockBkEnd <= blockBlEnd) {
                int currProteinID = m_E[blockBkStart].proteinIndex;

                while (blockBkEnd <= blockBlEnd &&
                       m_E[blockBkEnd].proteinIndex == currProteinID) {
                    blockBkEnd++;
                }

                if (blockBkEnd <= blockBlEnd) {
                    int BkLength = blockBkEnd - blockBkStart;
                    double J_Pi_Ga_Gb =
                        static_cast<double>(BkLength) /
                        static_cast<double>(c_T(currProteinID, currGenomeA) +
                                            c_T(currProteinID, currGenomeB) -
                                            BkLength);
                    S += J_Pi_Ga_Gb;
                    N += 1;

                    // Move on to the next Bk subblock
                    blockBkStart = blockBkEnd;
                } else {
                    // blockBkEnd > blockBlEnd
                    // Finish the last computation
                    int BkLength = blockBkEnd - blockBkStart;
                    double J_Pi_Ga_Gb =
                        static_cast<double>(BkLength) /
                        static_cast<double>(c_T(currProteinID, currGenomeA) +
                                            c_T(currProteinID, currGenomeB) -
                                            BkLength);
                    S += J_Pi_Ga_Gb;
                    N += 1;
                }
            }
            // DOUBLE CHECK THIS PART
            m_JAC[genomePair].S += S;
            m_JAC[genomePair].N += N;
        }
    }

  public:
    int computeJAC() {
        timer run_timer;
        // Prepare the JAC vector
        // PHASE 3: Compute the Jaccard Coefficient values
#pragma omp parallel default(none) shared(std::cout)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
#pragma omp single
            { prepJAC(nThreads); }
            // Chunks of JAC this thread is responsible for - INCLUSIVE bounds
            distributeGenomePairs(threadID, nThreads);
#pragma omp barrier
            //
            IdType currentLocalEIndex = m_threadEStarts[threadID];
            findEBlockExtents(threadID, currentLocalEIndex);
            // #pragma omp single
            // { print_e(); }
#pragma omp barrier
            computeEBlockJAC(threadID);
        }
        run_timer.elapsed();
        run_timer.print_elapsed("JAC Construction    : ", std::cout);
        return PFAAI_OK;
    }

    int computeAJI() {
        m_AJI.resize(m_nGenomePairs);
        // PHASE 4: Finalize output
#pragma omp parallel default(none) shared(std::cout)
        {
            int threadID = omp_get_thread_num();
            // AJI
            for (int genomePair = m_threadGPStarts[threadID];
                 genomePair <= m_threadGPEnds[threadID]; genomePair++) {
                m_AJI[genomePair] = m_JAC[genomePair].S / m_JAC[genomePair].N;
            }
        }
        return PFAAI_OK;
    }

    int run() {
        generateTetramerTuples();
        computeJAC();
        computeAJI();
        return PFAAI_OK;
    }
};

#endif  // !PFAAI_RUNNER_HPP
