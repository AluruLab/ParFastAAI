#ifndef PFAAI_RUNNER_HPP
#define PFAAI_RUNNER_HPP

#include <omp.h>
#include <utility>
#include <vector>

#include "pfaai/data.hpp"
#include "pfaai/psort.hpp"
#include "pfaai/utils.hpp"

template <typename IdType = int> struct ETriple {
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
        if (genomeB != genomeB) {
            return genomeB < genomeB;
        }
        return proteinIndex < proteinIndex;
    }
};

template <typename IdType = int, typename SType = double> struct JACTuple {
    IdType genomeA;
    IdType genomeB;
    SType S;
    IdType N;
};

template <typename IdType, typename IdPairType = std::pair<IdType, IdType>,
          typename SType = double>
class ParFAAIRunner {
    // Data references
    const ParFAAIData<IdType, IdPairType>& c_faaiDataRef;
    const std::vector<int>& c_Lc;
    const std::vector<int>& c_Lp;
    const std::vector<IdPairType>& c_F;
    const std::vector<std::vector<IdType>>& c_T;
    float m_slack;
    IdType m_nTetramers;
    IdType m_nGenomes;
    //
    int m_nGenomePairs;

    // tetramer tuples : Array E construction : Phase 2
    std::vector<IdType> m_tetramerStart, m_tetramerEnd;
    std::vector<ETriple<IdType>> m_E;
    std::vector<IdType> m_EStartIndex;
    std::vector<IdType> m_ESize;

    // Thread starts and end indices for Phase 3
    std::vector<IdType> m_genomePairEStartIndex;
    std::vector<IdType> m_genomePairEEndIndex;
    std::vector<IdType> m_genomePairStartIndex;
    std::vector<IdType> m_genomePairEndIndex;

    // Jaccard Step
    std::vector<JACTuple<IdType, SType>> m_JAC;

    // AJI step
    std::vector<SType> m_AJI;

    int genomePairToJACIndex(int genomeCount, int genomeA, int genomeB) {
        return (genomeCount * genomeA) + genomeB -
               static_cast<int>((genomeA + 2) * (genomeA + 1) / 2);
    }

    int genomePairToJACIndex(int genomeA, int genomeB) {
        // (37/2) * a - a^2 / 2 + b - 1
        // return (37 * genomeA - genomeA * genomeA) / 2 + genomeB - 1;
        return genomePairToJACIndex(m_nGenomes, genomeA, genomeB);
    }

  public:
    explicit ParFAAIRunner(const ParFAAIData<IdType, IdPairType>& fDataRef)
        : c_faaiDataRef(fDataRef), c_Lc(c_faaiDataRef.getLc()),
          c_Lp(c_faaiDataRef.getLp()), c_F(c_faaiDataRef.getF()),
          c_T(c_faaiDataRef.getT()),
          m_slack(c_faaiDataRef.getSlackPercentage()),
          m_nTetramers(c_faaiDataRef.getTetramerCount()),
          m_nGenomes(c_faaiDataRef.getGenomeCount()),
          m_nGenomePairs(m_nGenomes * (m_nGenomes - 1) / 2) {}

  private:
    void distributeTetramers(const int& nThreads) {
        // Ask 1 thread to carry out the distribution of tasks, i.e.
        // tetramer tuples for all threads
        m_tetramerStart.resize(nThreads, -1);
        m_tetramerEnd.resize(nThreads, -1);
        //
        float nTasks = static_cast<float>(c_F.size());
        std::vector<int> taskSumPerThread(nThreads, 0);
        int tasksPerProcessor =
            static_cast<int>((nTasks / nThreads) * (1 + m_slack));

        int tid = 0;
        for (IdType tetramer = 0; tetramer < IdType(c_Lc.size()); tetramer++) {
            IdType nTetraTasks = c_Lc[tetramer];
            if (taskSumPerThread[tid] + nTetraTasks <= tasksPerProcessor ||
                tid == nThreads - 1) {
                taskSumPerThread[tid] += nTetraTasks;
                if (m_tetramerStart[tid] == -1) {
                    m_tetramerStart[tid] = tetramer;
                }
                m_tetramerEnd[tid] = tetramer;
            } else {
                // Move to the next processor and assign the task there
                tid++;
                taskSumPerThread[tid] += nTetraTasks;
                m_tetramerStart[tid] = tetramer;
                m_tetramerEnd[tid] = tetramer;
            }
        }
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
        IdType currentLocalEIndex = m_EStartIndex[threadID];
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
                            ETriple newElement(currProteinID, genomeA_ID,
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
                    ETriple newElement(currProteinID, genomeA_ID, genomeB_ID);
                    m_E[currentLocalEIndex] = newElement;
                    currentLocalEIndex++;
                }
            }
        }
    }

    int generateTetramerTuples() {
        IdType totalESize(0);
        timer run_timer;
        // PHASE 2: Generate tetramer tuples
#pragma omp parallel default(none) shared(totalESize)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
#pragma omp single
            {
                distributeTetramers(nThreads);
                m_ESize.resize(nThreads);
                m_EStartIndex.resize(nThreads);
            }
            // Tetramers specific to this thread
            int tStart = m_tetramerStart[threadID], tEnd = m_tetramerEnd[threadID];
            // Size of tetramer tuples corresponding to this thread
            m_ESize[threadID] = countTetramerTuples(tStart, tEnd);

#pragma omp barrier
            // Get total length of E
#pragma omp for reduction(+ : totalESize)
            for (int i = 0; i < m_ESize.size(); i++) {
                totalESize += m_ESize[i];
            }

#pragma omp single
            { m_E.resize(totalESize); }
            // Get the start indices of each thread-partions of E by parallel
            // prefix sum (exclusive) on _ESize and store it in _EStartIndex
            int cumulativeSum = 0;
#pragma omp simd reduction(inscan, + : cumulativeSum)
            for (int i = 0; i < m_ESize.size(); i++) {
                m_EStartIndex[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
                cumulativeSum += m_ESize[i];
            }
            // Each thread construct its own part of E
            constructTetramerTuples(threadID, tStart, tEnd);
#pragma omp barrier
        }
        run_timer.print_elapsed("Time taken to construct E: ", std::cout);
        // Parallel Sort E
        run_timer.reset();
#pragma omp parallel default(none) shared(m_E)
        { parallelMergeSort(m_E, 0, m_E.size() - 1, 5); }
        run_timer.print_elapsed("Time taken to sort E: ", std::cout);
        return PFAAI_OK;
    }

    void prepJAC(const int& nThreads) {
        m_JAC.resize(m_nGenomePairs);
        m_genomePairEStartIndex.resize(m_nGenomePairs);
        m_genomePairEEndIndex.resize(m_nGenomePairs);
        for (int i = 0, gA = 0, gB = 1; i < m_JAC.size(); i++) {
            m_JAC[i].genomeA = gA;
            m_JAC[i].genomeB = gB;
            m_JAC[i].S = 0.0;
            m_JAC[i].N = 0;

            if (gB == m_nGenomes - 1) {
                gA += 1;
                gB = gA + 1;
            } else {
                gB += 1;
            }
        }
        m_genomePairStartIndex.resize(nThreads);
        m_genomePairEndIndex.resize(nThreads);
        std::cout << "Prepped JAC" << std::endl;
    }

    void distributeGP(const int& threadID, const int& nThreads) {
        if (threadID < m_nGenomePairs % nThreads) {
            // Each will be responsible for totalGenomePairs /
            // totalNumThreads + 1 pairs
            m_genomePairStartIndex[threadID] =
                threadID * (m_nGenomePairs / nThreads + 1);
            m_genomePairEndIndex[threadID] =
                (threadID + 1) * (m_nGenomePairs / nThreads + 1) - 1;
        } else {
            int buffer =
                (m_nGenomePairs % nThreads) * (m_nGenomePairs / nThreads + 1);
            m_genomePairStartIndex[threadID] =
                buffer + (threadID - m_nGenomePairs % nThreads) *
                             (m_nGenomePairs / nThreads);
            m_genomePairEndIndex[threadID] =
                buffer +
                (threadID - m_nGenomePairs % nThreads + 1) *
                    (m_nGenomePairs / nThreads) -
                1;
        }
    }

    void findEBlockStarts(const int& currentLocalEIndex) {
        // We need to know where in the sorted E does each genome pair start
        // and end Each thread can look though its local chunk of E and help
        // fill in the 2 arrays: genomePairEStartIndex, genomePairEEndIndex

        // Edge case - First genome pair seen in the local E chunk. It might
        // not start here
        if (currentLocalEIndex > 0) {
            ETriple firstElement = m_E[currentLocalEIndex];
            ETriple lastElementPrevThread = m_E[currentLocalEIndex - 1];
            if (firstElement.genomeA != lastElementPrevThread.genomeA ||
                firstElement.genomeB != lastElementPrevThread.genomeB) {
                // Then this means we hold the beginning of this genome pair
                // Safely write this information into the
                // genomePairEStartIndex array
                int genomePairIndexInJAC = genomePairToJACIndex(
                    firstElement.genomeA, firstElement.genomeB);
                m_genomePairEStartIndex[genomePairIndexInJAC] =
                    currentLocalEIndex;
            }  // Otherwise, let the previous processor fill in the starting
               // index
        } else {
            // currentLocalEIndex == 0
            // Safely write the beginning of the current genome pair index
            // into genomePairEStartIndex array
            ETriple firstElement = m_E[currentLocalEIndex];
            int genomePairIndexInJAC = genomePairToJACIndex(
                firstElement.genomeA, firstElement.genomeB);
            m_genomePairEStartIndex[genomePairIndexInJAC] = 0;
        }
    }

    void findEBlockEnds(const int& threadID, int& currentLocalEIndex) {
        while (currentLocalEIndex <
               m_EStartIndex[threadID] + m_ESize[threadID]) {
            int currentGenomeA = m_E[currentLocalEIndex].genomeA;
            int currentGenomeB = m_E[currentLocalEIndex].genomeB;

            while (currentLocalEIndex <
                       m_EStartIndex[threadID] + m_ESize[threadID] &&
                   m_E[currentLocalEIndex].genomeA == currentGenomeA &&
                   m_E[currentLocalEIndex].genomeB == currentGenomeB) {
                currentLocalEIndex++;
            }

            if (currentLocalEIndex <
                m_EStartIndex[threadID] + m_ESize[threadID]) {
                // Found the end index of the current genome pair
                // Safely write this information into the
                // genomePairEEndIndex array
                int genomePairIndexInJAC =
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

                if (currentLocalEIndex < m_E.size()) {
                    ETriple firstElementOfNextThread = m_E[currentLocalEIndex];

                    if (firstElementOfNextThread.genomeA != currentGenomeA ||
                        firstElementOfNextThread.genomeB != currentGenomeB) {
                        // We have seen the end of the current genome pair
                        // Safely write this information into
                        // genomePairEEndIndex array
                        int genomePairIndexInJAC = genomePairToJACIndex(
                            currentGenomeA, currentGenomeB);
                        m_genomePairEEndIndex[genomePairIndexInJAC] =
                            currentLocalEIndex - 1;
                    }  // Otherwise, this information will be filled by the
                       // next processor
                } else {
                    // currentLocalEIndex == E.size()
                    int genomePairIndexInJAC =
                        genomePairToJACIndex(currentGenomeA, currentGenomeB);
                    m_genomePairEEndIndex[genomePairIndexInJAC] =
                        m_E.size() - 1;
                }
            }
        }
    }

    void dbg_print() {
        std::cout << "Finished processing start and end index of pairs of "
                     "genomes in E "
                  << std::endl;
        std::cout << " First, E(Ga, Gb, ProteinID) is : ";
        for (int i = 0; i < m_E.size(); i++) {
            std::cout << "(" << m_E[i].genomeA << ", " << m_E[i].genomeB << ", "
                      << m_E[i].proteinIndex << ")   ";
        }
        std::cout << std::endl << std::endl;

        std::cout << "Each pair of genome STARTS at:" << std::endl;
        for (int i = 0; i < m_genomePairEStartIndex.size(); i++) {
            std::cout << "Genome (" << m_JAC[i].genomeA << ", "
                      << m_JAC[i].genomeB << "): " << m_genomePairEStartIndex[i]
                      << std::endl;
        }
        std::cout << std::endl << std::endl;

        std::cout << "Each pair of genome ENDS at:" << std::endl;
        for (int i = 0; i < m_genomePairEEndIndex.size(); i++) {
            std::cout << "Genome (" << m_JAC[i].genomeA << ", "
                      << m_JAC[i].genomeB << "): " << m_genomePairEEndIndex[i]
                      << std::endl;
        }
        std::cout << std::endl << std::endl;
    }

    void computeEBlockJAC(const int& threadID) {
        for (int genomePair = m_genomePairStartIndex[threadID];
             genomePair <= m_genomePairEndIndex[threadID]; genomePair++) {
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
                        (double)(BkLength) /
                        (double)(c_T[currProteinID][currGenomeA] +
                                 c_T[currProteinID][currGenomeB] + BkLength);
                    S += J_Pi_Ga_Gb;
                    N += 1;

                    // Move on to the next Bk subblock
                    blockBkStart = blockBkEnd;
                } else {
                    // blockBkEnd > blockBlEnd
                    // Finish the last computation
                    int BkLength = blockBkEnd - blockBkStart;
                    double J_Pi_Ga_Gb =
                        (double)(BkLength) /
                        (double)(c_T[currProteinID][currGenomeA] +
                                 c_T[currProteinID][currGenomeB] + BkLength);
                    S += J_Pi_Ga_Gb;
                    N += 1;
                }
            }
            // DOUBLE CHECK THIS PART
            m_JAC[genomePair].S += S;
            m_JAC[genomePair].N += N;
        }
    }

    int computeJAC() {
        timer run_timer;
        m_nGenomePairs = m_nGenomes * (m_nGenomes - 1) / 2;
        // Prepare the JAC vector
        // PHASE 3: Compute the Jaccard Coefficient values
#pragma omp parallel default(none)                                             
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
#pragma omp single
            { prepJAC(nThreads); }
            // Chunks of JAC this thread is responsible for - INCLUSIVE bounds
            distributeGP(threadID, nThreads);
#pragma omp barrier
            IdType currentLocalEIndex = m_EStartIndex[threadID];
            findEBlockStarts(currentLocalEIndex);
#pragma omp barrier
            findEBlockEnds(threadID, currentLocalEIndex);

#pragma omp barrier
            // #pragma omp single
            // {
            //  dbg_print();
            // }
            computeEBlockJAC(threadID);
        }
        run_timer.print_elapsed("Time taken to construct JAC: ", std::cout);
        return PFAAI_OK;
    }

    int computeAJI() {
        // PHASE 4: Finalize output
#pragma omp parallel default(none) shared( std::cout)
        {
            int threadID = omp_get_thread_num();
#pragma omp single
            {
                std::cout << "Finished JAC construction" << std::endl;
                m_AJI.resize(m_nGenomePairs);
                std::cout << "Prepped AJI" << std::endl;
            }
            // AJI
            for (int genomePair = m_genomePairStartIndex[threadID];
                 genomePair <= m_genomePairEndIndex[threadID]; genomePair++) {
                m_AJI[genomePair] = m_JAC[genomePair].S / m_JAC[genomePair].N;
            }
        }
        return PFAAI_OK;
    }


  public:
    int run() {
        generateTetramerTuples();
        computeJAC();
        computeAJI();
        return PFAAI_OK;
    }
};

#endif  // !PFAAI_RUNNER_HPP
