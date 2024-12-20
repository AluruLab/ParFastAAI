#ifndef PFAAI_RUNNER_HPP
#define PFAAI_RUNNER_HPP

#include <cmath>
#include <cstdlib>
#include <fmt/format.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "pfaai/interface.hpp"
#include "pfaai/utils.hpp"

template <typename IdType, typename ValueType> class ParFAAIImpl {
  public:
    using IdPairType = DPair<IdType, IdType>;
    using IdMatrixType = DMatrix<IdType>;
    using JACType = JACTuple<IdType, ValueType>;
    using DSIfx = DefaultDataStructInterface<IdType>;
    using IdTriple = ETriple<IdType>;

  private:
    // Data references
    const DSIfx& c_pfDSRef;
    const std::vector<int>& c_Lc;
    const std::vector<int>& c_Lp;
    const std::vector<IdPairType>& c_F;
    const IdMatrixType& c_T;
    const EParData<IdType>& c_pE;
    //
    float m_slack;
    IdType m_nTetramers;
    IdType m_nGenomePairs;

    // Starts and end indices for Phase 3
    // Thread distributions for genome pair ~(m_nGenomePairs/p)
    std::vector<IdType> m_threadGPStarts, m_threadGPEnds;  // of size nThreads
    // Start and end corresponding to genome pairs
    std::vector<IdType> m_genomePairEStartIndex,
        m_genomePairEEndIndex;  // of size m_nGenomePairs

    // Jaccard data
    std::vector<JACType> m_JAC;  // of size m_nGenomePairs

    // AJI data
    std::vector<ValueType> m_AJI;  // of size m_nGenomePairs

  public:
    explicit ParFAAIImpl(const DSIfx& fDataRef)
        : c_pfDSRef(fDataRef), c_Lc(c_pfDSRef.refLc()), c_Lp(c_pfDSRef.refLp()),
          c_F(c_pfDSRef.refF()), c_T(c_pfDSRef.refT()), c_pE(c_pfDSRef.refE()),
          m_slack(c_pfDSRef.slack()), m_nTetramers(c_pfDSRef.nTetramers()),
          m_nGenomePairs(c_pfDSRef.nGenomePairs()) {}

    const std::vector<JACTuple<IdType>>& getJAC() const { return m_JAC; }
    const std::vector<ValueType>& getAJI() const { return m_AJI; }

    inline IdType genomePairToIndex(IdType genomeA, IdType genomeB) const {
        return c_pfDSRef.genomePairToIndex(genomeA, genomeB);
    }

  private:
    void prepJAC(const int& nThreads) {
        m_genomePairEStartIndex.resize(m_nGenomePairs);
        m_genomePairEEndIndex.resize(m_nGenomePairs);
        m_JAC = c_pfDSRef.initJAC();
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
            IdTriple firstElement = c_pE.E[currentLocalEIndex];
            IdTriple lastElementPrevThread = c_pE.E[currentLocalEIndex - 1];
            if (firstElement.genomeA != lastElementPrevThread.genomeA ||
                firstElement.genomeB != lastElementPrevThread.genomeB) {
                // Then this means we hold the beginning of this genome pair
                // Safely write this information into the
                // genomePairEStartIndex array
                IdType genomePairIndexInJAC = genomePairToIndex(
                    firstElement.genomeA, firstElement.genomeB);
                m_genomePairEStartIndex[genomePairIndexInJAC] =
                    currentLocalEIndex;
            }  // Otherwise, let the previous processor fill in the starting
               // index
        } else {
            // currentLocalEIndex == 0
            // Safely write the beginning of the current genome pair index
            // into genomePairEStartIndex array
            IdTriple firstElement = c_pE.E[currentLocalEIndex];
            IdType genomePairIndexInJAC =
                genomePairToIndex(firstElement.genomeA, firstElement.genomeB);
            m_genomePairEStartIndex[genomePairIndexInJAC] = 0;
        }

        while (currentLocalEIndex <
               c_pE.threadEStarts[threadID] + c_pE.threadESize[threadID]) {
            IdType currentGenomeA = c_pE.E[currentLocalEIndex].genomeA;
            IdType currentGenomeB = c_pE.E[currentLocalEIndex].genomeB;

            while (currentLocalEIndex < c_pE.threadEStarts[threadID] +
                                            c_pE.threadESize[threadID] &&
                   c_pE.E[currentLocalEIndex].genomeA == currentGenomeA &&
                   c_pE.E[currentLocalEIndex].genomeB == currentGenomeB) {
                currentLocalEIndex++;
            }

            if (currentLocalEIndex <
                c_pE.threadEStarts[threadID] + c_pE.threadESize[threadID]) {
                // Found the end index of the current genome pair
                // Safely write this information into the
                // genomePairEEndIndex array
                IdType genomePairIndexInJAC =
                    genomePairToIndex(currentGenomeA, currentGenomeB);
                m_genomePairEEndIndex[genomePairIndexInJAC] =
                    currentLocalEIndex - 1;

                // Move on to the next pair
                // This means we are seeing the start index of a new genome
                // pair Safely write this information into the
                // genomePairEStartIndex array
                genomePairIndexInJAC =
                    genomePairToIndex(c_pE.E[currentLocalEIndex].genomeA,
                                      c_pE.E[currentLocalEIndex].genomeB);
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

                if (currentLocalEIndex < static_cast<IdType>(c_pE.E.size())) {
                    IdTriple firstElementOfNextThread =
                        c_pE.E[currentLocalEIndex];

                    if (firstElementOfNextThread.genomeA != currentGenomeA ||
                        firstElementOfNextThread.genomeB != currentGenomeB) {
                        // We have seen the end of the current genome pair
                        // Safely write this information into
                        // genomePairEEndIndex array
                        IdType genomePairIndexInJAC =
                            genomePairToIndex(currentGenomeA, currentGenomeB);
                        m_genomePairEEndIndex[genomePairIndexInJAC] =
                            currentLocalEIndex - 1;
                    }  // Otherwise, this information will be filled by the
                       // next processor
                } else {
                    // currentLocalEIndex == E.size()
                    IdType genomePairIndexInJAC =
                        genomePairToIndex(currentGenomeA, currentGenomeB);
                    m_genomePairEEndIndex[genomePairIndexInJAC] =
                        c_pE.E.size() - 1;
                }
            }
        }
    }

    void print_e() {
        std::cout << "E array : " << std::endl;
        for (int i = 0; i < c_pE.E.size(); i++) {
            fmt::print("{}", c_pE.E[i]);
        }
        std::cout << std::endl << std::endl;

        std::cout << "Genome Pair Extents:" << std::endl;
        for (int i = 0; i < m_genomePairEStartIndex.size(); i++) {
            fmt::print("GP {} : [{}, {}]", c_pE.E[i],
                       m_genomePairEStartIndex[i], m_genomePairEEndIndex[i]);
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
                int currProteinID = c_pE.E[blockBkStart].proteinIndex;

                while (blockBkEnd <= blockBlEnd &&
                       c_pE.E[blockBkEnd].proteinIndex == currProteinID) {
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
            {
                prepJAC(nThreads);
            }
            // Chunks of JAC this thread is responsible for - INCLUSIVE bounds
            distributeGenomePairs(threadID, nThreads);
#pragma omp barrier
            //
            IdType currentLocalEIndex = c_pE.threadEStarts[threadID];
            findEBlockExtents(threadID, currentLocalEIndex);
            // #pragma omp single
            // { print_e(); }
#pragma omp barrier
            computeEBlockJAC(threadID);
        }
        run_timer.elapsed().print("JAC Construction    : ", std::cout);
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
        computeJAC();
        computeAJI();
        return PFAAI_OK;
    }
};

#endif  // !PFAAI_RUNNER_HPP
