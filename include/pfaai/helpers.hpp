///
// @file helpers.hpp
// @brief The implementation classes and functions to construct the
//        data structures.
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

#ifndef PAR_FAST_AAI_DATA_HELPERS_H
#define PAR_FAST_AAI_DATA_HELPERS_H

#include <algorithm>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#include "pfaai/psort.hpp"
#include <pfaai/interface.hpp>
#include <pfaai/utils.hpp>

// Class with functions to construct data structures in parallel
template <typename IdType, typename IdPairType, typename IdMatrixType,
          typename JACType>
struct DataStructHelper {

    using IdTripleT = ETriple<IdType>;
    using DataStructIfx =
        DataStructInterface<IdType, IdPairType, IdMatrixType, JACType>;
    using DataBaseIfx =
        const DataBaseInterface<IdType, IdPairType, IdMatrixType>;

    // Parallel construction of the T matrix (nProteins x nGenomes)
    static PFAAI_ERROR_CODE constructT(const DataStructIfx& dataStructPtr,
                                       const DataBaseIfx& inDBIf,
                                       IdMatrixType& inT) {  // NOLINT
        // Assuming inT is a initialized matrix of size nProteins x nGenomes
        std::vector<timer> thTimers;
        std::vector<int> errCodes;
#pragma omp parallel default(none)                                             \
    shared(dataStructPtr, errCodes, thTimers, inDBIf, inT)
        {
            IdType nProteins = dataStructPtr.refProteinSet().size();
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            //
#pragma omp single
            {
                thTimers.resize(nThreads);
                errCodes.resize(nThreads);
            }
            //
            thTimers[threadID].reset();
            int proteinStart = BLOCK_LOW(threadID, nThreads, nProteins);
            int proteinEnd = BLOCK_HIGH(threadID, nThreads, nProteins);

            errCodes[threadID] = inDBIf.proteinTetramerCounts(
                IdPairType(proteinStart, proteinEnd), inT);
            thTimers[threadID].elapsed();
        }

        if (std::any_of(errCodes.begin(), errCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return PFAAI_ERR_CONSTRUCT;
        //
        // printThreadTimes("T construction : ", threadTimers);
        return PFAAI_OK;
    }

    // Parallel construction of the Lc array (size NTETRAMERS)
    static PFAAI_ERROR_CODE constructLc(const DataStructIfx& dataStruct,
                                        const DataBaseIfx& inDBIf,
                                        std::vector<IdType>& Lc) {  // NOLINT
        std::vector<int> errCodes;
        // Lc construction
#pragma omp parallel default(none) shared(errCodes, dataStruct, inDBIf, Lc)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            IdType tetraStart =
                BLOCK_LOW(threadID, nThreads, DataStructIfx::NTETRAMERS);
            IdType tetraEnd =
                BLOCK_HIGH(threadID, nThreads, DataStructIfx::NTETRAMERS);
#pragma omp single
            {
                errCodes.resize(nThreads, 0);
            }

            for (const std::string& protein : dataStruct.refProteinSet()) {
                int qryErrCode = inDBIf.tetramerOccCounts(
                    protein, IdPairType(tetraStart, tetraEnd), Lc);
                if (qryErrCode != SQLITE_OK) {
                    errCodes[threadID] = qryErrCode;
                }
            }
        }
        if (std::any_of(errCodes.begin(), errCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return PFAAI_ERR_CONSTRUCT;
        return PFAAI_OK;
    }

    // Parallel Prefix sum of the arrays
    static void parallelPrefixSum(const std::vector<IdType>& Lc,
                                  std::vector<IdType>& Lp) {  // NOLINT
        // Parallel prefix sum on Lc to construct Lp
        IdType cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan, + : cumulativeSum)
        for (std::size_t i = 0; i < Lc.size(); i++) {
            Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
            cumulativeSum += Lc[i];
        }
    }

    // Parallel construction of the (protein, genome)-pair F array
    //     (size: Total number of tetramer occurrences)
    static PFAAI_ERROR_CODE constructF(const DataStructIfx& dataStruct,
                                       const DataBaseIfx& dbIf,
                                       std::vector<IdPairType>& F) {  // NOLINT
        F.resize(dataStruct.refLp().back() + dataStruct.refLc().back(),
                 IdPairType(-1, -1));
        std::vector<IdType> tetramerStart, tetramerEnd;
        std::vector<int> errCodes;
        std::vector<timer> thTimers;
        //
#pragma omp parallel default(none) shared(                                     \
        tetramerStart, tetramerEnd, errCodes, thTimers, dataStruct, dbIf, F)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            IdType fCount = 0;
#pragma omp single
            {
                errCodes.resize(nThreads, 0);
                thTimers.resize(nThreads);
                distribute_bags_of_tasks(nThreads, IdType(F.size()),
                                         dataStruct.refLc(), dataStruct.slack(),
                                         tetramerStart, tetramerEnd);
            }
            thTimers[threadID].reset();
            IdType fBeginOffset = dataStruct.refLp()[tetramerStart[threadID]];
            errCodes[threadID] = dbIf.proteinSetGPPairs(
                IdPairType(tetramerStart[threadID], tetramerEnd[threadID]),
                F.begin() + fBeginOffset, &fCount);
            thTimers[threadID].elapsed();
        }
        if (std::any_of(errCodes.begin(), errCodes.end(),
                        [](int rc) { return rc != SQLITE_OK; }))
            return PFAAI_ERR_CONSTRUCT;
        //
        // printThreadTimes(" F construction : ", threadTimers);
        return PFAAI_OK;
    }

    // Distribute tasks across nThreads with total number of tasks being
    //   the sum of all the occurrences of the tetramers in ranges
    //    tetramerStart and tetramerEnd
    static void
    distributeTetramers(const DataStructIfx& dataStruct,
                        std::vector<IdType>& tetramerStart,  // NOLINT
                        std::vector<IdType>& tetramerEnd,    // NOLINT
                        const int& nThreads) {
        // Ask 1 thread to carry out the distribution of tasks, i.e.
        // tetramer tuples for all threads
        tetramerStart.resize(nThreads, -1);
        tetramerEnd.resize(nThreads, -1);
        //
        float nTasks = static_cast<float>(dataStruct.refF().size());
        std::vector<int> threadTaskSum(nThreads, 0);
        int nTasksPerThread =
            static_cast<int>((nTasks / nThreads) * (1 + dataStruct.slack()));

        int tid = 0;
        for (IdType tetramer = 0; tetramer < IdType(dataStruct.refLc().size());
             tetramer++) {
            IdType nTetraTasks = dataStruct.refLc()[tetramer];
            if (threadTaskSum[tid] + nTetraTasks <= nTasksPerThread ||
                tid == nThreads - 1) {
                threadTaskSum[tid] += nTetraTasks;
                if (tetramerStart[tid] == -1) {
                    tetramerStart[tid] = tetramer;
                }
                tetramerEnd[tid] = tetramer;
            } else {
                // Move to the next processor and assign the task there
                tid++;
                threadTaskSum[tid] += nTetraTasks;
                tetramerStart[tid] = tetramer;
                tetramerEnd[tid] = tetramer;
            }
        }
    }

    //
    // Count the number of tuples (protein, genome_A, genome_B) that should
    //   occur in the E array for tetramer range tStart and tEnd
    static IdType countTetramerTuples(const DataStructIfx& dsIfx,
                                      const IdType& tStart,
                                      const IdType& tEnd) {
        int nTetraTuples = 0;
        // Ask a thread to compute its own number of elements in its local
        // part of E
        for (IdType tetraID = tStart; tetraID <= tEnd; tetraID++) {
            // INCLUSIVE start and end index of the tetramer block in F
            IdType startIndexInF = dsIfx.refLp()[tetraID];
            IdType endIndexInF = 0;
            if (tetraID < tEnd) {
                endIndexInF = dsIfx.refLp()[tetraID + 1] - 1;
            } else {
                // tetramerID == tetramerEnd
                if (tetraID == dsIfx.nTetramers() - 1) {
                    endIndexInF = dsIfx.refF().size() - 1;
                } else {
                    endIndexInF = dsIfx.refLp()[tetraID + 1] - 1;
                }
            }

            int leftBoundary = startIndexInF;
            int currProteinID = dsIfx.refF()[startIndexInF].first;
            int rightBoundary = startIndexInF;
            int nQryOcc = 0, nTgtOcc = 0;

            while (rightBoundary <= endIndexInF) {
                int rProtienID = dsIfx.refF()[rightBoundary].first;
                int rGenomeID = dsIfx.refF()[rightBoundary].second;
                bool isQuery = dsIfx.isQryGenome(rGenomeID);
                //
                if (rProtienID == currProteinID) {
                    rightBoundary++;
                    //
                    nQryOcc += isQuery;
                    nTgtOcc += !isQuery;
                } else {
                    // int n = rightBoundary - leftBoundary;
                    // nTetraTuples += n * (n - 1) / 2;
                    // assert(rightBoundary - leftBoundary == nQryOcc);
                    int nbdTuples = dsIfx.countGenomePairs(nQryOcc, nTgtOcc);
                    nTetraTuples += nbdTuples;
                    //
                    leftBoundary = rightBoundary;
                    currProteinID = dsIfx.refF()[rightBoundary].first;
                    nQryOcc = nTgtOcc = 0;
                }
            }

            // Complete the last calculation for when the while loop
            // finishes because rightBoundary > endIndexInF, but the
            // calculation for currProteinID is not done yet
            // int n = rightBoundary - leftBoundary;
            // nTetraTuples += n * (n - 1) / 2;
            // assert(rightBoundary - leftBoundary == nQryOcc);
            int nbdTuples = dsIfx.countGenomePairs(nQryOcc, nTgtOcc);
            nTetraTuples += nbdTuples;
        }
        return nTetraTuples;
    }

    //
    // Construct the of tuples (protein, genome_A, genome_B) that
    //   occurs in the 'E' array for tetramer range tStart and tEnd
    static IdType constructTetramerTuples(const DataStructIfx& dsIfx,
                                          const int& threadID,
                                          const IdType& tStart,
                                          const IdType& tEnd,
                                          EParData<IdType>& pE) {  // NOLINT
        //
        IdType currentLocalEIndex = pE.threadEStarts[threadID];
        IdType nTuples = 0;
        const auto& c_Lp = dsIfx.refLp();
        const auto& c_F = dsIfx.refF();
        //
        for (IdType tetramerID = tStart; tetramerID <= tEnd; tetramerID++) {
            // INCLUSIVE start and end index of the tetramer block in F
            IdType startIndexInF = c_Lp[tetramerID];
            IdType endIndexInF = 0;
            if (tetramerID < tEnd) {
                endIndexInF = c_Lp[tetramerID + 1] - 1;
            } else {
                // tetramerID == tetramerEnd
                if (tetramerID == dsIfx.nTetramers() - 1) {
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
                    assert(std::is_sorted(
                        c_F.begin() + leftBoundary, c_F.begin() + rightBoundary,
                        [](const IdPairType& a, const IdPairType& b) {
                            return a.second < b.second;
                        }));
                    for (int i = leftBoundary; i < rightBoundary; i++) {
                        int genomeA_ID = c_F[i].second;
                        if (dsIfx.isQryGenome(genomeA_ID) == false)
                            continue;
                        for (int j = leftBoundary; j < rightBoundary; j++) {
                            IdType genomeB_ID = c_F[j].second;
                            // Skip if genomeA_ID is not one of the query genes
                            if (dsIfx.isValidPair(genomeA_ID, genomeB_ID) ==
                                false)
                                continue;
                            assert(currProteinID >= 0);
                            assert(currProteinID <
                                   IdType(dsIfx.refProteinSet().size()));
                            IdTripleT newElement(currProteinID, genomeA_ID,
                                                 genomeB_ID);
                            pE.E[currentLocalEIndex] = newElement;
                            currentLocalEIndex++;
                            nTuples++;
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
                if (dsIfx.isQryGenome(genomeA_ID) == false)
                    continue;
                for (IdType j = leftBoundary; j < rightBoundary; j++) {
                    IdType genomeB_ID = c_F[j].second;
                    // Skip if genomeA_ID is not one of the query genes
                    if (dsIfx.isValidPair(genomeA_ID, genomeB_ID) == false)
                        continue;
                    IdTripleT newElement(currProteinID, genomeA_ID, genomeB_ID);
                    pE.E[currentLocalEIndex] = newElement;
                    currentLocalEIndex++;
                    nTuples++;
                }
            }
        }
        return nTuples;
    }

    //
    // Construct 'E' array, the array of tuples (protein, genome_A, genome_B)
    //   share a tetramer
    static PFAAI_ERROR_CODE constructE(const DataStructIfx& dsIfx,
                                       EParData<IdType>& dsE) {  // NOLINT
        //
        IdType totalESize(0);
        timer run_timer;
        std::vector<IdType> tetramerStart, tetramerEnd;  // size nThreads
        //  Generate tetramer tuples
#pragma omp parallel default(none)                                             \
    shared(dsIfx, totalESize, tetramerStart, tetramerEnd, dsE)
        {
            int nThreads = omp_get_num_threads();
            int threadID = omp_get_thread_num();
            // Ask 1 thread to carry out the distribution of tasks, i.e.
#pragma omp single
            {
                distributeTetramers(dsIfx, tetramerStart, tetramerEnd,
                                    nThreads);
                dsE.threadESize.resize(nThreads);
                dsE.threadEStarts.resize(nThreads);
            }
            // Tetramers specific to this thread
            int tStart = tetramerStart[threadID], tEnd = tetramerEnd[threadID];
            // Size of tetramer tuples corresponding to this thread
            dsE.threadESize[threadID] =
                countTetramerTuples(dsIfx, tStart, tEnd);
#pragma omp barrier
            // Get total length of E
#pragma omp for reduction(+ : totalESize)
            for (std::size_t i = 0; i < dsE.threadESize.size(); i++) {
                totalESize += dsE.threadESize[i];
            }
#pragma omp single
            {
                dsE.E.resize(totalESize);
            }
            // Get the start indices of each thread-partions of E by parallel
            // prefix sum (exclusive) on _ESize and store it in _EStartIndex
            int cumulativeSum = 0;
#pragma omp simd reduction(inscan, + : cumulativeSum)
            for (std::size_t i = 0; i < dsE.threadESize.size(); i++) {
                dsE.threadEStarts[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
                cumulativeSum += dsE.threadESize[i];
            }
            // Each thread construct its own part of E
            IdType nTuples =
                constructTetramerTuples(dsIfx, threadID, tStart, tEnd, dsE);
            assert(nTuples == dsE.threadESize[threadID]);
#pragma omp barrier
        }
        run_timer.elapsed().print("E construction      : ", std::cout);
        // Parallel Sort E TODO(): Sorting speed is inconsistent, why ?
        timer srt_timer;
#pragma omp single
        {
            parallelMergeSort(dsE.E, 0, dsE.E.size() - 1, 5);
        }
        srt_timer.elapsed().print("E parallel sorting  : ", std::cout);
        return PFAAI_OK;
    }
};

template <typename IdType>
using DefaultDSHelper = DataStructHelper<IdType, DPair<IdType, IdType>,
                                         DMatrix<IdType>, JACTuple<IdType>>;

#endif  //! PAR_FAST_AAI_DATA_HELPERS_H
