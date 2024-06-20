/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "sqlite3.h"

#define SUCCESS 0
#define ERROR_SQLITE_DATABASE 1
#define ERROR_SQLITE_MEMORY_ALLOCATION 2
#define ERROR_IN_CONSTRUCTION 3
#define TETRAMER_COUNT 160000
#define SLACK_PERCENTAGE 0.0

struct ETriple {
	int proteinIndex;
	int genomeA;
	int genomeB;

	ETriple() : proteinIndex(-1), genomeA(-1), genomeB(-1) {}
	ETriple(const int proteinIndexVal, const int genomeAVal, const int genomeBVal) : proteinIndex(proteinIndexVal), genomeA(genomeAVal), genomeB(genomeBVal) {}
};

struct JACTuple {
	int genomeA;
	int genomeB;
	double S;
	int N;
};

int queryMetadataFromDatabase(sqlite3* db, int& GENOMECOUNT, std::vector <std::string>& proteinSet);
int constructLcandLp(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<int>& Lc, std::vector<int>& Lp);
int constructF(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::pair<int, int>>& F,
	std::vector<int>& Lc, std::vector<int>& Lp);
int constructT(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::vector<int>>& T, const int GENOMECOUNT);
bool customSortE(const ETriple& firstElement, const ETriple& secondElement);
void parallelMergeSort(std::vector<ETriple>& E, int left, int right, int serialThreshold);
int genomePairToJACIndex(int genomeA, int genomeB);

int parallelfastaai(const std::string pathToDatabase)
{
	sqlite3* database;
	int errorCode = sqlite3_open(pathToDatabase.c_str(), &database);

	if (errorCode != SQLITE_OK) {
		std::cerr << "Error in opening " << pathToDatabase << std::endl;
		std::cerr << "The error was: " << sqlite3_errstr(errorCode) << std::endl;
		sqlite3_close(database);
		return ERROR_SQLITE_DATABASE;
	}
	else if (database == nullptr) {
		std::cerr << "SQLite is unable to allocate memory for the database " << pathToDatabase << std::endl;
		sqlite3_close(database);
		return ERROR_SQLITE_MEMORY_ALLOCATION;
	}

	/** PHASE 1: Construction of the data structures **/

	int GENOMECOUNT = -1;
	std::vector<std::string> proteinSet;

	errorCode = queryMetadataFromDatabase(database, GENOMECOUNT, proteinSet);

	if (errorCode != SUCCESS) {
		std::cerr << "Error in querying metadata from database " << pathToDatabase << std::endl;
		std::cerr << "ERROR CODE = " << errorCode << std::endl;
		sqlite3_close(database);
		return ERROR_SQLITE_DATABASE;
	}

	// This is all the protein for modified_xantho_fastaai2.db
	/*std::vector<std::string> proteinSet = { "pf00411.19", "pf00237.19", "pf01016.19", "pf02033.18", "pf00347.23", "pf00119.20",
				"pf00297.22", "pf02601.15", "pf00318.20", "pf02367.17", "pf00825.18", "pf02410.15",
				"pf00406.22", "pf00380.19", "pf00213.18", "pf05221.17", "pf00252.18", "pf00177.21",
				"pf00709.21", "pf00312.22", "pf01192.22", "pf06026.14", "pf01649.18", "pf00572.18",
				"pf00338.22", "pf01142.18", "pf01746.21", "pf01632.19", "pf17136.4", "pf00164.25",
				"pf01808.18", "pf00750.19", "pf00889.19", "pf01196.19", "pf01250.17", "pf00162.19",
				"pf01725.16", "pf01668.18", "pf00749.21", "pf00238.19", "pf02565.15", "pf01715.17",
				"pf00203.21", "pf00828.19", "pf00573.22", "pf00121.18", "pf13393.6", "pf00276.20",
				"pf01139.17", "pf00886.19", "pf00189.20", "pf01176.19", "pf00687.21", "pf01025.19",
				"pf03948.14", "pf00829.21", "pf01195.19", "pf03840.14", "pf00231.19", "pf01351.18",
				"pf01264.21", "pf01193.24", "pf00410.19", "pf00584.20", "pf01245.20", "pf02130.17",
				"pf02699.15", "pf01765.19", "pf01783.23", "pf00281.19", "pf00416.22", "pf00366.20",
				"pf00344.20", "pf00831.23", "pf00334.19", "pf00830.19", "pf00861.22", "pf00453.18",
				"pf00181.23", "pf03652.15" };*/

	std::vector<int> Lc;
	std::vector<int> Lp;

	auto startTime = std::chrono::high_resolution_clock::now();
	errorCode = constructLcandLp(database, proteinSet, Lc, Lp);
	auto endTime = std::chrono::high_resolution_clock::now();

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing Lc and Lp, error code = " << errorCode << std::endl;
		sqlite3_close(database);
		return errorCode;
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
	auto minutes = duration_seconds.count() / 60;
	auto seconds = duration_seconds.count() % 60;
	std::cout << "Time taken to construct Lc, Lp: " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds)." << std::endl;

	std::vector<std::pair<int, int>> F(Lp[TETRAMER_COUNT - 1] + Lc[TETRAMER_COUNT - 1], std::make_pair(-1, -1));

	startTime = std::chrono::high_resolution_clock::now();
	errorCode = constructF(database, proteinSet, F, Lc, Lp);
	endTime = std::chrono::high_resolution_clock::now();

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing F, error code = " << errorCode << std::endl;
		sqlite3_close(database);
		return errorCode;
	}

	duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
	minutes = duration_seconds.count() / 60;
	seconds = duration_seconds.count() % 60;
	std::cout << "Time taken to construct F: " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds)." << std::endl;

	std::vector<std::vector<int>> T;
	startTime = std::chrono::high_resolution_clock::now();
	errorCode = constructT(database, proteinSet, T, GENOMECOUNT);
	endTime = std::chrono::high_resolution_clock::now();

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing T, error code = " << errorCode << std::endl;
		sqlite3_close(database);
		return errorCode;
	}

	duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
	duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
	minutes = duration_seconds.count() / 60;
	seconds = duration_seconds.count() % 60;
	std::cout << "Time taken to construct T: " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds)." << std::endl;

	sqlite3_close(database);

	std::vector<int> tetramerStartDistribution;
	std::vector<int> tetramerEndDistribution;
	int totalNumThreads = -1;
	float slack_percentage = SLACK_PERCENTAGE;

	std::vector<ETriple> E;
	std::vector<int> EChunkSize;
	std::vector<int> EChunkStartIndex;
	int ESize = 0;

	std::vector<JACTuple> JAC;
	std::vector<int> genomePairEStartIndex;
	std::vector<int> genomePairEEndIndex;

	std::vector<double> AJI;

#pragma omp parallel default(none) \
	shared(Lc, Lp, F, T, E, JAC, AJI, GENOMECOUNT, EChunkSize, EChunkStartIndex, genomePairEStartIndex, genomePairEEndIndex, tetramerStartDistribution, tetramerEndDistribution, slack_percentage, tota
	{
		/** PHASE 2: Generate tetramer tuples **/

		// Ask 1 thread to carry out the distribution of tasks, i.e. tetramer tuples for all threads
#pragma omp single
		{
			startTime = std::chrono::high_resolution_clock::now();
			totalNumThreads = omp_get_num_threads();
			tetramerStartDistribution.resize(totalNumThreads, -1);
			tetramerEndDistribution.resize(totalNumThreads, -1);

			int totalNumberOfTasks = F.size();
			std::vector<int> taskSumPerThread(totalNumThreads, 0);
			int targetTasksPerProcessor = ((float)(totalNumberOfTasks) / totalNumThreads) * (1 + slack_percentage);

			int currentThread = 0;
			for (int tetramer = 0; tetramer < Lc.size(); tetramer++) {
				int task = Lc[tetramer];
				if (taskSumPerThread[currentThread] + task <= targetTasksPerProcessor || currentThread == totalNumThreads - 1) {
					taskSumPerThread[currentThread] += task;
					if (tetramerStartDistribution[currentThread] == -1) {
						tetramerStartDistribution[currentThread] = tetramer;
					}
					tetramerEndDistribution[currentThread] = tetramer;
				}
				else {
					// Move to the next processor and assign the task there
					currentThread++;
					taskSumPerThread[currentThread] += task;
					tetramerStartDistribution[currentThread] = tetramer;
					tetramerEndDistribution[currentThread] = tetramer;
				}
			}

			EChunkSize.resize(totalNumThreads);
			EChunkStartIndex.resize(totalNumThreads);
		}

		// Tetramer tuples specific to this thread
		int threadID = omp_get_thread_num();
		int tetramerStart = tetramerStartDistribution[threadID];
		int tetramerEnd = tetramerEndDistribution[threadID];

		int countElementInE = 0;
		// Ask a thread to compute its own number of elements in its local part of E
		for (int tetramerID = tetramerStart; tetramerID <= tetramerEnd; tetramerID++) {
			// INCLUSIVE start and end index of the tetramer block in F
			int startIndexInF = Lp[tetramerID];
			int endIndexInF = 0;
			if (tetramerID < tetramerEnd) {
				endIndexInF = Lp[tetramerID + 1] - 1;
			}
			else {
				// tetramerID == tetramerEnd
				if (tetramerID == TETRAMER_COUNT - 1) {
					endIndexInF = F.size() - 1;
				}
				else {
					endIndexInF = Lp[tetramerID + 1] - 1;
				}
			}

			int currProteinID = F[startIndexInF].first;
			int leftBoundary = startIndexInF;
			int rightBoundary = startIndexInF;

			while (rightBoundary <= endIndexInF) {
				if (F[rightBoundary].first == currProteinID) {
					rightBoundary++;
				}
				else {
					int n = rightBoundary - leftBoundary;
					countElementInE += n * (n - 1) / 2;

					currProteinID = F[rightBoundary].first;
					leftBoundary = rightBoundary;
				}
			}

			// Complete the last calculation for when the while loop finishes because rightBoundary > endIndexInF,
			// but the calculation for currProteinID is not done yet
			int n = rightBoundary - leftBoundary;
			countElementInE += n * (n - 1) / 2;
		}

		EChunkSize[threadID] = countElementInE;

#pragma omp barrier

		// Get the length of E
#pragma omp for reduction(+: ESize)
		for (int i = 0; i < EChunkSize.size(); i++) {
			ESize += EChunkSize[i];
		}

#pragma omp single
		{
			E.resize(ESize);
		}

		// Get the beginning indices of each local chunks of E by parallel prefix sum (exclusive) on EChunkSize and store it in EChunkStartIndex
		int cumulativeSum = 0;
#pragma omp simd reduction(inscan, +:cumulativeSum)
		for (int i = 0; i < EChunkSize.size(); i++) {
			EChunkStartIndex[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
			cumulativeSum += EChunkSize[i];
		}

		int currentLocalEIndex = EChunkStartIndex[threadID];

		// Now we will let each thread go ahead and construct its own part of E
		for (int tetramerID = tetramerStart; tetramerID <= tetramerEnd; tetramerID++) {
			// INCLUSIVE start and end index of the tetramer block in F
			int startIndexInF = Lp[tetramerID];
			int endIndexInF = 0;
			if (tetramerID < tetramerEnd) {
				endIndexInF = Lp[tetramerID + 1] - 1;
			}
			else {
				// tetramerID == tetramerEnd
				if (tetramerID == TETRAMER_COUNT - 1) {
					endIndexInF = F.size() - 1;
				}
				else {
					endIndexInF = Lp[tetramerID + 1] - 1;
				}
			}

			int currProteinID = F[startIndexInF].first;

			// leftBoundary and rightBoundary defines the block B_l
			// leftBoundary is inclusive, rightBoundary is exclusive
			int leftBoundary = startIndexInF;
			int rightBoundary = startIndexInF;

			while (rightBoundary <= endIndexInF) {
				if (F[rightBoundary].first == currProteinID) {
					rightBoundary++;
				}
				else {
					for (int i = leftBoundary; i < rightBoundary; i++) {
						int genomeA_ID = F[i].second;
						for (int j = i + 1; j < rightBoundary; j++) {
							int genomeB_ID = F[j].second;
							ETriple newElement(currProteinID, genomeA_ID, genomeB_ID);
							E[currentLocalEIndex] = newElement;
							currentLocalEIndex++;
						}
					}

					currProteinID = F[rightBoundary].first;
					leftBoundary = rightBoundary;
				}
			}

			// Complete the last calculation for when the while loop finishes because rightBoundary > endIndexInF,
			// but the calculation for currProteinID is not done yet
			for (int i = leftBoundary; i < rightBoundary; i++) {
				int genomeA_ID = F[i].second;
				for (int j = i + 1; j < rightBoundary; j++) {
					int genomeB_ID = F[j].second;
					ETriple newElement(currProteinID, genomeA_ID, genomeB_ID);
					E[currentLocalEIndex] = newElement;
					currentLocalEIndex++;
				}
			}
		}

#pragma omp barrier

#pragma omp single
		{
			endTime = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
			auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime);
			auto minutes = duration_seconds.count() / 60;
			auto seconds = duration_seconds.count() % 60;
			std::cout << "Time taken to construct E: " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds)." << std::endl;
			// Sort E
			parallelMergeSort(E, 0, E.size() - 1, 5);

			std::cout << "Done with parallel sorting E" << std::endl;

			//// Print out E
			//std::cout << "E array ordered by thread chunks: " << std::endl;
			//for (int i = 0; i < totalNumThreads; i++) {
			//	std::cout << "Thread " << i << ": ";
			//	int startIndex = EChunkStartIndex[i];
			//	int chunkSize = EChunkSize[i];

			//	for (int j = 0; j < chunkSize; j++) {
			//		std::cout << "(" << E[startIndex + j].genomeA << ", " << E[startIndex + j].genomeB << ", " << E[startIndex + j].proteinIndex << ")   ";
			//	}
			//	std::cout << std::endl << std::endl;
			//}
		}

		/** PHASE 3: Compute the Jaccard Coefficient values **/

		// Prepare the JAC vector
		int totalGenomePairs = GENOMECOUNT * (GENOMECOUNT - 1) / 2;

#pragma omp single
		{
			JAC.resize(totalGenomePairs);
			genomePairEStartIndex.resize(totalGenomePairs);
			genomePairEEndIndex.resize(totalGenomePairs);

			int genomeA = 0;
			int genomeB = 1;

			for (int i = 0; i < JAC.size(); i++) {
				JAC[i].genomeA = genomeA;
				JAC[i].genomeB = genomeB;
				JAC[i].S = 0.0;
				JAC[i].N = 0;

				if (genomeB == GENOMECOUNT - 1) {
					genomeA += 1;
					genomeB = genomeA + 1;
				}
				else {
					genomeB += 1;
				}
			}
			std::cout << "Prepped JAC" << std::endl;
		}

		// Chunks of JAC this thread is responsible for - INCLUSIVE bounds
		int genomePairStartIndex = -1;
		int genomePairEndIndex = -1;

		if (threadID < totalGenomePairs % totalNumThreads) {
			// Each will be responsible for totalGenomePairs / totalNumThreads + 1 pairs
			genomePairStartIndex = threadID * (totalGenomePairs / totalNumThreads + 1);
			genomePairEndIndex = (threadID + 1) * (totalGenomePairs / totalNumThreads + 1) - 1;
		}
		else {
			int buffer = (totalGenomePairs % totalNumThreads) * (totalGenomePairs / totalNumThreads + 1);
			genomePairStartIndex = buffer + (threadID - totalGenomePairs % totalNumThreads) * (totalGenomePairs / totalNumThreads);
			genomePairEndIndex = buffer + (threadID - totalGenomePairs % totalNumThreads + 1) * (totalGenomePairs / totalNumThreads) - 1;
		}

#pragma omp barrier

		// We need to know where in the sorted E does each genome pair start and end
		// Each thread can look though its local chunk of E and help fill in the 2 arrays: genomePairEStartIndex, genomePairEEndIndex
		currentLocalEIndex = EChunkStartIndex[threadID];

		// Edge case - First genome pair seen in the local E chunk. It might not start here
		if (currentLocalEIndex > 0) {
			ETriple firstElement = E[currentLocalEIndex];
			ETriple lastElementPrevThread = E[currentLocalEIndex - 1];
			if (firstElement.genomeA != lastElementPrevThread.genomeA || firstElement.genomeB != lastElementPrevThread.genomeB) {
				// Then this means we hold the beginning of this genome pair
				// Safely write this information into the genomePairEStartIndex array
				int genomePairIndexInJAC = genomePairToJACIndex(firstElement.genomeA, firstElement.genomeB);
				genomePairEStartIndex[genomePairIndexInJAC] = currentLocalEIndex;
			} // Otherwise, let the previous processor fill in the starting index
		}
		else {
			// currentLocalEIndex == 0
			// Safely write the beginning of the current genome pair index into genomePairEStartIndex array
			ETriple firstElement = E[currentLocalEIndex];
			int genomePairIndexInJAC = genomePairToJACIndex(firstElement.genomeA, firstElement.genomeB);
			genomePairEStartIndex[genomePairIndexInJAC] = 0;
		}

#pragma omp barrier

		while (currentLocalEIndex < EChunkStartIndex[threadID] + EChunkSize[threadID]) {
			int currentGenomeA = E[currentLocalEIndex].genomeA;
			int currentGenomeB = E[currentLocalEIndex].genomeB;

			while (currentLocalEIndex < EChunkStartIndex[threadID] + EChunkSize[threadID] && E[currentLocalEIndex].genomeA == currentGenomeA
				&& E[currentLocalEIndex].genomeB == currentGenomeB) {
				currentLocalEIndex++;
			}

			if (currentLocalEIndex < EChunkStartIndex[threadID] + EChunkSize[threadID]) {
				// Found the end index of the current genome pair
				// Safely write this information into the genomePairEEndIndex array
				int genomePairIndexInJAC = genomePairToJACIndex(currentGenomeA, currentGenomeB);
				genomePairEEndIndex[genomePairIndexInJAC] = currentLocalEIndex - 1;

				// Move on to the next pair
				// This means we are seeing the start index of a new genome pair
				// Safely write this information into the genomePairEStartIndex array
				genomePairIndexInJAC = genomePairToJACIndex(E[currentLocalEIndex].genomeA, E[currentLocalEIndex].genomeB);
				genomePairEStartIndex[genomePairIndexInJAC] = currentLocalEIndex;
			}
			else {
				// currentLocalEIndex == EChunkStartIndex[threadID] + EChunkSize[threadID]
				// This means we've parsed our local E Chunk fully.
				// However, we still need to check if we have seen the end of the current genome pair. If yes, we must write this information
				// If our local E chunk doesn't store the end of the current genome pair, leave it to the next processor to fill this information

				if (currentLocalEIndex < E.size()) {
					ETriple firstElementOfNextThread = E[currentLocalEIndex];

					if (firstElementOfNextThread.genomeA != currentGenomeA || firstElementOfNextThread.genomeB != currentGenomeB) {
						// We have seen the end of the current genome pair
						// Safely write this information into genomePairEEndIndex array
						int genomePairIndexInJAC = genomePairToJACIndex(currentGenomeA, currentGenomeB);
						genomePairEEndIndex[genomePairIndexInJAC] = currentLocalEIndex - 1;
					} // Otherwise, this information will be filled by the next processor
				}
				else {
					// currentLocalEIndex == E.size()
					int genomePairIndexInJAC = genomePairToJACIndex(currentGenomeA, currentGenomeB);
					genomePairEEndIndex[genomePairIndexInJAC] = E.size() - 1;
				}
			}
		}

#pragma omp barrier

		/*#pragma omp single
		{
			std::cout << "Finished processing start and end index of pairs of genomes in E" << std::endl;
			std::cout << "First, E (Ga, Gb, ProteinID) is: ";
			for (int i = 0; i < E.size(); i++) {
				std::cout << "(" << E[i].genomeA << ", " << E[i].genomeB << ", " << E[i].proteinIndex << ")   ";
			}
			std::cout << std::endl << std::endl;


			std::cout << "Each pair of genome STARTS at:" << std::endl;
			for (int i = 0; i < genomePairEStartIndex.size(); i++) {
				std::cout << "Genome (" << JAC[i].genomeA << ", " << JAC[i].genomeB << "): " << genomePairEStartIndex[i] << std::endl;
			}
			std::cout << std::endl << std::endl;

			std::cout << "Each pair of genome ENDS at:" << std::endl;
			for (int i = 0; i < genomePairEEndIndex.size(); i++) {
				std::cout << "Genome (" << JAC[i].genomeA << ", " << JAC[i].genomeB << "): " << genomePairEEndIndex[i] << std::endl;
			}
			std::cout << std::endl << std::endl;
		}*/

		for (int genomePair = genomePairStartIndex; genomePair <= genomePairEndIndex; genomePair++) {
			int currGenomeA = JAC[genomePair].genomeA;
			int currGenomeB = JAC[genomePair].genomeB;

			// Block Bl of the same genome pair (Ga, Gb)
			// Note: Again, these are inclusive
			int blockBlStart = genomePairEStartIndex[genomePair];
			int blockBlEnd = genomePairEEndIndex[genomePair];

			double S = 0.0;
			int N = 0;

			// Subblocks Bk of the same protein Pi inside Bl
			int blockBkStart = blockBlStart;
			int blockBkEnd = blockBkStart;

			while (blockBkEnd <= blockBlEnd) {
				int currProteinID = E[blockBkStart].proteinIndex;

				while (blockBkEnd <= blockBlEnd && E[blockBkEnd].proteinIndex == currProteinID) {
					blockBkEnd++;
				}

				if (blockBkEnd <= blockBlEnd) {
					int BkLength = blockBkEnd - blockBkStart;
					double J_Pi_Ga_Gb = (double)(BkLength) / (double)(T[currProteinID][currGenomeA] + T[currProteinID][currGenomeB] - BkLength);
					S += J_Pi_Ga_Gb;
					N += 1;

					// Move on to the next Bk subblock
					blockBkStart = blockBkEnd;
				}
				else {
					// blockBkEnd > blockBlEnd
					// Finish the last computation
					int BkLength = blockBkEnd - blockBkStart;
					double J_Pi_Ga_Gb = (double)(BkLength) / (double)(T[currProteinID][currGenomeA] + T[currProteinID][currGenomeB] - BkLength);
					S += J_Pi_Ga_Gb;
					N += 1;
				}
			}

			// DOUBLE CHECK THIS PART
			JAC[genomePair].S += S;
			JAC[genomePair].N += N;
		}


		/** PHASE 4: Finalize output **/
#pragma omp single
		{
			std::cout << "Finished JAC construction" << std::endl;
			AJI.resize(totalGenomePairs);
			std::cout << "Prepped AJI" << std::endl;
		}

		for (int genomePair = genomePairStartIndex; genomePair <= genomePairEndIndex; genomePair++) {
			AJI[genomePair] = JAC[genomePair].S / JAC[genomePair].N;
		}
	}

	std::cout << "Computation Complete. Here is the final output." << std::endl;
	for (int i = 0; i < AJI.size(); i++) {
		int genomeA = JAC[i].genomeA;
		int genomeB = JAC[i].genomeB;

		std::cout << "Average Jaccard Index for (" << genomeA << ", " << genomeB << ") = " << AJI[i] << std::endl;
	}

	return SUCCESS;
}

int queryMetadataFromDatabase(sqlite3* db, int& GENOMECOUNT, std::vector <std::string>& proteinSet) {
	std::string sqlQuery = "SELECT count(*) as count_genome FROM genome_metadata";
	sqlite3_stmt* statement;
	int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

	if (errorCode != SQLITE_OK) {
		std::cerr << "Error in preparing sql statement " << sqlQuery << std::endl;
		std::cerr << "The error was: " << sqlite3_errmsg(db);
		return ERROR_SQLITE_DATABASE;
	}

	if (sqlite3_step(statement) == SQLITE_ROW) {
		GENOMECOUNT = sqlite3_column_int(statement, 0);
	}

	sqlQuery = "SELECT count(*) from scp_data";
	errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

	if (errorCode != SQLITE_OK) {
		std::cerr << "Error in preparing sql statement " << sqlQuery << std::endl;
		std::cerr << "The error was: " << sqlite3_errmsg(db);
		return ERROR_SQLITE_DATABASE;
	}

	if (sqlite3_step(statement) == SQLITE_ROW) {
		int scpCount = sqlite3_column_int(statement, 0);
		proteinSet.resize(scpCount);
	}

	sqlQuery = "SELECT scp_acc from scp_data";
	errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);
	if (errorCode != SQLITE_OK) {
		std::cerr << "Error in preparing sql statement " << sqlQuery << std::endl;
		std::cerr << "The error was: " << sqlite3_errmsg(db);
		return ERROR_SQLITE_DATABASE;
	}

	int index = 0;
	while (sqlite3_step(statement) == SQLITE_ROW) {
		const unsigned char* scp_acc = sqlite3_column_text(statement, 0);
		std::string proteinName(reinterpret_cast<const char*>(scp_acc));
		proteinSet[index] = proteinName;
		index++;
	}

	return SUCCESS;
}

int constructLcandLp(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<int>& Lc, std::vector<int>& Lp) {
	int totalNumThreads;
	int threadID;
	int tetramerStart = -1;
	int tetramerEnd = -1;

	Lc.resize(TETRAMER_COUNT, 0);
	Lp.resize(TETRAMER_COUNT, 0);
	int errorCodeFound = -1;

#pragma omp parallel default(none) \
	private(threadID, tetramerStart, tetramerEnd) \
	shared(totalNumThreads, Lc, proteinSet, db, errorCodeFound, std::cerr)
	{
		totalNumThreads = omp_get_num_threads();
		threadID = omp_get_thread_num();

		if (threadID < TETRAMER_COUNT % totalNumThreads) {
			// Each will be responsible for TETRAMER_COUNT / totalNumThreads + 1 tetramers
			tetramerStart = threadID * (TETRAMER_COUNT / totalNumThreads + 1);
			tetramerEnd = (threadID + 1) * (TETRAMER_COUNT / totalNumThreads + 1) - 1;
		}
		else {
			// Each will be responsible for TETRAMER_COUNT / totalNumThreads
			int buffer = (TETRAMER_COUNT % totalNumThreads) * (TETRAMER_COUNT / totalNumThreads + 1);
			tetramerStart = buffer + (threadID - TETRAMER_COUNT % totalNumThreads) * (TETRAMER_COUNT / totalNumThreads);
			tetramerEnd = buffer + (threadID - TETRAMER_COUNT % totalNumThreads + 1) * (TETRAMER_COUNT / totalNumThreads) - 1;
		}

		for (std::string protein : proteinSet) {
			std::string sqlQuery = "SELECT tetramer, genomes FROM `" + protein + "_tetras` WHERE tetramer BETWEEN ? AND ?";

			sqlite3_stmt* statement;
			int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

			sqlite3_bind_int(statement, 1, tetramerStart);
			sqlite3_bind_int(statement, 2, tetramerEnd);

			if (errorCode != SQLITE_OK) {
				errorCodeFound = ERROR_IN_CONSTRUCTION;
				std::cerr << "Error in preparing sql statement " << sqlQuery << std::endl;
				std::cerr << "The error was: " << sqlite3_errmsg(db);
			}

			if (errorCodeFound == -1) {
				while (sqlite3_step(statement) == SQLITE_ROW) {
					const int tetraID = sqlite3_column_int(statement, 0);
					size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
					int countGenomes = sizeOfBlobInBytes / sizeof(int);

					Lc[tetraID] += countGenomes;
				}
			}

			sqlite3_finalize(statement);
		}
	}

	if (errorCodeFound != -1) {
		return errorCodeFound;
	}

	// Parallel prefix sum on Lc to construct Lp
	int cumulativeSum = 0;
#pragma omp parallel for simd reduction(inscan,+: cumulativeSum)
	for (int i = 0; i < Lc.size(); i++) {
		Lp[i] = cumulativeSum;
#pragma omp scan exclusive(cumulativeSum)
		cumulativeSum += Lc[i];
	}

	return SUCCESS;
}

int constructF(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::pair<int, int>>& F, std::vector<int>& Lc, std::vector<int>& Lp) {
	int errorCodeFound = -1;
	std::vector<int> tetramerStartDistribution;
	std::vector<int> tetramerEndDistribution;
	int totalNumThreads = -1;
	float slack_percentage = SLACK_PERCENTAGE;

	std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> startTimes;
	std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> endTimes;

#pragma omp parallel default(none) \
	shared(F, db, proteinSet, Lc, Lp, tetramerStartDistribution, tetramerEndDistribution, totalNumThreads, errorCodeFound, slack_percentage, std::cerr, std::cout, startTimes, endTimes)
	{
		// Ask 1 thread to carry out the distribution of tasks
#pragma omp single
		{
			totalNumThreads = omp_get_num_threads();
			tetramerStartDistribution.resize(totalNumThreads, -1);
			tetramerEndDistribution.resize(totalNumThreads, -1);

			startTimes.resize(totalNumThreads);
			endTimes.resize(totalNumThreads);

			int totalNumberOfTasks = F.size();
			std::vector<int> taskSumPerThread(totalNumThreads, 0);
			int targetTasksPerProcessor = ((float)(totalNumberOfTasks) / totalNumThreads) * (1 + slack_percentage);

			int currentThread = 0;
			for (int tetramer = 0; tetramer < Lc.size(); tetramer++) {
				int task = Lc[tetramer];
				if (taskSumPerThread[currentThread] + task <= targetTasksPerProcessor || currentThread == totalNumThreads - 1) {
					taskSumPerThread[currentThread] += task;
					if (tetramerStartDistribution[currentThread] == -1) {
						tetramerStartDistribution[currentThread] = tetramer;
					}
					tetramerEndDistribution[currentThread] = tetramer;
				}
				else {
					// Move to the next processor and assign the task there
					currentThread++;
					taskSumPerThread[currentThread] += task;
					tetramerStartDistribution[currentThread] = tetramer;
					tetramerEndDistribution[currentThread] = tetramer;
				}
			}
		}

		int threadID = omp_get_thread_num();
		int tetramerStart = tetramerStartDistribution[threadID];
		int tetramerEnd = tetramerEndDistribution[threadID];

		auto startTime = std::chrono::high_resolution_clock::now();

		std::ostringstream oss;
		for (int proteinIndex = 0; proteinIndex < proteinSet.size(); proteinIndex++) {
			std::string protein = proteinSet[proteinIndex];
			oss << "SELECT tetramer, genomes, " << proteinIndex << " as source_table FROM `" + protein + "_tetras` WHERE tetramer BETWEEN " << tetramerStart << " AND " << tetramerEnd << " ";
			if (proteinIndex < proteinSet.size() - 1) {
				oss << " UNION ALL ";
			}
		}
		oss << " ORDER BY tetramer, source_table";

		std::string sqlQuery = oss.str();
		sqlite3_stmt* statement;
		int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

		if (errorCode != SQLITE_OK) {
			errorCodeFound = ERROR_IN_CONSTRUCTION;
			std::cerr << "Error in preparing sql statement " << std::endl;
			std::cerr << "The error was: " << sqlite3_errmsg(db);
		}

		// Starting location in F
		int indexInF = Lp[tetramerStart];

		if (errorCodeFound == -1) {
			while (sqlite3_step(statement) == SQLITE_ROW) {
				const int tetraID = sqlite3_column_int(statement, 0);
				const void* genomeBlob = sqlite3_column_blob(statement, 1);
				const int proteinIndex = sqlite3_column_int(statement, 2);

				const int* genomeArray = static_cast<const int*>(genomeBlob);
				size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);
				int countGenomes = sizeOfBlobInBytes / sizeof(int);

				for (int i = 0; i < countGenomes; i++) {
					int genomeID = genomeArray[i];

					F[indexInF] = std::make_pair(proteinIndex, genomeID);
					indexInF++;
				}
			}
		}
		sqlite3_finalize(statement);

		auto endTime = std::chrono::high_resolution_clock::now();
		startTimes[threadID] = startTime;
		endTimes[threadID] = endTime;
	}

	if (errorCodeFound != -1) {
		return errorCodeFound;
	}

	for (int threadID = 0; threadID < startTimes.size(); threadID++) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimes[threadID] - startTimes[threadID]);
		auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTimes[threadID] - startTimes[threadID]);
		auto minutes = duration_seconds.count() / 60;
		auto seconds = duration_seconds.count() % 60;
		std::cout << "ThreadID " << threadID << " took " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds) to construct local F" << std::endl;
	}

	return SUCCESS;
}

int constructT(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::vector<int>>& T, const int GENOMECOUNT) {
	int totalNumThreads;
	int proteinStart = -1;
	int proteinEnd = -1;
	int threadID = -1;
	int proteinCount = proteinSet.size();

	int errorCodeFound = -1;
	T.resize(proteinCount);
	for (auto& row : T) {
		row.resize(GENOMECOUNT, 0);
	}

	std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> startTimes;
	std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> endTimes;

#pragma omp parallel default(none) \
	shared(db, proteinSet, proteinCount, T, errorCodeFound, std::cerr, std::cout, startTimes, endTimes) \
	private(totalNumThreads, threadID, proteinStart, proteinEnd)
	{
		totalNumThreads = omp_get_num_threads();
		threadID = omp_get_thread_num();

#pragma omp single
		{
			startTimes.resize(totalNumThreads);
			endTimes.resize(totalNumThreads);
		}

		auto startTime = std::chrono::high_resolution_clock::now();
		if (threadID < proteinCount % totalNumThreads) {
			// Each will be responsible for proteinCount / totalNumThreads + 1 proteins
			proteinStart = threadID * (proteinCount / totalNumThreads + 1);
			proteinEnd = (threadID + 1) * (proteinCount / totalNumThreads + 1) - 1;
		}
		else {
			// Each will be responsible for proteinCount / totalNumThreads proteins
			int buffer = (proteinCount % totalNumThreads) * (proteinCount / totalNumThreads + 1);
			proteinStart = buffer + (threadID - proteinCount % totalNumThreads) * (proteinCount / totalNumThreads);
			proteinEnd = buffer + (threadID - proteinCount % totalNumThreads + 1) * (proteinCount / totalNumThreads) - 1;
		}

		std::ostringstream oss;
		for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd; proteinIndex++) {
			std::string protein = proteinSet[proteinIndex];
			oss << "SELECT genome_id, length(tetramers), " << proteinIndex << " as source_table from `" << protein << "_genomes` ";
			if (proteinIndex < proteinEnd) {
				oss << " UNION ALL ";
			}
		}
		std::string sqlQuery = oss.str();

		sqlite3_stmt* statement;
		int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

		if (errorCode != SQLITE_OK) {
			errorCodeFound = ERROR_IN_CONSTRUCTION;
			std::cerr << "Error in preparing sql statement " << sqlQuery << std::endl;
			std::cerr << "The error was: " << sqlite3_errmsg(db);
		}

		if (errorCodeFound == -1) {
			while (sqlite3_step(statement) == SQLITE_ROW) {
				const int genomeID = sqlite3_column_int(statement, 0);
				int sizeOfBlobInBytes = sqlite3_column_int(statement, 1);
				int countTetras = sizeOfBlobInBytes / sizeof(int);
				int proteinIndex = sqlite3_column_int(statement, 2);

				//// This code will read the tetramer values in the blob
				//const void* tetraBlob = sqlite3_column_blob(statement, 1);

				//auto* blobBytes = static_cast<const unsigned char*>(tetraBlob);
				//std::cout << "For genome " << genomeID << " the tetras blob include: ";
				//for (size_t i = 0; i < sizeOfBlobInBytes; i += sizeof(int)) {
				//	// From Kenji's python code, inferred that his sqlite database is storing in little endian and integer is 4 bytes
				//	int value = blobBytes[i] | (blobBytes[i + 1] << 8) | (blobBytes[i + 2] << 16) | (blobBytes[i + 3] << 24);
				//	std::cout << value << " ";
					//}
				//std::cout << std::endl;

				T[proteinIndex][genomeID] = countTetras;
			}
		}

		sqlite3_finalize(statement);
		auto endTime = std::chrono::high_resolution_clock::now();

		startTimes[threadID] = startTime;
		endTimes[threadID] = endTime;
	}


	if (errorCodeFound != -1) {
		return errorCodeFound;
	}

	for (int threadID = 0; threadID < startTimes.size(); threadID++) {
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTimes[threadID] - startTimes[threadID]);
		auto duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(endTimes[threadID] - startTimes[threadID]);
		auto minutes = duration_seconds.count() / 60;
		auto seconds = duration_seconds.count() % 60;
		std::cout << "ThreadID " << threadID << " took " << duration.count() << " milliseconds (i.e. " << minutes << " minutes " << seconds << " seconds) to construct local T" << std::endl;
	}

	return SUCCESS;
}

bool customSortE(const ETriple& firstElement, const ETriple& secondElement) {
	if (firstElement.genomeA != secondElement.genomeA) {
		return firstElement.genomeA < secondElement.genomeA;
	}

	if (firstElement.genomeB != secondElement.genomeB) {
		return firstElement.genomeB < secondElement.genomeB;
	}

	return firstElement.proteinIndex < secondElement.proteinIndex;
}

void parallelMergeSort(std::vector<ETriple>& E, int left, int right, int serialThreshold) {
	if (left < right) {
		if (right - left < serialThreshold) {
			std::sort(E.begin() + left, E.begin() + right + 1, customSortE);
		}
		else {
			int mid = left + (right - left) / 2;

#pragma omp task shared(E)
			parallelMergeSort(E, left, mid, serialThreshold);

#pragma omp task shared(E)
			parallelMergeSort(E, mid + 1, right, serialThreshold);

			// We must wait for both tasks to complete before we merge the sorted halves
#pragma omp taskwait

			std::vector<ETriple> temp;
			std::merge(E.begin() + left, E.begin() + mid + 1, E.begin() + mid + 1, E.begin() + right + 1, std::back_inserter(temp), customSortE);
			std::copy(temp.begin(), temp.end(), E.begin() + left);
		}
	}
}

int genomePairToJACIndex(int genomeA, int genomeB) {
	// (37/2) * a - a^2 / 2 + b - 1
	return (37 * genomeA - genomeA * genomeA) / 2 + genomeB - 1;
}



int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cout << "Missing Database location " << std::endl;
		std::cout << "Usage: " << argv[0] << " Database Location " << std::endl;
		return 1;
	}
	// const std::string pathToDatabse = "modified_xantho_fastaai2.db";
	const std::string pathToDatabase = argv[1];
	return parallelfastaai(pathToDatabase);
}
