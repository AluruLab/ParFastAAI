/* NOTE: This solution is written for a sqlite database with UTF-8 Encoding */
#include <algorithm>
#include <chrono>
#include <execution>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <set>
#include "sqlite-amalgamation-3450100/sqlite3.h"
#include <string>
#include <unordered_map>
#include <vector>

#define SUCCESS 0
#define ERROR_SQLITE_DATABASE 1
#define ERROR_SQLITE_MEMORY_ALLOCATION 2
#define ERROR_IN_CONSTRUCTION 3
#define TETRAMER_COUNT 160000
#define GENOME_COUNT 20
#define SLACK_PERCENTAGE 0.0

struct ETriple {
	int proteinIndex;
	int genomeA;
	int genomeB;

	ETriple() : proteinIndex(-1), genomeA(-1), genomeB(-1) {}
	ETriple(const int proteinIndexVal, const int genomeAVal, const int genomeBVal) : proteinIndex(proteinIndexVal), genomeA(genomeAVal), genomeB(genomeBVal) {}
};

int constructLcandLp(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<int>& Lc, std::vector<int>& Lp);
int constructF(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::pair<int, int>>& F,
	std::vector<int>& Lc, std::vector<int>& Lp);
int constructT(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::vector<int>>& T);
bool customSortE(const ETriple& firstElement, const ETriple& secondElement);
void parallelMergeSort(std::vector<ETriple>& E, int left, int right, int serialThreshold);

int main()
{
    sqlite3* database;
    const std::string pathToDatabse = "modified_xantho_fastaai2.db";
    int errorCode = sqlite3_open(pathToDatabse.c_str(), &database);

    if (errorCode != SQLITE_OK) {
        std::cerr << "Error in opening " << pathToDatabse << std::endl;
        std::cerr << "The error was: " << sqlite3_errstr(errorCode) << std::endl;
        sqlite3_close(database);
        return ERROR_SQLITE_DATABASE;
    }
    else if (database == nullptr) {
        std::cerr << "SQLite is unable to allocate memory for the database " << pathToDatabse << std::endl;
        sqlite3_close(database);
        return ERROR_SQLITE_MEMORY_ALLOCATION;
    }

	/** PHASE 1: Construction of the data structures **/

	//std::vector<std::string> proteinSet = { "pf00411.19", "pf00237.19", "pf01016.19", "pf02033.18", "pf00347.23", "pf00119.20",
	//			"pf00297.22", "pf02601.15", "pf00318.20", "pf02367.17", "pf00825.18", "pf02410.15",
	//			"pf00406.22", "pf00380.19", "pf00213.18", "pf05221.17", "pf00252.18", "pf00177.21",
	//			"pf00709.21", "pf00312.22", "pf01192.22", "pf06026.14", "pf01649.18", "pf00572.18",
	//			"pf00338.22", "pf01142.18", "pf01746.21", "pf01632.19", "pf17136.4", "pf00164.25",
	//			"pf01808.18", "pf00750.19", "pf00889.19", "pf01196.19", "pf01250.17", "pf00162.19",
	//			"pf01725.16", "pf01668.18", "pf00749.21", "pf00238.19", "pf02565.15", "pf01715.17",
	//			"pf00203.21", "pf00828.19", "pf00573.22", "pf00121.18", "pf13393.6", "pf00276.20",
	//			"pf01139.17", "pf00886.19", "pf00189.20", "pf01176.19", "pf00687.21", "pf01025.19",
	//			"pf03948.14", "pf00829.21", "pf01195.19", "pf03840.14", "pf00231.19", "pf01351.18",
	//			"pf01264.21", "pf01193.24", "pf00410.19", "pf00584.20", "pf01245.20", "pf02130.17",
	//			"pf02699.15", "pf01765.19", "pf01783.23", "pf00281.19", "pf00416.22", "pf00366.20",
	//			"pf00344.20", "pf00831.23", "pf00334.19", "pf00830.19", "pf00861.22", "pf00453.18",
	//			"pf00181.23", "pf03652.15" };

	std::vector<std::string> proteinSet = { "pf00121.18" }; // { "pf00347.23", "pf00411.19" };

    std::vector<int> Lc;
    std::vector<int> Lp;
    errorCode = constructLcandLp(database, proteinSet, Lc, Lp);

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing Lc and Lp, error code = " << errorCode << std::endl;
		return errorCode;
	}
	std::cout << "Done with constructing Lc and Lp" << std::endl;

	std::vector<std::pair<int, int>> F(Lp[TETRAMER_COUNT - 1] + Lc[TETRAMER_COUNT - 1], std::make_pair(-1, -1));
	errorCode = constructF(database, proteinSet, F, Lc, Lp);

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing F, error code = " << errorCode << std::endl;
		return errorCode;
	}
	std::cout << "Done with constructing F" << std::endl;

	std::vector<std::vector<int>> T;
	errorCode = constructT(database, proteinSet, T);

	if (errorCode != SUCCESS) {
		std::cerr << "Error in constructing T, error code = " << errorCode << std::endl;
		return errorCode;
	}

	std::cout << "Done with constructing T" << std::endl;

    sqlite3_close(database);

	std::vector<int> tetramerStartDistribution;
	std::vector<int> tetramerEndDistribution;
	int totalNumThreads = -1;
	float slack_percentage = SLACK_PERCENTAGE;

	std::vector<ETriple> E;
	std::vector<int> EChunkSize;
	std::vector<int> EChunkStartIndex;
	int ESize = 0;

	#pragma omp parallel default(none) \
	shared(Lc, Lp, F, T, E, EChunkSize, EChunkStartIndex, tetramerStartDistribution, tetramerEndDistribution, slack_percentage, totalNumThreads, ESize, std::cout)
	{
		/** PHASE 2: Generate tetramer tuples **/

		// Ask 1 thread to carry out the distribution of tasks, i.e. tetramer tuples for all threads
		#pragma omp single
		{
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
			} else {
				// tetramerID == tetramerEnd
				if (tetramerID == TETRAMER_COUNT - 1) {
					endIndexInF = F.size() - 1;
				} else {
					endIndexInF = Lp[tetramerID + 1] - 1;
				}
			}

			int currProteinID = F[startIndexInF].first;
			int leftBoundary = startIndexInF;
			int rightBoundary = startIndexInF;

			while (rightBoundary <= endIndexInF) {
				if (F[rightBoundary].first == currProteinID) {
					rightBoundary++;
				} else {
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
			} else {
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
			std::cout << "Done with constructing E" << std::endl;
			// Sort E
			parallelMergeSort(E, 0, E.size() - 1, 5);
		}
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
		} else {
			// Each will be responsible for TETRAMER_COUNT / totalNumThreads
			int buffer = (TETRAMER_COUNT % totalNumThreads) * (TETRAMER_COUNT / totalNumThreads + 1);
			tetramerStart = buffer + (threadID - TETRAMER_COUNT % totalNumThreads) * (TETRAMER_COUNT / totalNumThreads);
			tetramerEnd = buffer + (threadID - TETRAMER_COUNT % totalNumThreads + 1) * (TETRAMER_COUNT / totalNumThreads) - 1;
		}

		for (std::string protein : proteinSet) {
			std::string sqlQuery = "SELECT tetra, genomes FROM `" + protein + "_tetras` WHERE tetra BETWEEN ? AND ?";

			sqlite3_stmt* statement;
			int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

			sqlite3_bind_int(statement, 1, tetramerStart);
			sqlite3_bind_int(statement, 2, tetramerEnd);

			if (errorCode != SQLITE_OK) {
				errorCodeFound = ERROR_IN_CONSTRUCTION;
				std::cerr << "Error in preparing sql statement" << sqlQuery << std::endl;
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

	#pragma omp parallel default(none) \
	shared(F, db, proteinSet, Lc, Lp, tetramerStartDistribution, tetramerEndDistribution, totalNumThreads, errorCodeFound, slack_percentage, std::cerr)
	{
		// Ask 1 thread to carry out the distribution of tasks
		#pragma omp single
		{
			totalNumThreads = omp_get_num_threads();
			tetramerStartDistribution.resize(totalNumThreads, -1);
			tetramerEndDistribution.resize(totalNumThreads, -1);

			int totalNumberOfTasks = F.size();
			std::vector<int> taskSumPerThread(totalNumThreads, 0);
			int targetTasksPerProcessor = ((float) (totalNumberOfTasks) / totalNumThreads) * (1 + slack_percentage);

			int currentThread = 0;
			for (int tetramer = 0; tetramer < Lc.size(); tetramer++) {
				int task = Lc[tetramer];
				if (taskSumPerThread[currentThread] + task <= targetTasksPerProcessor || currentThread == totalNumThreads - 1) {
					taskSumPerThread[currentThread] += task;
					if (tetramerStartDistribution[currentThread] == -1) {
						tetramerStartDistribution[currentThread] = tetramer;
					}
					tetramerEndDistribution[currentThread] = tetramer;
				} else {
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

		for (int tetramerID = tetramerStart; tetramerID <= tetramerEnd; tetramerID++) {
			// Starting location in F
			int indexInF = Lp[tetramerID];

			for (int proteinIndex = 0; proteinIndex < proteinSet.size(); proteinIndex++) {
				std::string protein = proteinSet[proteinIndex];

				// Get the genome blob
				std::string sqlQuery = "SELECT genomes FROM `" + protein + "_tetras` WHERE tetra=?";

				sqlite3_stmt* statement;
				int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

				sqlite3_bind_int(statement, 1, tetramerID);

				if (errorCode != SQLITE_OK) {
					errorCodeFound = ERROR_IN_CONSTRUCTION;
					std::cerr << "Error in preparing sql statement" << sqlQuery << std::endl;
					std::cerr << "The error was: " << sqlite3_errmsg(db);
				}

				if (errorCodeFound == -1) {
					while (sqlite3_step(statement) == SQLITE_ROW) {
						const void* genomeBlob = sqlite3_column_blob(statement, 0);
						const int* genomeArray = static_cast<const int*>(genomeBlob);
						int countGenomes = Lc[tetramerID];

						for (int i = 0; i < countGenomes; i++) {
							int genomeID = genomeArray[i];
							F[indexInF] = std::make_pair(proteinIndex, genomeID);
							indexInF++;
						}
					}
				}
				sqlite3_finalize(statement);
			}
		}
	}

	if (errorCodeFound != -1) {
		return errorCodeFound;
	}

	return SUCCESS;
}

int constructT(sqlite3* db, std::vector<std::string>& proteinSet, std::vector<std::vector<int>>& T) {
	int totalNumThreads;
	int proteinStart = -1;
	int proteinEnd = -1;
	int threadID = -1;
	int proteinCount = proteinSet.size();

	int errorCodeFound = -1;

	T.resize(proteinCount);
	for (auto &row: T) {
		row.resize(GENOME_COUNT, 0);
	}

	#pragma omp parallel default(none) \
	shared(db, proteinSet, proteinCount, T, errorCodeFound, std::cerr) \
	private(totalNumThreads, threadID, proteinStart, proteinEnd)
	{
		totalNumThreads = omp_get_num_threads();
		threadID = omp_get_thread_num();


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

		for (int proteinIndex = proteinStart; proteinIndex <= proteinEnd; proteinIndex++) {
			std::string protein = proteinSet[proteinIndex];

			std::string sqlQuery = "SELECT genome, tetras from `" + protein + "_genomes`";

			sqlite3_stmt* statement;
			int errorCode = sqlite3_prepare_v2(db, sqlQuery.c_str(), -1, &statement, nullptr);

			if (errorCode != SQLITE_OK) {
				errorCodeFound = ERROR_IN_CONSTRUCTION;
				std::cerr << "Error in preparing sql statement" << sqlQuery << std::endl;
				std::cerr << "The error was: " << sqlite3_errmsg(db);
			}

			if (errorCodeFound == -1) {
				while (sqlite3_step(statement) == SQLITE_ROW) {
					const int genomeID = sqlite3_column_int(statement, 0);
					size_t sizeOfBlobInBytes = sqlite3_column_bytes(statement, 1);

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

					int countTetras = sizeOfBlobInBytes / sizeof(int);
					T[proteinIndex][genomeID] = countTetras;
				}
			}

			sqlite3_finalize(statement);
		}
	}

	if (errorCodeFound != -1) {
		return errorCodeFound;
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
		} else {
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
