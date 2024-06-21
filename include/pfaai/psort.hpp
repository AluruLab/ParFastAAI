#ifndef PAR_SORT_HPP
#define PAR_SORT_HPP

#include <algorithm>
#include <vector>

template <typename DataType>
inline bool customSortE(const DataType& firstElement,
                        const DataType& secondElement) {
    // if (firstElement.genomeA != secondElement.genomeA) {
    //	return firstElement.genomeA < secondElement.genomeA;
    // }

    // if (firstElement.genomeB != secondElement.genomeB) {
    //	return firstElement.genomeB < secondElement.genomeB;
    // }

    // return firstElement.proteinIndex < secondElement.proteinIndex;
    return firstElement < secondElement;
}

template <typename DataType>
void parallelMergeSort(std::vector<DataType>& E, int left, int right,
                       int serialThreshold) {
    if (left < right) {
        if (right - left < serialThreshold) {
            std::sort(E.begin() + left, E.begin() + right + 1,
                      customSortE<DataType>);
        } else {
            int mid = left + (right - left) / 2;

#pragma omp task shared(E)
            parallelMergeSort(E, left, mid, serialThreshold);

#pragma omp task shared(E)
            parallelMergeSort(E, mid + 1, right, serialThreshold);

            // We must wait for both tasks to complete before we merge the
            // sorted halves
#pragma omp taskwait

            std::vector<DataType> temp;
            std::merge(E.begin() + left, E.begin() + mid + 1,
                       E.begin() + mid + 1, E.begin() + right + 1,
                       std::back_inserter(temp), customSortE<DataType>);
            std::copy(temp.begin(), temp.end(), E.begin() + left);
        }
    }
}

#endif  // !PAR_SORT_HPP
