///
// @file psort.hpp
// @brief Implementation of a shared memory parallel merge sort.
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

#ifndef PAR_SORT_HPP
#define PAR_SORT_HPP

#include <algorithm>
#include <vector>

template <typename DataType>
void parallelMergeSort(std::vector<DataType>& E, int left, int right,  // NOLINT
                       int serialThreshold) {
    if (left < right) {
        if (right - left < serialThreshold) {
            std::sort(E.begin() + left, E.begin() + right + 1);
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
                       std::back_inserter(temp));
            std::copy(temp.begin(), temp.end(), E.begin() + left);
        }
    }
}

#endif  // !PAR_SORT_HPP
