///
// @file utils.hpp
// @brief Utility functions for timers/memory tracking.
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

#ifndef UTILS_HPP
#define UTILS_HPP
#include <chrono>  // NOLINT
#include <iomanip>
#include <iostream>
#include <limits>
#include <ratio>  // NOLINT
#include <string>
#include <vector>
#include <cassert>
#include <fmt/format.h>

/// macros for block decomposition
#define BLOCK_LOW(i, p, n) ((i * n) / p)
#define BLOCK_HIGH(i, p, n) ((((i + 1) * n) / p) - 1)
#define BLOCK_SIZE(i, p, n) (BLOCK_LOW((i + 1), p, n) - BLOCK_LOW(i, p, n))
#define BLOCK_OWNER(j, p, n) (((p) * ((j) + 1) - 1) / (n))

template <typename SizeType, typename T>
static inline SizeType block_low(const T& rank, const T& nproc,
                                 const SizeType& n) {
    return (rank * n) / nproc;
}

template <typename SizeType, typename T>
static inline SizeType block_high(const T& rank, const T& nproc,
                                  const SizeType& n) {
    return (((rank + 1) * n) / nproc) - 1;
}

template <typename SizeType, typename T>
static inline SizeType block_size(const T& rank, const T& nproc,
                                  const SizeType& n) {
    return block_low<SizeType, T>(rank + 1, nproc, n) -
           block_low<SizeType, T>(rank, nproc, n);
}

template <typename SizeType, typename T>
static inline T block_owner(const SizeType& j, const SizeType& n,
                            const T& nproc) {
    return (((nproc) * ((j) + 1) - 1) / (n));
}

//
template <typename IT>
inline std::vector<IT>
distribute_bags_of_tasks(int nproc, IT ntasks, const std::vector<IT>& bag_sizes,
                         float slack,
                         std::vector<IT>& bag_starts,  // NOLINT
                         std::vector<IT>& bag_ends) {  // NOLINT
    std::vector<IT> ntasks_dist(nproc, 0);
    int ntasks_per_proc = (static_cast<float>(ntasks) / nproc) * (1 + slack);
    bag_starts.resize(nproc, -1);
    bag_ends.resize(nproc, -1);
    IT n_bags = IT(bag_sizes.size());

    for (IT bag_id = 0, pid = 0; bag_id < n_bags; bag_id++) {
        IT ntasks_bag = bag_sizes[bag_id];
        if (ntasks_dist[pid] + ntasks_bag <= ntasks_per_proc ||
            pid == nproc - 1) {
            ntasks_dist[pid] += ntasks_bag;
            if (bag_starts[pid] == -1) {
                bag_starts[pid] = bag_id;
            }
            bag_ends[pid] = bag_id;
        } else {
            // Move to the next thread and assign the task there
            pid++;
            ntasks_dist[pid] += ntasks_bag;
            bag_starts[pid] = bag_id;
            bag_ends[pid] = bag_id;
        }
    }
    return ntasks_dist;
}

// timer definition
//

template <typename duration> class timer_impl {
  private:
    std::chrono::steady_clock::time_point start;
    typename duration::rep _total_elapsed;
    typename duration::rep _elapsed_time;

  public:
    const typename duration::rep& total_elapsed = _total_elapsed;

    timer_impl() : start(std::chrono::steady_clock::now()), _total_elapsed(0) {}

    inline timer_impl<duration>& accumulate() {
        _total_elapsed += _elapsed_time;
        return *this;
    }

    inline timer_impl<duration>& reset() {
        start = std::chrono::steady_clock::now();
        return *this;
    }

    inline timer_impl<duration>& elapsed() {
        std::chrono::steady_clock::time_point stop =
            std::chrono::steady_clock::now();
        _elapsed_time = duration(stop - start).count();
        return *this;
    }

    inline typename duration::rep elapsed_to_mins() const {
        return elapsed_to_seconds() / 60.0;
    }

    inline typename duration::rep elapsed_to_seconds() const {
        return _elapsed_time / 1000.0;
    }

    void print(const std::string prefix, std::ostream& ox,
               bool line_end = true) const {
        ox << prefix;
        ox.precision(3);
        ox << std::setw(10) << _elapsed_time << " ms (";
        ox.precision(4);
        ox << std::setw(10) << elapsed_to_seconds() << " sec/";
        ox.precision(4);
        ox << std::setw(10) << elapsed_to_mins() << " min).";
        if (line_end)
            ox << std::endl;
    }

    void measure_accumulate_print(const std::string prefix, std::ostream& ox,
                                  bool line_end = true) {
        elapsed().accumulate().print(prefix, ox, line_end);
    }
};

using timer = timer_impl<std::chrono::duration<double, std::milli>>;

template <typename TimerT>
void printThreadTimes(const std::string& prt_prefix,
                      const std::vector<TimerT>& threadTimers) {
    for (int threadID = 0; threadID < threadTimers.size(); threadID++) {
        std::cout << "Thread " << threadID;
        threadTimers[threadID].print_elapsed(prt_prefix, std::cout);
    }
}

//
// Compute Memory Usage
//
#include <sys/resource.h>
template <typename T> T get_mem_usage() {
    struct rusage usage;
    int ret;
    ret = getrusage(RUSAGE_SELF, &usage);
    if (ret != 0) {
        return std::numeric_limits<T>::max();
    }
    return T(usage.ru_maxrss);  // in KB
}

template <typename T>
void print_memory_usage(const std::string& prt_prefix, std::ostream& ox) {
    T used_mem = get_mem_usage<T>();
    ox << prt_prefix;
    if (used_mem < std::numeric_limits<T>::max()) {
        ox.precision(2);
        ox << std::setw(10) << used_mem << " KB (";
        ox.precision(2);
        ox << std::setw(10) << used_mem / 1024 << " MB/";
        ox.precision(2);
        ox << std::setw(10) << used_mem / (1024 * 1024) << " GB).";
    } else {
        ox << " Mem usage returned error ";
    }
    ox << std::endl;
}
#define PRINT_RUNTIME_MEMUSED(r_timer, sprefix, ostream)                       \
    {                                                                          \
        r_timer.measure_accumulate_print(sprefix, ostream, false);             \
        print_memory_usage<uint64_t>("; ", std::cout);                         \
    }

//
// Utility Class for Pair with first and second element
template <typename DT1, typename DT2> struct DPair {
    DT1 first;
    DT2 second;
    explicit DPair(DT1 a, DT2 b) : first(a), second(b) {}
    DPair() : first(DT1(-1)), second(DT2(-1)) {}
    DPair(std::initializer_list<DT1> l)
        : first(*l.begin()), second(*(l.begin() + 1)) {}
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


#endif  // !UTILS_HPP
