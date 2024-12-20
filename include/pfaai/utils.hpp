//
// Copyright [] <>
// TODO(srirampc)
//

#ifndef UTILS_HPP
#define UTILS_HPP
#include <chrono>  // NOLINT
#include <iomanip>
#include <iostream>
#include <limits>
#include <ratio>  // NOLINT
#include <string>
#include <vector>

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
inline void distribute_bags_of_tasks(int nproc, IT ntasks,
                                     const std::vector<IT>& bag_sizes,
                                     float slack,
                                     std::vector<IT>& proc_bag_start,  // NOLINT
                                     std::vector<IT>& proc_bag_end) {  // NOLINT
    std::vector<IT> ntasks_distr(nproc, 0);
    int ntasks_per_proc = (static_cast<float>(ntasks) / nproc) * (1 + slack);
    proc_bag_start.resize(nproc, -1);
    proc_bag_end.resize(nproc, -1);
    IT n_bags = IT(bag_sizes.size());

    for (IT bag_id = 0, pid = 0; bag_id < n_bags; bag_id++) {
        IT ntasks_bag = bag_sizes[bag_id];
        if (ntasks_distr[pid] + ntasks_bag <= ntasks_per_proc ||
            pid == nproc - 1) {
            ntasks_distr[pid] += ntasks_bag;
            if (proc_bag_start[pid] == -1) {
                proc_bag_start[pid] = bag_id;
            }
            proc_bag_end[pid] = bag_id;
        } else {
            // Move to the next thread and assign the task there
            pid++;
            ntasks_distr[pid] += ntasks_bag;
            proc_bag_start[pid] = bag_id;
            proc_bag_end[pid] = bag_id;
        }
    }
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

#endif  // !UTILS_HPP
