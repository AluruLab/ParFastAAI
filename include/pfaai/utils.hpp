//
// Copyright [] <>
// TODO(srirampc)
//

#ifndef UTILS_HPP
#define UTILS_HPP
#include <iomanip>
#include <string>
#include <iostream>
#include <chrono>   // NOLINT
#include <ratio>    // NOLINT

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

    inline void accumulate() { _total_elapsed += elapsed(); }

    inline void reset() { start = std::chrono::steady_clock::now(); }

    typename duration::rep elapsed() {
        std::chrono::steady_clock::time_point stop =
            std::chrono::steady_clock::now();
        _elapsed_time = duration(stop - start).count();
        return _elapsed_time;
    }

    inline typename duration::rep elapsed_to_mins() const {
      return elapsed_to_seconds()/60.0;
    }

    inline typename duration::rep elapsed_to_seconds() const {
      return _elapsed_time/1000.0;
    }

    void print_elapsed(const std::string prefix, std::ostream& ox) const {
	    ox << prefix;
        ox.precision(3);
        ox << std::setw(10) << _elapsed_time << " ms (" ;
        ox.precision(4);
        ox << std::setw(10) << elapsed_to_seconds() << " sec/" ;
        ox.precision(4);
        ox << std::setw(10) << elapsed_to_mins() << " min)." ;
        ox << std::endl;
    }
};

using timer = timer_impl<std::chrono::duration<double, std::milli>>;

#endif  // !UTILS_HPP
