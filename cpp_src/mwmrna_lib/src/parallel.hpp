#ifndef MWMRNA_LIB_PARALLEL_HPP
#define MWMRNA_LIB_PARALLEL_HPP

#include <functional>
#include <thread>
#include <vector>

namespace mwmrna::parallel {

/// @brief Simple parallel for each loop
void for_each(int st, int en, std::function<void(int)> f, int num_threads = -1) {
    if (num_threads == -1) {
        num_threads = std::thread::hardware_concurrency();
    }
    const int N = en - st;
    if (N < num_threads) {
        for (int i = st; i < en; ++i) {
            f(i);
        }
        return;
    }
    int chunk_size = N / num_threads;
    if (N % num_threads != 0) {
        ++chunk_size;
    }
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([i, st, N, chunk_size, &f]() {
            const int start = i * chunk_size;
            const int end = std::min((i + 1) * chunk_size, N);
            for (int j = start; j < end; ++j) {
                f(st + j);
            }
        });
    }
    for (auto &t : threads) {
        t.join();
    }
}
}  // namespace mwmrna::parallel

#endif  // MWMRNA_LIB_PARALLEL_HPP