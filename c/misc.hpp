#ifndef MISC_H
#define MISC_H

#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <cstdint>
#include <cassert>
#include <vector>
#include <thread>

uint32_t log2 (uint32_t value);

#ifdef _OPENMP

/**
 * This object is used to temporarily change the max number of omp threads.
 * When the object is destructed, the max threads is set to it's original value.
 */
class ThreadLimit {
public:
    ThreadLimit(uint32_t maxThreads):
        prev_max_threads(omp_get_max_threads())
    {
        omp_set_num_threads(maxThreads);
    }

    ~ThreadLimit() noexcept
    {
        omp_set_num_threads(prev_max_threads);
    }

private:
    uint32_t prev_max_threads;
};

#endif // _OPENMP

inline unsigned int threadCount() {
    unsigned int n = std::thread::hardware_concurrency();

    return n == 0 ? 1 : n;
}

inline std::vector<int> divideWork(int elementCount, int threadCount) {
    assert(elementCount > 0);
    assert(threadCount > 0);

    if (elementCount <= threadCount) {
        return std::vector<int>(elementCount, 1);
    }

    const int jobSize = elementCount / threadCount;

    std::vector<int> jobs(threadCount, jobSize);

    const int elementRest = elementCount % threadCount;

    for (int i = 0; i < elementRest; i++) {
        jobs[i] += 1;
    }

    return jobs;
}

template<typename Func>
void parallelFor(int begin, int end, int threadCount, Func&& func) {

    const int elementCount = end - begin;

    std::vector<std::thread> threads;
    auto jobs = divideWork(elementCount, threadCount);
    int  jobBegin = begin;
    int  k = 0;

    for (; k < jobs.size() - 1; jobBegin += jobs[k++]) {

        threads.push_back(std::thread(func, jobBegin, jobBegin + jobs[k], k));
    }

    func(jobBegin, jobBegin + jobs[k], k);

    for (k = 0; k < jobs.size() - 1; k++) {
        threads[k].join();
    }
}

template<typename Func>
void parallelBlock(int threadCount, Func&& func) {

    std::vector<std::thread> threads;
    int  k = 0;

    for (; k < threadCount - 1; k++) {
        threads.push_back(std::thread(func, k, threadCount));
    }

    func(k, threadCount);

    for (k = 0; k < threadCount - 1; k++) {
        threads[k].join();
    }
}


#endif // MISC_H
