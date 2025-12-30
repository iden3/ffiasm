#ifndef MISC_H
#define MISC_H

#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <cstdint>
#include <cassert>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <queue>

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

using ThreadJob = std::function<void(uint64_t)>;

class ThreadJobQueue
{
    std::queue<ThreadJob> queue;
    mutable std::mutex mutex;

public:
    template<typename Job>
    void unsafePush(Job&& value) {
        queue.push(std::forward<Job>(value));
    }

    bool tryPop(ThreadJob& value) {
        std::lock_guard<std::mutex> lock(mutex);

        if(queue.empty())
            return false;

        value = queue.front();
        queue.pop();

        return true;
    }
};

class ThreadWorker
{
    std::atomic<bool> stop;
    ThreadJob job;
    std::shared_ptr<ThreadJobQueue> queue;
    std::mutex mutex;
    std::condition_variable occupied;
    std::condition_variable finished;
    std::thread thread;
    uint64_t threadId = 0;

    void worker() {

        while (true) {
            {
                std::unique_lock<std::mutex> lock(mutex);

                occupied.wait(lock, [this]{return job || stop;});

                if (stop) {
                    return;
                }

                job(threadId);

                if (queue) {
                    ThreadJob queuedJob;

                    while (queue->tryPop(queuedJob)) {
                        queuedJob(threadId);
                    }
                }
                job = nullptr;
            }
            finished.notify_one();
        }
    }

public:
    ThreadWorker() :
        stop(false),
        thread(&ThreadWorker::worker, this)
    {}

    ~ThreadWorker() {
        stop = true;
        occupied.notify_one();

        thread.join();
    }

    void setQueue(std::shared_ptr<ThreadJobQueue>& newQueue) {
        queue = newQueue;
    }

    void setThreadId(uint64_t id) {
        threadId = id;
    }

    template<typename Job>
    void submit(Job&& newJob) {
        {
            std::lock_guard<std::mutex> lock(mutex);
            job = std::forward<Job>(newJob);
        }
        occupied.notify_one();
    }

    void wait() {
        std::unique_lock<std::mutex> lock(mutex);

        finished.wait(lock, [this]{return !job;});
    }
};

class ThreadPool {

    struct JobRange {
        int64_t begin;
        int64_t end;
    };

    static const uint64_t chunksPerThread = 4;

    uint64_t nThreads;
    std::vector<ThreadWorker> workers;
    std::shared_ptr<ThreadJobQueue> queue;
    std::mutex poolMutex;

public:
    ThreadPool(unsigned int _nThreads = 0) :
        nThreads(_nThreads==0 ? defaultThreadCount() : _nThreads),
        workers(nThreads-1),
        queue(std::make_shared<ThreadJobQueue>())
    {
        for (auto& worker: workers) {
            worker.setQueue(queue);
        }
    }

    static unsigned int defaultThreadCount() {
        unsigned int n = std::thread::hardware_concurrency();

        return n == 0 ? 1 : n;
    }

    static ThreadPool& defaultPool() {
        static ThreadPool pool;

        return pool;
    }

    unsigned int getThreadCount() const {
        return nThreads;
    }

    static std::vector<JobRange>
    divideWork(int64_t begin, int64_t end, uint64_t threadCount) {

        const uint64_t elementCount = end - begin;
        const uint64_t chunkCount   = std::min(threadCount * chunksPerThread, elementCount);
        const uint64_t minJobSize   = elementCount / chunkCount;
        const uint64_t elementRest  = elementCount % chunkCount;

        std::vector<JobRange> jobRanges;
        jobRanges.reserve(chunkCount);

        int64_t jobBegin = begin;

        for (int64_t i = 0; i < chunkCount; i++) {
            int64_t jobSize = minJobSize;

            if (i < elementRest) {
                jobSize += 1;
            }

            jobRanges.push_back({jobBegin, jobBegin += jobSize});
        }

        return jobRanges;
    }

    template<typename Func>
    void parallelFor(int begin, int end, Func&& func) {

        if (begin >= end) {
            return;
        }

        std::lock_guard<std::mutex> poolLock(poolMutex);

        const auto     jobs = divideWork(begin, end, nThreads);
        const uint64_t jobCount = jobs.size();
        const int64_t  threadCount = std::min(nThreads, jobCount);
        const uint64_t curThreadId = threadCount - 1;
        int64_t        k = 0;

        for (int i = 0; i < threadCount - 1; i++) {
            workers[i].setThreadId(i);
        }

        using namespace std::placeholders;

        for (k = threadCount; k < jobCount;  k++) {
            queue->unsafePush(std::bind(std::ref(func), jobs[k].begin, jobs[k].end, _1));
        }

        for (k = 0; k < threadCount - 1; k++) {
            workers[k].submit(std::bind(std::ref(func), jobs[k].begin, jobs[k].end, _1));
        }

        func(jobs[k].begin, jobs[k].end, curThreadId);

        ThreadJob queuedJob;

        while (queue->tryPop(queuedJob)) {
            queuedJob(curThreadId);
        }

        for (k = 0; k < threadCount - 1; k++) {
            workers[k].wait();
        }
    }

    template<typename Func>
    void parallelBlock(Func&& func) {

        std::lock_guard<std::mutex> poolLock(poolMutex);

        for (int i = 0; i < nThreads - 1; i++) {
            workers[i].setThreadId(i);
        }

        using namespace std::placeholders;

        int64_t k = 0;
        for (; k < nThreads - 1; k++) {
            workers[k].submit(std::bind(std::ref(func), nThreads, _1));
        }

        func(nThreads, k);

        for (k = 0; k < nThreads - 1; k++) {
            workers[k].wait();
        }
    }
};

#endif // MISC_H
