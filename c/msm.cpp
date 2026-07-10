#include <memory>
#include <cstring>
#include <algorithm>
#include "msm.hpp"
#include "misc.hpp"

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::preparePartition(Partition &p, uint64_t nThreads)
{
    ThreadPool &threadPool = ThreadPool::defaultPool();

#ifdef MSM_BITS_PER_CHUNK
    p.bitsPerChunk = MSM_BITS_PER_CHUNK;
    p.nSlices = 1;
#else
    calcChunkConfig(p.n, p.nBits, nThreads, p.bitsPerChunk, p.nSlices);
#endif

    p.nChunks = calcChunkCount(p.nBits, p.bitsPerChunk);
    p.nBuckets = calcBucketCount(p.bitsPerChunk);
    p.digits.reset(new int32_t[p.nChunks * p.n]);
    p.partials.reset(new typename Curve::Point[p.nSlices * p.nChunks]);

    // recode context for getBucketIndex
    scalars = p.scalars;
    scalarSize = p.scalarSize;
    bitsPerChunk = p.bitsPerChunk;

    const uint64_t nChunks = p.nChunks;
    const uint64_t nBuckets = p.nBuckets;
    const uint64_t nPoints = p.n;
    int32_t *digits = p.digits.get();

    threadPool.parallelFor(0, nPoints, [&, nChunks, nBuckets, nPoints] (int begin, int end, int numThread) {

        for (int i = begin; i < end; i++) {
            int carry = 0;

            for (uint64_t j = 0; j < nChunks; j++) {
                int bucketIndex = getBucketIndex(i, j) + carry;

                if (bucketIndex >= (int)nBuckets) {
                    bucketIndex -= nBuckets*2;
                    carry = 1;
                } else {
                    carry = 0;
                }

                digits[j*nPoints + i] = bucketIndex;
            }
        }
    });
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::prepare(typename Curve::PointAffine *_bases,
                                    uint8_t *_scalars,
                                    uint64_t _scalarSize,
                                    uint64_t _n,
                                    uint64_t parallelismShare)
{
    ThreadPool &threadPool = ThreadPool::defaultPool();

    const uint64_t nThreads = parallelismShare ? parallelismShare
                                               : threadPool.getThreadCount();

    partitions.clear();
    partitions.reserve(2);
    onesAcc.reset();
    nOnesBlocks = 0;
    trivial = false;
    prepared = true;

    if (_n == 0) {
        trivial = true;
        g.copy(trivialResult, g.zero());
        return;
    }
    if (_n == 1) {
        trivial = true;
        g.mulByScalar(trivialResult, _bases[0], _scalars, _scalarSize);
        return;
    }

    scalars = _scalars;
    scalarSize = _scalarSize;

    if (_scalarSize < 8) {
        partitions.emplace_back();
        Partition &p = partitions.back();
        p.bases = _bases;
        p.scalars = _scalars;
        p.scalarSize = _scalarSize;
        p.n = _n;
        p.nBits = _scalarSize*8;
        preparePartition(p, nThreads);
        return;
    }

    const uint64_t nBlocks = std::min<uint64_t>(threadPool.getThreadCount()*4, _n);
    const uint64_t blockSize = (_n + nBlocks - 1) / nBlocks;

    enum ScalarClass : uint8_t { CLS_ZERO = 0, CLS_ONE = 1, CLS_SMALL = 2, CLS_BIG = 3 };

    std::unique_ptr<uint8_t[]> classes(new uint8_t[_n]);
    std::unique_ptr<uint64_t[]> blockCounts(new uint64_t[nBlocks*3]);
    std::unique_ptr<uint64_t[]> blockMaxBits(new uint64_t[nBlocks*2]);

    threadPool.parallelFor(0, nBlocks, [&] (int begin, int end, int numThread) {
        for (int b = begin; b < end; b++) {
            const uint64_t i0 = (uint64_t)b*blockSize;
            const uint64_t i1 = std::min(i0 + blockSize, _n);
            uint64_t nSmall = 0, nBig = 0, nOnes = 0;
            uint64_t maxSmall = 0, maxBig = 0;

            for (uint64_t i = i0; i < i1; i++) {
                const uint64_t bits = significantBits(_scalars + i*_scalarSize);
                uint8_t cls;

                if (bits == 0) {
                    cls = CLS_ZERO;
                } else if (bits == 1) {
                    cls = CLS_ONE;
                    nOnes++;
                } else if (bits <= SMALL_SCALAR_BITS) {
                    cls = CLS_SMALL;
                    nSmall++;
                    if (bits > maxSmall) maxSmall = bits;
                } else {
                    cls = CLS_BIG;
                    nBig++;
                    if (bits > maxBig) maxBig = bits;
                }
                classes[i] = cls;
            }
            blockCounts[b*3]   = nSmall;
            blockCounts[b*3+1] = nBig;
            blockCounts[b*3+2] = nOnes;
            blockMaxBits[b*2]   = maxSmall;
            blockMaxBits[b*2+1] = maxBig;
        }
    });

    uint64_t nSmall = 0, nBig = 0, nOnes = 0;
    uint64_t maxSmallBits = 0, maxBigBits = 0;

    std::unique_ptr<uint64_t[]> blockOffsets(new uint64_t[nBlocks*2]);

    for (uint64_t b = 0; b < nBlocks; b++) {
        blockOffsets[b*2]   = nSmall;
        blockOffsets[b*2+1] = nBig;
        nSmall += blockCounts[b*3];
        nBig   += blockCounts[b*3+1];
        nOnes  += blockCounts[b*3+2];
        if (blockMaxBits[b*2]   > maxSmallBits) maxSmallBits = blockMaxBits[b*2];
        if (blockMaxBits[b*2+1] > maxBigBits)   maxBigBits   = blockMaxBits[b*2+1];
    }

    const uint64_t overallMaxBits = std::max(maxBigBits, std::max(maxSmallBits, (uint64_t)(nOnes ? 1 : 0)));

    // When almost every scalar is full width (e.g. the H MSM, whose scalars
    // are uniform field elements) partitioning saves nothing: run the whole
    // input in place instead of paying the gather.
    if (nBig >= _n - _n/16) {
        partitions.emplace_back();
        Partition &p = partitions.back();
        p.bases = _bases;
        p.scalars = _scalars;
        p.scalarSize = _scalarSize;
        p.n = _n;
        p.nBits = overallMaxBits + 2;
        preparePartition(p, nThreads);
        return;
    }

    // All scalars fit in 64 bits: gather only the scalars, bases stay in place.
    if (nSmall == _n) {
        partitions.emplace_back();
        Partition &p = partitions.back();
        p.ownScalars64.reset(new uint64_t[_n]);

        uint64_t *s64 = p.ownScalars64.get();

        threadPool.parallelFor(0, _n, [&, s64] (int begin, int end, int numThread) {
            for (int i = begin; i < end; i++) {
                std::memcpy(&s64[i], _scalars + (uint64_t)i*_scalarSize, sizeof(uint64_t));
            }
        });

        p.bases = _bases;
        p.scalars = (uint8_t *)s64;
        p.scalarSize = sizeof(uint64_t);
        p.n = _n;
        p.nBits = maxSmallBits + 2;
        preparePartition(p, nThreads);
        return;
    }

    Partition *small = NULL;
    Partition *big = NULL;

    if (nSmall > 0) {
        partitions.emplace_back();
        small = &partitions.back();
        small->ownScalars64.reset(new uint64_t[nSmall]);
        small->ownBases.reset(new typename Curve::PointAffine[nSmall]);
        small->bases = small->ownBases.get();
        small->scalars = (uint8_t *)small->ownScalars64.get();
        small->scalarSize = sizeof(uint64_t);
        small->n = nSmall;
        small->nBits = maxSmallBits + 2;
    }
    if (nBig > 0) {
        partitions.emplace_back();
        big = &partitions.back();
        big->ownScalars.reset(new uint8_t[nBig*_scalarSize]);
        big->ownBases.reset(new typename Curve::PointAffine[nBig]);
        big->bases = big->ownBases.get();
        big->scalars = big->ownScalars.get();
        big->scalarSize = _scalarSize;
        big->n = nBig;
        big->nBits = maxBigBits + 2;
    }

    nOnesBlocks = nBlocks;
    onesAcc.reset(new typename Curve::Point[nBlocks]);

    typename Curve::Point *ones = onesAcc.get();

    threadPool.parallelFor(0, nBlocks, [&, ones] (int begin, int end, int numThread) {
        for (int b = begin; b < end; b++) {
            const uint64_t i0 = (uint64_t)b*blockSize;
            const uint64_t i1 = std::min(i0 + blockSize, _n);
            uint64_t smallCur = blockOffsets[b*2];
            uint64_t bigCur   = blockOffsets[b*2+1];

            g.copy(ones[b], g.zero());

            for (uint64_t i = i0; i < i1; i++) {
                switch (classes[i]) {
                case CLS_ONE:
                    g.add(ones[b], ones[b], _bases[i]);
                    break;
                case CLS_SMALL:
                    std::memcpy(&small->ownScalars64[smallCur], _scalars + i*_scalarSize, sizeof(uint64_t));
                    small->ownBases[smallCur] = _bases[i];
                    smallCur++;
                    break;
                case CLS_BIG:
                    std::memcpy(&big->ownScalars[bigCur*_scalarSize], _scalars + i*_scalarSize, _scalarSize);
                    big->ownBases[bigCur] = _bases[i];
                    bigCur++;
                    break;
                default:
                    break;
                }
            }
        }
    });

    if (small) preparePartition(*small, nThreads);
    if (big)   preparePartition(*big, nThreads);
}

template <typename Curve, typename BaseField>
uint64_t MSM<Curve, BaseField>::maxBuckets() const
{
    uint64_t m = 0;

    for (const Partition &p : partitions) {
        if (p.nBuckets > m) m = p.nBuckets;
    }
    return m;
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::collectTasks(std::vector<Task> &tasks,
                                         typename Curve::Point *bucketArena,
                                         uint64_t bucketsPerThread)
{
    for (Partition &part : partitions) {
        Partition *p = &part;

        for (uint64_t s = 0; s < p->nSlices; s++) {
            const uint64_t i0 = p->n * s / p->nSlices;
            const uint64_t i1 = p->n * (s+1) / p->nSlices;

            for (uint64_t j = 0; j < p->nChunks; j++) {
                tasks.push_back([this, p, s, j, i0, i1, bucketArena, bucketsPerThread] (uint64_t threadId) {
                    typename Curve::Point *buckets = &bucketArena[threadId*bucketsPerThread];
                    const int32_t *digits = &p->digits[j*p->n];
                    typename Curve::PointAffine *bases = p->bases;
                    const uint64_t nBuckets = p->nBuckets;

                    for (uint64_t i = 0; i < nBuckets; i++) {
                        g.copy(buckets[i], g.zero());
                    }

                    for (uint64_t i = i0; i < i1; i++) {
                        const int32_t bucketIndex = digits[i];

                        if (bucketIndex > 0) {
                            g.add(buckets[bucketIndex-1], buckets[bucketIndex-1], bases[i]);

                        } else if (bucketIndex < 0) {
                            g.sub(buckets[-bucketIndex-1], buckets[-bucketIndex-1], bases[i]);
                        }
                    }

                    typename Curve::Point t, tmp;

                    g.copy(t, buckets[nBuckets - 1]);
                    g.copy(tmp, t);

                    for (int64_t i = nBuckets - 2; i >= 0 ; i--) {
                        g.add(tmp, tmp, buckets[i]);
                        g.add(t, t, tmp);
                    }

                    p->partials[s*p->nChunks + j] = t;
                });
            }
        }
    }
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::reducePartition(Partition &p, typename Curve::Point &r)
{
    typename Curve::Point chunkSum;

    for (int64_t j = p.nChunks - 1; j >= 0; j--) {
        g.copy(chunkSum, p.partials[j]);
        for (uint64_t s = 1; s < p.nSlices; s++) {
            g.add(chunkSum, chunkSum, p.partials[s*p.nChunks + j]);
        }

        if (j == (int64_t)p.nChunks - 1) {
            g.copy(r, chunkSum);
        } else {
            g.add(r, r, chunkSum);
        }

        if (j > 0) {
            for (uint64_t b = 0; b < p.bitsPerChunk; b++) {
                g.dbl(r, r);
            }
        }
    }
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::finish(typename Curve::Point &r)
{
    if (trivial) {
        g.copy(r, trivialResult);
        prepared = false;
        return;
    }

    typename Curve::Point acc, part;

    g.copy(acc, g.zero());

    for (Partition &p : partitions) {
        reducePartition(p, part);
        g.add(acc, acc, part);
    }

    for (uint64_t b = 0; b < nOnesBlocks; b++) {
        g.add(acc, acc, onesAcc[b]);
    }

    g.copy(r, acc);

    partitions.clear();
    onesAcc.reset();
    nOnesBlocks = 0;
    prepared = false;
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::run(typename Curve::Point &r,
                                typename Curve::PointAffine *_bases,
                                uint8_t* _scalars,
                                uint64_t _scalarSize,
                                uint64_t _n,
                                uint64_t _nThreads)
{
    ThreadPool &threadPool = ThreadPool::defaultPool();

    prepare(_bases, _scalars, _scalarSize, _n);

    if (!trivial) {
        const uint64_t nThreads = threadPool.getThreadCount();
        const uint64_t bucketsPerThread = maxBuckets();

        std::unique_ptr<typename Curve::Point[]> arena(
            new typename Curve::Point[nThreads * bucketsPerThread]);

        std::vector<Task> tasks;
        collectTasks(tasks, arena.get(), bucketsPerThread);

        if (!tasks.empty()) {
            threadPool.parallelFor(0, tasks.size(), [&] (int begin, int end, int numThread) {
                for (int t = begin; t < end; t++) {
                    tasks[t]((uint64_t)numThread);
                }
            });
        }
    }

    finish(r);
}
