#include <memory>
#include <cstring>
#include <algorithm>
#include "msm.hpp"
#include "misc.hpp"

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::run(typename Curve::Point &r,
                                typename Curve::PointAffine *_bases,
                                uint8_t* _scalars,
                                uint64_t _scalarSize,
                                uint64_t _n,
                                uint64_t _nThreads)
{
    if (_n == 0) {
        g.copy(r, g.zero());
        return;
    }
    if (_n == 1) {
        g.mulByScalar(r, _bases[0], _scalars, _scalarSize);
        return;
    }
    if (_scalarSize < 8) {
        runPartition(r, _bases, _scalars, _scalarSize, _scalarSize*8, _n);
        return;
    }

    ThreadPool &threadPool = ThreadPool::defaultPool();

    scalars = _scalars;
    scalarSize = _scalarSize;

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
        runPartition(r, _bases, _scalars, _scalarSize, overallMaxBits + 2, _n);
        return;
    }

    // All scalars fit in 64 bits: gather only the scalars, bases stay in place.
    if (nSmall == _n) {
        std::unique_ptr<uint64_t[]> smallScalars(new uint64_t[_n]);

        threadPool.parallelFor(0, _n, [&] (int begin, int end, int numThread) {
            for (int i = begin; i < end; i++) {
                std::memcpy(&smallScalars[i], _scalars + (uint64_t)i*_scalarSize, sizeof(uint64_t));
            }
        });

        runPartition(r, _bases, (uint8_t *)smallScalars.get(), sizeof(uint64_t), maxSmallBits + 2, _n);
        return;
    }

    std::unique_ptr<uint64_t[]> smallScalars(nSmall ? new uint64_t[nSmall] : nullptr);
    std::unique_ptr<typename Curve::PointAffine[]> smallBases(nSmall ? new typename Curve::PointAffine[nSmall] : nullptr);
    std::unique_ptr<uint8_t[]> bigScalars(nBig ? new uint8_t[nBig*_scalarSize] : nullptr);
    std::unique_ptr<typename Curve::PointAffine[]> bigBases(nBig ? new typename Curve::PointAffine[nBig] : nullptr);
    std::unique_ptr<typename Curve::Point[]> onesAcc(new typename Curve::Point[nBlocks]);

    threadPool.parallelFor(0, nBlocks, [&] (int begin, int end, int numThread) {
        for (int b = begin; b < end; b++) {
            const uint64_t i0 = (uint64_t)b*blockSize;
            const uint64_t i1 = std::min(i0 + blockSize, _n);
            uint64_t smallCur = blockOffsets[b*2];
            uint64_t bigCur   = blockOffsets[b*2+1];

            g.copy(onesAcc[b], g.zero());

            for (uint64_t i = i0; i < i1; i++) {
                switch (classes[i]) {
                case CLS_ONE:
                    g.add(onesAcc[b], onesAcc[b], _bases[i]);
                    break;
                case CLS_SMALL:
                    std::memcpy(&smallScalars[smallCur], _scalars + i*_scalarSize, sizeof(uint64_t));
                    smallBases[smallCur] = _bases[i];
                    smallCur++;
                    break;
                case CLS_BIG:
                    std::memcpy(&bigScalars[bigCur*_scalarSize], _scalars + i*_scalarSize, _scalarSize);
                    bigBases[bigCur] = _bases[i];
                    bigCur++;
                    break;
                default:
                    break;
                }
            }
        }
    });

    typename Curve::Point acc;

    g.copy(acc, onesAcc[0]);
    for (uint64_t b = 1; b < nBlocks; b++) {
        g.add(acc, acc, onesAcc[b]);
    }

    if (nSmall > 0) {
        typename Curve::Point rSmall;

        runPartition(rSmall, smallBases.get(), (uint8_t *)smallScalars.get(),
                     sizeof(uint64_t), maxSmallBits + 2, nSmall);
        g.add(acc, acc, rSmall);
    }
    if (nBig > 0) {
        typename Curve::Point rBig;

        runPartition(rBig, bigBases.get(), bigScalars.get(),
                     _scalarSize, maxBigBits + 2, nBig);
        g.add(acc, acc, rBig);
    }

    g.copy(r, acc);
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::runPartition(typename Curve::Point &r,
                                         typename Curve::PointAffine *_bases,
                                         uint8_t* _scalars,
                                         uint64_t _scalarSize,
                                         uint64_t _nBits,
                                         uint64_t _n)
{
    ThreadPool &threadPool = ThreadPool::defaultPool();

    const uint64_t nThreads = threadPool.getThreadCount();
    const uint64_t nPoints = _n;

    scalars = _scalars;
    scalarSize = _scalarSize;

#ifdef MSM_BITS_PER_CHUNK
    bitsPerChunk = MSM_BITS_PER_CHUNK;
#else
    bitsPerChunk = calcBitsPerChunk(nPoints, _nBits, nThreads);
#endif

    if (nPoints == 0) {
        g.copy(r, g.zero());
        return;
    }
    if (nPoints == 1) {
        g.mulByScalar(r, _bases[0], scalars, scalarSize);
        return;
    }

    const uint64_t nChunks = calcChunkCount(_nBits, bitsPerChunk);
    const uint64_t nBuckets = calcBucketCount(bitsPerChunk);
    const uint64_t matrixSize = nThreads * nBuckets;
    const uint64_t nSlices = nChunks*nPoints;

    std::unique_ptr<typename Curve::Point[]> bucketMatrix(new typename Curve::Point[matrixSize]);
    std::unique_ptr<typename Curve::Point[]> chunks(new typename Curve::Point[nChunks]);
    std::unique_ptr<int32_t[]> slicedScalars(new int32_t[nSlices]);

    threadPool.parallelFor(0, nPoints, [&] (int begin, int end, int numThread) {

        for (int i = begin; i < end; i++) {
            int carry = 0;

            for (int j = 0; j < nChunks; j++) {
                int bucketIndex = getBucketIndex(i, j) + carry;

                if (bucketIndex >= nBuckets) {
                    bucketIndex -= nBuckets*2;
                    carry = 1;
                } else {
                    carry = 0;
                }

                slicedScalars[j*nPoints + i] = bucketIndex;
            }
        }
    });

    threadPool.parallelFor(0, nChunks, [&] (int begin, int end, int numThread) {

        for (int j = begin; j < end; j++) {

            typename Curve::Point *buckets = &bucketMatrix[numThread*nBuckets];

            for (int i = 0; i < nBuckets; i++) {
                g.copy(buckets[i], g.zero());
            }

            for (int i = 0; i < nPoints; i++) {
                const int bucketIndex = slicedScalars[j*nPoints + i];

                if (bucketIndex > 0) {
                    g.add(buckets[bucketIndex-1], buckets[bucketIndex-1], _bases[i]);

                } else if (bucketIndex < 0) {
                    g.sub(buckets[-bucketIndex-1], buckets[-bucketIndex-1], _bases[i]);
                }
            }

            typename Curve::Point t, tmp;

            g.copy(t, buckets[nBuckets - 1]);
            g.copy(tmp, t);

            for (int i = nBuckets - 2; i >= 0 ; i--) {
                g.add(tmp, tmp, buckets[i]);
                g.add(t, t, tmp);
            }

            chunks[j] = t;
        }
    });

    g.copy(r, chunks[nChunks - 1]);

    for (int j = nChunks - 2; j >= 0; j--) {
        for (int i = 0; i < bitsPerChunk; i++) {
            g.dbl(r, r);
        }
        g.add(r, r, chunks[j]);
    }
}
