#include <memory>
#include "msm.hpp"
#include "misc.hpp"

template <typename Curve, typename BaseField>
uint64_t MSM<Curve, BaseField>::getBitsPerChunk(uint64_t n, uint64_t scalarSize) const
{
#ifdef MSM_BITS_PER_CHUNK
    return MSM_BITS_PER_CHUNK;
#else
    return calcBitsPerChunk(n, scalarSize);
#endif
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::run(typename Curve::Point &r,
                                typename Curve::PointAffine *_bases,
                                uint8_t* _scalars,
                                uint64_t _scalarSize,
                                uint64_t nPoints,
                                uint64_t _nThreads)
{
    if (nPoints == 0) {
        g.copy(r, g.zero());
        return;
    }
    if (nPoints == 1) {
        g.mulByScalar(r, _bases[0], _scalars, _scalarSize);
        return;
    }

    ThreadPool &threadPool = ThreadPool::defaultPool();

    scalars = _scalars;
    scalarSize = _scalarSize;
    bitsPerChunk = getBitsPerChunk(nPoints, scalarSize);

    const uint64_t nThreads = threadPool.getThreadCount();
    const uint64_t nChunks  = calcChunkCount(scalarSize, bitsPerChunk);
    const uint64_t nBuckets = calcBucketCount(bitsPerChunk);

    std::vector<typename Curve::Point> bucketMatrix(nThreads * nBuckets);
    std::vector<typename Curve::Point> chunks(nChunks);
    std::vector<int32_t>               slicedScalars(nChunks * nPoints);

    threadPool.parallelFor(0, nPoints, [&] (int64_t begin, int64_t end, uint64_t idThread) {

        for (int64_t i = begin; i < end; i++) {
            int32_t carry = 0;

            for (int64_t j = 0; j < nChunks; j++) {
                int32_t bucketIndex = getBucketIndex(i, j) + carry;

                if (bucketIndex >= nBuckets) {
                    bucketIndex -= nBuckets*2;
                    carry = 1;
                } else {
                    carry = 0;
                }

                slicedScalars[i*nChunks + j] = bucketIndex;
            }
        }
    });

    threadPool.parallelFor(0, nChunks, [&] (int64_t begin, int64_t end, uint64_t idThread) {

        for (int64_t j = begin; j < end; j++) {

            typename Curve::Point *buckets = &bucketMatrix[idThread*nBuckets];

            for (int64_t i = 0; i < nBuckets; i++) {
                g.copy(buckets[i], g.zero());
            }

            for (int64_t i = 0; i < nPoints; i++) {
                const int64_t bucketIndex = slicedScalars[i*nChunks + j];

                if (bucketIndex > 0) {
                    g.add(buckets[bucketIndex-1], buckets[bucketIndex-1], _bases[i]);

                } else if (bucketIndex < 0) {
                    g.sub(buckets[-bucketIndex-1], buckets[-bucketIndex-1], _bases[i]);
                }
            }

            typename Curve::Point t, tmp;

            g.copy(t, buckets[nBuckets - 1]);
            g.copy(tmp, t);

            for (int64_t i = nBuckets - 2; i >= 0 ; i--) {
                g.add(tmp, tmp, buckets[i]);
                g.add(t, t, tmp);
            }

            chunks[j] = t;
        }
    });

    g.copy(r, chunks[nChunks - 1]);

    for (int64_t j = nChunks - 2; j >= 0; j--) {
        for (int64_t i = 0; i < bitsPerChunk; i++) {
            g.dbl(r, r);
        }
        g.add(r, r, chunks[j]);
    }
}
