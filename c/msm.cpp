#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <memory>
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
    const uint64_t nPoints = _n;

#ifdef _OPENMP
    const uint64_t nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    ThreadLimit threadLimit(nThreads);
#else
    const uint64_t nThreads = 1;
#endif

    scalars = _scalars;
    scalarSize = _scalarSize;
    bitsPerChunk = calcBitsPerChunk(nPoints, scalarSize);

    if (nPoints == 0) {
        g.copy(r, g.zero());
        return;
    }
    if (nPoints == 1) {
        g.mulByScalar(r, _bases[0], scalars, scalarSize);
        return;
    }

    const uint64_t nChunks = calcChunkCount(scalarSize, bitsPerChunk);
    const uint64_t nBuckets = calcBucketCount(bitsPerChunk);
    const uint64_t matrixSize = nThreads * nBuckets;

    std::unique_ptr<typename Curve::Point[]> bucketMatrix(new typename Curve::Point[matrixSize]);
    std::unique_ptr<typename Curve::Point[]> chunks(new typename Curve::Point[nChunks]);

    #pragma omp parallel for
    for (int j = 0; j < nChunks; j++) {

#ifdef _OPENMP
        const int numThread = omp_get_thread_num();
#else
        const int numThread = 0;
#endif

        typename Curve::Point *buckets = &bucketMatrix[numThread*nBuckets];

        for (int i = 0; i < nBuckets; i++) {
            g.copy(buckets[i], g.zero());
        }

        for (int i = 0; i < nPoints; i++) {
            const int bucketIndex = getBucketIndex(i, j);

            if (bucketIndex > 0) {
                typename Curve::Point &bucket = buckets[bucketIndex];

                g.add(bucket, bucket, _bases[i]);
            }
        }

        typename Curve::Point t, tmp;

        g.copy(t, buckets[nBuckets - 1]);
        g.copy(tmp, t);

        for (int i = nBuckets - 2; i >= 1 ; i--) {
            g.add(tmp, tmp, buckets[i]);
            g.add(t, t, tmp);
        }

        chunks[j] = t;
    }

    g.copy(r, chunks[nChunks - 1]);

    for (int j = nChunks - 2; j >= 0; j--) {
        for (int i = 0; i < bitsPerChunk; i++) {
            g.dbl(r, r);
        }
        g.add(r, r, chunks[j]);
    }
}
