#ifndef MSM_HPP
#define MSM_HPP

#include <cstdint>

template <typename Curve, typename BaseField>
class MSM {
    const uint64_t MIN_CHUNK_SIZE_BITS = 3;
    const uint64_t MAX_CHUNK_SIZE_BITS = 16;

    // Scalars of at most this many significant bits go to the small
    // partition, which runs the bucket method over a single 64-bit word
    // and so pays ~4 windows instead of ~16. Scalars 0 and 1 are cheaper
    // still: they need no scalar multiplication at all.
    static const uint64_t SMALL_SCALAR_BITS = 64;

    Curve &g;
    uint8_t *scalars;
    uint64_t scalarSize;
    uint64_t bitsPerChunk;

private:
    // Estimated wall-clock cost in point additions. Chunks run in parallel,
    // so the elapsed time is the per-chunk cost times the number of waves of
    // nThreads chunks; the second term charges 1/8 of the total work so that
    // among near-equal wall costs the one burning fewer total additions wins.
    uint64_t calcCost(uint64_t nPoints, uint64_t nBits, uint64_t bitsPerChunk,
                      uint64_t nThreads) const {
        const uint64_t chunkCost = nPoints + ((uint64_t)1 << bitsPerChunk) + bitsPerChunk + 1;
        const uint64_t nChunks = calcChunkCount(nBits, bitsPerChunk);
        const uint64_t waves = (nChunks + nThreads - 1) / nThreads;

        return waves*chunkCost + nChunks*chunkCost/(8*nThreads);
    }

    uint64_t calcBitsPerChunk(uint64_t n, uint64_t nBits, uint64_t nThreads) const {
        uint64_t bitsPerChunk = MIN_CHUNK_SIZE_BITS;
        uint64_t minCost = calcCost(n, nBits, bitsPerChunk, nThreads);

        for (uint64_t k = MIN_CHUNK_SIZE_BITS + 1; k <= MAX_CHUNK_SIZE_BITS; k++) {
            const uint64_t curCost = calcCost(n, nBits, k, nThreads);

            if (curCost < minCost) {
                minCost = curCost;
                bitsPerChunk = k;
            }
        }
        return bitsPerChunk;
    }

    uint64_t calcChunkCount(uint64_t nBits, uint64_t bitsPerChunk) const {
        return ((nBits - 1) / bitsPerChunk) + 1;
    }

    uint64_t calcBucketCount(uint64_t bitsPerChunk) const {
        return ((uint64_t)1 << (bitsPerChunk-1));
    }

    uint64_t getBucketIndex(uint64_t scalarIdx, uint64_t chunkIdx) const {
        uint64_t bitStart = chunkIdx*bitsPerChunk;

        // Chunks past the scalar bytes exist only to absorb the signed-digit
        // carry; their digit is zero.
        if (bitStart >= scalarSize*8) return 0;

        uint64_t byteStart = bitStart/8;
        uint64_t efectiveBitsPerChunk = bitsPerChunk;

        if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
        if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;

        uint64_t shift = bitStart - byteStart*8;
        uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + byteStart);

        v = v >> shift;
        v = v & ( ((uint64_t)1 << efectiveBitsPerChunk) - 1);

        return uint64_t(v);
    }

    // Number of significant bits of a little-endian scalar; 0 for a zero scalar.
    uint64_t significantBits(const uint8_t *scalar) const {
        for (int64_t k = (int64_t)scalarSize - 1; k >= 0; k--) {
            if (scalar[k]) {
                return (uint64_t)k*8 + (32 - __builtin_clz((uint32_t)scalar[k]));
            }
        }
        return 0;
    }

public:
    MSM(Curve &_g): g(_g) {}

    // Partitions the scalars by significant bits (0, 1, up to 64 bits, wider)
    // and runs the bucket method separately per partition, so the mostly-0/1
    // scalars of a circom witness don't pay full-width window costs.
    void run(typename Curve::Point &r,
             typename Curve::PointAffine *_bases,
             uint8_t* _scalars,
             uint64_t _scalarSize,
             uint64_t _n,
             uint64_t _nThreads=0);

    // One bucket-method pass over all n points using only the lowest nBits
    // bits of each scalar. nBits must exceed the largest scalar's significant
    // bit count by at least 2, so a signed-digit carry can never propagate
    // out of the top chunk.
    void runPartition(typename Curve::Point &r,
                      typename Curve::PointAffine *bases,
                      uint8_t* scalars,
                      uint64_t scalarSize,
                      uint64_t nBits,
                      uint64_t n);
};

#include "msm.cpp"

#endif // MSM_HPP
