#ifndef MSM_HPP
#define MSM_HPP

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

// Pippenger bucket-method MSM with scalar-size partitioning and a
// task-based execution model.
//
// Usage (single MSM): run() — prepare, execute and reduce in one call.
//
// Usage (batched, e.g. the prover's A/B1/B2/C phase):
//     msm.prepare(bases, scalars, size, n, parallelismShare);
//     msm.collectTasks(tasks, arena);   // arena: nThreads*maxBuckets() Points
//     ... run all MSMs' tasks in one parallel region ...
//     msm.finish(r);
// Tasks from several MSM instances (even over different curves) can be
// mixed in one region; each task gets the executing thread id and uses
// that thread's row of its own curve's bucket arena.
template <typename Curve, typename BaseField>
class MSM {
public:
    typedef std::function<void(uint64_t threadId)> Task;

private:
    const uint64_t MIN_CHUNK_SIZE_BITS = 3;
    const uint64_t MAX_CHUNK_SIZE_BITS = 16;

    // Scalars of at most this many significant bits go to the small
    // partition, which runs the bucket method over a single 64-bit word
    // and so pays ~4 windows instead of ~16. Scalars 0 and 1 are cheaper
    // still: they need no scalar multiplication at all.
    static const uint64_t SMALL_SCALAR_BITS = 64;

    // Don't point-split below this many points per slice: the per-slice
    // running-sum cost (2^c bucket additions) would dominate.
    static const uint64_t MIN_POINTS_PER_SLICE = 4096;

    // One scalar-size class of the input, ready for bucket accumulation.
    struct Partition {
        typename Curve::PointAffine *bases;     // points (caller's or gathered)
        uint8_t *scalars;                       // scalars (caller's or gathered)
        uint64_t scalarSize;
        uint64_t n;
        uint64_t nBits;                         // significant bits + carry headroom
        uint64_t bitsPerChunk;
        uint64_t nChunks;
        uint64_t nBuckets;
        uint64_t nSlices;                       // point-split factor
        std::unique_ptr<int32_t[]> digits;      // chunk-major [nChunks][n]
        std::unique_ptr<typename Curve::Point[]> partials; // [nSlices][nChunks]

        // backing storage when the class was gathered
        std::unique_ptr<typename Curve::PointAffine[]> ownBases;
        std::unique_ptr<uint64_t[]> ownScalars64;
        std::unique_ptr<uint8_t[]> ownScalars;
    };

    Curve &g;

    // significantBits()/getBucketIndex() context for the classification and
    // recode passes; set before each pass.
    uint8_t *scalars;
    uint64_t scalarSize;
    uint64_t bitsPerChunk;

    std::vector<Partition> partitions;
    std::unique_ptr<typename Curve::Point[]> onesAcc;   // per-block partial sums of 1-scalar points
    uint64_t nOnesBlocks;
    bool prepared;

    // Set when prepare() resolved the whole MSM without bucket work
    // (n==0, n==1, or every scalar in {0,1}).
    bool trivial;
    typename Curve::Point trivialResult;

private:
    uint64_t calcChunkCount(uint64_t nBits, uint64_t bitsPerChunk) const {
        return ((nBits - 1) / bitsPerChunk) + 1;
    }

    uint64_t calcBucketCount(uint64_t bitsPerChunk) const {
        return ((uint64_t)1 << (bitsPerChunk-1));
    }

    // Estimated wall-clock cost in point additions of one partition executed
    // as nSlices*nChunks tasks on nThreads threads. Tasks run in waves; the
    // second term charges 1/8 of the total work so that among near-equal
    // wall costs the one burning fewer total additions wins.
    uint64_t calcCost(uint64_t n, uint64_t nBits, uint64_t bitsPerChunk,
                      uint64_t nSlices, uint64_t nThreads) const {
        const uint64_t sliceCost = n/nSlices + ((uint64_t)1 << bitsPerChunk) + bitsPerChunk + 1;
        const uint64_t nTasks = nSlices * calcChunkCount(nBits, bitsPerChunk);
        const uint64_t waves = (nTasks + nThreads - 1) / nThreads;

        return waves*sliceCost + nTasks*sliceCost/(8*nThreads);
    }

    // Pick window size and point-split factor minimizing estimated wall cost.
    void calcChunkConfig(uint64_t n, uint64_t nBits, uint64_t nThreads,
                         uint64_t &bestC, uint64_t &bestSlices) const {
        const uint64_t maxSlices = std::max<uint64_t>(1, std::min<uint64_t>(
            2*nThreads, n / MIN_POINTS_PER_SLICE));

        bestC = MIN_CHUNK_SIZE_BITS;
        bestSlices = 1;
        uint64_t minCost = calcCost(n, nBits, bestC, 1, nThreads);

        for (uint64_t k = MIN_CHUNK_SIZE_BITS; k <= MAX_CHUNK_SIZE_BITS; k++) {
            for (uint64_t s = 1; s <= maxSlices; s *= 2) {
                const uint64_t curCost = calcCost(n, nBits, k, s, nThreads);

                if (curCost < minCost) {
                    minCost = curCost;
                    bestC = k;
                    bestSlices = s;
                }
            }
        }
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

    // Recode a partition's scalars into signed digits and size its partials.
    void preparePartition(Partition &p, uint64_t nThreads);

    // Reduce one partition's task partials into a single point.
    void reducePartition(Partition &p, typename Curve::Point &r);

public:
    MSM(Curve &_g): g(_g), prepared(false), trivial(false) {}

    // Classify scalars by size, gather the classes and recode digits.
    // parallelismShare: number of threads this MSM should assume it has for
    // itself when sizing windows/slices — pass the pool size when the MSM
    // runs alone, or roughly poolSize/nMSMs when batched with others.
    void prepare(typename Curve::PointAffine *bases,
                 uint8_t *scalars,
                 uint64_t scalarSize,
                 uint64_t n,
                 uint64_t parallelismShare = 0);

    // Largest bucket row any of this MSM's tasks needs; the caller provides
    // an arena of nThreads*maxBuckets() Points to collectTasks().
    uint64_t maxBuckets() const;

    // Append one task per (partition, slice, chunk). Tasks only touch their
    // own partials and bucketArena[threadId*bucketsPerThread..]. When several
    // MSMs share an arena, bucketsPerThread is the max of their maxBuckets().
    void collectTasks(std::vector<Task> &tasks,
                      typename Curve::Point *bucketArena,
                      uint64_t bucketsPerThread);

    // Reduce all partials into the final result. Call after every task ran.
    void finish(typename Curve::Point &r);

    // Single-MSM convenience: prepare + run own tasks + finish.
    void run(typename Curve::Point &r,
             typename Curve::PointAffine *_bases,
             uint8_t* _scalars,
             uint64_t _scalarSize,
             uint64_t _n,
             uint64_t _nThreads=0);
};

#include "msm.cpp"

#endif // MSM_HPP
