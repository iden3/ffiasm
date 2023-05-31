#ifndef PAR_MULTIEXP3
#define PAR_MULTIEXP3

#include <memory>
#include "bitmap.h"
#include "batch_adder.h"
#include "../misc.hpp"

#define PME2_PACK_FACTOR 2
#define PME2_MAX_CHUNK_SIZE_BITS 16
#define PME2_MIN_CHUNK_SIZE_BITS 2

const size_t MAX_BATCH_SIZE = 1024;
const size_t MAX_COLLISION_SIZE = 2048;

template <typename Curve, typename BaseField>
class MSM {

    typename Curve::PointAffine *bases;
    uint8_t* scalars;
    uint32_t scalarSize;
    uint32_t n;
    uint32_t nThreads;
    uint32_t bitsPerSlice;
    Curve &g;
    uint32_t nSlices;
    uint32_t nBuckets;

    std::vector<typename Curve::PointAffine> curPoints;
//    std::vector<typename Curve::PointAffine> buckets;
    typename Curve::PointAffine* buckets;
//    Bitmap bitmap;
    std::unique_ptr<Bitmap> bitmap;

    //    BatchAdder<Curve, BaseField> batchAdder;
    std::unique_ptr<BatchAdder<Curve, BaseField>> batchAdder;



    std::vector<std::pair<uint32_t, uint32_t>> batchBucketsAndPoints;
    uint32_t maxBatchCnt;
    std::vector<std::pair<uint32_t, typename Curve::PointAffine>> collisionBucketsAndPoints;
    uint32_t maxCollisionCnt;

public:

    MSM(Curve &_g,typename Curve::PointAffine *_bases,uint8_t* _scalars,uint32_t _scalarSize,uint32_t _n) : g(_g){
        // calculate and stored all members of the MSM.
        // *******************************************

        // double check if to remove any constraints from here (all consts are from previous implementation)
        bitsPerSlice = log2(_n / PME2_PACK_FACTOR);
        if (bitsPerSlice > PME2_MAX_CHUNK_SIZE_BITS) bitsPerSlice = PME2_MAX_CHUNK_SIZE_BITS;
        if (bitsPerSlice < PME2_MIN_CHUNK_SIZE_BITS) bitsPerSlice = PME2_MIN_CHUNK_SIZE_BITS;
        nSlices = ((_scalarSize * 8 - 1 ) / bitsPerSlice) + 1;
        nBuckets = nSlices << bitsPerSlice;

//        g = _g;
        bases = _bases;
        scalars = _scalars;
        scalarSize = _scalarSize;
        n = _n;

        maxBatchCnt = MAX_BATCH_SIZE;
        maxCollisionCnt = MAX_COLLISION_SIZE;

//        batchAdder = BatchAdder<Curve, BaseField>(g, MAX_BATCH_SIZE);
        batchAdder =  std::unique_ptr<BatchAdder<Curve, BaseField>>(new BatchAdder<Curve, BaseField>(g, MAX_BATCH_SIZE));

        // we treat each bucket as one bit in the bit map, so every 32 buckets are represented by one uint32_t
//        bitmap = Bitmap(nBuckets / 32);
        bitmap =  std::unique_ptr<Bitmap>(new Bitmap(nBuckets / 32));

//        for(int i = 0; i < maxBatchCnt; i++) {
//            curPoints[i] = g.zeroAffine();
//        }

        buckets = new typename Curve::PointAffine[nBuckets]; // allocate memory
        // initialize with zeros
        for(int i = 0; i < nBuckets; i++) {
            buckets[i] = g.zeroAffine();
        }


    }

    ~MSM() {
        delete [] buckets; // deallocate
    }

    void run(typename Curve::Point &r, uint32_t _nThreads= 0);
    std::vector<uint32_t> slice(uint32_t scalarIdx, uint32_t scalarSize, uint32_t bitsPerChunk);
    void processPointAndSlices(uint32_t baseIdx, std::vector<uint32_t>& slices);
    void processSlices(uint32_t bucket_id, uint32_t currentPointIdx) ;
    void processBatch();
    void finalize();
    void reduce(typename Curve::Point &r);
    typename Curve::Point innerWindowReduce(size_t start, size_t end);
    typename Curve::Point intraWindowReduce(std::vector<typename Curve::Point> window_sums);

    };

#include "msm.cpp"

#endif // PAR_MULTIEXP3
