#ifdef USE_OPENMP
#include <omp.h>
#endif
#include <memory.h>
#include <iostream>
#include <vector>
#include <numeric>
//#include "../misc.hpp"
#include "msm.hpp"

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::run(typename Curve::Point &r, uint32_t _nThreads) {
#ifdef _OPENMP
    nThreads = _nThreads==0 ? omp_get_max_threads() : _nThreads;
    ThreadLimit threadLimit (nThreads);
#else
    nThreads = 1;
#endif

//    std::cout << "scalarSize: " << scalarSize << std::endl;
//    for (uint32_t i=0; i<n; i++) {
//        std::cout << "scalar " << i << ": " << (uint64_t)_scalars[i] << std::endl;
//    }

    for(uint32_t i=0; i<n; i++) {           // the for is iterating over the base points
        std::vector<uint32_t> slices = slice(i, scalarSize, bitsPerSlice);

         // print slices size and content
//         std::cout << "slices size: " << slices.size() << std::endl;
//         for (uint32_t j=0; j<slices.size(); j++) {
//             std::cout << "slice " << j << ": " << slices[j] << std::endl;
//         }

         // print the buckets size and elements
//        std::cout << "buckets size: " << buckets.size() << std::endl;
//        for (uint32_t j=0; j<buckets.size(); j++) {
//            typename Curve::Point point;
//            g.copy(point, buckets[j]); // Convert PointAffine to Point
//            std::string point_str = g.toString(point); // Now you can use the toString function
//            std::cout << "bucket " << j << ": " << point_str << std::endl;
//        }
        //curPoints
//        std::cout << "curPoints size: " << curPoints.size() << std::endl;
//        for (uint32_t j=0; j<curPoints.size(); j++) {
//            typename Curve::Point point;
//            g.copy(point, curPoints[j]); // Convert PointAffine to Point
//            std::string point_str = g.toString(point); // Now you can use the toString function
//            std::cout << "curPoints " << j << ": " << point_str << std::endl;
//        }


        processPointAndSlices(i, slices);
    }

    // finish processing remaining elements
    finalize();

    // combine all the results into a single point
    reduce(r);
}

template <typename Curve, typename BaseField>
std::vector<uint32_t> MSM<Curve, BaseField>::slice(uint32_t scalarIdx, uint32_t scalarSize, uint32_t bitsPerChunk) {
    std::vector<uint32_t> slices;
    uint32_t chunksCount = (scalarSize * 8 + bitsPerChunk - 1) / bitsPerChunk; // Calculate chunks count, same as slices in rust
    slices.resize(chunksCount);

    for (uint32_t chunkIdx = 0; chunkIdx < chunksCount; ++chunkIdx) {
        uint32_t bitStart = chunkIdx * bitsPerChunk;
        uint32_t byteStart = bitStart / 8;
        uint32_t efectiveBitsPerChunk = bitsPerChunk;
        if (byteStart > scalarSize-8) byteStart = scalarSize - 8;
        if (bitStart + bitsPerChunk > scalarSize*8) efectiveBitsPerChunk = scalarSize*8 - bitStart;
        uint32_t shift = bitStart - byteStart*8;
        uint64_t v = *(uint64_t *)(scalars + scalarIdx*scalarSize + byteStart);
        v = v >> shift;
        v = v & ( (1 << efectiveBitsPerChunk) - 1);
        slices[chunkIdx] = uint32_t(v);
    }
    return slices;
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::processPointAndSlices(uint32_t baseIdx, std::vector<uint32_t>& slices) {
    assert(nSlices == slices.size() && "slices size should equal num_windows");  // TODO Remove later for efficiency

    typename Curve::PointAffine point = bases[baseIdx];

    curPoints.push_back(point);

    for (int win = 0; win < slices.size(); win++) {
        if (slices[win] > 0) {
            uint32_t bucket_id = (win << bitsPerSlice) + slices[win] - 1; // skip slice == 0
            processSlices(bucket_id, curPoints.size() - 1);
        }
    }
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::processSlices(uint32_t bucket_id, uint32_t currentPointIdx) {
    // Check if the bucket_id bit is not set
    if (!bitmap->testAndSet(bucket_id)) {
        // If no collision found, add point to current batch
        batchBucketsAndPoints.push_back(std::make_pair(bucket_id, currentPointIdx));
    } else {
        // In case of collision, add it to the collision list
        collisionBucketsAndPoints.push_back(std::make_pair(bucket_id, curPoints[currentPointIdx]));
    }

    // If the count of collisions or the batch size reach their maximum limits, process the batch
    if (collisionBucketsAndPoints.size() >= maxCollisionCnt ||
            batchBucketsAndPoints.size() >= maxBatchCnt) {
        processBatch();
    }
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::processBatch(){
    if (batchBucketsAndPoints.empty()) {
        return;
    }

    // Separate bucket_ids and point_idxs from batch_buckets_and_points
    std::vector<uint32_t> bucket_ids, point_idxs;
    for (const auto &bp : batchBucketsAndPoints) {
        bucket_ids.push_back(bp.first);
        point_idxs.push_back(bp.second);
    }


    // Perform batch addition
    batchAdder->batchAddIndexed(buckets, nBuckets, bucket_ids, curPoints, point_idxs);

    // Clean up current batch
    bitmap->clear();
    batchBucketsAndPoints.clear();

    // Memorize the last point which is the current processing point
    auto slicing_point = curPoints.back();
    curPoints.pop_back();  // Remove the last point
    curPoints.clear();

    // Process collisions
    int next_pos = 0;
    for (int i = 0; i < collisionBucketsAndPoints.size(); i++) {
        auto bucket_and_point = collisionBucketsAndPoints[i];
        auto bucket_id = bucket_and_point.first; // or bucket_and_point.get<0>(), etc.
        auto point = bucket_and_point.second; // or bucket_and_point.get<1>(), etc.


        if (bitmap->testAndSet(bucket_id)) {
            // collision found
            std::swap(collisionBucketsAndPoints[next_pos], collisionBucketsAndPoints[i]);
            next_pos += 1;
        } else {
            batchBucketsAndPoints.push_back(std::make_pair(bucket_id, curPoints.size()));
            curPoints.push_back(point);
        }
    }


    // Remove processed collisions
    collisionBucketsAndPoints.resize(next_pos);

    // Push back the slicing_point to curPoints
    curPoints.push_back(slicing_point);
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::finalize() {
    processBatch();
    while ( !( (batchBucketsAndPoints.size() == 0)  && (collisionBucketsAndPoints.size() == 0) ) )  {
        processBatch();
    }
}

template <typename Curve, typename BaseField>
void MSM<Curve, BaseField>::reduce(typename Curve::Point &r) {
    std::vector<uint32_t> window_starts(nSlices);
    std::iota(window_starts.begin(), window_starts.end(), 0);

    std::vector<typename Curve::Point> window_sums;
    window_sums.reserve(window_starts.size());
    for (auto &w_start : window_starts) {
        size_t bucket_start = static_cast<size_t>(w_start << bitsPerSlice);
        size_t bucket_end = bucket_start + (1 << bitsPerSlice);
        window_sums.push_back(innerWindowReduce(bucket_start, bucket_end));
    }


    r = intraWindowReduce(window_sums);
}

template <typename Curve, typename BaseField>
typename Curve::Point MSM<Curve, BaseField>::intraWindowReduce(std::vector<typename Curve::Point> window_sums) {
    typename Curve::Point lowest = window_sums.front();
    typename Curve::Point total = g.zero();

    for (auto it = window_sums.rbegin(); it != window_sums.rend() - 1; ++it) {

        g.add(total, *it, total);
        for (int i = 0; i < bitsPerSlice; ++i) {
            g.dbl(total,total);
        }
    }

//    return lowest + total;
    g.add(total, lowest, total);
    return total;



//    for (auto rit = windowSums.rbegin() + 1; rit != windowSums.rend(); ++rit) {
//        total += *rit;
//        for (int i = 0; i < windowBits; ++i) {
//            total.double_in_place();
//        }
//    }
//
//    return lowest + total;
}



template <typename Curve, typename BaseField>
typename Curve::Point MSM<Curve, BaseField>::innerWindowReduce(size_t start, size_t end) {
    typename Curve::Point running_sum = g.zero();
    typename Curve::Point res = g.zero();

    for(int i = end - 1; i-- > start;) {
        g.add(running_sum, buckets[i], running_sum);
        g.add(res, running_sum, res);
    }

    return res;
}

