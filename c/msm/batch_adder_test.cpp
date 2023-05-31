#include <iostream>
#include <vector>
#include "batch_adder.h"
#include "../alt_bn128.hpp"

//#include "../../../../build/fq.hpp" //#include "fq.hpp"

using namespace std;

int main() {

    AltBn128::Engine engineInstance;

    const size_t max_batch_cnt = 10;
    const size_t array_size = 5;

    AltBn128::G1PointAffine* dest = new AltBn128::G1PointAffine[array_size];
    vector<AltBn128::G1PointAffine> src(array_size);

    vector<uint32_t> dest_indices = {0, 1, 2, 3, 4};
    vector<uint32_t> src_indices = {4, 3, 2, 1, 0};

    BatchAdder<AltBn128::Engine::G1, AltBn128::Engine::F1> batchAdder(
            engineInstance.g1,
            max_batch_cnt
    );

    batchAdder.batchAddIndexed(dest, array_size, dest_indices,src, src_indices);
    // Print the results
    for (size_t i = 0; i < array_size; ++i) {
        cout << "dest[" << i << "] = " << dest[i].x.v << endl;
    }

    return 0;
}
