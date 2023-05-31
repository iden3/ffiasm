//
// Created by Ron Kahat on 19/05/2023.
//

#include <iostream>
#include "bitmap.h"

int main() {
    Bitmap bitmap(1 << 4);

    std::cout << "Testing Bitmap class\n";

    // Test the bitmap functionality
    for (uint32_t i = 0; i < 2 << 4; ++i) {
        std::cout << "Setting bucket " << i << ": " << bitmap.testAndSet(i) << '\n';
    }

    for (uint32_t i = 0; i < 2 << 4; ++i) {
        std::cout << "Setting bucket " << i << ": " << bitmap.testAndSet(i) << '\n';
    }

    // Clear the bitmap
    bitmap.clear();

    std::cout << "Bitmap cleared\n";

    return 0;
}