//
// Created by yueyue on 2025/6/23.
//

#include "my_rand.h"

double my_rand(double range)
{
    static thread_local std::mt19937_64 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> distribution(0.0, range);
    return distribution(generator);
}

void randomlySwap(int& r1, int& r2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    if (dis(gen) == 1) {
        std::swap(r1, r2);
    }
}