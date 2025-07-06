#ifndef EFFPROPARECORDER_READPOLYMER_H
#define EFFPROPARECORDER_READPOLYMER_H

#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <algorithm>   // ← shuffle 需要
#include <random>      // ← 随机数引擎
#include "../include/DataStructureAndConstant.h"

void ReadPolymer(std::ifstream& inputFile) {
    std::cout << "start reading polymers " << std::endl;
    std::string line;
    std::getline(inputFile, line);                 // 跳过表头

    /* ------------ 1. 读入并分类 ------------- */
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        int chain1, chain2, num_M1, chain_feature;
        if (iss >> chain1 >> chain2 >> num_M1 >> chain_feature) {
            if (num_M1 > 1) {                      // 只保留 DP>1 的链
                M_initial[0] += num_M1;            // 总计单体数

                if (chain_feature == 0 || chain_feature == 00000) {
                    species[0]++;
                    unsaturated_Polymers.num_M1.push_back(num_M1);
                }
                else if (chain_feature == 1 || chain_feature == 11111) {
                    species[1]++;
                    saturated_Polymers.num_M1.push_back(num_M1);
                }
                else if (chain_feature == 2) {
                    species[2]++;
                    HH_Polymers.num_M1.push_back(num_M1);
                    HH_Polymers.chain1.push_back(chain1);
                    HH_Polymers.chain2.push_back(chain2);
                }
            }
        }
    }
    inputFile.close();

    /* ------------ 2. 随机“清洗”排序 ------------- */
    // 统一随机源，保证复现实验时可设种子：
    static std::mt19937 rng(std::random_device{}());

    // 2.1 单独一条向量的情况 —— 直接 shuffle
    std::shuffle(unsaturated_Polymers.num_M1.begin(),
                 unsaturated_Polymers.num_M1.end(), rng);

    std::shuffle(saturated_Polymers.num_M1.begin(),
                 saturated_Polymers.num_M1.end(),   rng);

    // 2.2 多条并行向量（HH）—— 先打乱索引，再重排
    const std::size_t nHH = HH_Polymers.num_M1.size();
    if (nHH > 1) {
        std::vector<std::size_t> idx(nHH);
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng);

        std::vector<int>   tmp_num;   tmp_num.reserve(nHH);
        std::vector<int>   tmp_c1;    tmp_c1.reserve(nHH);
        std::vector<int>   tmp_c2;    tmp_c2.reserve(nHH);

        for (std::size_t i : idx) {
            tmp_num.push_back(HH_Polymers.num_M1[i]);
            tmp_c1 .push_back(HH_Polymers.chain1[i]);
            tmp_c2 .push_back(HH_Polymers.chain2[i]);
        }
        HH_Polymers.num_M1.swap(tmp_num);
        HH_Polymers.chain1.swap(tmp_c1);
        HH_Polymers.chain2.swap(tmp_c2);
    }
}

#endif // EFFPROPARECORDER_READPOLYMER_H
