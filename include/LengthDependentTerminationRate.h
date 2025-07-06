#ifndef LENGTHDEPENDENTTERMINATIONRATE_H
#define LENGTHDEPENDENTTERMINATIONRATE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "DataStructureAndConstant.h"
#include <map>

// 关键修改：将全局变量改为extern声明
extern double isl;
extern double alpha_s;
extern double alpha_l;
extern double averageKt;

struct ChainData {
    int count;
    double terminationRate;
};

// 只保留函数声明
double CalculateKt_ii(double w_p, int i);
double CalculateKt_ij(double w_p, int i, int j);
double PopulationWeightedTerminationRate(const std::vector<std::vector<int>>& RadicalChains, double conversion);

#endif