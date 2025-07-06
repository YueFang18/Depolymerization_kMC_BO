//
// Created by yueyue on 2025/6/23.
//

#include "LengthDependentTerminationRate.h"

// 全局变量定义（之前在头文件中重复定义）
double isl = 100;
double alpha_s = 0.65;
double alpha_l = 0.15;
double averageKt;

double CalculateKt_ii(double w_p, int i) {
    double igel = 0.53 * pow(w_p, -1.8);
    double alphagel = 1.66 * w_p - 0.06;
    double terminationRate;

    if (i >= igel) {
        if (i < isl){
            terminationRate = k_11 * pow(igel, alphagel - alpha_s) * pow(i, -alpha_s);
        }
        if(i >= isl){
            terminationRate = k_11 * pow(isl, alpha_l-alpha_s)*pow(igel,alphagel-alpha_l) * pow(i,-alphagel);
        }
    }
    if(i < igel)
    {
        if (i < isl) {
            terminationRate = k_11 * pow(i, -alpha_s);
        }if (i >= isl) {
            terminationRate = k_11 * pow(isl, alpha_l - alpha_s) * pow(i, -alpha_l);
        }
    }
    return terminationRate;
}

double CalculateKt_ij(double w_p, int i, int j)
{
    double Kt_i = CalculateKt_ii(w_p,i);
    double Kt_j = CalculateKt_ii(w_p,j);
    return std::sqrt(Kt_i*Kt_j);
}

double PopulationWeightedTerminationRate(const std::vector<std::vector<int>>& RadicalChains, double conversion) {
    if (RadicalChains.empty())
    {
        return k_11;
    }
    double w_p = conversion;
    double AverageKt = 0.0;
    double totalRadicalNumber = RadicalChains.size();
    std::map<int, int> RadicalLengthCount;
    for (const auto& chain : RadicalChains)
    {
        RadicalLengthCount[chain[0]]++;
    }

    for (const auto& Radical_i : RadicalLengthCount) {
        double RadicalLength_i = Radical_i.first;
        double RadicalNumber_i = Radical_i.second;
        for (const auto& Radical_j : RadicalLengthCount) {
            double RadicalLength_j = Radical_j.first;
            double RadicalNumber_j = Radical_j.second;
            double NumberWeight = RadicalNumber_i * RadicalNumber_j / (totalRadicalNumber * totalRadicalNumber);
            AverageKt += NumberWeight * CalculateKt_ij(w_p, RadicalLength_i, RadicalLength_j);
        }
    }
    return AverageKt;
}