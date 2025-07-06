//
// Created by 13206 on 2023/8/30.
//
#include <map>
#include <fstream>
#include <iomanip> // for formatting output
#include "DataStructureAndConstant.h"

#ifndef EFFPROPARECORDER_MOLECULARWEIGHTDISTRIBUTION_H
#define EFFPROPARECORDER_MOLECULARWEIGHTDISTRIBUTION_H

std::map<int, int> chainLengthCount = {};
std::map<int, double> chainLengthNumberFraction = {};
std::map<int, double> chainLengthMassFraction = {};
std::map<int, double> logMFraction = {};

void PolymerLengthDistribution(std::ofstream& output,
                               const std::vector<int>& PolymerChains)
{
    // Clear previous data
    chainLengthCount.clear();
    chainLengthNumberFraction.clear();
    chainLengthMassFraction.clear();
    logMFraction.clear();

    for (int length : PolymerChains) {
        chainLengthCount[length]++;
    } // Count the number of polymer chains for each chain length

    double totalPolymerMass = 0;
    double totallogM = 0;
    for (const auto& entry : chainLengthCount) {
        double chainLength = entry.first;
        double chainNumber = entry.second;
        totalPolymerMass += chainNumber * chainLength * 100;
        totallogM  += log(chainNumber * chainLength * 100);
    }

    for (auto& entry : chainLengthCount) {
        double chainLength = entry.first;
        double chainNumber = entry.second;
        chainLengthNumberFraction[chainLength] = chainNumber / PolymerChains.size();
        chainLengthMassFraction[chainLength] = (chainLength * chainNumber / totalPolymerMass) ;
        logMFraction[chainLength] = log(chainLength * 100 * chainNumber) / totallogM ;
    }
    // Output the molecular weight distribution
    output << "Chain Length\tNumber Fraction\tMass Fraction\tLogM Fraction\n";
    for (auto& entry : chainLengthCount) {
        double chainLength = entry.first;
        output << chainLength << "\t"
               << chainLengthNumberFraction[chainLength] << "\t"
                << chainLengthMassFraction[chainLength]<< "\t"
               << logMFraction[chainLength] << "\n";
    }
}

void LogMolarMassDistribution(std::ofstream& spline_logMW_distribution,
                              double molarWeight)
{
    struct polymerChains{
        std::vector<double> log_MW;
        std::vector<double> numberFraction;
        std::vector<double> massFraction;
    };
    struct polymerChains polymers, spline_polymers;
    double chainLength,numberFraction,massFraction,logM;
    polymers = { {}, {} , {} };
    spline_polymers = { {}, {} , {} };

    for (auto& entry : chainLengthCount) {
        chainLength = entry.first;
        logM = log10((chainLength)*molarWeight);
        polymers.log_MW.push_back(logM);
        polymers.numberFraction.push_back(chainLengthNumberFraction[chainLength]);
        polymers.massFraction.push_back(chainLengthMassFraction[chainLength]);
    }

    double interval = 0.1; // set a interval for spline correction
    for (double logM = 0.0; logM < 8.0; logM += interval) {
        double sumNumberFraction = 0.0;  // 初始化聚合物摩尔分数之和为0
        double sumMassFraction = 0.0;
        // 遍历polymers中的数据
        for (size_t i = 0; i < polymers.log_MW.size(); i++) {
            // 判断数据点的分子量对数值是否在当前范围内
            if (polymers.log_MW[i] >= logM && polymers.log_MW[i] < logM + interval)
            {
                sumNumberFraction += polymers.numberFraction[i];  // 将摩尔分数累加到sumNumberFraction
                sumMassFraction += polymers.massFraction[i];
            }
        }

        spline_polymers.log_MW.push_back(logM + interval/2.0);  // 将当前范围的对数值中点存入spline_polymers的log_MW
        spline_polymers.numberFraction.push_back(sumNumberFraction);  // 存入对应的摩尔分数之和
        spline_polymers.massFraction.push_back(sumMassFraction);
    }

    spline_logMW_distribution << "log(M)" <<" "<< "number fraction" << "mass fraction" << std::endl;

    for (size_t i = 0; i < spline_polymers.log_MW.size(); i++) {
        spline_logMW_distribution   << spline_polymers.log_MW[i] << " "
                                    << spline_polymers.numberFraction[i] << " "
                                    << spline_polymers.massFraction[i] << std::endl;
    }
    spline_logMW_distribution.close();

    polymers.log_MW.clear();
    polymers.numberFraction.clear();
    polymers.massFraction.clear();
    spline_polymers.log_MW.clear();
    spline_polymers.numberFraction.clear();
    spline_polymers.massFraction.clear();

}



#endif // EFFPROPARECORDER_MOLECULARWEIGHTDISTRIBUTION_H