//
// Created by 13206 on 2023/8/30.
//

#ifndef EFFPROPARECORDER_MOLECULARWEIGHTOUTPUT_H
#define EFFPROPARECORDER_MOLECULARWEIGHTOUTPUT_H
#include "DataStructureAndConstant.h"

void MolecularWeightOutput(ofstream& output,
                           const vector<std::vector<int>>& Radicals,
                           const vector<int>& Polymers,
                           double conversion)
{
    double Mn = 0;
    long double Mw = 0;
    std::map<int, int> chainLengthCount;
    for (const auto& chain : Radicals) {
        chainLengthCount[chain[0]]++;
    }
    for (int length : Polymers) {
        chainLengthCount[length]++;
    }

    double totalPolymerCount = Radicals.size() + Polymers.size();
    long double totalPolymerMass = 0.0;
    for (auto &entry : chainLengthCount) {
        double chainLength = entry.first;
        double chainNumber = entry.second;
        totalPolymerMass += chainNumber * chainLength;
    }

    for (auto &entry : chainLengthCount) {
        double chainLength = entry.first;
        double chainNumber = entry.second;
        double NumberFraction = chainNumber / totalPolymerCount;
        double WeightFraction = (chainLength / totalPolymerMass) * chainNumber;
        Mn += NumberFraction * chainLength;
        Mw += WeightFraction * chainLength;
    }

    Mn = Mn * MonomerWeight;
    Mw = Mw * MonomerWeight;

    output << conversion << "\t" << Mn << "\t" << Mw << endl;
}

#endif //EFFPROPARECORDER_MOLECULARWEIGHTOUTPUT_H
