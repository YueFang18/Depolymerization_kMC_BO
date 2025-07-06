#ifndef INITIATOREFFICIENCY_H
#define INITIATOREFFICIENCY_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

double CalculateInitiatorEfficiency(const double T, double conversion)
{
    const double D_0I = 1.87e-8;
    const double E_I = 7.1e3; //J/mol
    const double V_m = 0.822e-6;
    const double V_p = 0.77e-6;
    const double Kmm = 1.49e-9;
    const double Kpm_Tg1 = -84.4;
    const double Kmp = 5.82e-10;
    const double Kpp_Tg1 = -327;
    const double Epsilon_cp = 0.354;
    const double Epsilon_mp = 0.59;
    double D_I, D_term, eff;
    double e = 2.718;
    double R = 8.314; //J/(mol.K)
    D_term = 5.3e-10;
    double w_p = conversion;
    double w_m = 1 - w_p;
    double VFH;
    VFH = Kmm*w_m*(Kpm_Tg1+T) + Kmp*w_p*(Kpp_Tg1+T);
    D_I = D_0I*pow(e,-E_I/(R*T))*pow(e,(-(w_m*V_m*Epsilon_cp/Epsilon_mp+w_p*V_p*Epsilon_cp)/VFH));
    eff = 1 - D_term/(D_I+D_term);
    return eff;
}

#endif

// =================================================================
// 修复后的 include/efficient_explicit_sequence_record_number.h
// =================================================================

#ifndef EFFICIENT_EXPLICIT_SEQUENCE_RECORD_NUMBER_H
#define EFFICIENT_EXPLICIT_SEQUENCE_RECORD_NUMBER_H

#include<iomanip>
#include<list>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>
#include "my_rand.h"
#include "LengthDependentTerminationRate.h"
#include "DataStructureAndConstant.h"
#include <vector>
#include <utility>
#include <random>
#include <iosfwd>
#include <unordered_map>
#include <limits>

using namespace std;

struct chain {
    vector<int> M1, M2;
    int dyad[4];
    int triad[8];
};

typedef vector<chain> chain_pool;

// 函数声明 - 移除inline和static
std::pair<int,int> ChooseChain(const std::vector<std::vector<int>>& Pr, double wp);

void efficient_explicit_sequence_record_number(int reaction_index,
                                               double *num_monomer,
                                               double *species,
                                               int &num_chain,
                                               vector<int> &M1M1,
                                               vector<int> &M1M2,
                                               vector<int> &M2M1,
                                               vector<int> &M2M2,
                                               poly_chain &poly,
                                               double w_p);

#endif // EFFICIENT_EXPLICIT_SEQUENCE_RECORD_NUMBER_H