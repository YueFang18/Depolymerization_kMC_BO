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

// 只有声明，没有实现
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