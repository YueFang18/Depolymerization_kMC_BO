#ifndef DEPOLYMERIZATION_EXPLICIT_SEQUENCE_RECORD_NUMBER_H
#define DEPOLYMERIZATION_EXPLICIT_SEQUENCE_RECORD_NUMBER_H

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
#include <cassert>
#include <vector>
#include <random>
#include <chrono>
#include <numeric>

using namespace std;

// 只有声明，没有实现
std::pair<int, int> generate_two_different_indices(int range);

void depolymerization_explicit_sequence_record_number(int reaction_index,
                                                      double *num_monomer,
                                                      int &num_chain,
                                                      vector<int> &M1M1,
                                                      vector<int> &M1M2,
                                                      vector<int> &M2M1,
                                                      vector<int> &M2M2,
                                                      poly_chain &poly,
                                                      double *tot_dyad,
                                                      double *tot_triad,
                                                      double monomerYield);

#endif // DEPOLYMERIZATION_EXPLICIT_SEQUENCE_RECORD_NUMBER_H