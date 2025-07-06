#ifndef EFFPROPARECORDER_DATASTRUCTUREANDCONSTANT_H
#define EFFPROPARECORDER_DATASTRUCTUREANDCONSTANT_H

#include<iomanip>
#include<list>
#include<cstdlib>
#include<vector>
#include<cmath>
#include<algorithm>
#include<iterator>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

struct poly_chain {
    poly_chain(){};
    std::vector<int> chain1, chain2;
    std::vector<int> num_M1;
    std::vector<int> chain_feature;
    std::vector<int> chain1_feature;
    std::vector<int> chain2_feature;
};

// 全局变量声明 - 使用extern避免重复定义
extern poly_chain poly;
extern poly_chain unsaturated_Polymers;
extern poly_chain saturated_Polymers;
extern poly_chain HH_Polymers;

extern std::vector<int> validIndices;
extern std::vector<std::vector<int>> Pr;

// 序列记录向量声明
extern std::vector<int> M1M1;
extern std::vector<int> M1M2;
extern std::vector<int> M2M1;
extern std::vector<int> M2M2;

// 常数定义
const double Na = 6.022e23;
const double R = 8.314*1e-3;    //kJ/mol/K
const double Pi = 3.1415926;
const double MonomerWeight = 100.0;
const double k_11 = pow(10, 9.1);
const double e = 2.71828;

// 其他全局变量声明
extern double DisproportionationCoefficient;
extern double species[6];
extern double num_monomer[2];
extern double M_initial[2];

#endif //EFFPROPARECORDER_DATASTRUCTUREANDCONSTANT_H