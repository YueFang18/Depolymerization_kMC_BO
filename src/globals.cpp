#include "DataStructureAndConstant.h"

// 全局变量定义
poly_chain poly;
poly_chain unsaturated_Polymers;
poly_chain saturated_Polymers;
poly_chain HH_Polymers;

std::vector<int> validIndices;
std::vector<std::vector<int>> Pr;

// 序列记录向量定义
std::vector<int> M1M1;
std::vector<int> M1M2;
std::vector<int> M2M1;
std::vector<int> M2M2;

// 其他全局变量定义
double DisproportionationCoefficient = 0.73;
double species[6] = {0};
double num_monomer[2] = {0.0};
double M_initial[2] = {0.0};