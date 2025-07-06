//
// Created by 13206 on 2024/1/26.
//

#ifndef GLASSEFFECTPROPAGATIONRATE_H
#define GLASSEFFECTPROPAGATIONRATE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include "../include/DataStructureAndConstant.h"
#include <map>

double CalculatePropagationRate (double conversion, double T, const double k_pchem)
{
    const double s_p = 5.9e-10; // s_p is the propagation reaction distance
//    std::cout << "initial_concentration =  " << initial_concentration  << std::endl;
//    std::cout << "Na =  " <<  Na << std::endl;
//    std::cout << "s_p =  " <<  s_p << std::endl;
    const double D0 = 1.27e-7;
    const double Ea = 12; // in kj/mol
    const double V_star = 0.872e-3;
    const double K_1m = 6.91e-7;
    const double K_1p = 2.92e-7;
    const double K2m_Tgm = 72.26;
    const double K2p_Tgp = -250.21;
    const double MMA_massfraction = 1-0.003*1.61;
    double k_app;
    const double k_chem = k_pchem;
    double k_diff;
    double Dm; // the translation diffusion coefficient
    double w_p = conversion*MMA_massfraction; // the mass fraction of the monomer;
    double w_m = (1-conversion)*MMA_massfraction; // the mass fraction of the polymer;
    double V_FH;
    V_FH = w_m *(K_1m*(K2m_Tgm + T)) + w_p *(K_1p*(K2p_Tgp + T));
    Dm = D0*pow(e,-Ea/(R*T))*pow(e,(-V_star*(MMA_massfraction)/V_FH));
    k_diff = 4*Pi*s_p*Na*Dm;

//    if(conversion > 0.8) {
//    std::cout << "conversion = " << conversion << std::endl;
//    std::cout << " V_FH = " <<  V_FH << std::endl;
//    std::cout << " Dm = " <<  Dm << std::endl;
//    std::cout << "k_diff = "<< k_diff << std::endl;
//    }

    k_app = pow((pow(k_chem,-1)+pow(k_diff,-1)),-1);

    return k_app;
}

#endif
