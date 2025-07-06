/*
 * ODE part
 */
#include <cstdio>
#include<iomanip>
#include<stdlib.h>
#include<vector>
#include<math.h>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<iterator>

#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;
typedef boost::numeric::ublas::vector<double> vector_type;
typedef boost::numeric::ublas::matrix<double> matrix_type;

 //ode for polymerization
double kt_re; // 全局变量
double kt_dis; // 全局变量
std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/Cooperation_project/P9_jinyang_depolymerization/J-H_MMA/PolymerizationOutput/20240624/SAA/";
double k[7] =
        {3.98805e-05 * 0.58,
         1035.49,
         1035.49,
         3.5595e+07,
         9.62385e+07,
         0.0143765,
         1035.49};

struct stiff_system {
    void operator()(const vector_type &x, vector_type &dxdt, double t
) {

        dxdt[0] = - k[0] * x[0];
        dxdt[1] =  2 * k[0] * x[0] - k[1] * x[1] * x[2];
        dxdt[2] =  - k[1] * x[1] * x[2] - k[2] * x[2] * x[3] - k[5] * x[2] * x[3] - k[6] * x[2] * x[5];
        dxdt[3] =  k[1] * x[1] * x[2] - k[3] * x[3] * x[3]- k[4] * x[3] * x[3] - k[5] * x[2] * x[3] + k[6] * x[5] * x[2];
        dxdt[4] =  k[3] * x[3] * x[3] + 2 * k[4] * x[3] * x[3] + k[5] * x[2] * x[3];
        dxdt[5] =  k[5] * x[2] * x[3] - k[6] * x[5] * x[2];
    }
};
struct stiff_system_jacobi {
    void operator()(const vector_type &x, matrix_type &J, const double &
 t
 , vector_type &dfdt) {
        //dxdt[0] = - k[0] * x[0];
        J(0,0) = -k[0];
        J(0,1) = 0;
        J(0,2) = 0;
        J(0,3) = 0;
        J(0,4) = 0;
        J(0,5) = 0;
        //dxdt[1] =  2 * k[0] * x[0] - k[1] * x[1] * x[2];
        J(1,0) = 2 * k[0];
        J(1,1) = -k[1] * x[2];
        J(1,2) = -k[1] * x[1];
        J(1,3) = 0;
        J(1,4) = 0;
        J(1,5) = 0;
        //dxdt[2] =  - k[1] * x[1] * x[2] - k[2] * x[2] * x[3] - k[5] * x[2] * x[3] - k[6] * x[2] * x[5];
        J(2,0) = 0;
        J(2,1) = -k[1] * x[2];
        J(2,2) = -k[1] * x[1] - k[2] * x[3] - k[5] * x[3] - k[6] * x[5];
        J(2,3) = -k[2] * x[2] - k[5] * x[2];
        J(2,4) = 0;
        J(2,5) = -k[6] * x[2];
        //dxdt[3] =  k[1] * x[1] * x[2] - k[3] * x[3] * x[3]- k[4] * x[3] * x[3] - k[5] * x[2] * x[3] + k[6] * x[5] * x[2];
        J(3,0) = 0;
        J(3,1) = k[1] * x[2];
        J(3,2) = k[1] * x[1] - k[5] * x[3] + k[6] * x[5];
        J(3,3) = - 2 * k[3] * x[3] - 2 * k[4] * x[3] - k[5] * x[2];
        J(3,4) = 0;
        J(3,5) = k[6] * x[2];
        //  dxdt[4] =  k[3] * x[3] * x[3] + 2 * k[4] * x[3] * x[3] + k[5] * x[2] * x[3];
        J(4,0) = 0;
        J(4,1) = 0;
        J(4,2) = k[5] * x[3];
        J(4,3) = 2 * k[3] * x[3] + 4 * k[4] * x[3] + k[5] * x[2];
        J(4,4) = 0;
        J(4,5) = 0;
        //dxdt[5] =  k[5] * x[2] * x[3] - k[6] * x[5] * x[2];
        J(5,0) = 0;
        J(5,1) = 0;
        J(5,2) = k[5] * x[3] - k[6] * x[5];
        J(5,3) = k[5] * x[2];
        J(5,4) = 0;
        J(5,5) = -k[6] * x[2];

        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
        dfdt[4] = 0.0;
        dfdt[5] = 0.0;
    }
};

void ode(stiff_system,stiff_system_jacobi,vector_type &r,vector_type &Rad,double end_t,double kt_re,double kt_dis,double kt_i,double kt_p ) {
    k[3] = kt_re;
    k[4] = kt_dis;
    k[0] = kt_i;
    k[2] = kt_p;
    vector_type x(6);
    x = r;
    ofstream t1(out_dir + "ode_data1.out", std::ios_base::app);
    size_t num_of_steps = integrate_const(make_dense_output<rosenbrock4<double> >(1.0e-6, 1.0e-6),
                                          make_pair(stiff_system(), stiff_system_jacobi()),
                                          x, 0.0, end_t, end_t/100,
                                          t1 << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " "
                                             << phoenix::arg_names::arg1[1] << " " << phoenix::arg_names::arg1[2] << " "
                                             << phoenix::arg_names::arg1[3] << " " << phoenix::arg_names::arg1[4] << " "
                                             << phoenix::arg_names::arg1[5] <<" "<< k[0]<< " "<< k[1]<<" "<< k[2]<< " "<< k[3]<<" "<< k[4]<< "\n");

    r = x;
    t1.close();
}

/*std::string out_dir = "/Users/yueyue/Dropbox/Mac/Documents/HKUST/Cooperation_project/P9_jinyang_depolymerization/J-H_MMA/DepolymerizationOutput/20241014ssa/";
double k[6] =
        {8.6e-4,
         2e-7,
         1.5e-2,
         5e3,
         1e4,
         1e4};

// 8.6e-4; 2e-7;1.5e-2; 5e3;1e4;1e4;

double phi_value;

struct stiff_system {
    void operator()(const vector_type &x, vector_type &dxdt, double *//* t *//*) {
        dxdt[0] = - k[0] * x[0] + (1 - phi_value) * k[5] * x[3] * x[3];
        dxdt[1] = - k[1] * x[1] +(1 - phi_value) *  k[5] * x[3] * x[3];
        dxdt[2] = - k[2] * x[2] + k[4] * x[3] * x[3];
        dxdt[3] =  2 * k[0] * x[0] +  2 * k[1] * x[1] + 2 * k[2] * x[2] - 2 *  k[4] * x[3] * x[3] -2 * k[5] * x[3] * x[3] ;
        dxdt[4] = (1 - phi_value) * k[3] * x[3] + 2 * phi_value * k[5] * x[3] * x[3] ;
    }
};
struct stiff_system_jacobi {
    void operator()(const vector_type &x, matrix_type &J, const double & *//* t *//* , vector_type &dfdt) {
        // dxdt[0] = - k[0] * x[0] + k[5] * x[3] * x[3];
        J(0,0) = -k[0];
        J(0,1) = 0;
        J(0,2) = 0;
        J(0,3) = 2 * k[5] * x[3] * (1 - phi_value);
        J(0,4) = 0;
        // dxdt[1] = - k[1] * x[1] + k[5] * x[3] * x[3];
        J(1,0) = 0;
        J(1,1) = -k[1];;
        J(1,2) = 0;
        J(1,3) =  2 * k[5] * x[3] * (1 - phi_value);
        J(1,4) = 0;
        //dxdt[2] =  - k[2] * x[2] + k[4] * x[3] * x[3];
        J(2,0) = 0;
        J(2,1) = 0;
        J(2,2) = -k[2];
        J(2,3) = 2 * k[4] * x[3] ;
        J(2,4) = 0;
        //dxdt[3] = 2 * k[0] * x[0] +  2 * k[1] * x[1] + 2 * k[2] * x[2] - 2 *  k[4] * x[3] * x[3] -2 * k[5] * x[3] * x[3];
        J(3,0) =  2 * k[0];
        J(3,1) =  2 * k[1];
        J(3,2) =  2 * k[2];
        J(3,3) = -4 *  k[4] * x[3]-4 * k[5] * x[3] ;
        J(3,4) = 0;
        // dxdt[4] = (1 - phi_value) * k[3] * x[3] + 2 * phi_value* k[5] * x[3] * x[3] ;
        J(4,0) = 0;
        J(4,1) = 0;
        J(4,2) = 0;
        J(4,3) = k[3] * (1 - phi_value) + 4 * phi_value* k[5] * x[3];
        J(4,4) = 0;

        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
        dfdt[4] = 0.0;
    }
};


void ode(stiff_system,stiff_system_jacobi,vector_type &r,vector_type &Rad,double end_t, double phi, double time) {

 //   std::cout << "end_t before calling ode: " << end_t << std::endl;
    end_t = end_t *1.0;

    phi_value = phi;
    vector_type x(5);
    x = r;
    ofstream t1(out_dir + "ode_data1.out", std::ios_base::app);
    size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-6 , 1.0e-6 ) ,
                                           make_pair( stiff_system() , stiff_system_jacobi() ) ,
                                           x , 0.0 , end_t ,end_t/10.0,
                                           t1 << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << " " << phoenix::arg_names::arg1[2] << " " << phoenix::arg_names::arg1[3] << " " << phoenix::arg_names::arg1[4] << " " << phi
                                              << " " <<time<<"\n" );

    r = x;
    t1.close();
}*/


