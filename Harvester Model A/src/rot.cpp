#include "rot.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

rot_optim::rot_optim()
{
    delta = 0.00005;
    st_st_delta = 0.001;
    Max_Counter = 200;

    magnet_a = 0.0025;
    magnet_b = 0.0020;

    mass = 9.83e-5;
    air_damp = 0.1;

    Inert = mass*(pow(magnet_b,2) + pow(magnet_a,2))/12.0;

    spring_k = Inert*pow(2.0*M_PI*510.0,2) + 0.5*pow(air_damp,2)*Inert;

    spring_k3= 0.0;
    spring_k5= 0.0;
    aext = 3.0;
    scale = 3.82;
    acc_scale = 5.0;


    cout<<Inert<<"\t"<<spring_k<<endl;

    experimental_freq = new double [10000];
    experimental_volt = new double [10000];
    experimental_volt_norm = new double [10000];
    experimental_freq_1 = new double [10000];
    experimental_volt_1 = new double [10000];
    experimental_volt_norm_1 = new double [10000];
    experimental_freq_2 = new double [10000];
    experimental_volt_2 = new double [10000];
    experimental_volt_norm_2 = new double [10000];
}

rot_optim::rot_optim(double A, double C, double K, double K3, double K5, double S, double as)
{
    delta = 0.00005;
    st_st_delta = 0.001;
    Max_Counter = 200;

    mass = 9.83e-5;
    air_damp = C;
    spring_k = K;
    spring_k3= K3;
    spring_k5= K5;
    aext = A;
    scale = S;
    acc_scale = as;

    magnet_a = 0.0025;
    magnet_b = 0.0020;
    Inert = 1/12.0*mass*(pow(magnet_b,2) + pow(magnet_a,2));

    spring_k = Inert*pow(2.0*M_PI*510.0,2) + 0.5*pow(air_damp,2)*Inert;

    cout<<Inert<<"\t"<<spring_k<<endl;

    experimental_freq = new double [10000];
    experimental_volt = new double [10000];
    experimental_volt_norm = new double [10000];
    experimental_freq_1 = new double [10000];
    experimental_volt_1 = new double [10000];
    experimental_volt_norm_1 = new double [10000];
    experimental_freq_2 = new double [10000];
    experimental_volt_2 = new double [10000];
    experimental_volt_norm_2 = new double [10000];
}

rot_optim::~rot_optim()
{
    delete experimental_freq;
    delete experimental_volt;
    delete experimental_volt_norm;
    delete experimental_freq_1;
    delete experimental_volt_1;
    delete experimental_volt_norm_1;
    delete experimental_freq_2;
    delete experimental_volt_2;
    delete experimental_volt_norm_2;
}

