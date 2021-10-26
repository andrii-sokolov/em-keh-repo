#include "transp.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include "gnuplot-iostream.h"

using namespace std;

progress_optim::progress_optim()
{
    delta = 0.00001;
    st_st_delta = 0.001;
    Max_Counter = 200;

    mass = 9.83e-5;
    air_damp = 0.15;
    spring_k = 697.3;
    spring_k3= 1.1e11;
    spring_k5= 0.0;
    spring_k7= 0.0;
    aext = 3.0;
    scale = 4.0;

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

progress_optim::progress_optim(double A, double C, double K, double K3, double K5, double K7, double S)
{
    delta = 0.00005;
    st_st_delta = 0.001;
    Max_Counter = 200;

    mass = 9.83e-5;
    air_damp = C;
    spring_k = K;
    spring_k3= K3;
    spring_k5= K5;
    spring_k7= K7;
    aext = A;
    scale = S;

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

progress_optim::~progress_optim()
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


