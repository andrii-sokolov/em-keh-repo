#include "optim.h"
#include <math.h>
#include <iostream>

double optimization::minsqr()
{
    int n_sol;
    double* sol = new double [10];
    double sol_1;
    double vrms;
    double Sum = 0.0;
    for(int i=0; i<imax; i++)
    {
        Solve(sol, n_sol, 2*M_PI*experimental_freq[i]);
        sol_1 = 0;
        for(int j=0; j<n_sol; j++)
        {
            if(sol[j]>sol_1)
                sol_1 = sol[j];
        }
        vrms = vrms_interp_progress(sol_1,2*M_PI*experimental_freq[i]);
        Sum += pow(experimental_volt[i] - vrms,2);
    }
    delete sol;
    return(Sum);

}

double optimization::minsqr   ( double* exp_freq,
                                double* exp_volt,
                                int length
                              )
{
    int n_sol;
    double* sol = new double [10];
    double vrms;
    double Sum = 0.0;
    for(int i=0; i<length; i++)
    {
        Solve(sol, n_sol, 2.0*M_PI*exp_freq[i]);
        vrms = vrms_interp_progress(sol[0],2.0*M_PI*exp_freq[i]);
        Sum += pow(exp_volt[i] - vrms,2);
    }
    delete sol;
    return(Sum);
}

double optimization::minsqrLin   (  double* exp_freq,
                                    double* exp_volt,
                                    int length
                                 )
{
    double Ampl;
    double vrms;
    double Sum = 0.0;
    for(int i=0; i<length; i++)
    {
        Ampl = Solve_Linear(exp_freq[i]);
        vrms = vrms_interp_progress(Ampl,2.0*M_PI*exp_freq[i]);
        Sum += pow(exp_volt[i] - vrms,2);
    }
    return(Sum);
}
