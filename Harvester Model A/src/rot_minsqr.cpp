#include "rot.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

double rot_optim::minsqr_norm()
{
    /** Definition of the maximal frequency of the maximum */
    double v_max = 0.0;
    double f_max = 0.0;
    for(int i=0; i<imax; i++)
    {
        if(v_max<experimental_volt[i])
        {
            v_max = experimental_volt[i];
            f_max = experimental_freq[i];
        }
    }
    /** Definition of the maximal theoretical value */
    int n_sol;
    double* sol = new double [10];
    Solve(sol, n_sol, 2*M_PI*f_max);
    double max_sol=0;
    for(int i=0; i<n_sol; i++)
    {
        if(sol[i]>max_sol)
            max_sol = sol[i];
    }
    double max_theor_vrms = vrms_interp_rotate(1.0,max_sol,2*M_PI*f_max*max_sol);

    /** Min square method */
    double sol_1;
    double vrms_1;
    double Sum = 0.0;
    for(int i=0; i<imax; i++)
    {
        Solve(sol, n_sol, 2*M_PI*experimental_freq[i]);
        sol_1 = 0;
        for(int j=0; j<n_sol; j++)
        {
            //cout<<sol[j]<<endl;
            if(sol[j]>sol_1)
                sol_1 = sol[j];
        }

        vrms_1 = vrms_interp_rotate(1.0,sol_1,2*M_PI*experimental_freq[i]*sol_1)/max_theor_vrms;
        Sum += pow(experimental_volt_norm[i] - vrms_1,2);
    }
    delete sol;
    return(Sum);

}

double rot_optim::minsqr()
{
    int n_sol;
    double* sol = new double [10];
    /** Min square method */
    double sol_1;
    double vrms_1;
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
        vrms_1 = vrms_interp_rotate(scale,sol_1,2*M_PI*experimental_freq[i]*sol_1);
        Sum += pow(experimental_volt[i] - vrms_1,2);
    }
    delete sol;
    return(Sum);

}

double rot_optim::minsqr_sets()
{
    int n_sol;
    double Sum = 0.0;
    double* sol = new double [10];

    /** Min square method for a= 3.0 m/s^2*/
    double vrms_1[10];
    double dif;
    aext = 3.0;
    for(int i=0; i<imax; i++)
    {
        Solve(sol, n_sol, 2*M_PI*experimental_freq[i]);

        for(int j=0; j<n_sol; j++)
        {
            vrms_1[j] = vrms_interp_rotate(scale,sol[j],2*M_PI*experimental_freq[i]*sol[j]);
        }
        dif = pow(experimental_volt[i] - vrms_1[0],2);
        for(int j=1; j<n_sol; j++)
        {
            if( pow(experimental_volt[i] - vrms_1[j],2)<dif )
            {
                dif = pow(experimental_volt[i] - vrms_1[j],2);
            }
        }

        Sum += dif;
    }
    /** Min square method for a= 4.0 m/s^2*/
    aext = 4.0;
    for(int i=0; i<imax1; i++)
    {
        Solve(sol, n_sol, 2*M_PI*experimental_freq_1[i]);

        for(int j=0; j<n_sol; j++)
        {
            vrms_1[j] = vrms_interp_rotate(scale,sol[j],2*M_PI*experimental_freq_1[i]*sol[j]);
        }
        dif = pow(experimental_volt_1[i] - vrms_1[0],2);
        for(int j=1; j<n_sol; j++)
        {
            if( pow(experimental_volt_1[i] - vrms_1[j],2)<dif )
            {
                dif = pow(experimental_volt_1[i] - vrms_1[j],2);
            }
        }

        Sum += dif;
    }
    /** Min square method for a= 5.0 m/s^2*/
    aext = 5.0;
    for(int i=0; i<imax2; i++)
    {
        Solve(sol, n_sol, 2*M_PI*experimental_freq_2[i]);

        for(int j=0; j<n_sol; j++)
        {
            vrms_1[j] = vrms_interp_rotate(scale,sol[j],2*M_PI*experimental_freq_2[i]*sol[j]);
        }
        dif = pow(experimental_volt_2[i] - vrms_1[0],2);
        for(int j=1; j<n_sol; j++)
        {
            if( pow(experimental_volt_2[i] - vrms_1[j],2)<dif )
            {
                dif = pow(experimental_volt_2[i] - vrms_1[j],2);
            }
        }

        Sum += dif;
    }


    delete sol;
    return(Sum);

}
