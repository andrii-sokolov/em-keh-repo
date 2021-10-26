#include "rot.h"
#include <math.h>
#include <fstream>

double rot_optim::AdvDuffEQN(double A0, double wext)
{
    return(pow((-Inert*A0*wext*wext + spring_k*A0 + 0.75*pow(A0,3)*spring_k3 + 0.625*pow(A0,5)*spring_k5),2) + pow(air_damp*A0*wext,2) - pow(Inert*acc_scale*aext,2));
    //return(pow((-Inert*A0*wext*wext + spring_k*A0),2) + pow(Inert*air_damp*A0*wext,2) - pow(Inert*aext,2));
}

void rot_optim::Solve(double* solutions, int& n_of_solutions, double wext)
{
    double A1, A2;
    n_of_solutions = 0;
    double dA = 1e-4;
    for(double A = 1e-30; A<=1.0; A+= dA)
    {
        if( AdvDuffEQN(A,wext)*AdvDuffEQN(A+dA,wext)<0)
        {

            A1 = A;
            A2 = A + dA;
            while( 2*abs(A2 - A1)/abs(A1 + A2)>0.00001 )
            {
                if( AdvDuffEQN(A1,wext)*AdvDuffEQN(0.5*(A1+A2),wext)<0 )
                {
                    A2 = 0.5*(A1 + A2);
                }
                else
                {
                    A1 = 0.5*(A1 + A2);
                }
            }
            solutions[n_of_solutions] = 0.5*(A1 + A2);
            n_of_solutions++;
        }
    }
}

double rot_optim::Solve(double wext)
{
    return( Inert*acc_scale*aext/sqrt(pow(spring_k-Inert*pow(wext,2),2)+pow(wext*air_damp,2)));
}

void rot_optim::Resonance_norm(double w_min, double w_max, double dw)
{
    ofstream ofs ("Norm_Resonance_aext_"+to_string(aext)+".dat");

    ofs<<"Aext = "<<to_string(aext)<<endl;
    ofs<<"AirD = "<<to_string(air_damp)<<endl;
    ofs<<"K = "<<to_string(spring_k)<<endl;
    ofs<<"k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"k5 = "<<to_string(spring_k5)<<endl;

    double* sol = new double[30];

    double freq[100000];
    double vrms[100000];
    double vrmsmax = 0.0;

    int imax = 0;

    for(double w= w_min; w<w_max; w+=dw)
    {
        int n_sol = 0;
        Solve(sol,n_sol,2*M_PI*w);
        for(int i=0; i<n_sol; i++)
        {
            if(vrmsmax<vrms_interp_rotate(1.0,sol[i],2*M_PI*w*sol[i]))
            {
                vrmsmax = vrms_interp_rotate(1.0,sol[i],2*M_PI*w*sol[i]);
            }
            freq[imax] = w;
            vrms[imax] = vrms_interp_rotate(1.0,sol[i],2*M_PI*w*sol[i]);
            imax++;
        }
    }

    for(int i=0;i<imax;i++)
    {
        ofs<<freq[i]<<"\t"<<vrms[i]/vrmsmax<<endl;
    }
    ofs.close();
}

void rot_optim::Resonance(double w_min, double w_max, double dw)
{
    ofstream ofs ("Resonance_rot_aext"+to_string(aext)+".dat");
    ofs<<"##Aext = "<<to_string(aext)<<endl;
    ofs<<"##AirD = "<<to_string(air_damp)<<endl;
    ofs<<"##K = "<<to_string(spring_k)<<endl;
    ofs<<"##k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"##k5 = "<<to_string(spring_k5)<<endl;
    ofs<<"##Scale = "<<to_string(scale)<<endl;

    double* sol = new double[30];

    double freq[10000];
    double vrms[10000];
    int i_max = 0;

    for(double w= w_min; w<w_max; w+=dw)
    {
        int n_sol = 0;
        Solve(sol,n_sol,2*M_PI*w);
        for(int i=0; i<n_sol; i++)
        {
            freq[i_max] = w;
            vrms[i_max] = vrms_interp_rotate(scale,sol[i],2*M_PI*w*sol[i]);
            i_max++;
        }
    }

    for(int i=0;i<i_max;i++)
    {
        ofs<<freq[i]<<"\t"<<vrms[i]<<endl;
    }

    ofs.close();
}
