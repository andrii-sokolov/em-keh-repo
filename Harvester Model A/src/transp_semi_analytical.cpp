#include "transp.h"
#include <math.h>
#include <fstream>

double progress_optim::AdvDuffEQN(double A0, double wext)
{
    return(pow((-mass*A0*wext*wext + spring_k*A0 + 0.75*pow(A0,3)*spring_k3 + 0.625*pow(A0,5)*spring_k5 + 0.546875*pow(A0,7)*spring_k7),2) + pow(air_damp*A0*wext,2) - pow(mass*aext,2));
}

void progress_optim::Solve(double* solutions, int& n_of_solutions, double wext)
{
    double A1, A2;
    n_of_solutions = 0;
    double dA = 5e-7;
    for(double A = 0.0; A<=8e-5; A+= dA)
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

double progress_optim::Solve_Linear(double fext)
{
    double wext=2.0*M_PI*fext;
    return(mass*aext/sqrt(pow(spring_k-mass*pow(wext,2),2)+pow(wext*air_damp,2)));
}


void progress_optim::Resonance_norm(double w_min, double w_max, double dw)
{
    ofstream ofs ("Norm_Resonance_aext_"+to_string(aext)+".dat");

    ofs<<"#Aext = "<<to_string(aext)<<endl;
    ofs<<"#AirD = "<<to_string(air_damp)<<endl;
    ofs<<"#K = "<<to_string(spring_k)<<endl;
    ofs<<"#k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"#k5 = "<<to_string(spring_k5)<<endl;

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
            if(vrmsmax<vrms_interp_progress(1.0,sol[i],2*M_PI*w*sol[i]))
            {
                vrmsmax = vrms_interp_progress(1.0,sol[i],2*M_PI*w*sol[i]);
            }
            freq[imax] = w;
            vrms[imax] = vrms_interp_progress(1.0,sol[i],2*M_PI*w*sol[i]);
            imax++;
        }
    }

    for(int i=0; i<imax; i++)
    {
        ofs<<freq[i]<<"\t"<<vrms[i]/vrmsmax<<endl;
    }
    ofs.close();
}

void progress_optim::Resonance(double w_min, double w_max, double dw)
{
    ofstream ofs ("Resonance_aext"+to_string(aext)+".dat");
    ofs<<"#Aext = "<<to_string(aext)<<endl;
    ofs<<"#AirD = "<<to_string(air_damp)<<endl;
    ofs<<"#K = "<<to_string(spring_k)<<endl;
    ofs<<"#k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"#k5 = "<<to_string(spring_k5)<<endl;
    ofs<<"#Scale = "<<to_string(scale)<<endl;

    double* sol = new double[30];

    bool start = true;

    double freq_start[100000];
    double vrms_start[100000];
    double freq_middle[100000];
    double vrms_middle[100000];
    double freq_min[100000];
    double vrms_min[100000];
    int imax_start = 0;
    int imax_min = 0;
    int imax_middle = 0;

    for(double w= w_min; w<w_max; w+=dw)
    {
        int n_sol = 0;
        Solve(sol,n_sol,2*M_PI*w);
        if ((n_sol == 1) & (start))
        {
            freq_start[imax_start] = w;
            vrms_start[imax_start] = vrms_interp_progress(scale,sol[0],2*M_PI*w*sol[0]);
            imax_start++;
        }
        else
        {
            start = false;
            if( n_sol == 1 )
            {
                freq_min[imax_min] = w;
                vrms_min[imax_min] = vrms_interp_progress(scale,sol[0],2*M_PI*w*sol[0]);
                imax_min++;
            }
            else
            {
                freq_min[imax_min] = w;
                vrms_min[imax_min] = vrms_interp_progress(scale,sol[0],2*M_PI*w*sol[0]);
                imax_min++;
                freq_middle[imax_middle] = w;
                vrms_middle[imax_middle] = vrms_interp_progress(scale,sol[1],2*M_PI*w*sol[1]);
                imax_middle++;
                freq_start[imax_start] = w;
                vrms_start[imax_start] = vrms_interp_progress(scale,sol[2],2*M_PI*w*sol[2]);
                imax_start++;
            }
        }
    }

    for(int i=0; i<imax_start; i++)
    {
        ofs<<freq_start[i]<<"\t"<<vrms_start[i]<<endl;
    }

    for(int i=imax_middle-1; i>=0; i--)
    {
        ofs<<freq_middle[i]<<"\t"<<vrms_middle[i]<<endl;
    }

    for(int i=0; i<imax_min; i++)
    {
        ofs<<freq_min[i]<<"\t"<<vrms_min[i]<<endl;
    }

    ofs.close();
}

void progress_optim::ResonanceLin(double w_min, double w_max, double dw)
{
    ofstream ofs ("Resonance_Lin_aext"+to_string(aext)+".dat");
    ofs<<"#Aext = "<<to_string(aext)<<endl;
    ofs<<"#AirD = "<<to_string(air_damp)<<endl;
    ofs<<"#K = "<<to_string(spring_k)<<endl;
    ofs<<"#k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"#k5 = "<<to_string(spring_k5)<<endl;
    ofs<<"#Scale = "<<to_string(scale)<<endl;

    for(double w= w_min; w<w_max; w+=dw)
    {
        ofs<<w<<"\t"<<vrms_interp_progress(scale,Solve_Linear(w),2*M_PI*w*Solve_Linear(w))<<endl;
    }

    ofs.close();
}
