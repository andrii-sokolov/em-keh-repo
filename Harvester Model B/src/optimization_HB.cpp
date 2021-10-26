#include "optim.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

void optimization::read_experiment()
{
    /** Reading the experimental data */
    std::ifstream infile(experimental_filename);    //trying to open file
    if(infile)
    {
        std::string file_f,file_v;                  //variables containing variables
        int i = 0;                                  //counter
        while (infile>>file_f>>file_v)
        {
            experimental_freq[i] = stod (file_f);   //reading frequency
            experimental_volt[i] = stod (file_v);   //reading RMS voltage
            i++;
        }
        imax = i;
    }
    else
    {
        std::cout<<"File \""<<experimental_filename<<"\" NOT FOUND!"<<std::endl;
    }
}

double optimization::AdvDuffEQN (   double A0,              // amplitude            [m]
                                    double wext             // external frequency   [Hz]
                                )
{
    double func = 0.0;
    if (tran)   //The semi-analytical equation with transducer force
        func = pow((-mass*A0*wext*wext + spring_k*A0 + 0.75*pow(A0,3)*spring_k3 + 0.625*pow(A0,5)*spring_k5 + mass*at0(A0,wext)),2) + pow(air_damp*A0*wext + mass*bt0(A0,wext),2) - pow(mass*aext,2);
    else        //The semi-analytical equation without transducer force
        func = pow((-mass*A0*wext*wext + spring_k*A0 + 0.75*pow(A0,3)*spring_k3 + 0.625*pow(A0,5)*spring_k5),2) + pow(air_damp*A0*wext,2) - pow(mass*aext,2);
    return( func );
}

double optimization::at0    (   double A0,              // amplitude            [m]
                                double wext             // external frequency   [Hz]
                            )
{
    return(0.0);
}

double optimization::at0_num    (   double A0,              // amplitude            [m]
                                    double wext             // external frequency   [Hz]
                                )
/* DEBUGGING function*/
{
    int n_of_integration_p = 100;
    double dt = 2*M_PI/wext/n_of_integration_p;
    double Sum = 0.0;
    for(double t = 0; t<2*M_PI/wext; t+= dt)
    {
        double x = A0*cos(wext*t);
        double dx = -A0*wext*sin(wext*t);
        Sum += wext/M_PI*force_interp_progress(dx,x)*cos(wext*t)*dt;
    }
    return(Sum);
}

double optimization::bt0    (   double A0,          // amplitude            [m]
                                double wext         // external frequency   [Hz]
                            )
{
    return( -(A0*wext*(262144*pow(emf00,2) + 4199*pow(A0,20)*pow(emf10,2) +
        9724*pow(A0,18)*emf10*(emf08 + 95*emf10*pow(zshift,2)) +
        262144*pow(zshift,4)*pow(emf02 + emf04*pow(zshift,2) +
           emf06*pow(zshift,4) + emf08*pow(zshift,6) + emf10*pow(zshift,8),
          2) + 5720*pow(A0,16)*
         (pow(emf08,2) + 306*emf08*emf10*pow(zshift,2) +
           emf10*(2*emf06 + 4845*emf10*pow(zshift,4))) +
        13728*pow(A0,14)*(emf04*emf10 +
           emf06*(emf08 + 120*emf10*pow(zshift,2)) +
           60*pow(zshift,2)*(pow(emf08,2) + 51*emf08*emf10*pow(zshift,2) +
              323*pow(emf10,2)*pow(zshift,4))) +
        14336*pow(A0,8)*(pow(emf04,2) +
           2*emf02*(emf06 + 45*emf08*pow(zshift,2) +
              495*emf10*pow(zshift,4)) +
           emf04*(90*emf06*pow(zshift,2) + 990*emf08*pow(zshift,4) +
              6006*emf10*pow(zshift,6)) +
           3*pow(zshift,4)*(165*pow(emf06,2) +
              2002*emf06*emf08*pow(zshift,2) +
              4290*pow(emf08,2)*pow(zshift,4) +
              8580*emf06*emf10*pow(zshift,4) +
              29172*emf08*emf10*pow(zshift,6) +
              41990*pow(emf10,2)*pow(zshift,8))) +
        1024*emf00*(21*pow(A0,10)*emf10 +
           28*pow(A0,8)*(emf08 + 45*emf10*pow(zshift,2)) +
           40*pow(A0,6)*(emf06 + 28*emf08*pow(zshift,2) +
              210*emf10*pow(zshift,4)) +
           512*pow(zshift,2)*(emf02 + emf04*pow(zshift,2) +
              emf06*pow(zshift,4) + emf08*pow(zshift,6) +
              emf10*pow(zshift,8)) +
           128*pow(A0,2)*(emf02 + 6*emf04*pow(zshift,2) +
              15*emf06*pow(zshift,4) + 28*emf08*pow(zshift,6) +
              45*emf10*pow(zshift,8)) +
           64*pow(A0,4)*(emf04 + 15*emf06*pow(zshift,2) +
              70*pow(zshift,4)*(emf08 + 3*emf10*pow(zshift,2)))) +
        8448*pow(A0,12)*(pow(emf06,2) +
           182*emf06*pow(zshift,2)*(emf08 + 20*emf10*pow(zshift,2)) +
           2*(emf02*emf10 + 910*pow(emf08,2)*pow(zshift,4) +
              18564*emf08*emf10*pow(zshift,6) +
              62985*pow(emf10,2)*pow(zshift,8) +
              emf04*(emf08 + 91*emf10*pow(zshift,2)))) +
        21504*pow(A0,10)*(emf02*(emf08 + 66*emf10*pow(zshift,2)) +
           emf04*(emf06 + 66*emf08*pow(zshift,2) +
              1001*emf10*pow(zshift,4)) +
           11*pow(zshift,2)*(3*pow(emf06,2) +
              91*emf06*pow(zshift,2)*(emf08 + 8*emf10*pow(zshift,2)) +
              26*pow(zshift,4)*
               (14*pow(emf08,2) + 153*emf08*emf10*pow(zshift,2) +
                 323*pow(emf10,2)*pow(zshift,4)))) +
        131072*pow(A0,2)*pow(zshift,2)*
         (3*pow(emf02,2) + emf02*
            (15*emf04*pow(zshift,2) + 28*emf06*pow(zshift,4) +
              45*emf08*pow(zshift,6) + 66*emf10*pow(zshift,8)) +
           pow(zshift,4)*(14*pow(emf04,2) +
              emf04*(45*emf06*pow(zshift,2) + 66*emf08*pow(zshift,4) +
                 91*emf10*pow(zshift,6)) +
              pow(zshift,4)*(33*pow(emf06,2) +
                 91*emf06*emf08*pow(zshift,2) +
                 60*pow(emf08,2)*pow(zshift,4) +
                 120*emf06*emf10*pow(zshift,4) +
                 153*emf08*emf10*pow(zshift,6) +
                 95*pow(emf10,2)*pow(zshift,8)))) +
        32768*pow(A0,4)*(pow(emf02,2) +
           10*emf02*(3*emf04*pow(zshift,2) + 14*emf06*pow(zshift,4) +
              42*emf08*pow(zshift,6) + 99*emf10*pow(zshift,8)) +
           pow(zshift,4)*(70*pow(emf04,2) +
              emf04*(420*emf06*pow(zshift,2) + 990*emf08*pow(zshift,4) +
                 2002*emf10*pow(zshift,6)) +
              pow(zshift,4)*(495*pow(emf06,2) +
                 2002*emf06*emf08*pow(zshift,2) +
                 1820*pow(emf08,2)*pow(zshift,4) +
                 3640*emf06*emf10*pow(zshift,4) +
                 6120*emf08*emf10*pow(zshift,6) +
                 4845*pow(emf10,2)*pow(zshift,8)))) +
        40960*pow(A0,6)*(emf02*
            (emf04 + 28*emf06*pow(zshift,2) + 210*emf08*pow(zshift,4) +
              924*emf10*pow(zshift,6)) +
           pow(zshift,2)*(14*pow(emf04,2) +
              21*emf04*(10*emf06*pow(zshift,2) + 44*emf08*pow(zshift,4) +
                 143*emf10*pow(zshift,6)) +
              pow(zshift,4)*(462*pow(emf06,2) +
                 3003*emf06*emf08*pow(zshift,2) +
                 4004*pow(emf08,2)*pow(zshift,4) +
                 8008*emf06*emf10*pow(zshift,4) +
                 18564*emf08*emf10*pow(zshift,6) +
                 19380*pow(emf10,2)*pow(zshift,8))))))/(262144.*resist))*scale*scale;
}

double optimization::bt0_num    (   double A0,          // amplitude            [m]
                                    double wext         // external frequency   [Hz]
                                )
/* DEBUGGING function*/
{
    int n_of_integration_p = 100;
    double dt = 2*M_PI/wext/n_of_integration_p;
    double Sum = 0.0;
    for(double t = 0; t<2*M_PI/wext; t+= dt)
    {
        double x = A0*cos(wext*t);
        double dx = -A0*wext*sin(wext*t);
        Sum += wext/M_PI*force_interp_progress(dx,x)*sin(wext*t)*dt;
    }
    return(Sum);
}

void optimization::Solve    (   double* solutions,              //pointer for the solutions
                                int& n_of_solutions,            //number of definite solutions
                                double wext                     //external circular frequency
                            )
{
    double A1 = 0.0;
    double A2 = A_max;
    while( 2*abs(A2 - A1)/abs(A1 + A2)>0.05 )                //second - fine search
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
    solutions[0] = 0.5*(A1+A2);
    n_of_solutions = 1;
}

/*{
    double A1, A2;                                      //variables for the points to the left and to the right from solutions
    n_of_solutions = 0;                                 //initialization of the number of defined solutions
    for(double A = 0.0; A<=A_max; A+= A_max/n_of_non_lin)
    {
        if( AdvDuffEQN(A,wext)*AdvDuffEQN(A+A_max/n_of_non_lin,wext)<0) //first - linear search
        {


            solutions[n_of_solutions] = 0.5*(A1 + A2);
            n_of_solutions++;
        }
    }
}*/

double optimization::Solve_Linear   (   double fext         //external frequency
                                    )
{
    double wext=2.0*M_PI*fext;
    return(mass*aext/sqrt(pow(spring_k-mass*pow(wext,2),2)+pow(wext*air_damp,2)));
}

void optimization::ResonanceHB  (   double f_min,           //minimal frequency
                                    double f_max,           //maximal frequency
                                    double df               //frequency step
                                )
{
    std::ofstream ofs ("ResonanceHB"+experimental_filename);
    ofs<<"#Aext = "<<std::to_string(aext)<<std::endl;
    ofs<<"#AirD = "<<std::to_string(air_damp)<<std::endl;
    ofs<<"#K = "<<std::to_string(spring_k)<<std::endl;
    ofs<<"#k3 = "<<std::to_string(spring_k3)<<std::endl;
    ofs<<"#k5 = "<<std::to_string(spring_k5)<<std::endl;
    ofs<<"#Scale = "<<std::to_string(scale)<<std::endl;

    double* sol = new double[30];
    int n_sol = 0;

    for(double f = f_min; f<f_max; f+=df)
    {
        n_sol = 0;
        Solve(sol,    n_sol,  2*M_PI*f);
        for(int i=0; i<n_sol; i++)
        {
            //std::cout<<f<<"\t"<<sol[i]<<"\t"<<vrms_interp_progress(sol[i],2.0*M_PI*f)<<std::endl;
            ofs<<f<<"\t"<<sol[i]<<"\t"<<vrms_interp_progress(sol[i],2.0*M_PI*f)<<std::endl;
        }
    }

    ofs.close();
}

void optimization::ResonanceHB()
{
    std::ofstream ofs ("ResonanceHB"+experimental_filename);
    ofs<<"#Aext = "<<std::to_string(aext)<<std::endl;
    ofs<<"#AirD = "<<std::to_string(air_damp)<<std::endl;
    ofs<<"#K = "<<std::to_string(spring_k)<<std::endl;
    ofs<<"#k3 = "<<std::to_string(spring_k3)<<std::endl;
    ofs<<"#k5 = "<<std::to_string(spring_k5)<<std::endl;
    ofs<<"#Scale = "<<std::to_string(scale)<<std::endl;

    double* sol = new double[30];
    int n_sol = 0;

    for(int i=0; i<imax; i++)
    {
        n_sol = 0;
        Solve(sol,    n_sol,  2*M_PI*experimental_freq[i]);
        for(int j=0; j<n_sol; j++)
        {
            ofs<<experimental_freq[i]<<"\t"<<sol[j]<<"\t"<<vrms_interp_progress(sol[j],2.0*M_PI*experimental_freq[i])<<std::endl;
            //std::cout<<experimental_freq[i]<<"\t"<<sol[j]<<"\t"<<vrms_interp_progress(sol[j],2.0*M_PI*experimental_freq[i])<<std::endl;
        }
    }

    ofs.close();
}

void optimization::ResonanceHBLinear()
{
    std::ofstream ofs ("ResonanceHB"+experimental_filename);

    double Ampl;

    for(int i=0; i<imax; i++)
    {
        Ampl = Solve_Linear(experimental_freq[i]);
        ofs<<experimental_freq[i]<<"\t"<<Ampl<<"\t"<<vrms_interp_progress(Ampl,2.0*M_PI*experimental_freq[i])<<std::endl;
            //std::cout<<experimental_freq[i]<<"\t"<<sol[j]<<"\t"<<vrms_interp_progress(sol[j],2.0*M_PI*experimental_freq[i])<<std::endl;
    }

    ofs.close();
}
