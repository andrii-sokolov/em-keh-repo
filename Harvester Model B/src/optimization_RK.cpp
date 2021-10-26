#include "optim.h"                      //Optimization module

#include <math.h>                       //Mathematical module
#include <fstream>                      //Writing to file module
#include <string>                       //Working with strings module
#include <iostream>                     //Input-output module

double optimization::acceleration   (   double x,       // displacement
                                        double v,       // velocity
                                        double t,       // time
                                        double wext     // external acceleration frequency
                                    )
{
    double acc;
    if( tran )
    {
        acc = aext*cos(wext*t)              //external force
            - air_damp/mass*v               //air damping force
            - spring_k/mass*x               //linear spring force
            - spring_k3/mass*pow(x,3)       //3-rd order spring force
            - spring_k5/mass*pow(x,5)       //5-th order spring force
            + force_interp_progress(v,x);   //transducer force
    }
    else
    {
        acc = aext*cos(wext*t)          //external force
            - air_damp/mass*v           //air damping force
            - spring_k/mass*x           //linear spring force
            - spring_k3/mass*pow(x,3)   //3-rd order spring force
            - spring_k5/mass*pow(x,5);  //5-th order spring force
    }
    return  (acc);
}

double optimization::velocity   (       double x,       // displacement
                                        double v,       // velocity
                                        double t,       // time
                                        double wext     // external acceleration frequency
                                )
{
    return(v);
}

void optimization::RK4  (   double fext,            // external frequency [Hz]
                            double& ampl_x,         // amplitude displacement result [m]
                            double& ampl_dx,        // amplitude velocity result [m/s]
                            double& rms_volt        // RMS voltage
                        )
{
    double wext = 2*M_PI*fext;              // setting up the external circular frequency
    double x = 0.0;                         // initial value of coordinate
    double dx = 0.0;                        // initial value of velocity
    double t = 0.0;                         // initial value of time
    double dt = delta/fext;                 // the time step
    double p1,p2,p3,p4,l1,l2,l3,l4 = 0.0;   // RK coefficients
    double x1,x2 = 0.0;                     // two values of coordinate to understand if there were half-period (if x1*x2<=0 then another half-period gone)

    double integr_ch = 0.0;                 // the even period path length
    double integr_ne = -1.0;                // the odd period path length
    double integr = 0.0;                    // the integration variable

    double rms_int = 0.0;                   // RMS integral
    double rms_time = 0.0;                  // time period

    double ampl_dx_min = 0.0;               // minimal value of the coordinate
    double ampl_dx_max = 0.0;               // maximal value of the coordinate
    ampl_dx = 0.0;                          // coordinate amplitude

    double ampl_x_min = 0.0;                // minimal value of the velocity
    double ampl_x_max = 0.0;                // maximal value of the velocity
    ampl_x = 0.0;                           // velocity amplitude

    int counter = 0;                        // counter of half-periods
    int counter_period = 0;                 // counter of periods
    int counter_tacts = 0;                  // counter of tacts

    while((fabs(1-fabs(integr_ch)/fabs(integr_ne))>st_st_delta)&(counter<Max_Counter))
    {
        x1 = x;     //seting up the "previous" coordinate

        /* RK coefficients calculating */
        p1 = dt*acceleration(x,dx,t,wext);
        l1 = dt*velocity(x,dx,t,wext);
        p2 = dt*acceleration(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        l2 = dt*velocity(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        p3 = dt*acceleration(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        l3 = dt*velocity(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        p4 = dt*acceleration(x + l3, dx + p3, t + dt, wext);
        l4 = dt*velocity(x + l3, dx + p3, t + dt, wext);

        t+=dt;                                  //iterating time
        x+=(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0;     //calculating coordinate
        dx+=(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;    //calculating velocity

        x2 = x;     //setting up the "next" coordinate

        counter_tacts++;    //iterating the global counter

        rms_int +=pow(emf_interp_progress(dx, x),2)*dt;
        rms_time+=dt;

        integr+=fabs(dx)*fabs(x2-x1);   //path integral calculating

        if(x>ampl_x_max)                //definition of the maximal coordinate
            ampl_x_max = x;
        if(x<ampl_x_min)                //definition of tht minimal coordinate
            ampl_x_min = x;

        if(dx>ampl_dx_max)              //definition of the maximal velocity
            ampl_dx_max = dx;
        if(dx<ampl_dx_min)              //definition of the minimal velocity
            ampl_dx_min = dx;

        if(x1*x2 <=0)                   //definition of the half periods iteration
            counter_period ++;

        if((x1*x2 <= 0)&(counter_period%5 == 0))    //if the number of half periods is divisible by 5
        {
            if(counter%2 == 0)              //for even periods
                integr_ch = integr;
            else                            //for odd periods
                integr_ne = integr;
            counter++;                      //periods counter
            integr = 0.0;                   //clearing of the integration variable

            ampl_x = (ampl_x_max - ampl_x_min)/2.0;     //setting the amplitude
            ampl_x_max = 0.0;   //clearification of the variables
            ampl_x_min = 0.0;

            ampl_dx = (ampl_dx_max - ampl_dx_min)/2.0;  //setting the velocity amplitude
            ampl_dx_max = 0.0;  //clearification of the variables
            ampl_dx_min = 0.0;

            rms_volt = sqrt(rms_int/rms_time);
            rms_time = 0.0;
            rms_int = 0.0;

        }
    }
}

void optimization::RK4  (   double fext,
                            double& ampl_x,
                            double& ampl_dx,
                            std::string filename
                        )
{
    double wext = 2*M_PI*fext;              // setting up the external circular frequency
    double x = 0.0;                         // initial value of coordinate
    double dx = 0.0;                        // initial value of velocity
    double t = 0.0;                         // initial value of time
    double dt = delta/fext;                 // the time step
    double p1,p2,p3,p4,l1,l2,l3,l4 = 0.0;   // RK coefficients
    double x1,x2 = 0.0;                     // two values of coordinate to understand if there were half-period (if x1*x2<=0 then another half-period gone)

    double integr_ch = 0.0;                 // the even period path length
    double integr_ne = -1.0;                // the odd period path length
    double integr = 0.0;                    // the integration variable

    double ampl_dx_min = 0.0;               // minimal value of the coordinate
    double ampl_dx_max = 0.0;               // maximal value of the coordinate
    ampl_dx = 0.0;                          // coordinate amplitude

    double ampl_x_min = 0.0;                // minimal value of the velocity
    double ampl_x_max = 0.0;                // maximal value of the velocity
    ampl_x = 0.0;                           // velocity amplitude

    int counter = 0;                        // counter of half-periods
    int counter_period = 0;                 // counter of periods
    int counter_tacts = 0;                  // counter of tacts

    std::ofstream ofs (filename);           //opening of the waveform file

    while((fabs(1-fabs(integr_ch)/fabs(integr_ne))>st_st_delta)&(counter<Max_Counter))
    {
        x1 = x;     //seting up the "previous" coordinate

        /* RK coefficients calculating */
        p1 = dt*acceleration(x,dx,t,wext);
        l1 = dt*velocity(x,dx,t,wext);
        p2 = dt*acceleration(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        l2 = dt*velocity(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        p3 = dt*acceleration(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        l3 = dt*velocity(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        p4 = dt*acceleration(x + l3, dx + p3, t + dt, wext);
        l4 = dt*velocity(x + l3, dx + p3, t + dt, wext);

        t+=dt;                                  //iterating time
        x+=(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0;     //calculating coordinate
        dx+=(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;    //calculating velocity

        x2 = x;     //setting up the "next" coordinate

        counter_tacts++;    //iterating the global counter

        integr+=fabs(dx)*fabs(x2-x1);   //path integral calculating

        if(counter_tacts%1000 == 0)
            ofs<<t<<"\t"<<x<<std::endl;     //selective writing in file

        if(x>ampl_x_max)                //definition of the maximal coordinate
            ampl_x_max = x;
        if(x<ampl_x_min)                //definition of tht minimal coordinate
            ampl_x_min = x;

        if(dx>ampl_dx_max)              //definition of the maximal velocity
            ampl_dx_max = dx;
        if(dx<ampl_dx_min)              //definition of the minimal velocity
            ampl_dx_min = dx;

        if(x1*x2 <=0)                   //definition of the half periods iteration
            counter_period ++;

        if((x1*x2 <= 0)&(counter_period%5 == 0))    //if the number of half periods is divisible by 5
        {
            if(counter%2 == 0)              //for even periods
                integr_ch = integr;
            else                            //for odd periods
                integr_ne = integr;
            counter++;                      //periods counter
            integr = 0.0;                   //clearing of the integration variable

            ampl_x = (ampl_x_max - ampl_x_min)/2.0;     //setting the amplitude
            ampl_x_max = 0.0;   //clearification of the variables
            ampl_x_min = 0.0;

            ampl_dx = (ampl_dx_max - ampl_dx_min)/2.0;  //setting the velocity amplitude
            ampl_dx_max = 0.0;  //clearification of the variables
            ampl_dx_min = 0.0;
        }
    }
    ofs.close();
}

void optimization::ResonanceRK  (   double f_min,   // start frequency [Hz]
                                    double f_max,   // stop frequency [Hz]
                                    double df       // step by frequency [Hz]
                                )
{
    std::ofstream ofs ("ResonanceRK"+experimental_filename);
    ofs<<"#Aext = "<<std::to_string(aext)<<std::endl;
    ofs<<"#AirD = "<<std::to_string(air_damp)<<std::endl;
    ofs<<"#K = "<<std::to_string(spring_k)<<std::endl;
    ofs<<"#k3 = "<<std::to_string(spring_k3)<<std::endl;
    ofs<<"#k5 = "<<std::to_string(spring_k5)<<std::endl;
    ofs<<"#Scale = "<<std::to_string(scale)<<std::endl;

    double Ampl_x;
    double Ampl_v;
    double RMS_volt;

    for(double f= f_min; f<f_max; f+=df)
    {
        RK4(f,Ampl_x,Ampl_v,RMS_volt);
        ofs<<f<<"\t"<<Ampl_x<<"\t"<<RMS_volt<<std::endl;
        std::cout<<f<<"\t"<<Ampl_x<<"\t"<<RMS_volt<<std::endl;
    }

    ofs.close();
}

void optimization::ResonanceRK  ()
{
    std::ofstream ofs ("ResonanceRK"+experimental_filename);
    ofs<<"#Aext = "<<std::to_string(aext)<<std::endl;
    ofs<<"#AirD = "<<std::to_string(air_damp)<<std::endl;
    ofs<<"#K = "<<std::to_string(spring_k)<<std::endl;
    ofs<<"#k3 = "<<std::to_string(spring_k3)<<std::endl;
    ofs<<"#k5 = "<<std::to_string(spring_k5)<<std::endl;
    ofs<<"#Scale = "<<std::to_string(scale)<<std::endl;

    double Ampl_x;
    double Ampl_v;
    double RMS_volt;

    for(int i = 0; i<imax; i++)
    {
        RK4(experimental_freq[i],Ampl_x,Ampl_v,RMS_volt);
        ofs<<experimental_freq[i]<<"\t"<<Ampl_x<<"\t"<<RMS_volt<<std::endl;
        std::cout<<experimental_freq[i]<<"\t"<<Ampl_x<<"\t"<<RMS_volt<<std::endl;
    }

    ofs.close();
}
