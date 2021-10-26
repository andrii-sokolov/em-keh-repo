#include "transp.h"
#include <math.h>
#include <fstream>

double progress_optim::acceleration(double x/** displacement */, double v/** velocity */, double t/** time */, double wext/** external acceleration frequency */)
{
    return(aext*cos(wext*t) - air_damp/mass*v - spring_k/mass*x - spring_k3/mass*pow(x,3) - spring_k5/mass*pow(x,5));
}

double progress_optim::velocity(double x/** displacement */, double v/** velocity */, double t/** time */, double wext/** external acceleration frequency */)
{
    return(v);
}

void progress_optim::RK4(double fext, double& ampl_x, double& ampl_dx)
{
    double wext = 2*M_PI*fext;
    double x = 0.0;
    double x1,x2;
    double dx = 0.0;
    double t = 0.0;
    double dt = delta/fext;
    double p1,p2,p3,p4,l1,l2,l3,l4;

    double integr_ch = 0.0;
    double integr_ne = -1.0;
    double integr = 0.0;

    double ampl_dx_min = 0.0;
    double ampl_dx_max = 0.0;
    ampl_dx = 0.0;

    double ampl_x_min = 0.0;
    double ampl_x_max = 0.0;
    ampl_x = 0.0;

    int counter = 0;
    int counter_period = 0;
    int counter_1 = 0;

    while((fabs(1-fabs(integr_ch)/fabs(integr_ne))>st_st_delta)&(counter<Max_Counter))
    {
        x1 = x;

        p1 = dt*acceleration(x,dx,t,wext);
        l1 = dt*velocity(x,dx,t,wext);

        p2 = dt*acceleration(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        l2 = dt*velocity(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);

        p3 = dt*acceleration(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        l3 = dt*velocity(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);

        p4 = dt*acceleration(x + l3, dx + p3, t + dt, wext);
        l4 = dt*velocity(x + l3, dx + p3, t + dt, wext);

        t+=dt;
        x+=(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0;
        dx+=(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;

        x2 = x;

        counter_1++;

        integr+=fabs(dx)*fabs(x2-x1);

        if(x>ampl_x_max)
            ampl_x_max = x;
        if(x<ampl_x_min)
            ampl_x_min = x;

        if(dx>ampl_dx_max)
            ampl_dx_max = dx;
        if(dx<ampl_dx_min)
            ampl_dx_min = dx;

        if(x1*x2 <=0)
            counter_period ++;

        if((x1*x2 <= 0)&(counter_period%5 == 0))
        {
            if(counter%2 == 0)
                integr_ch = integr;
            else
                integr_ne = integr;
            counter++;
            integr = 0.0;

            ampl_x = (ampl_x_max - ampl_x_min)/2.0;
            ampl_x_max = 0.0;
            ampl_x_min = 0.0;

            ampl_dx = (ampl_dx_max - ampl_dx_min)/2.0;
            ampl_dx_max = 0.0;
            ampl_dx_min = 0.0;
        }

    }
}

void progress_optim::RK4(double fext, double& ampl_x, double& ampl_dx, string filename)
{
    double wext = 2*M_PI*fext;
    double x = 0.0;
    double x1,x2;
    double dx = 0.0;
    double t = 0.0;
    double dt = delta/fext;
    double p1,p2,p3,p4,l1,l2,l3,l4;

    double integr_ch = 0.0;
    double integr_ne = -1.0;
    double integr = 0.0;

    double ampl_dx_min = 0.0;
    double ampl_dx_max = 0.0;
    ampl_dx = 0.0;

    double ampl_x_min = 0.0;
    double ampl_x_max = 0.0;
    ampl_x = 0.0;

    int counter = 0;
    int counter_period = 0;
    int counter_1 = 0;

    ofstream ofs (filename);

    while((fabs(1-fabs(integr_ch)/fabs(integr_ne))>st_st_delta)&(counter<Max_Counter))
    {
        x1 = x;

        p1 = dt*acceleration(x,dx,t,wext);
        l1 = dt*velocity(x,dx,t,wext);

        p2 = dt*acceleration(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);
        l2 = dt*velocity(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, wext);

        p3 = dt*acceleration(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);
        l3 = dt*velocity(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, wext);

        p4 = dt*acceleration(x + l3, dx + p3, t + dt, wext);
        l4 = dt*velocity(x + l3, dx + p3, t + dt, wext);

        t+=dt;
        x+=(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0;
        dx+=(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;

        x2 = x;

        counter_1++;

        if(counter_1%1000 == 0)
            ofs<<t<<"\t"<<x<<"\t"<<Solve_Linear(fext)*cos(2*M_PI*fext*t)<<endl;

        integr+=fabs(dx)*fabs(x2-x1);

        if(x>ampl_x_max)
            ampl_x_max = x;
        if(x<ampl_x_min)
            ampl_x_min = x;

        if(dx>ampl_dx_max)
            ampl_dx_max = dx;
        if(dx<ampl_dx_min)
            ampl_dx_min = dx;

        if(x1*x2 <=0)
            counter_period ++;

        if((x1*x2 <= 0)&(counter_period%5 == 0))
        {
            if(counter%2 == 0)
                integr_ch = integr;
            else
                integr_ne = integr;
            counter++;
            integr = 0.0;

            ampl_x = (ampl_x_max - ampl_x_min)/2.0;
            ampl_x_max = 0.0;
            ampl_x_min = 0.0;

            ampl_dx = (ampl_dx_max - ampl_dx_min)/2.0;
            ampl_dx_max = 0.0;
            ampl_dx_min = 0.0;
        }

    }
    ofs.close();
}

void progress_optim::ResonanceRK(double w_min, double w_max, double dw)
{
    ofstream ofs ("ResonanceRK_aext"+to_string(aext)+".dat");
    ofs<<"#Aext = "<<to_string(aext)<<endl;
    ofs<<"#AirD = "<<to_string(air_damp)<<endl;
    ofs<<"#K = "<<to_string(spring_k)<<endl;
    ofs<<"#k3 = "<<to_string(spring_k3)<<endl;
    ofs<<"#k5 = "<<to_string(spring_k5)<<endl;
    ofs<<"#Scale = "<<to_string(scale)<<endl;

    double Ampl_x;
    double Ampl_v;

    double freq[100000];
    double vrms[100000];
    int i_max = 0;

    for(double w= w_min; w<w_max; w+=dw)
    {
        RK4(w,Ampl_x,Ampl_v);
        freq[i_max] = w;
        vrms[i_max] = vrms_interp_progress(scale,Ampl_x,Ampl_v);
        i_max++;
    }

    for(int i=0; i<i_max; i++)
    {
        ofs<<freq[i]<<"\t"<<vrms[i]<<endl;
    }
    ofs.close();
}
