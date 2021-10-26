#include "transp.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

void progress_optim::coordinate_discent(double dc, double dk, double dk3, double dk5)
{
    double ms_prev = minsqr_norm();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double k3_prev= spring_k3;
    double k5_prev= spring_k5;

    double ms;
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        k3[0] = k3_prev - dk3;
        k5[0] = k5_prev - dk5;
        c0[1] = c_prev;
        k1[1] = k_prev;
        k3[1] = k3_prev;
        k5[1] = k5_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        k3[2] = k3_prev + dk3;
        k5[2] = k5_prev + dk5;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        air_damp = c0[i];
                        spring_k = k1[j];
                        spring_k3= k3[k];
                        spring_k5= k5[l];
                        ms = minsqr_norm();
                        if(ms<ms_next)
                        {
                            ms_next = ms;
                            c_prev= c0[i];
                            k_prev = k1[j];
                            k3_prev= k3[k];
                            k5_prev= k5[l];
                        }
                    }
                }
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }
        cout<<k_prev<<"\t"<<k3_prev<<"\t"<<k5_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

    }
}

void progress_optim::coordinate_discent(double ds)
{
    double ms_prev = minsqr();
    double sc_prev = scale;
    double ms_next = ms_prev;

    double ms;
    double sc[3];


    bool exit = true;

    while(exit)
    {
        sc[0] = sc_prev - ds;
        sc[1] = sc_prev;
        sc[2] = sc_prev + ds;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            scale = sc[i];
            ms = minsqr();
            if(ms == ms_prev)
            {
                exit = false;
            }
            if(ms<ms_next)
            {
                ms_next = ms;
                sc_prev= sc[i];
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }
        cout<<sc_prev<<"\t"<<ms_prev<<endl;
    }
}

void progress_optim::coordinate_discent(double ds, double dc, double dk, double dk3, double dk5)
{
    double ms_prev = minsqr();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double k3_prev= spring_k3;
    double k5_prev= spring_k5;
    double sc_prev= scale;

    double ms=0.0;
    double sc[3];
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        k3[0] = k3_prev - dk3;
        k5[0] = k5_prev - dk5;
        sc[0] = sc_prev - ds;
        c0[1] = c_prev;
        k1[1] = k_prev;
        k3[1] = k3_prev;
        k5[1] = k5_prev;
        sc[1] = sc_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        k3[2] = k3_prev + dk3;
        k5[2] = k5_prev + dk5;
        sc[2] = sc_prev + ds;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        for(int m=0; m<3; m++)
                        {
                            air_damp = c0[i];
                            spring_k = k1[j];
                            spring_k3= k3[k];
                            spring_k5= k5[l];
                            scale = sc[m];
                            ms = minsqr();
                            if(ms<ms_next)
                            {
                                ms_next = ms;
                                c_prev= c0[i];
                                k_prev = k1[j];
                                k3_prev= k3[k];
                                k5_prev= k5[l];
                                sc_prev= sc[m];
                            }
                        }
                    }
                }
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }
        cout<<scale<<"\t"<<k_prev<<"\t"<<k3_prev<<"\t"<<k5_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

    }
}

void progress_optim::coordinate_discent_sets(double ds, double dc, double dk, double dk3, double dk5)
{
    aext = 3.0;
    Resonance(300,600,1);
    aext = 4.0;
    Resonance(300,600,1);
    aext = 5.0;
    Resonance(300,600,1);

    Replot("experiment_3_progres.dat","experiment_4_progres.dat","experiment_5_progres.dat","Resonance_aext3.000000.dat","Resonance_aext4.000000.dat","Resonance_aext5.000000.dat");

    double ms_prev = minsqr_sets();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double k3_prev= spring_k3;
    double k5_prev= spring_k5;
    double sc_prev= scale;

    double ms=0.0;
    double sc[3];
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        k3[0] = k3_prev - dk3;
        k5[0] = k5_prev - dk5;
        sc[0] = sc_prev - ds;
        c0[1] = c_prev;
        k1[1] = k_prev;
        k3[1] = k3_prev;
        k5[1] = k5_prev;
        sc[1] = sc_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        k3[2] = k3_prev + dk3;
        k5[2] = k5_prev + dk5;
        sc[2] = sc_prev + ds;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int l=0; l<3; l++)
                    {
                        for(int m=0; m<3; m++)
                        {
                            if(c0[i]>0)
                                air_damp = c0[i];
                            else
                                dc = dc/10.0;
                            spring_k = k1[j];
                            spring_k3= k3[k];
                            spring_k5= k5[l];
                            scale = sc[m];
                            ms = minsqr_sets();
                            if(ms<ms_next)
                            {
                                ms_next = ms;
                                c_prev= c0[i];
                                k_prev = k1[j];
                                k3_prev= k3[k];
                                k5_prev= k5[l];
                                sc_prev= sc[m];
                            }
                        }
                    }
                }
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }
        cout<<scale<<"\t"<<k_prev<<"\t"<<k3_prev<<"\t"<<k5_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

        aext = 3.0;
        Resonance(300,600,1);
        aext = 4.0;
        Resonance(300,600,1);
        aext = 5.0;
        Resonance(300,600,1);

        Replot("experiment_3_progres.dat","experiment_4_progres.dat","experiment_5_progres.dat","Resonance_aext3.000000.dat","Resonance_aext4.000000.dat","Resonance_aext5.000000.dat");
    }
}

void progress_optim::coordinate_discent_linear(double ds, double dc, double dk)
{
    double ms_prev = minsqr();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double sc_prev= scale;

    spring_k3 = 0.0;
    spring_k5 = 0.0;

    double ms=0.0;
    double sc[3];
    double c0[3];
    double k1[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        sc[0] = sc_prev - ds;
        c0[1] = c_prev;
        k1[1] = k_prev;
        sc[1] = sc_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        sc[2] = sc_prev + ds;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int m=0; m<3; m++)
                {
                    air_damp = c0[i];
                    spring_k = k1[j];
                    scale = sc[m];
                    ms = minsqr();
                    if(ms<ms_next)
                    {
                        ms_next = ms;
                        c_prev= c0[i];
                        k_prev = k1[j];
                        sc_prev= sc[m];
                    }
                }
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }
        cout<<scale<<"\t"<<k_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

    }
}

void progress_optim::coordinate_discent_duffing(double ds, double dc, double dk, double dk3)
{
    spring_k5 = 0.0;


    double ms_prev = minsqr();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double k3_prev= spring_k3;
    double sc_prev= scale;

    Resonance(300,510,2);

    Gnuplot gp;
    gp << "plot \"Resonance_aext3.000000.dat\" title \" Theory \", "<< " \"experiment.dat\" title \"Experiment\" "<<endl;


    double ms=0.0;
    double sc[3];
    double c0[3];
    double k1[3];
    double k3[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        k3[0] = k3_prev - dk3;
        sc[0] = sc_prev - ds;
        c0[1] = c_prev;
        k1[1] = k_prev;
        k3[1] = k3_prev;
        sc[1] = sc_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        k3[2] = k3_prev + dk3;
        sc[2] = sc_prev + ds;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0; j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                    for(int m=0; m<3; m++)
                    {
                        if(c0[i]>0)
                        {
                            air_damp = c0[i];
                        }
                        spring_k = k1[j];
                        spring_k3= k3[k];
                        scale = sc[m];
                        ms = minsqr();
                        if(ms<ms_next)
                        {
                            ms_next = ms;
                            if(c0[i]>0)
                            {
                                c_prev= c0[i];
                            }
                            k_prev = k1[j];
                            k3_prev= k3[k];
                            sc_prev= sc[m];
                        }
                    }
                }
            }
        }

        if(ms_next == ms_prev)
        {
            exit = false;
            ms_prev = ms_next;
        }
        else
        {
            ms_prev = ms_next;
        }

        cout<<scale<<"\t"<<k_prev<<"\t"<<k3_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

        Resonance(300,510,2);
        gp << "plot \"Resonance_aext3.000000.dat\" title \" Theory \", "<< " \"experiment.dat\" title \"Experiment\"; replot "<<endl;
    }
}
