#include "rot.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

void rot_optim::coordinate_discent(double dc, double dk, double dk3, double dk5)
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
        c0[0] = c_prev - dc;    k1[0] = k_prev - dk;    k3[0] = k3_prev - dk3;   k5[0] = k5_prev - dk5;
        c0[1] = c_prev;         k1[1] = k_prev;         k3[1] = k3_prev;        k5[1] = k5_prev;
        c0[2] = c_prev + dc;    k1[2] = k_prev + dk;    k3[2] = k3_prev + dk3;   k5[2] = k5_prev + dk5;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0;j<3; j++)
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

void rot_optim::coordinate_discent(double dc, double dk, double dA)
{
    double ms_prev = minsqr();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double a_prev = aext;
    spring_k3 = 0.0;
    spring_k5 = 0.0;

    double ms;
    double c0[3];
    double k1[3];
    double a[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;    k1[0] = k_prev - dk;    a[0] = a_prev - dA;
        c0[1] = c_prev;         k1[1] = k_prev;         a[1] = a_prev;
        c0[2] = c_prev + dc;    k1[2] = k_prev + dk;    a[2] = a_prev + dA;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0;j<3; j++)
            {
                for(int k=0; k<3; k++)
                {
                        air_damp = c0[i];
                        spring_k = k1[j];
                        aext = a[k];
                        ms = minsqr();
                        if(ms<ms_next)
                        {
                            ms_next = ms;
                            c_prev= c0[i];
                            k_prev = k1[j];
                            a_prev = a[k];
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

        cout<<k_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

    }
}

void rot_optim::coordinate_discent(double da)
{
    double ms_prev = minsqr();
    double a_prev = aext;

    double ms;
    double a[3];


    bool exit = true;

    while(exit)
    {
        a[0] = a_prev - da;
        a[1] = a_prev;
        a[2] = a_prev + da;

        for(int i=0; i<3; i++)
        {
            aext = a[i];
            ms = minsqr();
            if(ms == ms_prev)
            {
                exit = false;
            }
            if(ms<ms_prev)
            {
                ms_prev = ms;
                a_prev= a[i];
            }
        }
        cout<<a_prev<<"\t"<<ms_prev<<endl;
    }
}

void rot_optim::coordinate_discent(double da, double dc, double dk, double dk3, double dk5)
{
    double ms_prev = minsqr();
    double ms_next = ms_prev;
    double c_prev = air_damp;
    double k_prev = spring_k;
    double k3_prev= spring_k3;
    double k5_prev= spring_k5;
    double a_prev= aext;

    double ms=0.0;
    double a[3];
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;    k1[0] = k_prev - dk;    k3[0] = k3_prev - dk;   k5[0] = k5_prev - dk5;  a[0] = a_prev - da;
        c0[1] = c_prev;         k1[1] = k_prev;         k3[1] = k3_prev;        k5[1] = k5_prev;        a[1] = a_prev;
        c0[2] = c_prev + dc;    k1[2] = k_prev + dk;    k3[2] = k3_prev + dk;   k5[2] = k5_prev + dk5;  a[2] = a_prev + da;

        ms_next = ms_prev;

        for(int i=0; i<3; i++)
        {
            for(int j=0;j<3; j++)
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
                            aext = a[m];
                            ms = minsqr();
                            if(ms<ms_next)
                            {
                                ms_next = ms;
                                c_prev= c0[i];
                                k_prev = k1[j];
                                k3_prev= k3[k];
                                k5_prev= k5[l];
                                a_prev= a[m];
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
        cout<<aext<<"\t"<<k_prev<<"\t"<<k3_prev<<"\t"<<k5_prev<<"\t"<<c_prev<<"\t"<<ms_prev<<endl;

    }
}
