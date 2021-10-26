#include "optim.h"
#include "gnuplot-iostream.h"
#include <vector>
#include <iostream>
#include <string>

void CoordinateDiscentMethod(   double initial_k,
                                double initial_k3,
                                double initial_k5,
                                double initial_air_damp,
                                double initial_scale,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dk3,
                                double dk5,
                                double dc,
                                double ds
                            )
{
    std::string* plot_strings = new std::string[num_files];

    int* exp_length = new int[num_files];
    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_air_damp,initial_k,initial_k3,initial_k5,initial_scale,trans[i],files_list[i]);
        opt.read_experiment();
        exp_length[i] = opt.GetExperimentLength();
    }

    double** exp_freq = new double* [num_files];
    double** exp_volt = new double* [num_files];
    for( int i=0; i<num_files; i++)
    {
        exp_freq[i] = new double[exp_length[i]+1];
        exp_volt[i] = new double[exp_length[i]+1];
    }

    double ms_prev = 0.0;

    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_air_damp,initial_k,initial_k3,initial_k5,initial_scale,trans[i],files_list[i]);
        opt.read_experiment();
        opt.ResonanceHB();
        plot_strings[i] = opt.GNUPlot_string();
        opt.GetExperiment(exp_freq[i],exp_volt[i]);
        ms_prev+=opt.minsqr(exp_freq[i],exp_volt[i],exp_length[i]);
    }
    std::cout<<ms_prev<<std::endl;

    Gnuplot gp1;
    gp1<<"set size 1.5,1.5"<<std::endl;
    gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
    gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
    gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;
    gp1<<"set multiplot layout 3,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
    gp1<<"set format y \"%.2f\""<<std::endl;
    gp1<<"set key box opaque"<<std::endl;
    //gp1<<"set yrange [0.0: 3.0]"<<std::endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;


    for(int i=0; i<num_files; i++)
    {
        gp1<<plot_strings[i]<<std::endl;
    }
    gp1<<"unset multiplot"<<std::endl;


    double ms_next = ms_prev;
    double c_prev = initial_air_damp;
    double k_prev = initial_k;
    double k3_prev= initial_k3;
    double k5_prev= initial_k5;
    double sc_prev= initial_scale;

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
                                initial_air_damp = c0[i];
                            else
                                dc = dc/10.0;
                            initial_k = k1[j];
                            initial_k3= k3[k];
                            initial_k5= k5[l];
                            initial_scale = sc[m];

                            ms = 0.0;
                            for(int n=0; n<num_files; n++)
                            {
                                optimization opt(acc_list[n],initial_air_damp,initial_k,initial_k3,initial_k5,initial_scale,trans[n],files_list[n]);
                                ms+=opt.minsqr(exp_freq[n],exp_volt[n],exp_length[n]);
                            }
                            //std::cout<<ms<<std::endl;
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
            for(int n=0; n<num_files; n++)
            {
                optimization opt(acc_list[n],c_prev,k_prev,k3_prev,k5_prev,sc_prev,trans[n],files_list[n]);
                opt.read_experiment();
                opt.ResonanceHB();
            }
        }
        std::cout<<sc_prev<<"\t"<<k_prev<<"\t"<<k3_prev<<"\t"<<k5_prev<<"\t"<<c_prev<<"\t"<<ms_next<<"\t"<<ms_prev<<std::endl;

        gp1<<"set size 1.5,1.5"<<std::endl;

        gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
        gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
        gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;

        gp1<<"set multiplot layout 3,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
        gp1<<"set format y \"%.2f\""<<std::endl;
        gp1<<"set key box opaque"<<std::endl;
        //gp1<<"set yrange [0.0: 3.0]"<<std::endl;
        gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;

        for(int i=0; i<num_files; i++)
        {
            gp1<<plot_strings[i]<<std::endl;
        }
        gp1<<"unset multiplot"<<std::endl;

    }

    delete[] plot_strings;
    for(int i=0; i<num_files; i++)
    {
        delete exp_freq[i];
        delete exp_volt[i];
    }
    delete[] exp_freq;
    delete[] exp_volt;
    delete[] exp_length;

}


void CoordinateDiscentMethod(   double initial_k,
                                double initial_k3,
                                double initial_k5,
                                double initial_c,
                                double initial_e00,
                                double initial_e02,
                                double initial_e04,
                                double initial_e06,
                                double initial_e08,
                                double initial_e10,
                                double initial_shift,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dk3,
                                double dk5,
                                double dc,
                                double de00,
                                double de02,
                                double de04,
                                double de06,
                                double de08,
                                double de10,
                                double dshift
                            )
{
    std::string* plot_strings = new std::string[num_files];

    int* exp_length = new int[num_files];
    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_c,initial_k,initial_k3,initial_k5,1.0,trans[i],files_list[i]);
        opt.read_experiment();
        exp_length[i] = opt.GetExperimentLength();
    }

    double** exp_freq = new double* [num_files];
    double** exp_volt = new double* [num_files];
    for( int i=0; i<num_files; i++)
    {
        exp_freq[i] = new double[exp_length[i]+1];
        exp_volt[i] = new double[exp_length[i]+1];
    }

    double ms_prev = 0.0;

    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_c,initial_k,initial_k3,initial_k5,1.0,trans[i],files_list[i]);
        opt.set_emf_c(initial_e00,initial_e02,initial_e04,initial_e06,initial_e08,initial_e10);
        opt.set_z_shift(initial_shift);
        opt.read_experiment();
        opt.ResonanceHB();
        plot_strings[i] = opt.GNUPlot_string();
        opt.GetExperiment(exp_freq[i],exp_volt[i]);
        ms_prev+=opt.minsqr(exp_freq[i],exp_volt[i],exp_length[i]);
    }
    std::cout<<ms_prev<<std::endl;

    Gnuplot gp1;
    gp1<<"set size 1.5,1.5"<<std::endl;
    gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
    gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
    gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;
    gp1<<"set multiplot layout 1,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
    gp1<<"set format y \"%.2f\""<<std::endl;
    gp1<<"set key box opaque"<<std::endl;
    gp1<<"set autoscale y"<<std::endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;


    for(int i=0; i<num_files; i++)
    {
        gp1<<plot_strings[i]<<std::endl;
    }
    gp1<<"unset multiplot"<<std::endl;


    double ms_next = ms_prev;
    double c_prev = initial_c;
    double k_prev = initial_k;
    double k3_prev= initial_k3;
    double k5_prev= initial_k5;
    double e00_prev=initial_e00;
    double e02_prev=initial_e02;
    double e04_prev=initial_e04;
    double e06_prev=initial_e06;
    double e08_prev=initial_e08;
    double e10_prev=initial_e10;
    double zs_prev =initial_shift;

    double ms=0.0;
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];
    double e00[3];
    double e02[3];
    double e04[3];
    double e06[3];
    double e08[3];
    double e10[3];
    double zs[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        k3[0] = k3_prev - dk3;
        k5[0] = k5_prev - dk5;
        e00[0]= e00_prev - de00;
        e02[0]= e02_prev - de02;
        e04[0]= e04_prev - de04;
        e06[0]= e06_prev - de06;
        e08[0]= e08_prev - de08;
        e10[0]= e10_prev - de10;
        zs[0] = zs_prev - dshift;
        c0[1] = c_prev;
        k1[1] = k_prev;
        k3[1] = k3_prev;
        k5[1] = k5_prev;
        e00[1]= e00_prev;
        e02[1]= e02_prev;
        e04[1]= e04_prev;
        e06[1]= e06_prev;
        e08[1]= e08_prev;
        e10[1]= e10_prev;
        zs[1] = zs_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        k3[2] = k3_prev + dk3;
        k5[2] = k5_prev + dk5;
        e00[2]= e00_prev + de00;
        e02[2]= e02_prev + de02;
        e04[2]= e04_prev + de04;
        e06[2]= e06_prev + de06;
        e08[2]= e08_prev + de08;
        e10[2]= e10_prev + de10;
        zs[2] = zs_prev + dshift;

        ms_next = ms_prev;
        int c = 0;
        for(int i1=0; i1<3; i1++)
        {
            for(int i2=0; i2<3; i2++)
            {
                for(int i3=0; i3<3; i3++)
                {
                    for(int i4=0; i4<3; i4++)
                    {
                        for(int i5=0; i5<3; i5++)
                        {
                            for(int i6=0; i6<3; i6++)
                            {
                                if(c0[i1]>0)
                                    initial_c = c0[i1];
                                else
                                    dc = dc/10.0;
                                initial_k = k1[i2];
                                initial_shift = zs[i3];
                                initial_e00 = e00[i4];
                                initial_e02 = e02[i5];
                                initial_e04 = e04[i6];

                                ms = 0.0;
                                for(int n=0; n<num_files; n++)
                                {
                                    optimization opt(acc_list[n],initial_c,initial_k,initial_k3,initial_k5,1.0,trans[n],files_list[n]);
                                    opt.set_emf_c(initial_e00,initial_e02,initial_e04,initial_e06,initial_e08,initial_e10);
                                    opt.set_z_shift(initial_shift);
                                    ms+=opt.minsqr(exp_freq[n],exp_volt[n],exp_length[n]);
                                }
                                //
                                if(ms<ms_next)
                                {
                                    ms_next = ms;
                                    c_prev = c0[i1];
                                    k_prev = k1[i2];
                                    zs_prev= zs[i3];
                                    e00_prev = e00[i4];
                                    e02_prev = e02[i5];
                                    e04_prev = e04[i6];
                                }
                                //std::cout<<zs_prev<<std::endl;
                                c++;
                            }
                        }
                    }
                }
            }
        }
        for(int n=0; n<num_files; n++)
        {
            optimization opt(acc_list[n],c_prev,k_prev,0.0,0.0,1.0,trans[n],files_list[n]);
            opt.read_experiment();
            opt.set_emf_c(e00_prev,e02_prev,e04_prev,0.0,0.0,0.0);
            opt.set_z_shift(0.0);
            opt.ResonanceHBLinear();
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
        std::cout<<"k= "<<k_prev<<"\tk3= "<<k3_prev<<"\tc= "<<c_prev<<"\te00= "<<e00_prev<<"\te02= "<<e02_prev<<"\te04= "<<e04_prev<<"\tMSQR= "<<ms_next<<std::endl;

        gp1<<"set size 1.5,1.5"<<std::endl;

        gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
        gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
        gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;

        gp1<<"set multiplot layout 1,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
        gp1<<"set format y \"%.2f\""<<std::endl;
        gp1<<"set key box opaque"<<std::endl;
        gp1<<"set yrange [0.0: 3.0]"<<std::endl;
        gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;

        for(int i=0; i<num_files; i++)
        {
            gp1<<plot_strings[i]<<std::endl;
        }
        gp1<<"unset multiplot"<<std::endl;

    }

    delete[] plot_strings;
    for(int i=0; i<num_files; i++)
    {
        delete exp_freq[i];
        delete exp_volt[i];
    }
    delete[] exp_freq;
    delete[] exp_volt;
    delete[] exp_length;
}

void CoordinateDiscentMethodLinear(   double initial_k,
                                double initial_c,
                                double initial_e00,
                                double initial_e02,
                                double initial_e04,
                                double initial_e06,
                                double initial_e08,
                                double initial_e10,
                                double initial_shift,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dc,
                                double de00,
                                double de02,
                                double de04,
                                double de06,
                                double de08,
                                double de10,
                                double dshift
                            )
{
    std::string* plot_strings = new std::string[num_files];

    int* exp_length = new int[num_files];
    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_c,initial_k,0.0,0.0,1.0,trans[i],files_list[i]);
        opt.read_experiment();
        exp_length[i] = opt.GetExperimentLength();
    }

    double** exp_freq = new double* [num_files];
    double** exp_volt = new double* [num_files];
    for( int i=0; i<num_files; i++)
    {
        exp_freq[i] = new double[exp_length[i]+1];
        exp_volt[i] = new double[exp_length[i]+1];
    }

    double ms_prev = 0.0;

    for(int i=0; i<num_files; i++)
    {
        optimization opt(acc_list[i],initial_c,initial_k,0.0,0.0,1.0,trans[i],files_list[i]);
        opt.set_emf_c(initial_e00,initial_e02,initial_e04,initial_e06,initial_e08,initial_e10);
        opt.set_z_shift(initial_shift);
        opt.read_experiment();
        opt.ResonanceHBLinear();
        plot_strings[i] = opt.GNUPlot_string();
        opt.GetExperiment(exp_freq[i],exp_volt[i]);
        ms_prev+=opt.minsqrLin(exp_freq[i],exp_volt[i],exp_length[i]);
    }
    std::cout<<ms_prev<<std::endl;

    Gnuplot gp1;
    gp1<<"set size 1.5,1.5"<<std::endl;
    gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
    gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
    gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
    gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;
    gp1<<"set multiplot layout 4,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
    gp1<<"set format y \"%.2f\""<<std::endl;
    gp1<<"set key box opaque"<<std::endl;
    gp1<<"set yrange [0.0: 3.0]"<<std::endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;


    for(int i=0; i<num_files; i++)
    {
        gp1<<plot_strings[i]<<std::endl;
    }
    gp1<<"unset multiplot"<<std::endl;


    double ms_next = ms_prev;
    double c_prev = initial_c;
    double k_prev = initial_k;
    double e00_prev=initial_e00;
    double e02_prev=initial_e02;
    double e04_prev=initial_e04;
    double e06_prev=initial_e06;
    double e08_prev=initial_e08;
    double e10_prev=initial_e10;
    double zs_prev =initial_shift;

    double ms=0.0;
    double c0[3];
    double k1[3];
    double k3[3];
    double k5[3];
    double e00[3];
    double e02[3];
    double e04[3];
    double e06[3];
    double e08[3];
    double e10[3];
    double zs[3];

    bool exit = true;

    while(exit)
    {
        c0[0] = c_prev - dc;
        k1[0] = k_prev - dk;
        e00[0]= e00_prev - de00;
        e02[0]= e02_prev - de02;
        e04[0]= e04_prev - de04;
        e06[0]= e06_prev - de06;
        e08[0]= e08_prev - de08;
        e10[0]= e10_prev - de10;
        zs[0] = zs_prev - dshift;
        c0[1] = c_prev;
        k1[1] = k_prev;
        e00[1]= e00_prev;
        e02[1]= e02_prev;
        e04[1]= e04_prev;
        e06[1]= e06_prev;
        e08[1]= e08_prev;
        e10[1]= e10_prev;
        zs[1] = zs_prev;
        c0[2] = c_prev + dc;
        k1[2] = k_prev + dk;
        e00[2]= e00_prev + de00;
        e02[2]= e02_prev + de02;
        e04[2]= e04_prev + de04;
        e06[2]= e06_prev + de06;
        e08[2]= e08_prev + de08;
        e10[2]= e10_prev + de10;
        zs[2] = zs_prev + dshift;

        ms_next = ms_prev;
        int c = 0;
        for(int i1=0; i1<3; i1++)
        {
            for(int i2=0; i2<3; i2++)
            {
                for(int i3=0; i3<3; i3++)
                {
                    for(int i4=0; i4<3; i4++)
                    {
                        for(int i5=0; i5<3; i5++)
                        {
                            for(int i6=0; i6<3; i6++)
                            {
                                if(c0[i1]>0)
                                    initial_c = c0[i1];
                                else
                                    dc = dc/10.0;
                                initial_k = k1[i2];
                                initial_e00 = e00[i4];
                                initial_e02 = e02[i5];
                                initial_e04 = e04[i6];

                                ms = 0.0;
                                for(int n=0; n<num_files; n++)
                                {
                                    optimization opt(acc_list[n],initial_c,initial_k,0.0,0.0,1.0,trans[n],files_list[n]);
                                    opt.set_emf_c(initial_e00,initial_e02,initial_e04,initial_e06,initial_e08,initial_e10);
                                    opt.set_z_shift(initial_shift);
                                    ms+=opt.minsqrLin(exp_freq[n],exp_volt[n],exp_length[n]);
                                }
                                //
                                if(ms<ms_next)
                                {
                                    ms_next = ms;
                                    c_prev = c0[i1];
                                    k_prev = k1[i2];
                                    e00_prev = e00[i4];
                                    e02_prev = e02[i5];
                                    e04_prev = e04[i6];
                                }
                                //std::cout<<zs_prev<<std::endl;
                                c++;
                            }
                        }
                    }
                }
            }
        }
        for(int n=0; n<num_files; n++)
        {
            optimization opt(acc_list[n],c_prev,k_prev,0.0,0.0,1.0,trans[n],files_list[n]);
            opt.read_experiment();
            opt.set_emf_c(e00_prev,e02_prev,e04_prev,0.0,0.0,0.0);
            opt.set_z_shift(0.0);
            opt.ResonanceHBLinear();
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

        std::cout<<"k= "<<k_prev<<"\tc= "<<c_prev<<"\te00= "<<e00_prev<<"\te02= "<<e02_prev<<"\te04= "<<e04_prev<<"\tMSQR= "<<ms_next<<std::endl;

        gp1<<"set size 1.5,1.5"<<std::endl;

        gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<std::endl;
        gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<std::endl;
        gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<std::endl;
        gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<std::endl;

        gp1<<"set multiplot layout 4,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<std::endl;
        gp1<<"set format y \"%.2f\""<<std::endl;
        gp1<<"set key box opaque"<<std::endl;
        gp1<<"set yrange [0.0: 3.0]"<<std::endl;
        gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;

        for(int i=0; i<num_files; i++)
        {
            gp1<<plot_strings[i]<<std::endl;
        }
        gp1<<"unset multiplot"<<std::endl;

    }

    delete[] plot_strings;
    for(int i=0; i<num_files; i++)
    {
        delete exp_freq[i];
        delete exp_volt[i];
    }
    delete[] exp_freq;
    delete[] exp_volt;
    delete[] exp_length;
}

