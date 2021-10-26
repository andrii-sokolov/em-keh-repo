#include <iostream>
#include <string>

#include <boost/tuple/tuple.hpp>        //Boost module
#include "gnuplot-iostream.h"           //GNU Plot module header

#include "optim.h"                      //Optimization module

using namespace std;

optimization::optimization()
{
    delta = 0.00001;
    st_st_delta = 0.001;
    Max_Counter = 200;
    A_max = 0.01;
    n_of_non_lin = 10;

    mass = 0.003;
    air_damp = 0.15;
    spring_k = 697.3;
    spring_k3= 1.1e11;
    spring_k5= 0.0;
    aext = 3.0;
    scale = 1.0;
    zshift = 1e-6;

    experimental_size= 100000;

    set_emf_c(5.11189, -1.85556e06, 1.14521e11, 1.40326e16, -3.36402e21, 2.1979e26);
    set_force_c(-0.00838619, 6097.03, -1.51023e09, 1.24942e14, 5.12612e18, - 1.13488e24);

    experimental_freq = new double [experimental_size];
    experimental_volt = new double [experimental_size];
}

optimization::optimization      (   double A,       // external acceleration                        [m/s^2]
                                    double C,       // air damping coefficient                      [kg/s]
                                    double K,       // linear spring constant                       [N/m]
                                    double K3,      // 3-rd order spring constant                   [N/m^3]
                                    double K5,      // 5-th order spring constant                   [N/m^5]
                                    double S,       // scale factor                                 [%]
                                    bool t,         // EM coupling
                                    std::string filename // the filename of the experimental data file   <example.dat>
                                )
{
    delta = 0.0001;
    st_st_delta = 0.001;
    Max_Counter = 20;
    A_max = 0.01;
    n_of_non_lin = 1000;

    mass = 0.003;
    air_damp    = C;
    spring_k    = K;
    spring_k3   = K3;
    spring_k5   = K5;
    aext        = A;
    scale       = S;
    resist      = 3170.0;
    zshift     = 1e-6;
    experimental_filename = filename;

    tran = t;

    experimental_size= 100000;

    set_emf_c(5.11189, -1.85556e06, 1.14521e11, 0.0, 0.0, 0.0);
    set_force_c(-0.00838619, 6097.03, -1.51023e09, 1.24942e14, 5.12612e18, - 1.13488e24);

    experimental_freq = new double [experimental_size];
    experimental_volt = new double [experimental_size];
}

optimization::~optimization()
{
    delete experimental_freq;
    delete experimental_volt;
}

void FullOptimization()
/** Function for the implementation of the Coordinate discent method for the given set of files */
{
    int num_of_f = 4;                                  // The number of used files
    string* names = new string[num_of_f];               // The array of file names
        names[0] = "RMS_P1_2400ohm_0.05g.txt";          // Loaded RMS voltage vs frequency file for Aext = 0.5 m/s^2
        names[1] = "RMS_P1_2400ohm_0.1g.txt";           // Loaded RMS voltage vs frequency file for Aext = 1 m/s^2
        names[2] = "RMS_P1_2400ohm_0.3g.txt";           // Loaded RMS voltage vs frequency file for Aext = 3 m/s^2
        names[3] = "RMS_P1_2400ohm_0.5g.txt";           // Loaded RMS voltage vs frequency file for Aext = 5 m/s^2
    bool* trans = new bool[num_of_f];                   //The array of switches (If the circuit is opened or closed)
        trans[0] = true;
        trans[1] = true;
        trans[2] = true;
        trans[3] = true;
    double* acc = new double[num_of_f];                 //The array of accelerations
        acc[0] = 0.5;
        acc[1] = 1.0;
        acc[2] = 3.0;
        acc[3] = 5.0;

    CoordinateDiscentMethod(    404.05,                 // Initial linear spring constant       [N/m]
                                1.768e06,                // Initial 3-rd order spring constant   [N/m^3]
                                0.0,                // Initial 5-th order spring constant   [N/m^5]
                                0.0258,                 // Initial air damping coefficient      [kg/s]
                                1.05,                    // Initial scale factor                 [%]
                                num_of_f,               // The number of files
                                names,                  // Array of names
                                acc,                    // Array of acceleration
                                trans,                  // Array of switches
                                0.01,                   // Step by the linear spring constant   [N/m]
                                1.0e02,                 // Step by the 3-rd order spring const  [N/m^3]
                                0.0,                 // Step by the 5-th order spring const  [N/m^5]
                                0.0001,                 // Step by the air damping coefficient  [kg/s]
                                0.001                    // Step by the scale factor             [%]
                            );                          //  Loading of the CD method

    delete[] names;
    delete[] trans;
    delete[] acc;
}

void FullOptimization_coef()
{
    int num_of_f = 2;                                  // The number of used files
    string* names = new string[num_of_f];               // The array of file names
        names[0] = "RMS_prototype1_0.05g_Voc.txt";      // Unloaded RMS voltage vs frequency file for Aext = 0.5 m/s^2
        names[1] = "RMS_P1_2400ohm_0.05g.txt";          // Loaded RMS voltage vs frequency file for Aext = 0.5 m/s^2
    bool* trans = new bool[num_of_f];                   //The array of switches (If the circuit is opened or closed)
        trans[0] = false;
        trans[1] = true;
    double* acc = new double[num_of_f];                 //The array of accelerations
        acc[0] = 0.5;
        acc[1] = 0.5;

    CoordinateDiscentMethod(    405,                 // Initial linear spring constant       [N/m]
                                0.0,                // Initial 3-rd order spring constant   [N/m^3]
                                0.0,//8.09e08,                // Initial 5-th order spring constant   [N/m^5]
                                0.01,                 // Initial air damping coefficient      [kg/s]
                                0.0,
                                1.0e8,
                                -1.0e16,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                num_of_f,               // The number of files
                                names,                  // Array of names
                                acc,                    // Array of acceleration
                                trans,                  // Array of switches
                                1.0,                    // Step by the linear spring constant   [N/m]
                                0.0,                 // Step by the 3-rd order spring const  [N/m^3]
                                0.0,                 // Step by the 5-th order spring const  [N/m^5]
                                0.001,                 // Step by the air damping coefficient  [kg/s]
                                1.0,
                                1.0e06,
                                1.0e12,
                                0.0,
                                0.0,
                                0.0,
                                1.0e-8
                                                    // Step by the scale factor             [%]
                            );                          //  Loading of the CD method

    delete[] names;
    delete[] trans;
    delete[] acc;
}

void FullOptimization_coef_lin()
{
    int num_of_f = 8;                                  // The number of used files
    string* names = new string[num_of_f];               // The array of file names
        names[0] = "RMS_prototype1_0.05g_Voc.txt";      // Unloaded RMS voltage vs frequency file for Aext = 0.5 m/s^2
        names[1] = "RMS_prototype1_0.1g_Voc.txt";       // Unloaded RMS voltage vs frequency file for Aext = 1 m/s^2
        names[2] = "RMS_prototype1_0.3g_Voc.txt";       // Unloaded RMS voltage vs frequency file for Aext = 3 m/s^2
        names[3] = "RMS_prototype1_0.5g_Voc.txt";       // Unloaded RMS voltage vs frequency file for Aext = 5 m/s^2
        names[4] = "RMS_P1_2400ohm_0.05g.txt";          // Loaded RMS voltage vs frequency file for Aext = 0.5 m/s^2
        names[5] = "RMS_P1_2400ohm_0.1g.txt";           // Loaded RMS voltage vs frequency file for Aext = 1 m/s^2
        names[6] = "RMS_P1_2400ohm_0.3g.txt";           // Loaded RMS voltage vs frequency file for Aext = 3 m/s^2
        names[7] = "RMS_P1_2400ohm_0.5g.txt";           // Loaded RMS voltage vs frequency file for Aext = 5 m/s^2
    bool* trans = new bool[num_of_f];                   //The array of switches (If the circuit is opened or closed)
        trans[0] = false;
        trans[1] = false;
        trans[2] = false;
        trans[3] = false;
        trans[4] = true;
        trans[5] = true;
        trans[6] = true;
        trans[7] = true;
    double* acc = new double[num_of_f];                 //The array of accelerations
        acc[0] = 0.5;
        acc[1] = 1.0;
        acc[2] = 3.0;
        acc[3] = 5.0;
        acc[4] = 0.5;
        acc[5] = 1.0;
        acc[6] = 3.0;
        acc[7] = 5.0;

    CoordinateDiscentMethodLinear(    405,                 // Initial linear spring constant       [N/m]
                                0.0227,                 // Initial air damping coefficient      [kg/s]
                                50.11189,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                0.0,
                                num_of_f,               // The number of files
                                names,                  // Array of names
                                acc,                    // Array of acceleration
                                trans,                  // Array of switches
                                1.0,                    // Step by the linear spring constant   [N/m]
                                0.001,                 // Step by the air damping coefficient  [kg/s]
                                0.01*50.11189,
                                0.5*1.85556e06,
                                0.5*1.14521e11,
                                0.01*1.40326e16,
                                -0.01*3.36402e21,
                                0.01*2.1979e26,
                                1.0e-8
                                                    // Step by the scale factor             [%]
                            );                          //  Loading of the CD method

    delete[] names;
    delete[] trans;
    delete[] acc;
}

void Show_RK_HB_Graphs()
{
    int num_of_f = 6;
    string* names = new string[num_of_f];
        names[0] = "RMSFreq_prototype1_0.1g_Voc.dat";
        names[1] = "RMSFreq_prototype1_0.3g_Voc.dat";
        names[2] = "RMSFreq_prototype1_0.5g_Voc.dat";
        names[3] = "RMSFreq_P1_2400ohm_0.1g.dat";
        names[4] = "RMSFreq_P1_2400ohm_0.3g.dat";
        names[5] = "RMSFreq_P1_2400ohm_0.5g.dat";
    bool* trans = new bool[num_of_f];
        trans[0] = false;
        trans[1] = false;
        trans[2] = false;
        trans[3] = true;
        trans[4] = true;
        trans[5] = true;
    double* acc = new double[num_of_f];
        acc[0] = 1.0;
        acc[1] = 3.0;
        acc[2] = 5.0;
        acc[3] = 1.0;
        acc[4] = 3.0;
        acc[5] = 5.0;

        std::string* plot_strings = new std::string[num_of_f];

    for(int i=0; i<num_of_f; i++)
    {
        optimization opt(acc[i], 0.0194, 403.57, 2.35e06, 8.43e08, 1.17, trans[i],names[i]);
        opt.read_experiment();
        opt.ResonanceHB();
        opt.ResonanceRK();
        plot_strings[i] = opt.GNUPlot_string();
    }

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
    gp1<<"set yrange [0.0: 3.0]"<<std::endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<std::endl;


    for(int i=0; i<num_of_f; i++)
    {
        gp1<<"plot \"ResonanceHB"+names[i]+"\" using 1:3 title \" Theory HB \" lt 1, \"ResonanceRK"+names[i]+"\" using 1:3 title \" Theory HB \" lt 2, \""+names[i]+"\" title \"Experiment\" lt 3"<<std::endl;
    }
    gp1<<"unset multiplot"<<std::endl;


    delete[] names;
    delete[] trans;
    delete[] acc;

}

int main()
{
    FullOptimization();
    //Show_RK_HB_Graphs();
    //optimization opt(1.0,0.0194,403.57,2.35e6,8.43e8,117.0, false, "RMSFreq_prototype1_0.1g_Voc.dat");
    //for(double x=0; x<1e-2; x+=1e-4)
    //{
    //    cout<<x<<"\t"<<pow(opt.emf_interp_progress(1.0, x),2)/2400.0<<"\t"<<opt.force_interp_progress(1.0,x)<<endl;
    //}
    /*optimization opt(1.0,0.0194,403.57,2.35e6,8.43e8,117.0, false, "RMSFreq_prototype1_0.1g_Voc.dat");
    opt.read_experiment();
    opt.ResonanceHB();
    cout << opt.GNUPlot_string() << endl;
    */


    return 0;
}
