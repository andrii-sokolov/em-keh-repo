#include "rot.h"
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>

void rot_optim::read_experiment(string filename)
{
    /** Reading the experimental data */
    ifstream infile(filename);
    string file_f,file_v,file_df,file_dv;
    int i = 0;
    double v_max = 0.0;
    while (infile>>file_f>>file_v>>file_df>>file_dv)
    {
        experimental_freq[i] = stod (file_f);
        experimental_volt[i] = stod (file_v)/1000;
        if(experimental_volt[i] > v_max)
        {
            v_max = experimental_volt[i];
        }
        i++;
    }

    imax = i;
    /** Normalizing of the experimental data */
    for(i = 0; i<imax; i++)
    {
        experimental_volt_norm[i] = experimental_volt[i]/v_max;
    }
    //cout<<experimental_volt_norm[10]<<endl;
}

void rot_optim::read_experiment(string file1, string file2, string file3)
{
    /** Reading the experimental data */
    ifstream infile(file1);
    string file_f,file_v,file_df,file_dv;
    int i = 0;
    double v_max = 0.0;
    while (infile>>file_f>>file_v>>file_df>>file_dv)
    {
        experimental_freq[i] = stod (file_f);
        experimental_volt[i] = stod (file_v);
        if(experimental_volt[i] > v_max)
        {
            v_max = experimental_volt[i];
        }
        i++;
    }
    imax = i;

    ifstream infile1(file2);
    i = 0;
    double v_max1 = 0.0;
    while (infile1>>file_f>>file_v>>file_df>>file_dv)
    {
        experimental_freq_1[i] = stod (file_f);
        experimental_volt_1[i] = stod (file_v);
        if(experimental_volt_1[i] > v_max1)
        {
            v_max1 = experimental_volt_1[i];
        }
        i++;
    }
    imax1 = i;

    ifstream infile2(file3);
    i = 0;
    double v_max2 = 0.0;
    while (infile2>>file_f>>file_v>>file_df>>file_dv)
    {
        experimental_freq_2[i] = stod (file_f);
        experimental_volt_2[i] = stod (file_v);
        if(experimental_volt_2[i] > v_max2)
        {
            v_max2 = experimental_volt_2[i];
        }
        i++;
    }
    imax2 = i;

    /** Normalizing of the experimental data */
    for(i = 0; i<imax; i++)
        experimental_volt_norm[i] = experimental_volt[i]/v_max;
    for(i = 0; i<imax1; i++)
        experimental_volt_norm_1[i] = experimental_volt_1[i]/v_max1;
    for(i = 0; i<imax2; i++)
        experimental_volt_norm_2[i] = experimental_volt_2[i]/v_max2;

}
