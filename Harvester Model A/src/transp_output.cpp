#include "transp.h"
#include <fstream>
#include "gnuplot-iostream.h"

void progress_optim::generate_experiment_file()
{
    ofstream ofs("experiment.dat");
    for(int i=0; i<imax; i++)
    {
        ofs<<experimental_freq[i]<<"\t"<<experimental_volt[i]<<endl;
    }
    ofs.close();
}

void progress_optim::Replot(string exp_name_1, string exp_name_2, string exp_name_3, string theor_name_1, string theor_name_2, string theor_name_3)
{
    gp1<<"set size 1.5,1.5"<<endl;

    gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<endl;
    gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<endl;
    gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<endl;
    gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<endl;
    gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<endl;

    gp1<<"set multiplot layout 2,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<endl;
    gp1<<"set format y \"%.2f\""<<endl;
    gp1<<"set key box opaque"<<endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<endl;

    gp1<<"plot \""<<theor_name_1<<"\" title \" Theory \" with lines, "<<" \""<<exp_name_1<<"\" title \"Experiment\" lt 2"<<endl;
    gp1<<"set xlabel 'Frequency, Hz'"<<endl;
    gp1<<"plot \""<<theor_name_2<<"\" title \" Theory \" with lines, "<<" \""<<exp_name_2<<"\" title \"Experiment\" lt 4"<<endl;

    gp1<<"plot \""<<theor_name_3<<"\" title \" Theory \" with lines, "<<" \""<<exp_name_3<<"\" title \"Experiment\" lt 6"<<endl;
    gp1<<"set xlabel 'Frequency, Hz'"<<endl;
    gp1<<"unset multiplot"<<endl;

}

void progress_optim::Replot(string f_name_11, string f_name_12, string f_name_13, string f_name_21, string f_name_22, string f_name_23, string f_name_31, string f_name_32, string f_name_33)
{
    gp1<<"set size 1.5,1.5"<<endl;

    gp1<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<endl;
    gp1<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<endl;
    gp1<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<endl;
    gp1<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<endl;
    gp1<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<endl;

    gp1<<"set multiplot layout 2,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<endl;
    gp1<<"set format y \"%.4f\""<<endl;
    gp1<<"set key box opaque"<<endl;
    gp1<<"set ylabel \"RMS Voltage, V\""<<endl;

    gp1<<"plot \""<<f_name_11<<"\" title \""<<f_name_11<<"\" lt 1, "<<" \""<<f_name_12<<"\" title \""<<f_name_12<<"\" lt 2,"<<" \""<<f_name_13<<"\" title \""<<f_name_13<<"\" lt 3"<<endl;
    gp1<<"set xlabel 'Frequency, Hz'"<<endl;
    gp1<<"plot \""<<f_name_21<<"\" title \""<<f_name_21<<"\" lt 1, "<<" \""<<f_name_22<<"\" title \""<<f_name_22<<"\" lt 2,"<<" \""<<f_name_23<<"\" title \""<<f_name_23<<"\" lt 3"<<endl;

    gp1<<"plot \""<<f_name_31<<"\" title \""<<f_name_31<<"\" lt 1, "<<" \""<<f_name_32<<"\" title \""<<f_name_32<<"\" lt 2,"<<" \""<<f_name_33<<"\" title \""<<f_name_33<<"\" lt 3"<<endl;
    gp1<<"set xlabel 'Frequency, Hz'"<<endl;
    gp1<<"unset multiplot"<<endl;

}
