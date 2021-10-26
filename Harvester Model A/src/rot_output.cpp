#include "rot.h"
#include "gnuplot-iostream.h"

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>


/** Replot graphs */
void rot_optim::Replot(string exp_name_1, string exp_name_2, string exp_name_3, string theor_name_1, string theor_name_2, string theor_name_3)
{
    gp<<"set size 1.5,1.5"<<endl;

    gp<<"if (!exists(\"MP_LEFT\"))   MP_LEFT = .1"<<endl;
    gp<<"if (!exists(\"MP_RIGHT\"))  MP_RIGHT = .95"<<endl;
    gp<<"if (!exists(\"MP_BOTTOM\")) MP_BOTTOM = .1"<<endl;
    gp<<"if (!exists(\"MP_TOP\"))    MP_TOP = .9"<<endl;
    gp<<"if (!exists(\"MP_GAP\"))    MP_GAP = 0.05"<<endl;

    gp<<"set multiplot layout 2,2 columnsfirst title \"{/:Bold=Theory and experiment}\" margins screen MP_LEFT, MP_RIGHT, MP_BOTTOM, MP_TOP spacing screen MP_GAP"<<endl;
    gp<<"set format y \"%.2f\""<<endl;
    gp<<"set key box opaque"<<endl;
    gp<<"set ylabel \"RMS Voltage, V\""<<endl;

    gp<<"plot \""<<theor_name_1<<"\" title \" Theory \" lt 1, "<<" \""<<exp_name_1<<"\" title \"Experiment\" lt 2"<<endl;
    gp<<"set xlabel 'Frequency, Hz'"<<endl;
    gp<<"plot \""<<theor_name_2<<"\" title \" Theory \" lt 3, "<<" \""<<exp_name_2<<"\" title \"Experiment\" lt 4"<<endl;

    gp<<"plot \""<<theor_name_3<<"\" title \" Theory \" lt 5, "<<" \""<<exp_name_3<<"\" title \"Experiment\" lt 6"<<endl;
    gp<<"set xlabel 'Frequency, Hz'"<<endl;
    gp<<"unset multiplot"<<endl;
}
