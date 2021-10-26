/**
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

    Authors: Andrii Sokolov, Elena Blokhina
    2019.
*/


#include <iostream>
#include <boost/tuple/tuple.hpp>    /** Library used for the GUI */

#include "gnuplot-iostream.h"       /** Gnuplot external library */
#include "transp.h"                 /** Progressive mode header file */
#include "rot.h"                    /** Rotational mode header file */

using namespace std;

int main()
{

    progress_optim opt(3.0,0.00339,556.59,-8.1e10,1190000000000000000000.0,-2.9e29,4.16);
    //progress_optim opt(3.0,4.98,556.31,0.0,0.0,1.0);
    opt.read_experiment("experiment_3_progres.dat","experiment_4_progres.dat","experiment_5_progres.dat");
    opt.coordinate_discent_sets(0.01,0.00001,0.01,1e8,1e19);

    //opt.set_aext(3.0);
   // opt.ResonanceRK(300.0,500.0,2.0);
    //opt.Resonance(300.0,600.0,2.0);
    //opt.ResonanceLin(300.0,600.0,2.0);
    //opt.set_aext(4.0);
   // opt.ResonanceRK(300.0,500.0,2.0);
   // opt.Resonance(300.0,600.0,2.0);
   // opt.ResonanceLin(300.0,600.0,2.0);
   // opt.set_aext(5.0);
  //  opt.ResonanceRK(300.0,500.0,2.0);
   // opt.Resonance(300.0,600.0,2.0);
   // opt.ResonanceLin(300.0,600.0,2.0);
   // opt.Replot("Resonance_aext3.000000.dat","Resonance_Lin_aext3.000000.dat","ResonanceRK_aext3.000000.dat","Resonance_aext4.000000.dat","Resonance_Lin_aext4.000000.dat","ResonanceRK_aext4.000000.dat","Resonance_aext5.000000.dat","Resonance_Lin_aext5.000000.dat","ResonanceRK_aext5.000000.dat");
/*
    rot_optim ropt(3.0,0.001,556.59,0.0,0.0,4.16,1.0);
    ropt.read_experiment("experiment_3_rot.dat","experiment_4_rot.dat","experiment_5_rot.dat");
    ropt.SetAext(3.0);
    ropt.ResonanceRK(460,590,0.05);
    ropt.SetAext(4.0);
    ropt.ResonanceRK(460,590,0.05);
    ropt.SetAext(5.0);
    ropt.ResonanceRK(460,590,0.05);
    ropt.Replot("experiment_3_rot.dat","experiment_4_rot.dat","experiment_5_rot.dat","Resonance_rot_RK_aext3.000000.dat","Resonance_rot_RK_aext4.000000.dat","Resonance_rot_RK_aext5.000000.dat");
*/
    //ropt.coordinate_discent(0.1,0.0001,1,1e11,1e22);
    return 0;
}
