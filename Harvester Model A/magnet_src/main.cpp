#include <iostream>
#include <fstream>
#include <math.h>
#include "forces.h"

using namespace std;

void GenerateForceRot(double a_min, double a_max, double da, double omega)
{
    magnet m;
    ofstream file("force_rot.dat");
    for(double a = a_min; a<a_max; a+=da)
    {
        file<<a<<"\t"<<m.Force_rot(omega,a)<<endl;
    }
}


void GenerateForcePerspect(double z_min, double z_max, double dz, double v)
{
    magnet m;
    ofstream file("force_perspect.dat");
    for(double z = z_min; z<z_max; z+=dz)
    {
        file<<z<<"\t"<<m.Force_transp(v,z)<<endl;
    }
}


int main()
{
    magnet m;
    /*m.GenerateDerivFlux_transp(0.002,0.004,0.0001);
    m.GenerateEMF_transp(1.0);*/
    //m.Force_transp(1.0,0.004);
    m.GenerateDerivFlux_rot(0.0,0.035,0.001);
    m.GenerateEMF_rot(1.0);

    //GenerateForcePerspect(0.002,0.004,0.0001,1.0);
    //GenerateForceRot(-0.035,0.035,0.001,1.0);
    //magnet m;
    for(double z=0.002;z<=0.004;z+=0.0001)
    {
        cout<<z<<"\t"<<m.Force_progres_interp(1.0,z)*1.0<<"\t"<<pow(m.emf_interp_transp(1.0,z),2)/m.GetResist()<<endl;
    }
    for(double a=-0.03;a<=0.03;a+=0.001)
    {
        cout<<a<<"\t"<<m.Force_rot_interp(1.0,a)*1.0<<"\t"<<pow(m.emf_interp_rot(1.0,a),2)/m.GetResist()<<endl;
    }
    return 0;
}
