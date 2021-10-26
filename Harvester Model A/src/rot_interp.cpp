#include "rot.h"
#include <math.h>

double rot_optim::emf_interp_rotate(double scale, double omega, double phi)
{
    return(scale*omega*(-0.00012035207598931708*phi - 4.098063691078718e-11*pow(phi,2) - 0.0009464994355293829*pow(phi,3) + 0.00407300408196066*pow(phi,5) - 0.007116892648305839*pow(phi,7)));
}

double rot_optim::vrms_interp_rotate(double scale, double phimax, double omegamax)
{
    return(scale*sqrt(2.0153058715948363e-6*pow(phimax,2)*(0.0008984133872117744 + 0.007065501420544434*pow(phimax,2) - 0.0016383952959970593*pow(phimax,4) -
     0.08136925048172404*pow(phimax,6) + 0.3059080614991014*pow(phimax,8) - 0.46353045474682364*pow(phimax,10) + 0.32903887900147005*pow(phimax,12))*pow(omegamax,2)));
}
