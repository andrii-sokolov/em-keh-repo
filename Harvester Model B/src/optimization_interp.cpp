#include "optim.h"                      //Optimization module

#include <math.h>                       //Mathematical module

double optimization::emf_interp_progress(double v, double z)
{
    return(v*scale*(emf00 + emf02*pow(z-zshift,2) + emf04*pow(z-zshift,4) + emf06*pow(z-zshift,6) + emf08*pow(z-zshift,8) + emf10*pow(z-zshift,10)));
}

double optimization::vrms_interp_progress(double A0, double wext)
{

    return( scale*sqrt(pow(A0,2)*pow(wext,2)*(262144*pow(emf00,2) +
       4199*pow(A0,20)*pow(emf10,2) +
       9724*pow(A0,18)*emf10*(emf08 + 95*emf10*pow(zshift,2)) +
       262144*pow(zshift,4)*pow(emf02 + emf04*pow(zshift,2) +
          emf06*pow(zshift,4) + emf08*pow(zshift,6) + emf10*pow(zshift,8),2)
         + 5720*pow(A0,16)*(pow(emf08,2) + 306*emf08*emf10*pow(zshift,2) +
          emf10*(2*emf06 + 4845*emf10*pow(zshift,4))) +
       13728*pow(A0,14)*(emf04*emf10 +
          emf06*(emf08 + 120*emf10*pow(zshift,2)) +
          60*pow(zshift,2)*(pow(emf08,2) + 51*emf08*emf10*pow(zshift,2) +
             323*pow(emf10,2)*pow(zshift,4))) +
       14336*pow(A0,8)*(pow(emf04,2) +
          2*emf02*(emf06 + 45*emf08*pow(zshift,2) +
             495*emf10*pow(zshift,4)) +
          emf04*(90*emf06*pow(zshift,2) + 990*emf08*pow(zshift,4) +
             6006*emf10*pow(zshift,6)) +
          3*pow(zshift,4)*(165*pow(emf06,2) +
             2002*emf06*emf08*pow(zshift,2) +
             4290*pow(emf08,2)*pow(zshift,4) +
             8580*emf06*emf10*pow(zshift,4) +
             29172*emf08*emf10*pow(zshift,6) +
             41990*pow(emf10,2)*pow(zshift,8))) +
       1024*emf00*(21*pow(A0,10)*emf10 +
          28*pow(A0,8)*(emf08 + 45*emf10*pow(zshift,2)) +
          40*pow(A0,6)*(emf06 + 28*emf08*pow(zshift,2) +
             210*emf10*pow(zshift,4)) +
          512*pow(zshift,2)*(emf02 + emf04*pow(zshift,2) +
             emf06*pow(zshift,4) + emf08*pow(zshift,6) + emf10*pow(zshift,8)
             ) + 128*pow(A0,2)*
           (emf02 + 6*emf04*pow(zshift,2) + 15*emf06*pow(zshift,4) +
             28*emf08*pow(zshift,6) + 45*emf10*pow(zshift,8)) +
          64*pow(A0,4)*(emf04 + 15*emf06*pow(zshift,2) +
             70*pow(zshift,4)*(emf08 + 3*emf10*pow(zshift,2)))) +
       8448*pow(A0,12)*(pow(emf06,2) +
          182*emf06*pow(zshift,2)*(emf08 + 20*emf10*pow(zshift,2)) +
          2*(emf02*emf10 + 910*pow(emf08,2)*pow(zshift,4) +
             18564*emf08*emf10*pow(zshift,6) +
             62985*pow(emf10,2)*pow(zshift,8) +
             emf04*(emf08 + 91*emf10*pow(zshift,2)))) +
       21504*pow(A0,10)*(emf02*(emf08 + 66*emf10*pow(zshift,2)) +
          emf04*(emf06 + 66*emf08*pow(zshift,2) + 1001*emf10*pow(zshift,4)) +
          11*pow(zshift,2)*(3*pow(emf06,2) +
             91*emf06*pow(zshift,2)*(emf08 + 8*emf10*pow(zshift,2)) +
             26*pow(zshift,4)*(14*pow(emf08,2) +
                153*emf08*emf10*pow(zshift,2) +
                323*pow(emf10,2)*pow(zshift,4)))) +
       131072*pow(A0,2)*pow(zshift,2)*
        (3*pow(emf02,2) + emf02*
           (15*emf04*pow(zshift,2) + 28*emf06*pow(zshift,4) +
             45*emf08*pow(zshift,6) + 66*emf10*pow(zshift,8)) +
          pow(zshift,4)*(14*pow(emf04,2) +
             emf04*(45*emf06*pow(zshift,2) + 66*emf08*pow(zshift,4) +
                91*emf10*pow(zshift,6)) +
             pow(zshift,4)*(33*pow(emf06,2) +
                91*emf06*emf08*pow(zshift,2) +
                60*pow(emf08,2)*pow(zshift,4) +
                120*emf06*emf10*pow(zshift,4) +
                153*emf08*emf10*pow(zshift,6) +
                95*pow(emf10,2)*pow(zshift,8)))) +
       32768*pow(A0,4)*(pow(emf02,2) +
          10*emf02*(3*emf04*pow(zshift,2) + 14*emf06*pow(zshift,4) +
             42*emf08*pow(zshift,6) + 99*emf10*pow(zshift,8)) +
          pow(zshift,4)*(70*pow(emf04,2) +
             emf04*(420*emf06*pow(zshift,2) + 990*emf08*pow(zshift,4) +
                2002*emf10*pow(zshift,6)) +
             pow(zshift,4)*(495*pow(emf06,2) +
                2002*emf06*emf08*pow(zshift,2) +
                1820*pow(emf08,2)*pow(zshift,4) +
                3640*emf06*emf10*pow(zshift,4) +
                6120*emf08*emf10*pow(zshift,6) +
                4845*pow(emf10,2)*pow(zshift,8)))) +
       40960*pow(A0,6)*(emf02*(emf04 + 28*emf06*pow(zshift,2) +
             210*emf08*pow(zshift,4) + 924*emf10*pow(zshift,6)) +
          pow(zshift,2)*(14*pow(emf04,2) +
             21*emf04*(10*emf06*pow(zshift,2) + 44*emf08*pow(zshift,4) +
                143*emf10*pow(zshift,6)) +
             pow(zshift,4)*(462*pow(emf06,2) +
                3003*emf06*emf08*pow(zshift,2) +
                4004*pow(emf08,2)*pow(zshift,4) +
                8008*emf06*emf10*pow(zshift,4) +
                18564*emf08*emf10*pow(zshift,6) +
                19380*pow(emf10,2)*pow(zshift,8))))))/(512.*sqrt(2)));
}

double optimization::force_interp_progress(double v/** velocity */, double z/** displacement */)
{
    return(pow(emf_interp_progress(v,z),2)/(resist*v));
}

