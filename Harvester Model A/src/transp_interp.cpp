#include "transp.h"
#include <math.h>

double progress_optim::emf_interp_progress(double scale, double v, double z)
{
    double h = 0.003;
    return(scale*v*(-287.3081201368885 +
                    721237.8199818128 * (z+h) +
                    -7.933384389830811e8 * pow(z+h,2) +
                    4.980186190109384e11 * pow(z+h,3) +
                    -1.9488231667351844e14 * pow(z+h,4) +
                    4.8644226666230424e16 * pow(z+h,5) +
                    -7.560515532157858e18 * pow(z+h,6) +
                    6.688351807224581e20 * pow(z+h,7) +
                    -2.578102130059701e22 * pow(z+h,8)));
}

double progress_optim::vrms_interp_progress(double scale, double amax, double vmax)
{
    double h = 0.003;
    return(sqrt((259325.76862828282 +
                4.556236038929291e43*pow(amax,16) +
                -1.3019858397420669e9*h +
                3.0663488204751416e12*pow(h,2) +
                -4.4941769672094945e15*pow(h,3) +
                4.585934061123207e18*pow(h,4) +
                -3.453425207301887e21*pow(h,5) +
                1.984702156902755e24*pow(h,6) +
                -8.877603537927839e26*pow(h,7) +
                3.1229399610364277e29*pow(h,8) +
                -8.667267084308521e31*pow(h,9) +
                1.8912932051605642e34*pow(h,10) +
                -3.210450140286099e36*pow(h,11) +
                4.155694542801182e38*pow(h,12) +
                -3.965215590253493e40*pow(h,13) +
                -1.0141204801825835e31*pow(amax,3)*pow(h,13) +
                2.6300663575518796e42*pow(h,14) +
                -1.0834256043708406e44*pow(h,15) +
                2.0880943010298603e45*pow(h,16) +
                pow(amax,14)*(6.886587325376931e40 - 4.2552757349404045e43*h + 6.560979896058179e45*pow(h,2)) +
                pow(amax,12)*(1.3392374991449123e37 - 1.6612084845886215e40*h + 7.712977804422165e42*pow(h,2) - 1.5886362743777512e45*pow(h,3) + 1.2247162472641934e47*pow(h,4)) +
                pow(amax,10)*(7.757257286791377e32 - 1.4484648093868923e36*h + 1.1249594992817264e39*pow(h,2) - 4.6513837568481405e41*pow(h,3) + 1.0798168926191031e44*pow(h,4) - 1.334454470477311e46*pow(h,5) +
                6.858410984679483e47*pow(h,6)) + pow(amax,8)*(1.7078577911917962e28 - 4.2659205180581e31*h + 4.654354372074826e34*pow(h,2) - 2.8969296187737847e37*pow(h,3) + 1.1249594992817265e40*pow(h,4) -
                2.7908302541088845e42*pow(h,5) + 4.3192675704764124e44*pow(h,6) - 3.8127270585066034e46*pow(h,7) + 1.4696594967170322e48*pow(h,8)) +
                pow(amax,6)*(1.5505485600802774e23 - 4.8549394348042877e26*h + 6.831431164767184e29*pow(h,2) - 5.687894024077466e32*pow(h,3) + 3.1029029147165506e35*pow(h,4) - 1.1587718475095139e38*pow(h,5) +
                2.999891998084604e40*pow(h,6) - 5.31586715068359e42*pow(h,7) + 6.170382243537731e44*pow(h,8) - 4.2363633983406695e46*pow(h,9) + 1.3063639970818063e48*pow(h,10)) +
                pow(amax,4)*(5.732417576404009e17 - 2.1583907545636796e21*h + 3.7213165441926656e24*pow(h,2) - 3.883951547843429e27*pow(h,3) + 2.732572465906874e30*pow(h,4) - 1.365094565778592e33*pow(h,5) +
                4.96464466354648e35*pow(h,6) - 1.324310682868016e38*pow(h,7) + 2.571335998358232e40*pow(h,8) - 3.543911433789059e42*pow(h,9) + 3.2908705298867898e44*pow(h,10) -
                1.8485949374577467e46*pow(h,11) + 4.7504145348429316e47*pow(h,12)) + pow(amax,2)*
                (7.665872051187854e11 - 3.3706327254071205e15*h + 6.878901091684811e18*pow(h,2) - 8.633563018254718e21*pow(h,3) + 7.442633088385331e24*pow(h,4) - 4.6607418574121154e27*pow(h,5) +
                2.186057972725499e30*pow(h,6) - 7.80054037587767e32*pow(h,7) + 2.1277048558056346e35*pow(h,8) - 4.414368942893386e37*pow(h,9) + 6.856895995621951e39*pow(h,10) -
                7.732170400994312e41*pow(h,11) + 5.983400963430528e43*pow(h,12) - 2.8439922114734567e45*pow(h,13) + 6.264282903089581e46*pow(h,14)))*pow(scale,2)*pow(vmax,2))/sqrt(2*M_PI));
}

