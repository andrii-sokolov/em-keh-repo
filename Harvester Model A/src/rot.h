#ifndef ROT_H_INCLUDED
#define ROT_H_INCLUDED

#include <string>
#include "gnuplot-iostream.h"

using namespace std;

class rot_optim
{
    public:
        /** Default constructor */
        rot_optim();
        /** Constructro with the initial value */
        rot_optim(double A, double C, double K, double K3, double K5, double S, double as);
        /** Destructor */
        ~rot_optim();

        /** Generate fine with the normalized resonance curve for the given parameters*/
        void Resonance_norm(double w_min/** starting frequency*/, double w_max/** ending frequency*/, double dw /** frequency step*/);
        /** Generate fine with the normalized resonance curve for the given parameters*/
        void Resonance(double w_min/** starting frequency*/, double w_max/** ending frequency*/, double dw /** frequency step*/);
        /** Generate fine with the normalized resonance curve for the given parameters*/
        void ResonanceRK(double w_min/** starting frequency*/, double w_max/** ending frequency*/, double dw /** frequency step*/);

        /** Optimization with coordinate discent */
        void coordinate_discent(double dc, double dk, double dk3, double dk5);
        void coordinate_discent(double ds);
        void coordinate_discent(double ds, double dc, double dk, double dk3, double dk5);
        void coordinate_discent(double dc, double dk, double dA);

        /** Reading experimental data */
        void read_experiment(string filename);
        void read_experiment(string filename1, string filename2, string filename3);

        /** Replot graphs */
        void Replot(string exp_name_1, string exp_name_2, string exp_name_3, string theor_name_1, string theor_name_2, string theor_name_3);

        inline void SetAext(double ae) { aext = ae; };

    private:
        /** Time step for the Runge-Kutta method */
        double delta;
        /** Steady state condition */
        double st_st_delta;
        /** The maximal number of iterations in RK */
        int Max_Counter;

        /** Mass of the magnet */
        double mass;
        /** Air damping coefficient */
        double air_damp;
        /** Linear spring constant */
        double spring_k;
        /** Cubic spring constant */
        double spring_k3;
        /** Five order spring constant */
        double spring_k5;
        /** External acceleration */
        double aext;
        /** Scale factor of e.m.f */
        double scale;
        /** Acceleration scale */
        double acc_scale;


        /** Length and Width of the magnet */
        double magnet_a;
        /** Height of the magnet */
        double magnet_b;
        /** Distance between the magnet and the coil */
        double magnet_h;

        /** Inertia moment of the magnet */
        double Inert;

        /** Experimental data set */
        /** Read first dataset */
        double* experimental_freq;
        double* experimental_volt;
        double* experimental_volt_norm;
        /** Read second dataset */
        double* experimental_freq_1;
        double* experimental_volt_1;
        double* experimental_volt_norm_1;
        /** Read third dataset */
        double* experimental_freq_2;
        double* experimental_volt_2;
        double* experimental_volt_norm_2;

        /** The counters for the length of files */
        int imax,imax1,imax2;


        /** The expression for the \ddot{x} */
        /** a = aext*cos(wext*t) - air_damp/mass*v - spring_k/mass*x - spring_k3/mass*x^3 - spring_k5/mass*x^5*/
        double acceleration(double phi/** rotational displacement */, double omega/** angular velocity */, double t/** time */, double wext/** external acceleration frequency */);
        /** The expression for the velocity */
        double velocity(double phi/** rotational displacement */, double omega/** angular velocity */, double t/** time */, double wext/** external acceleration frequency */);

        /** The EMF interpolation formuls for the progressive mode of vibrations */
        double emf_interp_rotate(double scale/** scale factor */, double omega/** angular velocity */, double phi/** rotational displacement */);
        /** The root means square interpolation formula for the progressive mode */
        double vrms_interp_rotate(double scale/** scale factor */, double phimax/** amplitude of rotational displacement */, double omegamax /** amplitude of angular velocity */);

        /** 4-th order Runge-Kutta method return Root Mean Square of the ampletude Voltage in steady state*/
        void RK4(double fext /** external frequency */, double& ampl_x /** return amplitude of displacement */, double& ampl_dx /** return amplitude of velocity */);
        void RK4(double fext /** external frequency */, double& ampl_x /** return amplitude of displacement */, double& ampl_dx /** return amplitude of velocity */, string file_name);

        /** Non-linear semi-analythical equation for the 3-rd and 5-th order non-linear spring */
        double AdvDuffEQN(double A0 /** amplitude of displacement */, double wext /** external frequency */);
        /** Solution of non-linear equation */
        void Solve(double* solutions /** array of solutions */, int& n_of_solutions /** number of solutions */, double wext/** external frequency */);
        /** Solution of non-linear equation */
        double Solve(double wext/** external frequency */);

        /** Minimal square methods */
        /** Simple minimal square fumction for the single file*/
        double minsqr();
        /** Multiple minimal square method for 3 files at the same time */
        double minsqr_sets();
        /** Minimal square method for normalized functions */
        double minsqr_norm();

        Gnuplot gp;
};

#endif // ROT_H_INCLUDED
