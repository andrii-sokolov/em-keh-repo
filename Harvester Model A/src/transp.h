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

#ifndef TRANSP_H_INCLUDED
#define TRANSP_H_INCLUDED

#include <string>
#include "gnuplot-iostream.h"

using namespace std;

class progress_optim
{
    public:
        /** Default constructor */
        progress_optim();
        /** Constructro with the initial values */
        progress_optim( double A /** External acceleration ampliude [m/s^2]*/,
                        double C /** Air damping coefficient [kg/s]*/,
                        double K /** Linear spring constant [N/m]*/,
                        double K3/** 3-rd order spring constant [N/m^3]*/,
                        double K5/** 5-th order spring constant [N/m^5]*/,
                        double K7/** 7-th order spring constant [N/m^7]*/,
                        double S /** Magnetization scale-factor [1] */);

        /** Destructor (Free all variables) */
        ~progress_optim();

        /** Generate the file with the normalized resonance curve for given parameters (Harmonic balance method) */
        void Resonance_norm( double w_min /** starting frequency*/,
                             double w_max /** ending frequency*/,
                             double dw    /** frequency step*/);
        /** Generate file with unnormalized resonance curve for given parameters (Harmonic balance method) */
        void Resonance( double w_min/** starting frequency*/,
                        double w_max/** ending frequency*/,
                        double dw /** frequency step*/);
        /** Generate file with unnormalized resonance curve obtained by the Runge-Kutta method */
        void ResonanceRK(   double w_min/** starting frequency*/,
                            double w_max/** ending frequency*/,
                            double dw /** frequency step*/);
        /** Generate file with unnormalized linear resonance curve for given parameters */
        void ResonanceLin(  double w_min/** starting frequency*/,
                            double w_max/** ending frequency*/,
                            double dw /** frequency step*/);

        /*****************************************/
        /** Optimization with coordinate discent */
        /*****************************************/

        /** Coordinate discent optimization of SINGLE curve with given Aext (in our case either 3.0 m/s^2, 4.0 m/s^2 or 5.0 m/s^2)
            by [air_damp], [spring_k], [spring_k3], [spring_k5] */
        void coordinate_discent(    double dc   /** Air damping coefficient step [kg/s]    RECOMMENDED: 0.00001 */,
                                    double dk   /** Linear spring constant step [kg/s]     RECOMMENDED: 0.01 */,
                                    double dk3  /** 3-rd order spring constant [kg/s^3]    RECOMMENDED: 1.0e08 */,
                                    double dk5  /** 5-th order spring constant [kg/s^5]    RECOMMENDED: 1.0e19 */);
        /** Coordinate discent optimization of SINGLE curve with given Aext (in our case either 3.0 m/s^2, 4.0 m/s^2 or 5.0 m/s^2)
            by [scale] */
        void coordinate_discent(    double ds   /** Magnetization scale factor step [1]    RECOMMENDED: 0.01*/);
        /** Coordinate discent optimization of SINGLE curve with given Aext (in our case either 3.0 m/s^2, 4.0 m/s^2 or 5.0 m/s^2)
            by [scale], [air_damp], [spring_k], [spring_k3], [spring_k5] */
        void coordinate_discent(    double ds   /** Magnetization scale factor step [1]    RECOMMENDED: 0.01 */,
                                    double dc   /** Air damping coefficient step [kg/s]    RECOMMENDED: 0.00001 */,
                                    double dk   /** Linear spring constant step [kg/s]     RECOMMENDED: 0.01 */,
                                    double dk3  /** 3-rd order spring constant [kg/s^3]    RECOMMENDED: 1.0e08 */,
                                    double dk5  /** 5-th order spring constant [kg/s^5]    RECOMMENDED: 1.0e19 */);
        /** Coordinate discent optimization of SINGLE curve with given Aext (in our case either 3.0 m/s^2, 4.0 m/s^2 or 5.0 m/s^2) in linear regine
            by [scale], [air_damp], [spring_k] */
        void coordinate_discent_linear( double ds/** Magnetization scale factor step [1]    RECOMMENDED: 0.01 */,
                                        double dc/** Air damping coefficient step [kg/s]    RECOMMENDED: 0.00001 */,
                                        double dk/** Linear spring constant step [kg/s]     RECOMMENDED: 0.01 */);
        /** Coordinate discent optimization of SINGLE curve with given Aext (in our case either 3.0 m/s^2, 4.0 m/s^2 or 5.0 m/s^2) in 3-rd orfer non-linear regine
            by [scale], [air_damp], [spring_k], [spring_k3] */
        void coordinate_discent_duffing(double ds/** Magnetization scale factor step [1]    RECOMMENDED: 0.01 */,
                                        double dc/** Air damping coefficient step [kg/s]    RECOMMENDED: 0.00001 */,
                                        double dk/** Linear spring constant step [kg/s]     RECOMMENDED: 0.01 */,
                                        double dk3/** 3-rd order spring constant [kg/s^3]    RECOMMENDED: 1.0e08 */);
        /** Coordinate discent optimization of MULTIPLE curves with given Aext (in our case 3.0 m/s^2, 4.0 m/s^2 and 5.0 m/s^2)
            by [scale], [air_damp], [spring_k], [spring_k3], [spring_k5] */
        void coordinate_discent_sets(   double ds   /** Magnetization scale factor step [1]    RECOMMENDED: 0.01 */,
                                        double dc   /** Air damping coefficient step [kg/s]    RECOMMENDED: 0.00001 */,
                                        double dk   /** Linear spring constant step [kg/s]     RECOMMENDED: 0.01 */,
                                        double dk3  /** 3-rd order spring constant [kg/s^3]    RECOMMENDED: 1.0e08 */,
                                        double dk5  /** 5-th order spring constant [kg/s^5]    RECOMMENDED: 1.0e19 */);

        /***************************************/
        /** Reading of experimental data files */
        /***************************************/

        /** Reading of experimental file data
        EXAMPLES: read_experiment("experiment_3_progres.dat");
                  read_experiment("experiment_4_progres.dat");
                  read_experiment("experiment_5_progres.dat");
        DATA STRUCTURE: <frequency, Hz> "\t" <RMS voltage, V> "\t" <frequency uncertinty, Hz> "\t" <RMS voltage uncertinty, V> */
        void read_experiment(string filename);
        /** Reading of experimental file data for all three external amplitudes
        EXAMPLES: read_experiment("experiment_3_progres.dat","experiment_4_progres.dat","experiment_5_progres.dat");
        DATA STRUCTURE: <frequency, Hz> "\t" <RMS voltage, V> "\t" <frequency uncertinty, Hz> "\t" <RMS voltage uncertinty, V> */
        void read_experiment(string file1, string file2, string file3);

        /** DEBUG: Generation of the experimental sets to check the memory*/
        void generate_experiment_file();

        /********************************/
        /** Parameters settings methods**/
        /********************************/

        /** Set delta_t step for the RK-4 method */
        inline void set_delta(double delt   /** delta_t [s] DEFAULT VALUE: 0.00001 s */) { delta = delt; };
        /** Set steady_state condition (comparisment of the energy per neighboring cycles) */
        inline void set_st_st(double st     /** st_st_delta [1] DEFAULT VALUE: 0.001 */  ) { st_st_delta = st; };
        /** Set the maximal number of cycles (Not to get stuck) */
        inline void set_Max_Counter(int mc  /** max_counter [1] DEFAULT VALUE: 200 */) { Max_Counter = mc; };
        /** Set resonators mass */
        inline void set_mass (double m /** Resonators mass [kg] DEFAULT VALUE: 9.83e-5 [kg]*/) {mass = m;};

        /** Set air damping coefficient */
        inline void set_air_damp (double ad /** Air damping coefficient [kg/s] */) {air_damp = ad;};
        /** Set linear spring constant */
        inline void set_k (double k /** Linear spring constant */){ spring_k = k; };
        /** Set 3-rd order spring constant */
        inline void set_k3(double k3/** 3-rd order spring constant */){ spring_k3 = k3; };
        /** Set 5-th order spring constant */
        inline void set_k5(double k5/** 5-th order spring constant */){ spring_k5 = k5; };
        /** Set 7-th order spring constant */
        inline void set_k7(double k7/** 7-th order spring constant */){ spring_k7 = k7; };
        /** Set external acceleration */
        inline void set_aext (double ae/** external acceleration */) {aext = ae; };

        /********************/
        /**GNUplot interface*/
        /********************/

        /** Replot graphs of experimental data and theoretical data
            USAGE: Replot("experiment_3_progres.dat","experiment_4_progres.dat","experiment_5_progres.dat","Resonance_aext3.000000.dat","Resonance_aext4.000000.dat","Resonance_aext5.000000.dat");*/
        void Replot(string exp_name_1, string exp_name_2, string exp_name_3, string theor_name_1, string theor_name_2, string theor_name_3);
        /** DEBUG Replot graphs */
        void Replot(string f_name_11, string f_name_12, string f_name_13, string f_name_21, string f_name_22, string f_name_23, string f_name_31, string f_name_32, string f_name_33);

        /**********************/
        /**Numerical solutions*/
        /**********************/

        /** 4-th order Runge-Kutta method return Root Mean Square of the ampletude Voltage in steady state*/
        void RK4(   double fext     /** external frequency */,
                    double& ampl_x  /** return amplitude of displacement */,
                    double& ampl_dx /** return amplitude of velocity */);
        /** 4-th order Runge-Kutta method return Root Mean Square of the ampletude Voltage in steady state with export of hte wave-form */
        void RK4(   double fext /** external frequency */,
                    double& ampl_x /** return amplitude of displacement */,
                    double& ampl_dx /** return amplitude of velocity */,
                    string file_name);
        /** Solution of non-linear equation (Harmonic balance method)*/
        void Solve( double* solutions /** array of solutions */,
                    int& n_of_solutions /** number of solutions */,
                    double wext/** external frequency */);

        /*******************/
        /**Linear resonance*/
        /*******************/

        /** Function for linear resonance */
        double Solve_Linear(double fext /** external frequency */);

    private:

        /**********************/
        /**Numerical variables*/
        /**********************/

        double delta;       /** Time step for the Runge-Kutta method */
        double st_st_delta; /** Steady state condition */
        int Max_Counter;    /** The maximal number of iterations in RK */

        /************************************/
        /**Physical parameters of the system*/
        /************************************/

        double mass;        /** Mass of the magnet [kg] DEFAULT: 9.83e-5 kg*/
        double air_damp;    /** Air damping coefficient [kg/s] */
        double spring_k;    /** Linear spring constant [N/m] */
        double spring_k3;   /** 3-rd order spring constant [N/m^3] */
        double spring_k5;   /** 5-th order spring constant [N/m^5] */
        double spring_k7;   /** 7-th order spring constant [N/m^7] */
        double aext;        /** External acceleration [m/s^2]      */
        double scale;       /** Scale factor of e.m.f */

        /*******************************/
        /**1st experimental data set */
        /*******************************/

        double* experimental_freq;      /** Frequency data */
        double* experimental_volt;      /** RMS voltage data */
        double* experimental_volt_norm; /** Normalized RMS voltage data*/

        /********************************/
        /**2nd experimental data set */
        /********************************/

        double* experimental_freq_1;      /** Frequency data */
        double* experimental_volt_1;      /** RMS voltage data */
        double* experimental_volt_norm_1; /** Normalized RMS voltage data*/

        /******************************/
        /**3rd experimental data set*/
        /******************************/

        double* experimental_freq_2;      /** Frequency data */
        double* experimental_volt_2;      /** RMS voltage data */
        double* experimental_volt_norm_2; /** Normalized RMS voltage data*/

        int imax,imax1,imax2;   /** Sizes of each set */

        /*****************************************/
        /** Functions for the Runge-Kutta method */
        /*****************************************/

        /** The expression for the \ddot{x} */
        /** a = aext*cos(wext*t) - air_damp/mass*v - spring_k/mass*x - spring_k3/mass*x^3 - spring_k5/mass*x^5*/
        double acceleration(    double x/** displacement */,
                                double v/** velocity */,
                                double t/** time */,
                                double wext/** external acceleration frequency */);
        /** The expression for the velocity */
        double velocity(    double x/** displacement */,
                            double v/** velocity */,
                            double t/** time */,
                            double wext/** external acceleration frequency */);

        /****************************/
        /** Interpolation functions */
        /****************************/

        /** The EMF interpolation formuls for the progressive mode of vibrations */
        double emf_interp_progress( double scale/** scale factor */,
                                    double v/** velocity */,
                                    double z/** displacement */);
        /** The root means square interpolation formula for the progressive mode */
        double vrms_interp_progress(    double scale/** scale factor */,
                                        double amax/** amplitude of displacement */,
                                        double vmax /** amplitude of velocity */);

        /****************************/
        /**Harmonic Balance method***/
        /****************************/

        /**Non-linear semi-analythical equation for the 3-rd and 5-th order non-linear spring */
        double AdvDuffEQN(  double A0 /** amplitude of displacement */,
                            double wext /** external frequency */);

        /**Minimal square method for the single acceleration */
        double minsqr();
        /**Minimal square method for the set of acceleration */
        double minsqr_sets();
        /**Minimal square method for normalized data sets */
        double minsqr_norm();

        /** Gnuplot initialization of method */
        Gnuplot gp1;

};

#endif // TRANSP_H_INCLUDED
