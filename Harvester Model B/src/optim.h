#ifndef OPTIM_H_INCLUDED
#define OPTIM_H_INCLUDED

#include <string>

class optimization
{
    private:

        int experimental_size;              //size of the experimental stack

        std::string experimental_filename;   // The filename of experimental data EXAMPLE: "example.dat"

        /**********************/
        /**Numerical constants*/
        /**********************/
        double delta;               // Time step for the Runge-Kutta method
        double st_st_delta;         // Steady state condition
        int Max_Counter;            // The maximal number of iterations in RK
        double A_max;               // Maximal value of the amplitude
        int n_of_non_lin;           // Number of points in non-linear equation solution

        /*********************/
        /**Physical constants*/
        /*********************/
        double mass;                // Mass of the magnet [kg]
        double air_damp;            // Air damping coefficient [kg/s]
        double spring_k;            // Linear spring constant [N/m]
        double spring_k3;           // 3-rd order spring constant [N/m^3]
        double spring_k5;           // 5-th order spring constant [N/m^5]
        double aext;                // External acceleration [m/s^2]
        double scale;               // Scale factor of e.m.f [%]
        double resist;              // Resistance
        bool tran;                  // true if there is electromagnetic coupling, false if not
        double zshift;             // Shift of the antisimmetry

        /***************************/
        /** Experimental data sets */
        /***************************/
        double* experimental_freq;  // Container for frequency
        double* experimental_volt;  // Container for RMS voltage
        int imax;                   // The number of points in the data set

        /******************************************************************************************************/
        /** The expression for the \ddot{x} *******************************************************************/
        /******************************************************************************************************/
        double acceleration (   double x,       // displacement                             [m]
                                double v,       // velocity                                 [m/s]
                                double t,       // time                                     [s]
                                double wext     // external acceleration circular frequency [1/s]
                            );
        /* a = aext*cos(wext*t) - air_damp/mass*v - spring_k/mass*x - spring_k3/mass*x^3 - spring_k5/mass*x^5*/

        double velocity(        double x,       // displacement                             [m]
                                double v,       // velocity                                 [m/s]
                                double t,       // time                                     [s]
                                double wext);   // external acceleration circular frequency [1/s]
        /* The expression for the velocity \dot{x} = v */

        double AdvDuffEQN(      double A0,      // amplitude of displacement        [m]
                                double wext);   // external circular frequency      [1/s]
        /* Non-linear semi-analythical equation for the 3-rd and 5-th order non-linear spring */

        double at0  (   double A0,          // amplitude of displacement        [m]
                        double wext         // external circular frequency      [1/s]
                    );
        /* Transducer force fundamental harmonic coefficient (cos) */

        double bt0  (   double A0,          // amplitude of displacement        [m]
                        double wext         // external circular frequency      [1/s]
                    );
        /* Transducer force fundamental harmonic coefficient (sin) */

        double at0_num  (   double A0,          // amplitude of displacement        [m]
                            double wext         // external circular frequency      [1/s]
                        );
        /*DEBUG: Numerical transducer force fundamental harmonic coefficient (cos) */

        double bt0_num  (   double A0,          // amplitude of displacement        [m]
                            double wext         // external circular frequency      [1/s]
                        );
        /*DEBUG: Numerical transducer force fundamental harmonic coefficient (sin) */
        /****************************/
        /** Coefficients for e.m.f. */
        /****************************/
        double emf00, emf02, emf04, emf06, emf08, emf10;

        /************************************/
        /** Coefficients for magnetic force */
        /************************************/
        double force_00, force_02, force_04, force_06, force_08, force_10;

    public:

        double minsqr();
        /* Function returns minimal square difference between the experimental data and theory */

        double minsqr   (   double* exp_freq,
                            double* exp_volt,
                            int length
                        );
        /* Function returns minimal square difference between the experimental data and theory */

        double minsqrLin   (   double* exp_freq,
                            double* exp_volt,
                            int length
                        );
        /* Function returns minimal square difference between the experimental data and theory */

        std::string GNUPlot_string();

        optimization();
        /* Default constructor */

        optimization    (   double A,       // external acceleration                        [m/s^2]
                            double C,       // air damping coefficient                      [kg/s]
                            double K,       // linear spring constant                       [N/m]
                            double K3,      // 3-rd order spring constant                   [N/m^3]
                            double K5,      // 5-th order spring constant                   [N/m^5]
                            double S,       // scale factor                                 [%]
                            bool t,         // is the electromagnetic coupling true
                            std::string filename // the filename of the experimental data file   <example.dat>
                        );
        /* Constructro with the initial values */

        optimization    (   double A,       // external acceleration                        [m/s^2]
                            double C,       // air damping coefficient                      [kg/s]
                            double K,       // linear spring constant                       [N/m]
                            double K3,      // 3-rd order spring constant                   [N/m^3]
                            double K5,      // 5-th order spring constant                   [N/m^5]
                            double S,       // scale factor                                 [%]
                            double e00,
                            double e02,
                            double e04,
                            double e06,
                            double e08,
                            double e10,
                            double z_shift,
                            bool t,         // is the electromagnetic coupling true
                            std::string filename // the filename of the experimental data file   <example.dat>
                        );
        /* Constructro with the initial values */

        ~optimization();
        /* Destructor */

        void ResonanceHB(   double f_min,       // starting frequency   [Hz]
                            double f_max,       // ending frequency     [Hz]
                            double df           // frequency step       [Hz]
                        );
        /* Generating of the resonance curve with harmonic balance method */

        void ResonanceHB();
        /* Generating of the resonance curve with harmonic balance method */

        void ResonanceHBLinear();
        /* Generating of the resonance curve with harmonic balance method */

        void ResonanceRK    (   double f_min,   // starting frequency   [Hz]
                                double f_max,   // ending frequency     [Hz]
                                double df       // frequency step       [Hz]
                            );
        /* Generating of the resonance curve using Runge-Kutta method */

        void ResonanceRK    ();
        /* Generating of the resonance curve using Runge-Kutta method */

        void ResonanceLin   (   double f_min,   // starting frequency   [Hz]
                                double f_max,   // ending frequency     [Hz]
                                double df       // frequency step       [Hz]
                            );
        /* Linear resonance curve for given parameters */

        void read_experiment();
        /* Read experimental data file */

        /**************************/
        /** Setting of parameters */
        /**************************/

        inline void set_delta(double delt) { delta = delt; };
        inline void set_st_st(double st  ) { st_st_delta = st; };
        inline void set_Max_Counter(int mc) { Max_Counter = mc; };

        inline void set_mass (double m) {mass = m;};
        inline void set_air_damp (double ad) {air_damp = ad;};

        inline void set_k (double k ){ spring_k = k; };
        inline void set_k3(double k3){ spring_k3 = k3; };
        inline void set_k5(double k5){ spring_k5 = k5; };

        inline void set_aext (double ae) {aext = ae; };

        inline void set_emf_c( double e00, double e02, double e04, double e06, double e08, double e10)
                    {emf00 = e00;  emf02 = e02; emf04 = e04; emf06 = e06; emf08 = e08; emf10 = e10; };

        inline void set_force_c( double f00, double f02, double f04, double f06, double f08, double f10)
                    {force_00 = f00;  force_02 = f02; force_04 = f04; force_06 = f06; force_08 = f08; force_10 = f10; };

        inline void set_z_shift( double zs ){  zshift = zs; };

        /**************************************************************************************************/
        /** 4-th order Runge-Kutta method return Root Mean Square of the ampletude Voltage in steady state*/
        /**************************************************************************************************/

        void RK4(   double fext,        // external frequency [Hz]
                    double& ampl_x,     // return amplitude of displacement
                    double& ampl_dx,    // return amplitude of velocity
                    double& rms_volt    // return RMS voltage
                );
        /*  Method withot waveforms generation */

        void RK4(   double fext,            // external frequency [Hz]
                    double& ampl_x,         // return amplitude of displacement
                    double& ampl_dx,        // return amplitude of velocity
                    std::string file_name   // the name of the file with waveform info
                );
        /* Method with waveform file generation */

        /*********************************************/
        /** Harmonic balance equation solving method */
        /*********************************************/

        void Solve  (   double* solutions,      // array of solutions
                        int& n_of_solutions,    // number of solutions      [%]
                        double wext             // external frequency       [Hz]
                    );
        /* Solution method for the non-linear equation */

        double Solve_Linear (    double fext    // external frequency   [Hz]
                            );
        /* The same thing for the linear resonance */

        double emf_interp_progress( double v,           // velocity     [m/s]
                                    double z);          // displacement [m]
        /* The EMF interpolation formula for the progressive mode of vibrations */

        double force_interp_progress(   double v,          // velocity      [m/s]
                                        double z);         // displacement  [m]
        /* Interpolation formula for the electromagnetic force */

        double vrms_interp_progress(    double amax,        // amplitude of displacement    [m]
                                        double vmax);       // amplitude of velocity        [m/s]
        /* The root means square interpolation formula for the progressive mode */

        inline void    GetExperiment    (   double* freq,
                                            double* volt
                                        ){
                                            for(int i=0; i<imax; i++)
                                            {
                                                freq[i] = experimental_freq[i];
                                                volt[i] = experimental_volt[i];
                                            }
                                        };

        inline int GetExperimentLength(){return(imax);};
};

void CoordinateDiscentMethod(   double initial_k,
                                double initial_k3,
                                double initial_k5,
                                double initial_air_damp,
                                double initial_scale,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dk3,
                                double dk5,
                                double dc,
                                double ds
                            );

void CoordinateDiscentMethod(   double initial_k,
                                double initial_k3,
                                double initial_k5,
                                double initial_c,
                                double initial_e00,
                                double initial_e02,
                                double initial_e04,
                                double initial_e06,
                                double initial_e08,
                                double initial_e10,
                                double initial_shift,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dk3,
                                double dk5,
                                double dc,
                                double de00,
                                double de02,
                                double de04,
                                double de06,
                                double de08,
                                double de10,
                                double dshift
                            );

void CoordinateDiscentMethodLinear(   double initial_k,
                                double initial_c,
                                double initial_e00,
                                double initial_e02,
                                double initial_e04,
                                double initial_e06,
                                double initial_e08,
                                double initial_e10,
                                double initial_shift,
                                int num_files,
                                std::string* files_list,
                                double* acc_list,
                                bool* trans,
                                double dk,
                                double dc,
                                double de00,
                                double de02,
                                double de04,
                                double de06,
                                double de08,
                                double de10,
                                double dshift
                            );
#endif // OPTIM_H_INCLUDED
