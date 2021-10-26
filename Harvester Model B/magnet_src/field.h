#ifndef FIELD_H_INCLUDED
#define FIELD_H_INCLUDED

class magnetic
{
    private:
        /** Function returnes the closest distance (in meters) between the "wire" BC and the observation point A */
        double dist(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
        /** Function returnes the cos(alpha1) the angle between vectors \vec{AB} and \vec{CB} */
        double cosine1(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
        /** Function returnes the cos(alpha1) the angle between vectors \vec{AC} and \vec{CB} */
        double cosine2(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
        /** Function findes the coordinates of the closest point of the line BC to the point A */
        void Crosspoint(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/);
        /** Function findes the direction of the magnet flux density in point A created by the current in line-segment BC */
        void nB(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/);
        /** Function calculates the magnet flux density produced by the block-magnet at the point */
        void B(double x_coord/**x-coordinate of the observation point*/, double y_coord/**y-coordinate of the observation point*/, double z_coord/**z-coordinate of the observation point*/, double magnet_current_dens/**current density of the equivalent coil*/, double* ans/**pointer on the array with answer[3] ( ans[0] = Bx, ans[1] = By, ans[2] = Bz )*/);
        /** Function returnes the flux through rount surface placed on the xc,yc,zc coordinate with radius*/
        double Flux(double xc/**x-coordinate of the center of the plate*/, double yc /**y-coordinate of the center of the plate*/, double zc /**y-coordinate of the center of the plate*/, double radius/**radius of the plate*/);

        /** Return the file with derivative dPhi/dz vs shift for z-axis for the loop with radius r and shift x*/
            /** File has the following structure:       */
            /** Flux_<radius>_<x-shift>.dat             */
            /** It contains the following information:  */
            /** <z_shift>   <dPhi/dz>                   */
        void GenerateFluxFileSingle(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step */, double r /** Radius of the loop */, double x /** x-shift */);

        /** Method use the file <fit.dat> to generate the electromotive force                                                           */
            /** <fit.dat> contains the dPhi/dz(z) = v*(b1 + b1*z^2 + b2*z^4 + b3*z^6 + b4*z^8 + b5*z^10) parameters for the each loop   */
            /** <fit.dat> has the following structure:                                                                                  */
            /** <x_shift>   <radius>    <b1>    <b2>    <b3>    <b4>    <b5>                                                            */
        double EMFFromFile(double z /** z-shift of the coil */, double v /** speed of the coil */);

        /** Method use the file <fit.dat> to generate the electromotive force                                                           */
            /** <fit.dat> contains the dPhi/dz(z) = v*(b1 + b1*z^2 + b2*z^4 + b3*z^6 + b4*z^8 + b5*z^10) parameters for the each loop   */
            /** <fit.dat> has the following structure:                                                                                  */
            /** <x_shift>   <radius>    <b1>    <b2>    <b3>    <b4>    <b5>                                                            */
        double EMForce(double zc /** z-shift of the coil */, double vc /** speed of the coil */);

        /** Interpolation function for the electromagnetic force:                                                                       */
            /** F(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                                 */
        double EMForceInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */);

        /** Interpolation function for the EMF:                                                                                         */
            /** EMF(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                               */
        double EMFInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */);

        /** Dimention right part of equation:           */
            /** if magn_force = True - closed circuit       */
            /** if magn_force = False - open circuit        */
        double RightPart_dv(double x /** displacement */, double dx /** velocity */, double t/** time */, double freq/** frequency */, bool magn_force/** open circuit or closed*/);

        /** Dimention right part of equation dx/dt:     */
        double RightPart_dx(double dx /** velocity */);

        double magnet_a1, magnet_a2, magnet_b1, magnet_b2, magnet_b3, magnet_i;
        double mu0;
        double coil_r_min, coil_r_max, coil_d, coil_resist;
        double emf_a, emf_b, emf_c, emf_d, emf_e, emf_f;
        double force_a, force_b, force_c, force_d, force_e, force_f;

        double RK45_step, RK45_exit;
        /** Spring constant (N/m)*/
        double springK;
        /** Nonlinear spring constant (N/m^2)*/
        double springK2;
        /** Nonlinear spring constant (N/m^2)*/
        double springK3;
        /** Nonlinear spring constant (N/m^2)*/
        double springK4;
        /** Nonlinear spring constant (N/m^2)*/
        double springK5;

        /** Oscillator's mass (kg)*/
        double mass;
        /** Dumping coefficient (kg/s)*/
        double DampingB;
        /** Amplitude of shaker*/
        double Aext;

        int magnet_n, flux_n, coil_n1, coil_n2, coil_force_n;

        double export_v_rms;
        double export_p_rms;
        double export_ampl;

    public:
        /** The constructor of the method   */
        magnetic();

        /** Interface method generates the set of files for each loop of the coil described in the paper            */
            /** usage is following:                                                                                 */
            /** magnetic mag();
                mag.GenerateFlux(0.0, 0.0021, 0.0001);
            */
            /** The method is designed to work with 8-thread processor.                                             */
        void GenerateFlux(double Zmin /** Minimal z-shift*/, double Zmax /** Maximal z-shift*/, double dZ /** z-shift step*/);

        /** Interface generates <Force.dat> file with self-containment controle                                                            */
            /** <z-shift>   <force>                                                                                                        */
            /** Usage:                                                                                                                     */
            /** magnetic mag;
                mag.GenerateForce(-0.002,0.0021,0.0001,0.01);
            */
        void GenerateForce(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step*/, double v /** velocity*/);

        /** Interface generates <EMF.dat> file with self-containment controle                                                            */
            /** <z-shift>   <EMF>                                                                                                        */
            /** Usage:                                                                                                                   */
            /** magnetic mag;
                mag.GenerateForce(-0.002,0.0021,0.0001,0.01);
            */
        void GenerateEMF(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step*/, double v /** velocity*/);

        /** Runge-Kutta method*/
        void RK4(double fext/** external frequency*/, bool magn_force/** switcher of EM force*/);

        /** RMS Voltage */
        double GetRMSV();
        /** RMS Power */
        double GetRMSP();
        double GetAmpl();
};

void GenerateResonance(double f_min, double f_max, double df);
#endif // FIELD_H_INCLUDED
