#ifndef HARMONIC_BALLANCE_H_INCLUDED
#define HARMONIC_BALLANCE_H_INCLUDED

class harmonic
{
    private:
        /** emf interpolation parameters*/
        double emf_a, emf_b, emf_c, emf_d, emf_e, emf_f;
        /** force interpolation parameters*/
        double force_a, force_b, force_c, force_d, force_e, force_f;
        /** Spring constant (N/m)*/
        double springK;
        /** Nonlinear spring constant (N/m^2)*/
        double springK2;
        /** Nonlinear spring constant (N/m^3)*/
        double springK3;
        /** Nonlinear spring constant (N/m^4)*/
        double springK4;
        /** Nonlinear spring constant (N/m^5)*/
        double springK5;
        /** Oscillator's mass (kg)*/
        double mass;
        /** Dumping coefficient (kg/s)*/
        double DampingB;
        /** Amplitude of shaker*/
        double Aext;

        double dA;
        double Amax;
        double delta;

        /** EMF interpolation function*/
        double EMFInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */);
        /** Electromagnetic Force interpolation function*/
        double EMForceInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */);
        /** Modelling function for the displacement*/
        double x(double A0/** amplitude*/, double wext/** external force frequency*/, double t/** time */);

        /** FFT coefficients*/
        double a1_sqr(double A0, double wext);
        double b1_sqr(double A0, double wext);

        double a1_cube(double A0, double wext);
        double b1_cube(double A0, double wext);

        double a1_quad(double A0, double wext);
        double b1_quad(double A0, double wext);

        double a1_five(double A0, double wext);
        double b1_five(double A0, double wext);

        double eqn(double A0, double wext);
    public:
        /** Empty constructor*/
        harmonic();
        /** Constructor for fitting*/
        harmonic(double m, double k2/** Nonlinear spring constant (N/m^2)*/, double k3/** Nonlinear spring constant (N/m^3)*/, double k4/** Nonlinear spring constant (N/m^4)*/, double k5/** Oscillator's mass (kg)*/, double B/** Dumping coefficient (kg/s)*/, double Ae);
        void Solve(double wext, int & n_of_roots, double* roots);
        double RMS_V(double A0, double wext);
        double Linear_A(double wext);
};

void GenerateHArmonicResonance(double fmin, double fmax, double df);
double Compare_Resonance_Full(double m, double k2, double k3, double k4, double k5, double B);
double Compare_Resonance_Linear(double m, double B);
double Compare_Resonance_SQR(double k2, double k, double B);
double Compare_Resonance_Cube(double k3, double k2, double k, double B);
#endif // HARMONIC_BALLANCE_H_INCLUDED
