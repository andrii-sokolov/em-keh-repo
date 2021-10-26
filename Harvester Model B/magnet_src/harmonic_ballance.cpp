#include "harmonic_ballance.h"
#include "INIReader.h"
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

harmonic::harmonic()
{
    INIReader reader("param.ini");
    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load 'param.ini'\n";
    }

    emf_a = reader.GetReal("Intermolation parameters", "emf_a", -1);
    emf_b = reader.GetReal("Intermolation parameters", "emf_b", -1);
    emf_c = reader.GetReal("Intermolation parameters", "emf_c", -1);
    emf_d = reader.GetReal("Intermolation parameters", "emf_d", -1);
    emf_e = reader.GetReal("Intermolation parameters", "emf_e", -1);
    emf_f = reader.GetReal("Intermolation parameters", "emf_f", -1);

    force_a = reader.GetReal("Intermolation parameters", "force_a", -1);
    force_b = reader.GetReal("Intermolation parameters", "force_b", -1);
    force_c = reader.GetReal("Intermolation parameters", "force_c", -1);
    force_d = reader.GetReal("Intermolation parameters", "force_d", -1);
    force_e = reader.GetReal("Intermolation parameters", "force_e", -1);
    force_f = reader.GetReal("Intermolation parameters", "force_f", -1);

    springK = reader.GetReal("Dynamics","k",-1);
    springK2= reader.GetReal("Dynamics","k2",-1);
    springK3= reader.GetReal("Dynamics","k3",-1);
    springK4= reader.GetReal("Dynamics","k4",-1);
    springK5= reader.GetReal("Dynamics","k5",-1);
    mass = reader.GetReal("Dynamics","mass",-1);
    DampingB = reader.GetReal("Dynamics","c",-1);
    Aext = reader.GetReal("External force","Aext",-1);

    dA = 1e-5;
    Amax = 1e-2;
    delta = 1e-15;
}

harmonic::harmonic(double m, double k2, double k3, double k4, double k5, double B, double Ae)
{
    INIReader reader("param.ini");
    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load 'param.ini'\n";
    }

    emf_a = reader.GetReal("Intermolation parameters", "emf_a", -1);
    emf_b = reader.GetReal("Intermolation parameters", "emf_b", -1);
    emf_c = reader.GetReal("Intermolation parameters", "emf_c", -1);
    emf_d = reader.GetReal("Intermolation parameters", "emf_d", -1);
    emf_e = reader.GetReal("Intermolation parameters", "emf_e", -1);
    emf_f = reader.GetReal("Intermolation parameters", "emf_f", -1);

    force_a = reader.GetReal("Intermolation parameters", "force_a", -1);
    force_b = reader.GetReal("Intermolation parameters", "force_b", -1);
    force_c = reader.GetReal("Intermolation parameters", "force_c", -1);
    force_d = reader.GetReal("Intermolation parameters", "force_d", -1);
    force_e = reader.GetReal("Intermolation parameters", "force_e", -1);
    force_f = reader.GetReal("Intermolation parameters", "force_f", -1);

    springK = reader.GetReal("Dynamics","k",-1);
    springK2= k2;
    springK3= k3;
    springK4= k4;
    springK5= k5;
    mass = m;
    DampingB = B;
    Aext = Ae;

    dA = 1e-5;
    Amax = 1e-2;
    delta = 1e-15;
}

double harmonic::EMFInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */)
/** Interpolation function for the EMF:                                                                                       */
/** EMF(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                                 */
{
    return(vc*(emf_a + emf_b*zc*zc + emf_c*pow(zc,4) + emf_d*pow(zc,6) + emf_e*pow(zc,8) + emf_f*pow(zc,10)));
}

double harmonic::EMForceInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */)
/** Interpolation function for the electromagnetic force:                                                                   */
/** F(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                                 */
{
    return(vc*(force_a + force_b*zc*zc + force_c*pow(zc,4) + force_d*pow(zc,6) + force_e*pow(zc,8) + force_f*pow(zc,10)));
}

double harmonic::x(double A0/** amplitude*/,  double wext/** external force frequency*/, double t/** time */)
{
    return(A0*cos(wext*t));
}

double harmonic::a1_cube(double A0, double wext)
{
    double t1 = 0.0;
    double t3 = 2.0*M_PI/wext;

    double Int = 2.0*springK3*(pow(A0,3)*(12*(-t1 + t3)*wext - 8*sin(2*t1*wext) - sin(4*t1*wext) + 8*sin(2*t3*wext) + sin(4*t3*wext)))/(32.*wext)/(t3-t1);
    return (Int);
}

double harmonic::b1_cube(double A0, double wext)
{
    double t1 = 0.0;
    double t3 = 2.0*M_PI/wext;

    double Int = 2.0*springK3*(pow(A0,3)*(pow(cos(t1*wext),4) - pow(cos(t3*wext),4)))/(4.*wext)/(t3-t1);
    return (Int);
}

double harmonic::a1_five(double A0, double wext)
{
    double t1 = 0.0;
    double t3 = 2.0*M_PI/wext;

    double Int = 2.0*springK5*((pow(A0,5)*(60*(-t1 + t3)*wext - 45*sin(2*t1*wext) - 9*sin(4*t1*wext) - sin(6*t1*wext) + 45*sin(2*t3*wext) + 9*sin(4*t3*wext) +
       sin(6*t3*wext)))/(192.*wext))/(t3-t1);
    return (Int);
}

double harmonic::b1_five(double A0, double wext)
{
    double t1 = 0.0;
    double t3 = 2.0*M_PI/wext;

    double Int = 2.0*springK5*((pow(A0,5)*(pow(cos(t1*wext),6) - pow(cos(t3*wext),6)))/(6.*wext))/(t3-t1);
    return (Int);
}

double harmonic::a1_sqr(double A0, double wext)
{
    double Int = 0;

    double t1 = -0.5*M_PI/wext;
    double t2 = 0.5*M_PI/wext;
    double t3 = 1.5*M_PI/wext;

    Int+= (pow(A0,2)*(-9*sin(t1*wext) - sin(3*t1*wext) + 9*sin(t2*wext) + sin(3*t2*wext)))/(12.*wext);
    Int+= -(pow(A0,2)*(-9*sin(t2*wext) - sin(3*t2*wext) + 9*sin(t3*wext) + sin(3*t3*wext)))/(12.*wext);

    Int = 2.0*springK2*Int/(t3-t1);

    return(Int);
}

double harmonic::b1_sqr(double A0, double wext)
{
    double Int = 0;

    double t1 = -0.5*M_PI/wext;
    double t2 = 0.5*M_PI/wext;
    double t3 = 1.5*M_PI/wext;

    Int+= (pow(A0,2)*(pow(cos(t1*wext),3) - pow(cos(t2*wext),3)))/(3.*wext);
    Int+= -(pow(A0,2)*(pow(cos(t2*wext),3) - pow(cos(t3*wext),3)))/(3.*wext);

    Int = 2.0*springK2*Int/(t3-t1);

    return(Int);
}

double harmonic::a1_quad(double A0, double wext)
{
    double Int = 0;

    double t1 = -0.5*M_PI/wext;
    double t2 = 0.5*M_PI/wext;
    double t3 = 1.5*M_PI/wext;

    Int+= (pow(A0,4)*(-150*sin(t1*wext) - 25*sin(3*t1*wext) - 3*sin(5*t1*wext) + 25*(6*sin(t2*wext) + sin(3*t2*wext)) + 3*sin(5*t2*wext)))/(240.*wext);
    Int+= -(pow(A0,4)*(-150*sin(t2*wext) - 25*sin(3*t2*wext) - 3*sin(5*t2*wext) + 25*(6*sin(t3*wext) + sin(3*t3*wext)) + 3*sin(5*t3*wext)))/(240.*wext);

    Int = 2.0*springK4*Int/(t3-t1);

    return(Int);
}

double harmonic::b1_quad(double A0, double wext)
{
    double Int = 0;

    double t1 = -0.5*M_PI/wext;
    double t2 = 0.5*M_PI/wext;
    double t3 = 1.5*M_PI/wext;

    Int+= (pow(A0,4)*(pow(cos(t1*wext),5) - pow(cos(t2*wext),5)))/(5.*wext);
    Int+= -(pow(A0,4)*(pow(cos(t2*wext),5) - pow(cos(t3*wext),5)))/(5.*wext);

    Int = 2.0*springK4*Int/(t3-t1);

    return(Int);
}

double harmonic::eqn(double A0, double wext)
{
    double a1sqr = a1_sqr(A0,wext);
    double a1cube = a1_cube(A0,wext);
    double a1quad = a1_quad(A0,wext);
    double a1five = a1_five(A0,wext);
    double b1sqr = b1_sqr(A0,wext);
    double b1cube = b1_cube(A0,wext);
    double b1quad = b1_quad(A0,wext);
    double b1five = b1_five(A0,wext);
    return(-(pow(Aext,2)*pow(mass,2)) + pow(b1cube +  b1five + b1quad + b1sqr - A0*DampingB*wext,2) + pow(a1cube +  a1five + a1quad + a1sqr + A0*springK - A0*mass*pow(wext,2),2));
}

double harmonic::RMS_V(double A0, double wext)
{
    double Int = pow(A0,2)*pow(emf_a,2)*M_PI*wext + (pow(A0,4)*emf_a*
            (4096*emf_b + 2048*pow(A0,2)*emf_c + 429*pow(A0,12)*emf_d*emf_e + 672*pow(A0,8)*emf_f)*M_PI*wext)/8192. +
            (pow(A0,6)*(4194304*pow(emf_b,2) + 1835008*pow(A0,4)*pow(emf_c,2) +
            73216*pow(A0,10)*emf_c*(17*pow(A0,4)*emf_d*emf_e + 24*emf_f) +
            2048*emf_b*(2560*pow(A0,2)*emf_c + 715*pow(A0,12)*emf_d*emf_e + 1056*pow(A0,8)*emf_f) +
            323*pow(A0,16)*(1035*pow(A0,8)*pow(emf_d,2)*pow(emf_e,2) + 2576*pow(A0,4)*emf_d*emf_e*emf_f + 1664*pow(emf_f,2)))*
            M_PI*wext)/3.3554432e7;
    return(sqrt(Int*wext/(2*M_PI)));
}

void harmonic::Solve(double fext, int & n_of_roots, double* roots)
{
    n_of_roots =0;
    for(double A = dA; A <  Amax; A+=dA)
    {
        if(eqn(A,2*M_PI*fext)*eqn(A+dA,2*M_PI*fext)<=0)
        {
            double A1 = A;
            double A2 = A+dA;
            while(abs(A2-A1)>delta)
            {
                if((eqn(A1,2*M_PI*fext)*eqn(0.5*(A1+A2),2*M_PI*fext))<0)
                    A2 = 0.5*(A1+A2);
                else
                    A1 = 0.5*(A1+A2);
            }
            //cout<<0.5*(A1+A2)<<"\t"<<RMS_V(0.5*(A1+A2),2*M_PI*fext)<<endl;
            roots[n_of_roots] = 0.5*(A1+A2);
            n_of_roots++;
        }
    }

}

double harmonic::Linear_A(double wext)
{
    return(Aext/sqrt(pow(springK/mass - wext*wext,2) + pow(DampingB*wext/mass,2)));
}

void GenerateHArmonicResonance(double fmin, double fmax, double df)
{

    double volts[100000];
    double freqs[100000];
    int counter = 0;

    double max_v = 0.0;

    for(double f = fmin; f<fmax; f+=df)
    {
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms;
        harms.Solve(f, n_of_roots, roots);
        for(int i=0; i<n_of_roots; i++)
        {
            volts[counter] = harms.RMS_V(roots[i],2*M_PI*f);
            freqs[counter] = f;

            if(volts[counter] > max_v)
            {
                max_v = volts[counter];
            }
            counter++;
        }
        delete roots;
    }
    cout<<freqs[0]<<"\t"<<freqs[counter-1]<<endl;
    ofstream ofs("harmonic_resonance.dat");
    for(int i=0; i<counter; i++)
        ofs<<freqs[i]<<"\t"<<volts[i]<<endl;
    ofs.close();
}

double Compare_Resonance_Full(double m, double k2, double k3, double k4, double k5, double B)
{
    double ans = 0.0;
    /** Import A_ext = 1 m/s^2 */
    ifstream infile("RMSData/Freq_prototype1_0.1g_RMS.dat");
    string a0,a1;
    double freq,volt;
    while (infile>>a0>>a1)
    {
        freq = stod(a0);
        volt = stod(a1);
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms(m,k2,k3,k4,k5,B,1.0);
        harms.Solve(freq,n_of_roots,roots);
        if(n_of_roots == 1)
            ans += pow((volt - harms.RMS_V(roots[0],2*M_PI*freq)),2);
        else
        {
            double minim = pow(volt -  harms.RMS_V(roots[0],2*M_PI*freq),2);
            for(int j=1; j<n_of_roots; j++)
            {
                if(pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2)<minim)
                    minim = pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2);
            }
            ans += minim;
        }
        delete roots;
    }

    ifstream infile2("RMSData/Freq_prototype1_0.3g_RMS.dat");
    string a2,a3;
    while (infile2>>a2>>a3)
    {
        freq = stod(a2);
        volt = stod(a3);
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms(m,k2,k3,k4,k5,B,3.0);
        harms.Solve(freq,n_of_roots,roots);
        if(n_of_roots == 1)
            ans += pow((volt - harms.RMS_V(roots[0],2*M_PI*freq)),2);
        else
        {
            double minim = pow(volt - harms.RMS_V(roots[0],2*M_PI*freq),2);
            for(int j=1; j<n_of_roots; j++)
            {
                if(pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2)<minim)
                    minim = pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2);
            }
            ans += minim;
        }
        delete roots;
    }

    ifstream infile3("RMSData/Freq_prototype1_0.5g_RMS.dat");
    string a4,a5;
    while (infile3>>a4>>a5)
    {
        freq = stod(a4);
        volt = stod(a5);
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms(m,k2,k3,k4,k5,B,5.0);
        harms.Solve(freq,n_of_roots,roots);
        if(n_of_roots == 1)
            ans += pow((volt - harms.RMS_V(roots[0],2*M_PI*freq)),2);
        else
        {
            double minim = pow(volt - harms.RMS_V(roots[0],2*M_PI*freq),2);
            for(int j=1; j<n_of_roots; j++)
            {
                if(pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2)<minim)
                    minim = pow(volt - harms.RMS_V(roots[j],freq),2);
            }
            ans += minim;
        }
        delete roots;
    }

    ifstream infile4("RMSData/Freq_prototype1_0.8g_RMS.dat");
    string a6,a7;
    while (infile4>>a6>>a7)
    {
        freq = stod(a6);
        volt = stod(a7);
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms(m,k2,k3,k4,k5,B,8.0);
        harms.Solve(freq,n_of_roots,roots);
        if(n_of_roots == 1)
            ans += pow((volt - harms.RMS_V(roots[0],2*M_PI*freq)),2);
        else
        {
            double minim = pow(volt - harms.RMS_V(roots[0],2*M_PI*freq),2);
            for(int j=1; j<n_of_roots; j++)
            {
                if(pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2)<minim)
                    minim = pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2);
            }
            ans += minim;
        }
        delete roots;
    }
    ifstream infile5("RMSData/Freq_prototype1_1g_RMS.dat");
    string a8,a9;
    while (infile5>>a8>>a9)
    {
        freq = stod(a8);
        volt = stod(a9);
        int n_of_roots;
        double* roots = new double[50];
        harmonic harms(m,k2,k3,k4,k5,B,10.0);
        harms.Solve(freq,n_of_roots,roots);
        if(n_of_roots == 1)
            ans += pow((volt - harms.RMS_V(roots[0],2*M_PI*freq)),2);
        else
        {
            double minim = pow(volt - harms.RMS_V(roots[0],2*M_PI*freq),2);
            for(int j=1; j<n_of_roots; j++)
            {
                if(pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2)<minim)
                    minim = pow(volt - harms.RMS_V(roots[j],2*M_PI*freq),2);
            }
            ans += minim;
        }
        delete roots;
    }


    return(ans);

}

double Compare_Resonance_Linear(double m, double B)
{
    double ans = 0.0;
    double freq,volt;

    ifstream infile2("RMSLinearData/Freq_prototype1_0.3g_RMS.dat");
    string a2,a3;
    while (infile2>>a2>>a3)
    {
        freq = stod(a2);
        volt = stod(a3);
        harmonic harms(m,0.0,0.0,0.0,0.0,B,3.0);
        double root = harms.Linear_A(2*M_PI*freq);
        ans += pow((volt - harms.RMS_V(root,2*M_PI*freq)),2);
    }

    ifstream infile3("RMSLinearData/Freq_prototype1_0.5g_RMS.dat");
    string a4,a5;
    while (infile3>>a4>>a5)
    {
        freq = stod(a4);
        volt = stod(a5);
        harmonic harms(m,0.0,0.0,0.0,0.0,B,5.0);
        double root = harms.Linear_A(2*M_PI*freq);
        ans += pow((volt - harms.RMS_V(root,2*M_PI*freq)),2);
    }

    ifstream infile4("RMSLinearData/Freq_prototype1_0.8g_RMS.dat");
    string a6,a7;
    while (infile4>>a6>>a7)
    {
        freq = stod(a6);
        volt = stod(a7);
        harmonic harms(m,0.0,0.0,0.0,0.0,B,8.0);
        double root = harms.Linear_A(2*M_PI*freq);
        ans += pow((volt - harms.RMS_V(root,2*M_PI*freq)),2);
    }
    ifstream infile5("RMSLinearData/Freq_prototype1_1g_RMS.dat");
    string a8,a9;
    while (infile5>>a8>>a9)
    {
        freq = stod(a8);
        volt = stod(a9);
        harmonic harms(m,0.0,0.0,0.0,0.0,B,10.0);
        double root = harms.Linear_A(2*M_PI*freq);
        ans += pow((volt - harms.RMS_V(root,2*M_PI*freq)),2);
    }


    return(ans);

}

double Compare_Resonance_SQR(double k2, double k, double B)
{
    double ans = 0.0;
    /** Import A_ext = 1 m/s^2 */
    ifstream infile("RMSData/Freq_prototype1_0.3g_RMS.dat");
    string a0,a1;
    const int Size = 1000;
    double freq[Size];
    double frequency, voltage;
    double volt[Size];
    double theor[Size];
    int i = 0;
    double max_volt = 0.0;
    double max_theor = 0.0;
    while (infile>>a0>>a1)
    {
        harmonic harms(k,k2,0.0,0.0,0.0,B,3.0);
        frequency = stod(a0);
        voltage = stod(a1);
        double roots[50];
        int nroots;
        harms.Solve(frequency,nroots,roots);
        for(int j=0; j<nroots; j++)
        {
            freq[i] = frequency;
            volt[i] = voltage;
            theor[i] = harms.RMS_V(roots[j],2*M_PI*frequency);
            if(volt[i]>max_volt)
                max_volt = volt[i];
            if(theor[i]>max_theor)
                max_theor = theor[i];
            i++;
        }
    }
    double i_max = i;
    for(i =0; i<i_max; i++)
        ans += pow((volt[i] - theor[i]),2);
    return(ans);
}

double Compare_Resonance_Cube(double k3, double k2, double k, double B)
{
    double ans = 0.0;
    /** Import A_ext = 1 m/s^2 */
    ifstream infile("RMSData/Freq_prototype1_0.3g_RMS.dat");
    string a0,a1;
    const int Size = 1000;
    double freq[Size];
    double frequency, voltage;
    double volt[Size];
    double theor[Size];
    int i = 0;
    double max_volt = 0.0;
    double max_theor = 0.0;
    while (infile>>a0>>a1)
    {
        harmonic harms(k,k2,k3,0.0,0.0,B,3.0);
        frequency = stod(a0);
        voltage = stod(a1);
        double roots[50];
        int nroots;
        harms.Solve(frequency,nroots,roots);
        for(int j=0; j<nroots; j++)
        {
            freq[i] = frequency;
            volt[i] = voltage;
            theor[i] = harms.RMS_V(roots[j],2*M_PI*frequency);
            if(volt[i]>max_volt)
                max_volt = volt[i];
            if(theor[i]>max_theor)
                max_theor = theor[i];
            i++;
        }
    }
    double i_max = i;
    for(i =0; i<i_max; i++)
        ans += pow((volt[i] - theor[i]),2);
    return(ans);
}
