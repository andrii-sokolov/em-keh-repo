#include <iostream>
#include <fstream>
#include "INIReader.h"
#include "field.h"
#include "harmonic_ballance.h"
#include <thread>
#include <ctime>
#include <math.h>


using namespace std;

void mkrmf(string name)
{
    ifstream infile(name);
    ofstream ofs("rms_export.dat");
    string a0,a1,a2,a3,a4,a5,a6,a7;
    int i=0;
    const int window = 1000;
    double vrms, frms;
    while (infile>>a0>>a1>>a2)
    {
        frms += stod (a1);
        vrms += pow(stod (a2),2);
        i++;
        if (i>window)
        {
            ofs<<frms/window<<"\t"<<sqrt(vrms/window)<<endl;
            i=0;
            frms = 0;
            vrms = 0;
        }

    }
    ofs.close();
}

void mknorm(string name)
{
    ifstream infile(name);
    string a0,a1,a2,a3,a4,a5,a6,a7;
    int i=0;
    const int Size = 100000;
    double vrms[Size];
    double frms[Size];
    double max_v = 0;
    while (infile>>a0>>a1)
    {
        frms[i] = stod (a0);
        vrms[i] = stod (a1);
        if (vrms[i] > max_v)
        {
            max_v = vrms[i];
        }
        i++;
    }
    ofstream ofs("norm_rms_export.dat");
    for(int j = 0; j<i; j++)
    {
        ofs<<frms[j]<<"\t"<<vrms[j]/max_v<<endl;
    }
    ofs.close();
}

void LinearOptimization(double Mmin, double Mmax, double dM, double Bmin, double Bmax, double dB)
{
    ofstream ofs("Fitting.dat");
    double minsqr;
    double minimum = Compare_Resonance_Linear(Mmin, Bmin);
    for(double M = Mmin; M<Mmax; M+=dM)
    {
        for(double B = Bmin; B<Bmax; B+=dB)
        {
            minsqr = Compare_Resonance_Linear(M, B);
            ofs<<M<<"\t"<<B<<"\t"<<minsqr<<endl;
            //cout<<K<<"\t"<<B<<"\t"<<minsqr<<endl;
            if(minimum > minsqr)
            {
                minimum = minsqr;
                cout<<M<<"\t"<<B<<"\t"<<minsqr<<endl;
            }
        }
    }
    ofs.close();
}

void SqrOptimization(double K2min, double K2max, double dK2, double K, double B)
{
    ofstream ofs("FittingSQR.dat");
    double minsqr;
    double minimum = Compare_Resonance_SQR(K2min,K,B);
    for(double K2 = K2min; K2<K2max; K2+=dK2)
    {
        minsqr = Compare_Resonance_SQR(K2,K,B);
        ofs<<K2<<"\t"<<minsqr<<endl;
            //cout<<K<<"\t"<<B<<"\t"<<minsqr<<endl;
            if(minimum > minsqr)
            {
                minimum = minsqr;
                cout<<K2<<"\t"<<minsqr<<endl;
            }
    }
    ofs.close();
}

void CubeOptimization(double K3min, double K3max, double dK3, double K2min, double K2max, double dK2, double K, double B)
{
    ofstream ofs("FittingCube.dat");
    double minsqr;
    double minimum = Compare_Resonance_Cube(K3min,K2min,K,B);
    for(double K2 = K2min; K2<K2max; K2+= dK2)
    {
        for(double K3 = K3min; K3<K3max; K3+=dK3)
        {
            minsqr = Compare_Resonance_Cube(K3,K2,K,B);
            ofs<<K2<<"\t"<<K3<<"\t"<<minsqr<<endl;
                //cout<<K<<"\t"<<B<<"\t"<<minsqr<<endl;
                if(minimum > minsqr)
                {
                    minimum = minsqr;
                    cout<<K2<<"\t"<<K3<<"\t"<<minsqr<<endl;
                }
                //cout<<K3<<"\t"<<minsqr<<endl;
        }
    }
    ofs.close();
}

void FullOptimization(double K3min, double K3max, double dK3, double K, double B)
{
    ofstream ofs("FittingCube.dat");
    double minsqr;
    double minimum = Compare_Resonance_Full(K,0.0,K3min,0.0,0.0,B);
    for(double K3 = K3min; K3<K3max; K3+= dK3)
    {
            minsqr = Compare_Resonance_Full(K,0.0,K3,0.0,0.0,B);
            ofs<<K3<<"\t"<<minsqr<<endl;
                //cout<<K<<"\t"<<B<<"\t"<<minsqr<<endl;
                if(minimum > minsqr)
                {
                    minimum = minsqr;
                    cout<<K3<<"\t"<<minsqr<<"min"<<endl;
                }
                cout<<K3<<"\t"<<minsqr<<endl;
    }
    ofs.close();
}


int main()
{
//    cout<<Compare_Resonance(384.1, 1.0e6, 1.95e6, 3.0e8, 1.84e10, 0.00864744)<<endl;
    //LinearOptimization(0.0025, 0.0035, 0.00001, 0.0001, 0.008, 0.0001);
    //SqrOptimization(6950, 7050, 1.0, 384.682, 0.007857);
    //FullOptimization(0.0, 1000000, 1000, 384.682, 0.007857);
    //CubeOptimization(0.0, 1000000, 10000, 0.0, 100000,1000,  384.682, 0.007857);
    //GenerateResonance(55.0,65.01,0.01);
    //GenerateHArmonicResonance(52.0, 66.01, 0.01);
    //mknorm("RMSData/Freq_prototype1_0.3g_RMS.dat");
    //mkrmf("Data/Freq_prototype1_0.8g_Voc.dat");


    magnetic mag;
    mag.GenerateEMF(-0.002,0.0021,0.0001,1.0);
    mag.GenerateForce(-0.002,0.0021,0.0001,1.0);
    return 0;
}

