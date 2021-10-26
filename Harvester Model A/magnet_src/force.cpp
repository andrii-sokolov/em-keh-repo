#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "forces.h"
#include "INIReader.h"

using namespace std;

double dist(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
/** Function returnes the closest distance (in meters) between the "wire" BC and the observation point A
*/
{
    return(sqrt((x0 * x0 * y1 * y1 - 2.0 * x0 * x0 * y1 * y2 + x0 * x0 * y2 * y2 + y1 * y1 * z0 * z0 - 2.0 * y1 * y2 * z0 * z0 + y2 * y2 * z0 * z0 + x2 * x2 * (y0 * y0 - 2.0 * y0 * y1 + y1 * y1 + pow((z0 - z1),2.0)) - 2.0 * y0 * y1 * z0 * z1 + 2.0 * y0 * y2 * z0 * z1 + 2.0 * y1 * y2 * z0 * z1 - 2.0 * y2 * y2 * z0 * z1 + x0 * x0 * z1 * z1 + y0 * y0 * z1 * z1 - 2.0 * y0 * y2 * z1 * z1 + y2 * y2 * z1 * z1 + x1 * x1 * (y0 * y0 - 2.0 * y0 * y2 + y2 * y2 + pow((z0 - z2),2.0)) - 2.0 * x0 * x2 * (y1 * y1 - y1 * y2 + y0 * (-y1 + y2) - (z0 - z1) * (z1 - z2)) + 2.0 * y0 * y1 * z0 * z2 - 2.0 * y1 * y1 * z0 * z2 - 2.0 * y0 * y2 * z0 * z2 + 2.0 * y1 * y2 * z0 * z2 - 2.0 * x0 * x0 * z1 * z2 - 2.0 * y0 * y0 * z1 * z2 + 2.0 * y0 * y1 * z1 * z2 + 2.0 * y0 * y2 * z1 * z2 - 2.0 * y1 * y2 * z1 * z2 + x0 * x0 * z2 * z2 + y0 * y0 * z2 * z2 - 2.0 * y0 * y1 * z2 * z2 + y1 * y1 * z2 * z2 - 2.0 * x1 * (x2 * (y0 * y0 + y1 * y2 - y0 * (y1 + y2) + (z0 - z1) * (z0 - z2)) + x0 * (y0 * (y1 - y2) - y1 * y2 + y2 * y2 + z0 * z1 - z0 * z2 - z1 * z2 + z2 * z2)))/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2)));
}

double cosine1(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
/** Function returnes the cos(alpha1) the angle between vectors \vec{AB} and \vec{CB}
*/
{
    double vec1[3] = {x1-x0,y1-y0,z1-z0};
    double vec2[3] = {x1-x2,y1-y2,z1-z2};
    double vecmult = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    double modv1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
    double modv2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
    return(vecmult/(modv1*modv2));
}

double cosine2(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
/** Function returnes the cos(alpha1) the angle between vectors \vec{AC} and \vec{CB}
*/
{
    double vec1[3] = {x2-x0,y2-y0,z2-z0};
    double vec2[3] = {x1-x2,y1-y2,z1-z2};
    double vecmult = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
    double modv1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
    double modv2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);
    return(vecmult/(modv1*modv2));
}

void Crosspoint(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/)
/** Function findes the coordinates of the closest point of the line BC to the point A
*/
{
    ans[0] = -((-x0 * x1 * x1 + 2.0 * x0 * x1 * x2 - x0 * x2 * x2 - x1 * y0 * y1 + x2 * y0 * y1 - x2 * y1 * y1 + x1 * y0 * y2 - x2 * y0 * y2 + x1 * y1 * y2 + x2 * y1 * y2 - x1 * y2 * y2 - x1 * z0 * z1 + x2 * z0 * z1 - x2 * z1 * z1 + x1 * z0 * z2 - x2 * z0 * z2 + x1 * z1 * z2 + x2 * z1 * z2 - x1 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
    ans[1] = -((-x0 * x1 * y1 + x0 * x2 * y1 + x1 * x2 * y1 - x2 * x2 * y1 - y0 * y1 * y1 + x0 * x1 * y2 - x1 * x1 * y2 - x0 * x2 * y2 + x1 * x2 * y2 + 2.0 * y0 * y1 * y2 - y0 * y2 * y2 - y1 * z0 * z1 + y2 * z0 * z1 - y2 * z1 * z1 + y1 * z0 * z2 - y2 * z0 * z2 + y1 * z1 * z2 + y2 * z1 * z2 - y1 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
    ans[2] = -((-x0 * x1 * z1 + x0 * x2 * z1 + x1 * x2 * z1 - x2 * x2 * z1 - y0 * y1 * z1 + y0 * y2 * z1 + y1 * y2 * z1 - y2 * y2 * z1 - z0 * z1 * z1 + x0 * x1 * z2 - x1 * x1 * z2 - x0 * x2 * z2 + x1 * x2 * z2 + y0 * y1 * z2 - y1 * y1 * z2 - y0 * y2 * z2 + y1 * y2 * z2 + 2.0 * z0 * z1 * z2 - z0 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
}

void nB(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/)
/** Function findes the direction of the magnet flux density in point A created by the current in line-segment BC
*/
{
    double* D_point = new double (3);
    Crosspoint(x0,y0,z0,x1,y1,z1,x2,y2,z2,D_point);
    double x3 = D_point[0];
    double y3 = D_point[1];
    double z3 = D_point[2];
    if((sqrt(pow((x2-x1),2.0)+pow((y2-y1),2.0)+pow((z2-z1),2.0))!=0.0) & (sqrt(pow((x0-x3),2.0)+pow((y0-y3),2.0)+pow((z0-z3),2.0))!=0.0))
    {
        double vecBC [3] = {(x2-x1)/sqrt(pow((x2-x1),2.0)+pow((y2-y1),2.0)+pow((z2-z1),2.0)), (y2-y1)/sqrt(pow((x2-x1),2.0)+pow((y2-y1),2.0)+pow((z2-z1),2.0)), (z2-z1)/sqrt(pow((x2-x1),2.0)+pow((y2-y1),2.0)+pow((z2-z1),2.0))};
        double vecDA [3] = {(x0-x3)/sqrt(pow((x0-x3),2.0)+pow((y0-y3),2.0)+pow((z0-z3),2.0)), (y0-y3)/sqrt(pow((x0-x3),2.0)+pow((y0-y3),2.0)+pow((z0-z3),2.0)), (z0-z3)/sqrt(pow((x0-x3),2.0)+pow((y0-y3),2.0)+pow((z0-z3),2.0))};
        ans[0] = vecBC[1]*vecDA[2] - vecBC[2]*vecDA[1];
        ans[1] = vecBC[2]*vecDA[0] - vecBC[0]*vecDA[2];
        ans[2] = vecBC[0]*vecDA[1] - vecBC[1]*vecDA[0];
    }
    else
    {
        ans[0] = 0.0;
        ans[1] = 0.0;
        ans[2] = 0.0;
    }
    delete D_point;
}

void magnet::B(double x_coord/**x-coordinate of the observation point*/, double y_coord/**y-coordinate of the observation point*/, double z_coord/**z-coordinate of the observation point*/, double* ans/**pointer on the array with answer[3] ( ans[0] = Bx, ans[1] = By, ans[2] = Bz )*/)
/** Function calculates the magnet flux density produced by the block-magnet at the point */
{
    double Bind [3] = {0.0,0.0,0.0};
    double* norm = new double [3];
    double val_B;
    double z, x1,x2,y1,y2;
    for(int i=0; i<magnet_n + 1; i++)
    {
        z = b*i/magnet_n;
        x1 =  a*0.5;
        y1 = -a*0.5;
        x2 =  a*0.5;
        y2 =  a*0.5;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z, norm);
        val_B = (mu0*(Curr*b/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        x1 = a*0.5;
        y1 = a*0.5;
        x2 = -a*0.5;
        y2 = a*0.5;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z, norm);
        val_B = (mu0*(Curr*b/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        x1 = -a*0.5;
        y1 =  a*0.5;
        x2 = -a*0.5;
        y2 = -a*0.5;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z, norm);
        val_B = (mu0*(Curr*b/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        x1 = -a*0.5;
        y1 = -a*0.5;
        x2 =  a*0.5;
        y2 = -a*0.5;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z, norm);
        val_B = (mu0*(Curr*b/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z , x2 , y2 , z));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
    }
    ans[0] = Bind[0];
    ans[1] = Bind[1];
    ans[2] = Bind[2];
    delete norm;
}

double magnet::Flux_single(double c, double z, double alpha_x)
{
    vector<vector<double> > points(flux_n*flux_n);
    for ( int i = 0 ; i < (flux_n*flux_n) ; i++ )
        points[i].resize(3);

    double dc = c/flux_n;
    int k=0;
    double temp_x, temp_y,temp_z;
    for(int i = 0; i <flux_n; i++)
    {
        for(int j = 0; j<flux_n; j++ )
        {
            temp_x = -0.5*c + dc*(0.5 + i);
            temp_y = -0.5*c + dc*(0.5 + j);
            temp_z = z;
            points[k][0] = temp_x;
            points[k][1] = temp_y*cos(alpha_x) - temp_z*sin(alpha_x);
            points[k][2] = temp_y*sin(alpha_x) + temp_z*cos(alpha_x);
            k++;
        }
    }
    //for ( int i = 0 ; i < (flux_n*flux_n) ; i++ )
    //    cout<<points[i][0]<<"\t"<<points[i][1]<<"\t"<<points[i][2]<<endl;
    double n_x = 0.0;
    double n_y = -sin(alpha_x);
    double n_z = cos(alpha_x);
    double Flux = 0.0;
    double* Induction = new double [3];
    for ( int i = 0 ; i < (flux_n*flux_n) ; i++ )
    {
        B(points[i][0],points[i][1],points[i][2],Induction);
        Flux += (Induction[0]*n_x + Induction[1]*n_y + Induction[2]*n_z)*pow(dc,2);
    }
    delete Induction;
    return(Flux);
}

void magnet::GenerateDerivFlux_rot(double alpha_0, double alpha_1, double d_alpha)
{
    ofstream file("magnet_src/flux_derivative_rot.dat");
    double da = 0.5*M_PI*1e-3;
    double F1,F2,dF,alpha,c;
    int i,j;
    int imax = (int)((alpha_1-alpha_0)/d_alpha);
    for( i=0; i<imax; i++)
    {
        alpha = alpha_0 + i*d_alpha;
        for( j=0; j<coil_n + 1; j++)
        {
            c = c1 + (c2-c1)/coil_n*j;
            F2 = Flux_single(c, d, alpha+0.5*da);
            F1 = Flux_single(c, d, alpha-0.5*da);
            dF = (F2 - F1)/da;
            file<<alpha<<"\t"<<c<<"\t"<<d<<"\t"<<dF<<endl;
            cout<<alpha<<"\t"<<c<<"\t"<<d<<"\t"<<dF<<endl;
            F2 = Flux_single(c, d + Dl, alpha+0.5*da);
            F1 = Flux_single(c, d + Dl, alpha-0.5*da);
            dF = (F2-F1)/da;
            file<<alpha<<"\t"<<c<<"\t"<<d + Dl<<"\t"<<dF<<endl;
            cout<<alpha<<"\t"<<c<<"\t"<<d + Dl<<"\t"<<dF<<endl;
        }
    }
    file.close();
}

void magnet::GenerateDerivFlux_transp(double z_1, double z_2, double dz)
{
    ofstream file("magnet_src/flux_derivative.dat");
    double dc = (c2-c1)/coil_n;
    double F;
    for(double z = z_1; z<z_2; z+=dz)
    {
        for(double c = c1; c<c2; c+=dc)
        {
            F = (Flux_single(c, z + 0.5*dz, 0) - Flux_single(c, z - 0.5*dz, 0))/dz;
            file<<c<<"\t"<<z<<"\t"<<F<<endl;
            cout<<c<<"\t"<<z<<"\t"<<F<<endl;
            F = (Flux_single(c, z + Dl + 0.5*dz, 0) - Flux_single(c, z + Dl - 0.5*dz, 0))/dz;
            file<<c<<"\t"<<z + Dl<<"\t"<<F<<endl;
            cout<<c<<"\t"<<z + Dl<<"\t"<<F<<endl;
        }
    }
    file.close();
}

void magnet::GenerateEMF_transp(double v)
{
    ifstream infile("magnet_src/flux_derivative.dat");
    string file_a,file_b,file_c;
    int i = 0;
    double dPhi[100000];
//    double radius[100000];
    double z_i[100000];
    while (infile>>file_a>>file_b>>file_c)
    {
//        radius[i] = stod (file_a);
        z_i[i] = stod (file_b);
        dPhi[i] = stod (file_c);
        i++;
    }
    int i_max = i;
    ofstream file("magnet_src/emf.dat");
    for( int j = 145; j<i_max; j+=146)
    {
        double emf = 0.0;
        for(int k = j-145; k<j+1; k++)
        {
            emf+=dPhi[k]*v;
        }
        file<<z_i[j-1]<<"\t"<<emf<<endl;
    }
    file.close();
}

void magnet::GenerateEMF_rot(double omega)
{
    ifstream infile("magnet_src/flux_derivative_rot.dat");
    string file_a,file_b,file_c,file_d;
    int i = 0;
    double alpha[100000];
    double dPhi[100000];
    double radius[100000];
    double z_i[100000];
    while (infile>>file_a>>file_b>>file_c>>file_d)
    {
        alpha[i] = stod (file_a);
        radius[i] = stod (file_b);
        z_i[i] = stod (file_c);
        dPhi[i] = stod (file_d);
        i++;
    }
    int i_max = i;
    ofstream file("magnet_src/emf_rot.dat");
    for( int j = 145; j<i_max; j+=146)
    {
        double emf = 0.0;
        for(int k = j-145; k<j+1; k++)
        {
            emf+=dPhi[k]*omega;
        }
        file<<alpha[j-1]<<"\t"<<emf<<"\t"<<emf_interp_rot(omega, alpha[j-1])<<endl;
    }
    file.close();
}


double magnet::emf_interp_transp(double v, double z)
{
    return(v*(-0.8540921267394467 + 1150.9278844217931*z - 483435.9688260379*pow(z,2) + 6.2248686795744255e7*pow(z,3) - 3.8224484114907404e11*pow(z,5))/
   (4.4838175041083 - 8473.454048648566*z + 6.604677180107996e6*pow(z,2) - 2.5308414452052097e9*pow(z,3) + 3.93246741422158e11*pow(z,4)));
}

double magnet::emf_interp_rot(double omega, double alpha)
{
    return(omega*(-0.00012035207598931708*alpha - 4.098063691078718e-11*pow(alpha,2) - 0.0009464994355293829*pow(alpha,3) + 0.00407300408196066*pow(alpha,5) - 0.007116892648305839*pow(alpha,7)));
}

double magnet::Force_progres_interp(double omega, double alpha)
{
    return(omega*((-0.001463271490329916 + 2.885580641532994*alpha - 1985.9480976001541*pow(alpha,2) + 612504.5449322374*pow(alpha,3) - 8.910000738487323e7*pow(alpha,4) +
     5.02599104322372e9*pow(alpha,5))/
   (2.739666627358009 - 6.010763843017788e6*pow(alpha,2) + 6.02102325648272e9*pow(alpha,3) - 2.3905294829962466e12*pow(alpha,4) + 3.55783018480984e14*pow(alpha,5))));
}

double magnet::Force_rot_interp(double v, double z)
{
    return(v*(-1.8461832348701756e-17 - 7.541357584924355e-11*pow(z,2) - 1.1921442114522662e-9*pow(z,4) + 8.700445322549457e-10*pow(z,6) + 1.7777384605376496e-8*pow(z,8)));
}

double magnet::Force_transp(double v, double z)
{
    double dc = (c2-c1)/coil_n;
    double e = emf_interp_transp(v,z);
    double* Induction = new double[3];
    double dlx,dly,dlz,dFx,dFy,dFz, Fx,Fy,Fz,x1,x2,y1,y2,z1,z2;
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
    for(double c = c1; c<c2; c+=dc)
    {
        double dx = c/force_n;
        /** The first wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = -0.5*c + dx*i;
            x2 = -0.5*c + dx*(i+1);
            y1 = 0.5*c;
            y2 = 0.5*c;
            z1 = z;
            z2 = z;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The second wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = 0.5*c - dx*i;
            x2 = 0.5*c - dx*(i+1);
            y1 = -0.5*c;
            y2 = -0.5*c;
            z1 = z;
            z2 = z;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The third wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = 0.5*c;
            x2 = 0.5*c;
            y1 = 0.5*c - dx*i;
            y2 = 0.5*c - dx*(i+1);
            z1 = z;
            z2 = z;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The fourth wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = -0.5*c;
            x2 = -0.5*c;
            y1 = -0.5*c + dx*i;
            y2 = -0.5*c + dx*(i+1);
            z1 = z;
            z2 = z;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The fifth wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = -0.5*c + dx*i;
            x2 = -0.5*c + dx*(i+1);
            y1 = 0.5*c;
            y2 = 0.5*c;
            z1 = z + Dl;
            z2 = z + Dl;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The sixth wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = 0.5*c - dx*i;
            x2 = 0.5*c - dx*(i+1);
            y1 = -0.5*c;
            y2 = -0.5*c;
            z1 = z + Dl;
            z2 = z + Dl;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The seventh wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = 0.5*c;
            x2 = 0.5*c;
            y1 = 0.5*c - dx*i;
            y2 = 0.5*c - dx*(i+1);
            z1 = z + Dl;
            z2 = z + Dl;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
        /** The eigth wire */
        for(int i=0; i<force_n; i++)
        {
            x1 = -0.5*c;
            x2 = -0.5*c;
            y1 = -0.5*c + dx*i;
            y2 = -0.5*c + dx*(i+1);
            z1 = z + Dl;
            z2 = z + Dl;
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);
            Fx += dFx;
            Fy += dFy;
            Fz += dFz;
        }
    }
    cout<<Fx<<"\t"<<Fy<<"\t"<<Fz<<endl;
    cout<<Fz*v<<"\t"<<e*e/R<<endl;
    delete Induction;
    return(Fz);
}

double magnet::Force_rot(double omega, double alpha)
{
    double dc = (c2-c1)/coil_n;
    double z = d;
    double e = emf_interp_rot(omega,alpha);
    double* Induction = new double[3];
    double dlx, dly, dlz, dFx, dFy, dFz, Mx, My, Mz, dMx, dMy, dMz, x1, x2, y1, y2, z1, z2;
    double temp_x1, temp_x2;
    double temp_y1, temp_y2;
    double temp_z1, temp_z2;
    Mx = 0;
    My = 0;
    Mz = 0;
    for(double c = c1; c<c2; c+=dc)
    {
        double dx = c/force_n;
        /** The first wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = -0.5*c + dx*i;
            temp_x2 = -0.5*c + dx*(i+1);
            temp_y1 = 0.5*c;
            temp_y2 = 0.5*c;
            temp_z1 = z;
            temp_z2 = z;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The second wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = 0.5*c - dx*i;
            temp_x2 = 0.5*c - dx*(i+1);
            temp_y1 = -0.5*c;
            temp_y2 = -0.5*c;
            temp_z1 = z;
            temp_z2 = z;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The third wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = 0.5*c;
            temp_x2 = 0.5*c;
            temp_y1 = 0.5*c - dx*i;
            temp_y2 = 0.5*c - dx*(i+1);
            temp_z1 = z;
            temp_z2 = z;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The fourth wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = -0.5*c;
            temp_x2 = -0.5*c;
            temp_y1 = -0.5*c + dx*i;
            temp_y2 = -0.5*c + dx*(i+1);
            temp_z1 = z;
            temp_z2 = z;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The fifth wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = -0.5*c + dx*i;
            temp_x2 = -0.5*c + dx*(i+1);
            temp_y1 = 0.5*c;
            temp_y2 = 0.5*c;
            temp_z1 = z + Dl;
            temp_z2 = z + Dl;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The sixth wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = 0.5*c - dx*i;
            temp_x2 = 0.5*c - dx*(i+1);
            temp_y1 = -0.5*c;
            temp_y2 = -0.5*c;
            temp_z1 = z + Dl;
            temp_z2 = z + Dl;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The seventh wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = 0.5*c;
            temp_x2 = 0.5*c;
            temp_y1 = 0.5*c - dx*i;
            temp_y2 = 0.5*c - dx*(i+1);
            temp_z1 = z + Dl;
            temp_z2 = z + Dl;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
        /** The eigth wire */
        for(int i=0; i<force_n; i++)
        {
            temp_x1 = -0.5*c;
            temp_x2 = -0.5*c;
            temp_y1 = -0.5*c + dx*i;
            temp_y2 = -0.5*c + dx*(i+1);
            temp_z1 = z + Dl;
            temp_z2 = z + Dl;
            x1 = temp_x1;
            x2 = temp_x2;
            y1 = temp_y1*cos(alpha) - temp_z1*sin(alpha);
            y2 = temp_y2*cos(alpha) - temp_z2*sin(alpha);
            z1 = temp_y1*sin(alpha) + temp_z1*cos(alpha);
            z2 = temp_y2*sin(alpha) + temp_z2*cos(alpha);
            dlx = x2-x1;
            dly = y2-y1;
            dlz = z2-z1;
            B(0.5*(x1+x2),0.5*(y1+y2),0.5*(z1+z2),Induction);
            dFx = e/R*(dly*Induction[2] - dlz*Induction[1]);
            dFy = e/R*(dlz*Induction[0] - dlx*Induction[2]);
            dFz = e/R*(dlx*Induction[1] - dly*Induction[0]);

            dMx = 0.5*(y1+y2)*dFz - 0.5*(z1+z2)*dFy;
            dMy = 0.5*(z1+z2)*dFx - 0.5*(x1+x2)*dFz;
            dMz = 0.5*(x1+x2)*dFy - 0.5*(y1+y2)*dFx;

            Mx+=dMx;
            My+=dMy;
            Mz+=dMz;
        }
    }
    delete Induction;
    //cout<<Mx<<"\t"<<My<<"\t"<<Mz<<endl;
    cout<<Mx*omega<<"\t"<<e*e/R<<endl;


    return(Mx);
}



double magnet::GetResist()
{
    return(R);
}

magnet::magnet()
{
    INIReader reader("magnet_src/param.ini");
    /** Initialization of the magnet's parameters */
    a = reader.GetReal("Magnet","a",0.0025);
    b = reader.GetReal("Magnet","b",0.002);
    Curr = reader.GetReal("Magnet","i",1114249.0);
    magnet_n = reader.GetInteger("Numerical","magnet_n",100);
    /** Initialization of the coil's parameters */
    c1 = reader.GetReal("Coil","c1",0.00009);
    c2 = reader.GetReal("Coil","c2",0.0028 );
    d = reader.GetReal("Coil","d",0.0008);
    Dl = reader.GetReal("Coil","Dl",0.0002);
    coil_n = reader.GetInteger("Coil","n",72);
    R = reader.GetReal("Coil","R",192.0);
    /** Initialization of the constants */
    mu0 = reader.GetReal("Constants","mu0",1.25663706e-06);
    /** Initialization of numerical constants */
    flux_n = reader.GetInteger("Numerical","flux_n",100);
    force_n = reader.GetInteger("Numerical","force_n",100);
}
