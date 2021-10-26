#include <math.h>
#include <iostream>
#include "field.h"
#include "INIReader.h"
#include <fstream>
#include "progressbar.h"
#include <iomanip>
#include <cmath>
#include <thread>

using namespace std;

double magnetic::dist(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
/** Function returnes the closest distance (in meters) between the "wire" BC and the observation point A
*/
{
    return(sqrt((x0 * x0 * y1 * y1 - 2.0 * x0 * x0 * y1 * y2 + x0 * x0 * y2 * y2 + y1 * y1 * z0 * z0 - 2.0 * y1 * y2 * z0 * z0 + y2 * y2 * z0 * z0 + x2 * x2 * (y0 * y0 - 2.0 * y0 * y1 + y1 * y1 + pow((z0 - z1),2.0)) - 2.0 * y0 * y1 * z0 * z1 + 2.0 * y0 * y2 * z0 * z1 + 2.0 * y1 * y2 * z0 * z1 - 2.0 * y2 * y2 * z0 * z1 + x0 * x0 * z1 * z1 + y0 * y0 * z1 * z1 - 2.0 * y0 * y2 * z1 * z1 + y2 * y2 * z1 * z1 + x1 * x1 * (y0 * y0 - 2.0 * y0 * y2 + y2 * y2 + pow((z0 - z2),2.0)) - 2.0 * x0 * x2 * (y1 * y1 - y1 * y2 + y0 * (-y1 + y2) - (z0 - z1) * (z1 - z2)) + 2.0 * y0 * y1 * z0 * z2 - 2.0 * y1 * y1 * z0 * z2 - 2.0 * y0 * y2 * z0 * z2 + 2.0 * y1 * y2 * z0 * z2 - 2.0 * x0 * x0 * z1 * z2 - 2.0 * y0 * y0 * z1 * z2 + 2.0 * y0 * y1 * z1 * z2 + 2.0 * y0 * y2 * z1 * z2 - 2.0 * y1 * y2 * z1 * z2 + x0 * x0 * z2 * z2 + y0 * y0 * z2 * z2 - 2.0 * y0 * y1 * z2 * z2 + y1 * y1 * z2 * z2 - 2.0 * x1 * (x2 * (y0 * y0 + y1 * y2 - y0 * (y1 + y2) + (z0 - z1) * (z0 - z2)) + x0 * (y0 * (y1 - y2) - y1 * y2 + y2 * y2 + z0 * z1 - z0 * z2 - z1 * z2 + z2 * z2)))/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2)));
}

double magnetic::cosine1(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
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

double magnetic::cosine2(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/)
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

void magnetic::Crosspoint(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/)
/** Function findes the coordinates of the closest point of the line BC to the point A
*/
{
    ans[0] = -((-x0 * x1 * x1 + 2.0 * x0 * x1 * x2 - x0 * x2 * x2 - x1 * y0 * y1 + x2 * y0 * y1 - x2 * y1 * y1 + x1 * y0 * y2 - x2 * y0 * y2 + x1 * y1 * y2 + x2 * y1 * y2 - x1 * y2 * y2 - x1 * z0 * z1 + x2 * z0 * z1 - x2 * z1 * z1 + x1 * z0 * z2 - x2 * z0 * z2 + x1 * z1 * z2 + x2 * z1 * z2 - x1 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
    ans[1] = -((-x0 * x1 * y1 + x0 * x2 * y1 + x1 * x2 * y1 - x2 * x2 * y1 - y0 * y1 * y1 + x0 * x1 * y2 - x1 * x1 * y2 - x0 * x2 * y2 + x1 * x2 * y2 + 2.0 * y0 * y1 * y2 - y0 * y2 * y2 - y1 * z0 * z1 + y2 * z0 * z1 - y2 * z1 * z1 + y1 * z0 * z2 - y2 * z0 * z2 + y1 * z1 * z2 + y2 * z1 * z2 - y1 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
    ans[2] = -((-x0 * x1 * z1 + x0 * x2 * z1 + x1 * x2 * z1 - x2 * x2 * z1 - y0 * y1 * z1 + y0 * y2 * z1 + y1 * y2 * z1 - y2 * y2 * z1 - z0 * z1 * z1 + x0 * x1 * z2 - x1 * x1 * z2 - x0 * x2 * z2 + x1 * x2 * z2 + y0 * y1 * z2 - y1 * y1 * z2 - y0 * y2 * z2 + y1 * y2 * z2 + 2.0 * z0 * z1 * z2 - z0 * z2 * z2)/(x1 * x1 - 2.0 * x1 * x2 + x2 * x2 + y1 * y1 - 2.0 * y1 * y2 + y2 * y2 + z1 * z1 - 2.0 * z1 * z2 + z2 * z2));
}

void magnetic::nB(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/)
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

void magnetic::B(double x_coord/**x-coordinate of the observation point*/, double y_coord/**y-coordinate of the observation point*/, double z_coord/**z-coordinate of the observation point*/, double magnet_current_dens/**current density of the equivalent coil*/, double* ans/**pointer on the array with answer[3] ( ans[0] = Bx, ans[1] = By, ans[2] = Bz )*/)
/** Function calculates the magnet flux density produced by the block-magnet at the point */
{
    double Bind [3] = {0.0,0.0,0.0};
    double* norm = new double [3];
    double val_B;
    double x1,x2,y1,y2,z1,z2;
    /** Magnet #1*/
    for(int i=0; i<magnet_n; i++)
    {
        x1 = 0.5*magnet_a2 + magnet_b3*i/(magnet_n-1);
        x2 = 0.5*magnet_a2 + magnet_b3*i/(magnet_n-1);
        y1 = 0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = 0.5*magnet_a1 + magnet_b1;
        z2 = 0.5*magnet_a1 + magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = 0.5*magnet_a1 + magnet_b1;
        z2 = 0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = 0.5*magnet_a1;
        z2 = 0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = 0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = 0.5*magnet_a1;
        z2 = 0.5*magnet_a1 + magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
    }
    /** Magnet #2*/
    for(int i=0; i<magnet_n; i++)
    {
        x1 = -0.5*magnet_a2 - magnet_b3*i/(magnet_n-1);
        x2 = -0.5*magnet_a2 - magnet_b3*i/(magnet_n-1);
        y1 = 0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = 0.5*magnet_a1 + magnet_b1;
        z2 = 0.5*magnet_a1 + magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = 0.5*magnet_a1 + magnet_b1;
        z2 = 0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = 0.5*magnet_a1;
        z2 = 0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = 0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = 0.5*magnet_a1;
        z2 = 0.5*magnet_a1 + magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
    }
    /** Magnet #3*/
    for(int i=0; i<magnet_n; i++)
    {
        x1 = 0.5*magnet_a2 + magnet_b3*i/(magnet_n-1);
        x2 = 0.5*magnet_a2 + magnet_b3*i/(magnet_n-1);
        y1 = 0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = -0.5*magnet_a1 - magnet_b1;
        z2 = -0.5*magnet_a1 - magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = -0.5*magnet_a1 - magnet_b1;
        z2 = -0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = -0.5*magnet_a1;
        z2 = -0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = 0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = -0.5*magnet_a1;
        z2 = -0.5*magnet_a1 - magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
    }
    /** Magnet #4*/
    for(int i=0; i<magnet_n; i++)
    {
        x1 = -0.5*magnet_a2 - magnet_b3*i/(magnet_n-1);
        x2 = -0.5*magnet_a2 - magnet_b3*i/(magnet_n-1);
        y1 = 0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = -0.5*magnet_a1 - magnet_b1;
        z2 = -0.5*magnet_a1 - magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = -0.5*magnet_b2;
        z1 = -0.5*magnet_a1 - magnet_b1;
        z2 = -0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = -0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = -0.5*magnet_a1;
        z2 = -0.5*magnet_a1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
        y1 = 0.5*magnet_b2;
        y2 = 0.5*magnet_b2;
        z1 = -0.5*magnet_a1;
        z2 = -0.5*magnet_a1 - magnet_b1;
        nB(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2, norm);
        val_B = (mu0*(magnet_current_dens*magnet_b3/magnet_n)/(4.0*M_PI*dist(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)))*(cosine1(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2)-cosine2(x_coord, y_coord, z_coord, x1 , y1 , z1 , x2 , y2 , z2));
        Bind[0]+=val_B*norm[0];
        Bind[1]+=val_B*norm[1];
        Bind[2]+=val_B*norm[2];
    }
    ans[0] = Bind[0];
    ans[1] = Bind[1];
    ans[2] = Bind[2];
    delete norm;
}

double magnetic::Flux(double xc/**x-coordinate of the center of the wire loop*/, double yc/**y-coordinate of the center of the wire loop*/, double zc/**z-coordinate of the center of the wire loop*/, double radius/**radius of the loop*/)
/** Function returnes the flux through the round plane with the center in xc,yc,zc and the known radius*/
{
    double F=0.0;
    double *res = new double [3];
    double dy = 2.0*radius/flux_n;
    double dz = 2.0*radius/flux_n;
    double dS = dy*dz;
    for(double y=yc - radius + 0.5*dy; y<yc + radius; y+=dy)
    {
        for(double z=zc - radius + 0.5*dz; z< zc + radius; z+=dz)
        {
            if(sqrt((y-yc)*(y-yc) + (z-zc)*(z-zc))<=radius)
            {
                B(xc,y,z,magnet_i,res);
                F+=dS*res[0];

            }
        }
    }
    delete res;
    return(F);
}

void magnetic::GenerateFluxFileSingle(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step */, double r /** Radius of the loop */, double x /** x-shift */)
/** Return the file with derivative dPhi/dz vs shift for z-axis for the loop with radius r and shift x */
    /** File has the following structure:       */
    /** Flux_<radius>_<x-shift>.dat             */
    /** It contains the following information:  */
    /** <z_shift>   <dPhi/dz>                   */
{
    std::ofstream file("Flux_"+std::to_string(r)+"_"+std::to_string(x)+".dat");
    cout<<"r = \t"<<r<<"\tx = \t"<<x<<endl;
    for(double z = z_min; z<z_max+dz; z+=dz)
    {
        double F2 = Flux(x, 0.0, z+dz/100.0, r);
        double F1 = Flux(x, 0.0, z-dz/100.0, r);
        file<<z<<"\t"<<100.0*(F2-F1)/(2.0*dz)<<"\n";
        std::cout<<z<<"\t"<<100.0*(F2-F1)/(2.0*dz)<<"\n";
    }
    file.close();
}

void magnetic::GenerateFlux(double Zmin /** Minimal z-shift*/, double Zmax /** Maximal z-shift*/, double dZ /** z-shift step*/)
/** Interface method generates the set of files for each loop of the coil described in the paper        */
/** usage is following:                                                                                 */
/** magnetic mag();
    mag.GenerateFlux(0.0, 0.0021, 0.0001);
*/
/** The method is designed to work with 8-thread processor.                                             */
{
    double dx = coil_d/coil_n1;
    double dr = (coil_r_max-coil_r_min)/coil_n2;
    for(double x = 1.6129e-05; x<coil_d*0.5+dx; x+=dx)
    {
        for( double r = coil_r_min; r< coil_r_max + dr; r+= dr)
        {
            thread th1(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th2(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th3(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th4(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th5(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th6(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th7(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            r+=dr;
            thread th8(&magnetic::GenerateFluxFileSingle,this,Zmin,Zmax,dZ,r,x);
            th1.join();
            th2.join();
            th3.join();
            th4.join();
            th5.join();
            th6.join();
            th7.join();
            th8.join();
        }
    }
}


double magnetic::EMFFromFile(double z /** z-shift of the coil */, double v /** speed of the coil */)
/** Method use the file <fit.dat> to generate the electromotive force                                                       */
/** <fit.dat> contains the dPhi/dz(z) = v*(b1 + b1*z^2 + b2*z^4 + b3*z^6 + b4*z^8 + b5*z^10) parameters for the each loop   */
/** <fit.dat> has the following structure:                                                                                  */
/** <x_shift>   <radius>    <b1>    <b2>    <b3>    <b4>    <b5>                                                            */
{
    std::ifstream infile("magnet_src/fit_deriv.dat");
    string a0,a1,a2,a3,a4,a5,a6,a7;
    int i=0;
    double data[9][5000];
    while (infile>>a0>>a1>>a2>>a3>>a4>>a5>>a6>>a7)
    {
        data[0][i] = stod (a0);
        data[1][i] = stod (a1);
        data[2][i] = stod (a2);
        data[3][i] = stod (a3);
        data[4][i] = stod (a4);
        data[5][i] = stod (a5);
        data[6][i] = stod (a6);
        data[7][i] = stod (a7);
        i++;
    }
    int n_of_p = i;
    double EMF = 0.0;
    for(i = 0; i<n_of_p; i++)
    {
        EMF+=(data[2][i] + data[3][i]*z*z + data[4][i]*std::pow(z,4) + data[5][i]*std::pow(z,6) + data[6][i]*std::pow(z,8) + data[7][i]*std::pow(z,10))*v;
    }
    return(EMF);
}

double magnetic::EMForce(double zc /** z-shift of the coil */, double vc /** speed of the coil */)
/** Method use the file <fit.dat> to generate the electromotive force                                                       */
/** <fit.dat> contains the dPhi/dz(z) = v*(b1 + b1*z^2 + b2*z^4 + b3*z^6 + b4*z^8 + b5*z^10) parameters for the each loop   */
/** <fit.dat> has the following structure:                                                                                  */
/** <x_shift>   <radius>    <b1>    <b2>    <b3>    <b4>    <b5>                                                            */
{
    double curr = EMFFromFile(zc,vc)/coil_resist;
    double dalpha = 2*M_PI/coil_force_n;
    double y,z;
    double dlx,dly,dlz;
    double* Ind = new double[3];
    double ForceX = 0;
    double ForceY = 0;
    double ForceZ = 0;

    std::ifstream infile("magnet_src/fit_deriv.dat");
    string a0,a1,a2,a3,a4,a5,a6,a7;
    int i=0;
    double data[9][5000];
    while (infile>>a0>>a1>>a2>>a3>>a4>>a5>>a6>>a7)
    {
        data[0][i] = stod (a0);
        data[1][i] = stod (a1);
        data[2][i] = stod (a2);
        data[3][i] = stod (a3);
        data[4][i] = stod (a4);
        data[5][i] = stod (a5);
        data[6][i] = stod (a6);
        data[7][i] = stod (a7);
        i++;
    }
    int n_of_p = i;

    for(i=0; i<n_of_p; i++)
    {
    double x = data[1][i];
    double r = data[0][i];

        for (double alpha = 0; alpha< 2*M_PI; alpha+=dalpha)
            {
                z = zc+r*cos(alpha+0.5*dalpha);
                y = r*sin(alpha+0.5*dalpha);
                B(x,y,z,magnet_i,Ind);
                dlx = 0.0;
                dly = r*sin(alpha+dalpha) - r*sin(alpha);
                dlz = r*cos(alpha+dalpha) - r*cos(alpha);
                ForceX += curr*(dly*Ind[2]-dlz*Ind[1]);
                ForceY += curr*(dlz*Ind[0]-dlx*Ind[2]);
                ForceZ += curr*(dlx*Ind[1]-dly*Ind[0]);
            }
    }
    delete Ind;
    return(ForceZ);

}

double magnetic::EMFInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */)
/** Interpolation function for the EMF:                                                                                       */
/** EMF(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                                 */
{
    return(vc*(emf_a + emf_b*zc*zc + emf_c*pow(zc,4) + emf_d*pow(zc,6) + emf_e*pow(zc,8) + emf_f*pow(zc,10)));
}

double magnetic::EMForceInterp(double zc /** z-shift of the coil */, double vc /** speed of the coil */)
/** Interpolation function for the electromagnetic force:                                                                   */
/** F(z,v) = v*(a + b*z^2 + c*z^4 + d*z^6 + e*z^8 + f*z^10)                                                                 */
{
    return(vc*(force_a + force_b*zc*zc + force_c*pow(zc,4) + force_d*pow(zc,6) + force_e*pow(zc,8) + force_f*pow(zc,10)));
}

void magnetic::GenerateForce(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step*/, double v /** velocity*/)
/** Interface generates <Force.dat> file with self-containment controle                                                     */
/** <z-shift>   <force>                                                                                                     */
/** Usage:                                                                                                                  */
/** magnetic mag;
    mag.GenerateForce(-0.002,0.0021,0.0001,0.01);
*/
{
    std::ofstream out("magnet_src/Force.dat");
    for(double zc = z_min; zc<z_max; zc+=dz)
    {
        double F = EMForce(zc,v);
        std::cout<<"\t"<<zc<<"\t"<<F<<"\t"<<F*v<<"\t"<<EMFFromFile(zc,v)*EMFFromFile(zc,v)/coil_resist<<std::endl;
        out<<zc<<"\t"<<F<<std::endl;
    }
    out.close();
}

void magnetic::GenerateEMF(double z_min /** Minimal z-shift */, double z_max /** Maximal z-shift */, double dz /** z-shift step*/, double v /** velocity*/)
/** Interface generates <EMF.dat> file with self-containment controle                                                     */
/** <z-shift>   <EMF>                                                                                                     */
/** Usage:                                                                                                                  */
/** magnetic mag;
    mag.GenerateForce(-0.002,0.0021,0.0001,0.01);
*/
{
    std::ofstream out("magnet_src/EMF.dat");
    for(double zc = z_min; zc<z_max; zc+=dz)
    {
        double EMF = EMFFromFile(zc,v);
        std::cout<<zc<<"\t"<<EMF<<std::endl;
        out<<zc<<"\t"<<EMF<<std::endl;
    }
    out.close();
}

double magnetic::RightPart_dv(double x /** displacement */, double dx /** velocity */, double t/** time */, double freq/** frequency */, bool magn_force/** open circuit or closed*/)
/** Dimention right part of equation dv/dt:     */
/** if magn_force = True - closed circuit       */
/** if magn_force = False - open circuit        */
{
    double rp;
    if (magn_force)
        rp = Aext*cos(2.0*M_PI*freq*t) - DampingB*dx/mass - springK*x/mass - springK2*x*abs(x)/mass - springK3*pow(x,3)/mass - springK4*abs(x)*pow(x,3)/mass - springK5*pow(x,5)/mass - EMForceInterp(x,dx)/mass;
    else
        rp = Aext*cos(2.0*M_PI*freq*t) - DampingB*dx/mass - springK*x/mass - springK2*x*abs(x)/mass - springK3*pow(x,3)/mass - springK4*abs(x)*pow(x,3)/mass - springK5*pow(x,5)/mass;
    return(rp);
}

double magnetic::RightPart_dx(double dx /** velocity */)
/** Dimention right part of equation dx/dt:     */
{
    return(dx);
}

void magnetic::RK4(double fext/** frequency */, bool magn_force /** open circuit or closed */)
/** Runge-Kutta method*/
{
    double x = 0.0;
    double dx = 0.0;
    double t = 0.0;
    double x1,x2;
    double dt = RK45_step;
    double p1,p2,p3,p4,l1,l2,l3,l4;

    double integr_ch = 0.0;
    double integr_ne = -1.0;
    double integr = 0.0;

    int counter = 0;
    int counter_period = 0;
    int counter_1 = 0;

    double v_rms = 0.0;
    int n_rms = 0;

    double ampl_max = 0.0;
    double ampl_min = 0.0;

    //ofstream ofs("waveform"+to_string(fext)+".dat");

    while((fabs(1-fabs(integr_ch)/fabs(integr_ne))>RK45_exit)&(counter<200)&(t<8.0))
    {
        x1 = x;

        p1 = dt*RightPart_dv(x,dx,t,fext,magn_force);
        l1 = dt*RightPart_dx(dx);

        p2 = dt*RightPart_dv(x + 0.5*l1, dx + 0.5*p1, t + 0.5*dt, fext, magn_force);
        l2 = dt*RightPart_dx(dx + 0.5*p1);

        p3 = dt*RightPart_dv(x + 0.5*l2, dx + 0.5*p2, t + 0.5*dt, fext, magn_force);
        l3 = dt*RightPart_dx(dx + 0.5*p2);

        p4 = dt*RightPart_dv(x + l3, dx + p3, t + dt, fext, magn_force);
        l4 = dt*RightPart_dx(dx + p3);

        t+=dt;
        x+=(l1 + 2.0*l2 + 2.0*l3 + l4)/6.0;
        dx+=(p1 + 2.0*p2 + 2.0*p3 + p4)/6.0;

        x2 = x;

        counter_1++;

        //if(counter_1%1000 == 0)
        //    ofs<<t<<"\t"<<x<<"\t"<<EMFInterp(x,dx)<<endl;

        integr+=fabs(dx)*fabs(x2-x1);

        v_rms += EMFInterp(x,dx)*EMFInterp(x,dx);
        n_rms ++;

        if(x>ampl_max)
            ampl_max = x;

        if (x<ampl_min)
            ampl_min = x;

        if(x1*x2 <=0)
            counter_period ++;

        if((x1*x2 <= 0)&(counter_period%5 == 0))
        {
            //cout<<counter<<"\t"<<counter_period<<"\t"<<integr_ch<<"\t"<<integr_ne<<"\t"<<fabs(1-(fabs(integr_ch)/fabs(integr_ne)))<<"\n";
            if(counter%2 == 0)
                integr_ch = integr;
            else
                integr_ne = integr;
            counter++;
            integr = 0.0;

            export_v_rms = sqrt(v_rms/n_rms);
            export_p_rms = export_v_rms*export_v_rms/coil_resist;
            export_ampl = 0.5*(ampl_max-ampl_min);

            v_rms = 0.0;
            n_rms = 0;
        }

    }
    //ofs.close();
}

double magnetic::GetRMSV()
{
    return(export_v_rms);
}

double magnetic::GetRMSP()
{
    return(export_p_rms);
}

double magnetic::GetAmpl()
{
    return(export_ampl);
}

magnetic::magnetic()
{
    INIReader reader("magnet_src/param.ini");
    if (reader.ParseError() < 0)
    {
        std::cout << "Can't load 'param.ini'\n";
    }
    magnet_a1   = reader.GetReal("Magnet", "a1", -1);
    magnet_a2   = reader.GetReal("Magnet", "a2", -1);
    magnet_b1   = reader.GetReal("Magnet", "b1", -1);
    magnet_b2   = reader.GetReal("Magnet", "b2", -1);
    magnet_b3   = reader.GetReal("Magnet", "b3", -1);
    magnet_i    = reader.GetReal("Magnet", "i", -1);
    magnet_n    = reader.GetInteger("Numerical", "magnet_n", -1);

    mu0         = reader.GetReal("Constants","mu0",-1);
    flux_n      = reader.GetInteger("Numerical", "flux_n", -1);

    coil_n1     = reader.GetInteger("Coil", "n1", -1);
    coil_n2     = reader.GetInteger("Coil", "n2", -1);
    coil_d      = reader.GetReal("Coil", "d", -1);
    coil_r_min  = reader.GetReal("Coil", "r_min", -1);
    coil_r_max  = reader.GetReal("Coil", "r_max", -1);
    coil_resist = reader.GetReal("Coil", "resist", -1);
    coil_force_n     = reader.GetInteger("Numerical","coil_force_n",-1);

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
    Aext = reader.GetReal("Dynamics","Aext",-1);

    RK45_step = reader.GetReal("Numerical", "RK45_step", -1);
    RK45_exit = reader.GetReal("Numerical", "RK45_exit", -1);

    //std::cout<<flux_n<<std::endl;
}

void GenerateResonance(double f_min, double f_max, double df)
{
    magnetic mag;
    ofstream ofs("resonance_curve.dat");
    for(double freq = f_min; freq<f_max; freq+=df)
    {
        mag.RK4(freq,false);
        cout<<freq<<"\t"<<mag.GetAmpl()<<"\t"<<mag.GetRMSV()<<"\t"<<mag.GetRMSP()<<endl;
        ofs<<freq<<"\t"<<mag.GetAmpl()<<"\t"<<mag.GetRMSV()<<"\t"<<mag.GetRMSP()<<endl;
    }
    ofs.close();
}
