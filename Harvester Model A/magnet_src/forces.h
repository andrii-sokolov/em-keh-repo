#ifndef FORCES_H_INCLUDED
#define FORCES_H_INCLUDED

double dist(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
/** Function returnes the closest distance (in meters) between the "wire" BC and the observation point A*/

double cosine1(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
/** Function returnes the cos(alpha1) the angle between vectors \vec{AB} and \vec{CB}*/

double cosine2(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/);
/** Function returnes the cos(alpha1) the angle between vectors \vec{AC} and \vec{CB}*/

void Crosspoint(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/);
/** Function findes the coordinates of the closest point of the line BC to the point A*/

void nB(double x0/*! x-coordinate of the observation point A*/, double y0 /*! y-coordinate of the observation point A*/, double z0 /*! z-coordinate of the observation point A*/, double x1 /*!x-coordinate of the point B*/, double y1 /*!y-coordinate of the point B*/, double z1 /*!z-coordinate of the point B*/, double x2 /*!x-coordinate of the point C*/, double y2 /*!y-coordinate of the point C*/, double z2/*!z-coordinate of the point C*/, double* ans /*! pointer on array with answer [3]*/);
/** Function findes the direction of the magnet flux density in point A created by the current in line-segment BC*/

class magnet
{
    private:
        double a;           //Width and length of the magnet
        double b;           //Height of the magnet
        double d;           //The distance between the coil and the magnet
        double c1;          //Internal width of the coil
        double c2;          //External width of the coil
        double Curr;           //Current in equivalent coil
        double Dl;          //The distance between layers
        double mu0;         //Magnetic constant
        double R;           //Coil resistance
        int coil_n;         //The number of loops per layer
        int magnet_n;       //The number of layers for the equivalent coil
        int flux_n;         //Number of mesh cells in flux calculation
        int force_n;
    public:
        magnet();           //Constructor
        void B(double x_coord/**x-coordinate of the observation point*/, double y_coord/**y-coordinate of the observation point*/, double z_coord/**z-coordinate of the observation point*/, double* ans/**pointer on the array with answer[3] ( ans[0] = Bx, ans[1] = By, ans[2] = Bz )*/);
        double Flux_single(double c, double z, double alpha_x);
        void GenerateDerivFlux_rot(double alpha_0, double alpha_1,double d_alpha);
        void GenerateDerivFlux_transp(double z_1, double z_2,double dz);
        void GenerateEMF_transp(double v);
        double emf_interp_transp(double v, double z);
        double emf_interp_rot(double omega, double alpha);
        double Force_transp(double v, double z);
        double Force_rot(double omega, double alpha);
        double Force_rot_interp(double omega, double alpha);
        double Force_progres_interp(double omega, double alpha);
        void GenerateEMF_rot(double omega);
        double GetResist();

};

#endif // FORCES_H_INCLUDED
