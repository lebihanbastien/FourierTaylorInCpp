#ifndef ORBIT_H_INCLUDED
#define ORBIT_H_INCLUDED

#include "Invman.h"
#include "ftc_errno.h"



class Orbit
{
    private:

        //Parameterization
        int order;            //order of the pm (can be <= invman.ofts_order)
        int ofs_order;        //order of the Fourier coeffs (can be <= invman.ofs_order)
        int reduced_nv;       //reduced number of variables (== invman.reduced_nv)
        int fwrk;             //framework (== invman.fwrk)
        Ofsc ofs;             //Auxiliary Ofs object

        //Time and projection
        double    tfx;         //final time after computation
        double    t0x;         //initial time
        double    tprojx;      //default time between each projection
        double    tprojminx;   //minimum time between each projection
        double    ePmaxx;      //maximum projection distance allowed

        //Initial, current and final state
        double   *z0x;         //Initial position in NC coordinates dim = 6
        double   *six;         //Initial RCM configuration dim = REDUCED_NV
        double   *s0dx;        //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
        cdouble  *s0x;         //Initial position in CCM8 coordinates (real+imag part) dim = 4
        double   *xfx;         //Final position dim = 6


        //Macro structures and objects
        Invman const *invman; //The associated manifold (cannot be null)
        QBCP_L const *qbcp_l; //QBCP around a given Li point (parent)
        OdeStruct *driver;    //NC ode struct


    public:
         Orbit(Invman const *invman, QBCP_L const *qbcp_l, OdeStruct *driver, int ofts_order,
                int ofs_order, double t0, double tf);
        ~Orbit();

        //--------------------------------------------------------------------------------
        //Getters
        //--------------------------------------------------------------------------------
        const double  getT0() const;
        const double  getN()  const;
        const double  getTf() const;
        const double* getZ0() const;
        const double* getSi() const;
        const double  getSi(int dim) const;
        const Invman* getInvman() const;

        //--------------------------------------------------------------------------------
        //Setters
        //--------------------------------------------------------------------------------
        void setTf(double tf);
        void setT0(double t0);
        void setSi(double value, int i);
        void addSi(double value, int i);

        //--------------------------------------------------------------------------------
        //Update
        //--------------------------------------------------------------------------------
        void update_ic(const double st0[], double t0);
        void update_ic(const double st0[]);
        void ccm8torcm(const double s0d[]);
        void evalRCMtoNC(double const t, double z1[]) const;

        //--------------------------------------------------------------------------------
        //Integrates
        //--------------------------------------------------------------------------------
        int traj_int_grid(double tf, double **yNCE, double *tNCE, int N, int isResetOn);
        int traj_int_grid(double **yNCE, double *tgridNCE, int N, int isResetOn);
        int traj_int_var_grid(double tf, double **yNCE, double *tNCE, int N, int isResetOn);
        int gslc_proj_step(double yv[], double *t, double t0, double t1, double *projdist, int *nreset, int isResetOn);
        int gslc_proj_evolve(double yv[], double *t, double t0, double t1, double *projdist, int *nreset, int isResetOn);

        //--------------------------------------------------------------------------------
        //Projection on (un)stable manifold
        //--------------------------------------------------------------------------------
        void NCprojCCMtoCM(double *yv, double tv, double sti[5]);

};

//========================================================================================
//          Orbit on a grid
//========================================================================================
int oo_gridOrbit(double st0[], double t0, double tf, double dt);

#endif // ORBIT_H_INCLUDED
