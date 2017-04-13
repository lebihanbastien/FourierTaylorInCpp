#ifndef SINGLE_ORBIT_H_INCLUDED
#define SINGLE_ORBIT_H_INCLUDED

#include "env.h"
#include "ode.h"
#include "pmcoc.h"
#include "vf.h"
#include "gslc.h"
#include "eminsem.h"
#include "timec.h"
#include "ephemerides.h"

#include <list>
#include <algorithm>    // std::sort
#include <vector>       // std::vector

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#define TYPE_STROB         -1
#define TYPE_ORBIT         0
#define TYPE_CU            1
#define TYPE_CS            2
#define TYPE_MAN           3
//#define TYPE_MAN_SORT_DR   4
//#define TYPE_MAN_SORT_DH   5
#define TYPE_MAN_PROJ      6
#define TYPE_MAN_SORT      7
#define TYPE_MAN_SORT_IN   8
#define TYPE_CONT_ATF      9
#define TYPE_CONT_ATF_TRAJ 10
//3D
#define TYPE_CU_3D         11
//For JPL
#define TYPE_COMP_FOR_JPL  12
#define TYPE_CONT_JPL_TRAJ 13

#define TYPE_MAN_PROJ_3D    61
#define TYPE_MAN_PROJ_ORBIT 62

#define GSIZE 50;
#define ORBIT_SEM_UNSTABLE_MIN 1e-5

extern "C" {
    #include "nrutil.h"
    #include "gnuplot_i.h"
}


extern int COMPLETION;


//========================================================================================
//
//          SingleOrbit structure - deprecated
//
//========================================================================================
/**
 *  \struct SingleOrbit
 *  \brief  Defines a given orbit with proper arrays to store results.
 **/
typedef struct SingleOrbit SingleOrbit;
struct SingleOrbit
{
    //-----------
    //Parent
    //-----------
    //-----------
    QBCP_L   *qbcp_l;              //QBCP around a given Li point (parent)

    //-----------
    //Parameterization (common to all orbits)
    //-----------
    vector<Oftsc>*  W;             //z(t) = W(s(t), t)
    vector<Oftsc>*  Wh;            //zh(t) = Wh(s(t), t)
    Ofsc* ofs;                     //Auxiliary Ofs object
    double  n;                     //Pulsation of the QBCP
    int order;                     //order of the pm
    int ofs_order;                 //order of the Fourier coefficients
    int reduced_nv;                //reduced number of variables
    bool isGS;                     //was the pm obtained through graph style?

    //-----------
    //COC (common to all orbits)
    //-----------
    matrix<Ofsc>* PC;               //COC matrix
    matrix<Ofsc>* CQ;               //COC matrix
    vector<Ofsc>* V;                //COC vector

    //Characteristics
    //-----------
    double   *z0;                    //Initial position in NC coordinates dim = 6
    double   *si;                    //Initial RCM configuration dim = REDUCED_NV
    double   *s0d;                   //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
    cdouble  *s0;                    //Initial position in CCM8 coordinates (real+imag part) dim = 4
    double   *xf;                    //Final position dim = 6
    double    tf;                    //final time after computation
    double    t0;                    //initial time
    double    tproj;                 //default time between each projection
    double    tprojmin;              //minimum time between each projection
    double    ePmax;                 //maximum projection distance allowed

    //-----------
    //ODE integration
    //-----------
    OdeStruct *driver;              //NC ode struct
};



//========================================================================================
//
//          Integration - main routines
//
//========================================================================================
/**
 * \brief Integrates a generic vector field from any input type to any output type.
 **/
int ode78(double **yv, double *tv, OdeEvent *odeEvent,
          double t0NC, double tfNC, const double *y0NC,
          int nvar, int nGrid, int dcs,
          int inputType, int outputType);

/**
 * \brief Integrates a generic vector field from any input type to any output type.
 **/
int ode78(double **yv, double *tv, int *ode78coll,
          double t0NC, double tfNC, const double *y0NC,
          int nvar, int nGrid, int dcs,
          int inputType, int outputType);


//========================================================================================
//
//          Integration - with event
//
//========================================================================================
/**
 * \brief Integrates the QBCP vector field from any input type to any output type.
 *        The output are states and time corresponding to a particular event stored
 *        in the structure val_par.
 **/
int ode78_qbcp_event(double **ye, double *te, int *ode78coll,
                     double t0NC, double tfNC, const double *y0NC,
                     int nvar, int dcs,
                     int inputType, int outputType,
                     value_params* val_par);

//========================================================================================
//
//          Integration - variable and given time grids (not really used)
//
//========================================================================================
/**
 * \brief Integrates the QBCP vector field from any input type to any output type. Return the last position that is filled on the grid.
 **/
int ode78_qbcp_vg(double **yv, double *tv, int *ode78coll,
                  double t0NC, double tfNC, const double *y0NC,
                  int nvar, int nGridmax, int dcs,
                  int inputType, int outputType, int dli);

/**
 * \brief Integrates the QBCP vector field from any input type to any output type.
 *        The time grid is given in inputs (tvi).
 **/
int ode78_qbcp_gg(double **yv, double *tvf, int *ode78coll,
                  const double *tvi,
                  const double *y0NC,
                  int nvar, int nGrid, int dcs,
                  int inputType, int outputType);


//========================================================================================
//
//          Integration - JPL (deprecated)
//
//========================================================================================
/**
 * \brief Integrates the Ephemerides vector field. Only VECLI, VSEM and VEM coordinates are used for now.
 **/
int ode78_jpl(double **yv, double *tv, int *ode78coll,
              double t0NC, double tfNC, const double *y0NC,
              int nvar, int nGrid, int dcs,
              int inputType, int outputType);


//========================================================================================
//
//          Integration - subroutines
//
//========================================================================================
/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int ode78_grid(OdeStruct *odestruct, double t0, double tf, double *y0, double **yNCE, double *tNCE, int N);

/**
 *  \brief Stepper that includes detections of event as described in OdeParams.
 **/
int gslc_event_step(OdeStruct* odestruct, OdeParams* odP, double yv[], double* t, double t1);

/**
 *   \brief Integrates the current state yv up to t = t1, includes detections of event as described in OdeParams.
 **/
int gslc_event_evolve(OdeStruct* odestruct, OdeParams* odP, double yv[], double* t, double t1);

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid. Return the last position that is filled on the grid.
 **/
int ode78_variable_grid(OdeStruct *odestruct, double t0, double tf, double *y0, double **yNCE, double *tNCE, int N);

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid. The time grid is given as input
 **/
int ode78_grid_gg(OdeStruct *odestruct, double *y0, double **yNCE, const double *tNCEi, int N);

//========================================================================================
//
//          Integration - SingleOrbit (deprecated)
//
//========================================================================================
/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N. Return the last position that is filled on the grid.
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);

/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(SingleOrbit &orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *proj_dist_SEM,
                   int *nreset,
                   int isResetOn);


/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(SingleOrbit &orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *proj_dist_SEM,
                     int *nreset,
                     int isResetOn);

//========================================================================================
//
//          Projection on (un)stable manifold - deprecated
//
//========================================================================================
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CM_TFC, MIcoc, Vcoc). Then the CCM state is extended by adding a non-null direction
 *        along the hyperbolic direction (sti[4]).
 **/
void NCprojCCMtoCUS(double *yv, double tv, double n, double sti[5], vector<Oftsc> &CM_TFC, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc);


/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CM_TFC, MIcoc, Vcoc).
 **/
void NCprojCCMtoCM(double *yv, double tv, double n, double sti[5], vector<Oftsc> &CM_TFC, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc);


//========================================================================================
//
//                  Energy on vectors
//
//========================================================================================
/**
 *  \brief Hamiltonian along one orbit in SEM units and SEM coordinates
 **/
void HSEM_vec(double *tSEM, double **ySEM, double *Hvec, int N, QBCP_L *qbcp_l);

/**
 *  \brief Hamiltonian along the dyneq of SEMLi for a given time vector in SEM units and SEM coordinates
 **/
void HSEMLi_vec(double *tSEM, double *Hvec, int N, QBCP_L *qbcp_l);


/**
 *  \brief Hamiltonian along the dyneq of EMLi for a given time vector in SEM units and SEM coordinates
 **/
void HEMLi_in_SEM_vec(double *tEM, double *Hvec, int N, QBCP_L *qbcp_l);


//========================================================================================
//
//         Update points
//
//========================================================================================
/**
 *  \brief Update some key positions of notable points in the EM system
 **/
void emPoints(double t, double **emP);

/**
 *  \brief Update some key positions of notable points in the SEM system
 **/
void semPoints(double t, double **semP);

/**
 *  \brief Update some key positions of notable points in the EM system (Normalized version)
 **/
void emNCPoints(double t, double **emP);

/**
 *  \brief Update some key positions of notable points in the SEM system (Normalized version)
 **/
void semNCPoints(double t, double **semP);

//========================================================================================
//
//        Plots
//
//========================================================================================
/**
 *  \brief Sets notable points in EM system on the gnuplot window ctrl h1
 **/
void emPlot(gnuplot_ctrl *h1, double **emP, int isPlot);

/**
 *  \brief Sets notable points in SEM system on the gnuplot window ctrl h2
 **/
void semPlot(gnuplot_ctrl *h2, double **semP, int isPlot);


//========================================================================================
//
//          I/O orbit
//
//========================================================================================
/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, size, and type of the data.
 **/
string filenameOrbit(int ofts_order, int sizeOrbit, int type);

/**
 * \brief Store the orbit (tNCE, yNCE) in the  data file filenameOrbit(ofts_order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, int ofts_order, int sizeOrbit, int N, int type);

/**
 * \brief Store the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(ofts_order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, double **sNCE, int ofts_order, int sizeOrbit, int N, int type);

/**
 * \brief Get the length of the data file filenameOrbit(ofts_order, sizeOrbit, type).
 **/
int getLineNumber(int ofts_order, int sizeOrbit, int type);

/**
 * \brief Read the orbit (tNCE, yNCE) in the  data file filenameOrbit(ofts_order, size, type)
 **/
int readOrbit(double *tNCE, double **yNCE, int ofts_order, int sizeOrbit, int N, int type);

/**
 * \brief Read the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(ofts_order, sizeOrbit, type)
 **/
int readOrbit(double *tNCE, double **yNCE, double **sNCE, int ofts_order, int sizeOrbit, int N, int type);

//========================================================================================
//
//          QBCP test
//
//========================================================================================
/**
 *  \brief Derivatives of the QBCP in EM inertial coordinates (primaries + fourth body motion)
 *          - First 8 variable is the primaries' motion (z, Z in real form).
 *          - Last  6 variables is the state (Xin, Vin).
 */
int qbcp_derivatives_em_in(double t, const double y[], double f[], void *params);

/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbcp_test(QBCP_L &qbcp_l);

//========================================================================================
//
//          Orbit on a grid - deprecated
//
//========================================================================================
/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initial conditions st0, on the center manifold CM.
 **/
int gridOrbit(double st0[],
              double t0,
              double tf,
              double dt,
              vector<Oftsc> &CM,
              vector<Oftsc> &CM_TFC,
              matrix<Ofsc>  &Mcoc,
              matrix<Ofsc>  &MIcoc,
              vector<Ofsc>  &Vcoc);


//========================================================================================
//
//          SingleOrbit structure - deprecated
//
//========================================================================================
/**
 *   \brief Initialize one SingleOrbit structure
 **/
void init_orbit(SingleOrbit &orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Ofsc>*  PC,
                matrix<Ofsc>*  CQ,
                vector<Ofsc>*  V,
                Ofsc* orbit_ofs,
                int ofts_order,
                int ofs_order,
                int isGS,
                double t0,
                double tf,
                double tproj,
                OdeStruct *driver,
                QBCP_L *qbcp_l);

/**
 *   \brief Free one orbit
 **/
void free_orbit(SingleOrbit *orbit);

//========================================================================================
//
//          Initial conditions - deprecated
//
//========================================================================================
/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0);

/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 *          and an array of initial NC conditions z0.
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], const double z0[]);


#endif // SINGLE_ORBIT_H_INCLUDED
