#ifndef OOLENCONREF_H_INCLUDED
#define OOLENCONREF_H_INCLUDED

#include "lencon_io.h"
#include "ephemerides.h"
#include "Orbit.h"
#include "ftc_errno.h"


//----------------------------------------------------------------------------------------
// Limits for the domain of practical convergence
//----------------------------------------------------------------------------------------
#define SI_NORM_EM_MAX  56.0
#define SI_NORM_SEM_MAX 1.0

//----------------------------------------------------------------------------------------
//Parameters for the type of refinements (used in structure RefSt)
//----------------------------------------------------------------------------------------
#define REF_PLANAR     0
#define REF_3D         1
#define REF_MIXED      101

#define REF_SINGLE     2

#define REF_CONT             30
#define REF_CONT_D           31
#define REF_CONT_D_HARD_CASE 32

#define REF_COMP       4

#define REF_FIXED_TIME 5
#define REF_VAR_TN     6
#define REF_VAR_TIME   7

#define REF_FIXED_GRID 8
#define REF_VAR_GRID   9
#define REF_GIVEN_GRID 10

#define REF_COND_S5    11
#define REF_COND_T     12

//----------------------------------------------------------------------------------------
// RefSt structure
//----------------------------------------------------------------------------------------
/**
 *  \struct RefSt
 *  \brief  Define a given refinement structures, with a set of parameters
 **/
typedef struct RefSt RefSt;
struct RefSt
{
    //------------------------------------------------------------------------------------
    // Parameters that change often
    //------------------------------------------------------------------------------------
    int type;             //single solution or continuation procedure
    int dim;              //planar or 3d
    double t0_des;        //desired  initial time
    double t0xT_des;      //desired  initial time as a percent

    // Limits for domain of research of the first guess
    double s1_CMU_EM_MIN;
    double s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN;
    double s3_CMU_EM_MAX;

    double s2_CMU_EM_MIN;
    double s2_CMU_EM_MAX;
    double s4_CMU_EM_MIN;
    double s4_CMU_EM_MAX;

    // The domain of search for first guess fixed by the user if true
    int isLimUD;

    // Direction of the continuation procedure
    int isDirUD;          //the direction of refinement is fixed by the user if true
    int Dir;              //the direction of refinement if isDirUD = false

    // Limits for the time of flight during transfers - not used if negative
    double tof_MIN;
    double tof_MAX;

    // Maximum number of steps in the continuation procedure
    int cont_step_max;    //with fixed time
    int cont_step_max_vt; //with variable time

    // Initial step in the continuation procedure
    double ds0;           //with fixed time
    double ds0_vt;        //with variable time

    // Desired number of iterations in Newton's method in the continuation procedure
    int nu0;              //with fixed time
    int nu0_vt;           //with variable time

    // User parameters
    int isFlagOn;         //do we have steps in the procedure - asking the user to press enter to go on?
    int isPlotted;        //do we plot the results during the computation?
    int isSaved;          //do we save the results in data files?
    int isFromServer;     //does the raw data comes from server files?
    int isPar;            //is parallel computation allowed?

    // Maximum angle around SEMLi if REF_COND_T is used (in degrees)
    double thetaMax;      //should be a multiple of 90°

    //------------------------------------------------------------------------------------
    // Parameters that are stable
    //------------------------------------------------------------------------------------
    int isDebug;          //if yes, additionnal tests are performed
    int gridSize;         //number of points on the refinement grid
    int mplot;            //number of points per plot between to pach points (e.g. total plot points is gridSize*mplot)

    int time;             //type of constraints on the times in REF_CONT
    int grid;             //type of grid
    int termination;      //termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
    int coord_type;       //coordinates system in the refinement procedure (usually NCSEM)

    // Maximum/Minimum step in the continuation procedure
    double dsmin;           //with fixed time
    double dsmin_vt;        //with variable time
    double dsmax;           //with fixed time
    double dsmax_vt;        //with variable time

    double xps;           //position of the poincaré section in NCSEM coordinates
    int isJPL;            //is the JPL refinement performed when possible?
    int djplcoord;        //coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with $NJ2000
    int sidim;            //0 or 2 - component of s0 that stays constant when t0 is free.

    //Sampling frequencies in REF_COMP (complete trajectory) in days
    int sf_eml2;          // orbit at EML2
    int sf_man;           // transfer leg
    int sf_seml2;         // orbit at SEML2

    // Integration window for each orbit
    double tspan_EM;
    double tspan_SEM;

    // Storing the orbits at each step?
    int isSaved_EM;    //0: don't save, 1: save using projection method
    int isSaved_SEM;   //0: don't save, 1: save using projection method, 2: save using integration in reduced coordinates


    //Check if the type is of continuation type
    bool isCont(){return (type == REF_CONT || type == REF_CONT_D || type == REF_CONT_D_HARD_CASE );}
    // Check if the trajectory are 3D
    bool is3D(){return (dim == REF_3D || dim == REF_MIXED);}
};


//----------------------------------------------------------------------------------------
// Pointer structures
//----------------------------------------------------------------------------------------
/**
 *  \brief  Pointer to a given differential corrector. Examples: msft3d, msvt3d, etc.
 **/
typedef int (*diffcorrptr)(double**, double*, double**, double*, double*,
                           int, int, int, double, int, Orbit&, Orbit&,
                           gnuplot_ctrl*, RefSt&, int*);

/**
 *  \brief  Pointer to a given predictor. Examples: ufvarft3d, ufvarvt3d, etc.
 **/
typedef int (*predictorptr)(double**, double*, double*, double, double*, Orbit&, Orbit&, int, int, RefSt&);



//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS with new implementation
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msft3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           double precision, int isFirst,
           Orbit &orbit_EM, Orbit &orbit_SEM,
           gnuplot_ctrl *h1, RefSt &refSt, int *niter);

int msftmixed(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
               int nov, int mgs, int coord_type, double precision, int isFirst,
               Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt, int *niter);


int msvt3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           double precision, int isFirst,
           Orbit &orbit_EM, Orbit &orbit_SEM,
           gnuplot_ctrl *h1, RefSt &refSt, int *niter);

int msvtmixed(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
               int nov, int mgs, int coord_type,
               double precision, int isFirst,
               Orbit &orbit_EM, Orbit &orbit_SEM,
               gnuplot_ctrl *h1, RefSt &refSt, int *niter);


int msftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type, double precision, int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refSt, int *niter);

int msvtplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type, double precision, int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refSt, int *niter);

int msvltplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type, double precision, int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refSt, int *niter);

int msvftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
              int nov, int mgs, int coord_type, double precision, int isFirst,
              Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter);


;//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: JACOBIAN MATRICES
//
//========================================================================================
int jacvftplan(int k, gsl_matrix* DF, int dim, int mgs, gsl_matrix **Ji, gsl_matrix *Phi0, gsl_matrix *PhiN, gsl_vector *K4, gsl_matrix *Id);
int jacftplan(int k, gsl_matrix* DF, int mgs, gsl_matrix** Ji, gsl_matrix* Phi0, gsl_matrix* PhiN, gsl_matrix *Id);

int nullvectorjac(double *nullvector, gsl_matrix* DF, int ncs, int nfv);

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: UPDATE FREE VARIABLES with NEW IMPLEMENTATION
//
//========================================================================================

int ufvarft3d(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
              double *nullvector,
              Orbit &orbit_EM, Orbit &orbit_SEM,
              int mgs, int coord_type,  RefSt &refSt);

int ufvarvt3d(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt);

int ufvarftmixed(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt);

int ufvarvtmixed(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt);


int ufvarftplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt);

int ufvarvtplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int man_grid_size, int coord_type,  RefSt &refSt);

int ufvarvltplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int man_grid_size, int coord_type,  RefSt &refSt);

int ufvarvftplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt);

//========================================================================================
//
//          SUBROUTINES:
//
//========================================================================================
/**
 *  \brief Computation of the matrix Phi0. This matrix gives at first order the linear
 *         relationship between the parameterization at EML2 (4 dimensions) at t = t0 and
 *         the final state at the first patch point at t = t1, in coord_type coordinates
 *         (6 dimensions).
 *         It is computed thanks to the following equation:
 *
 *              Phi0 =  J1 * COORD_J_RCM
 *
 *         Where:   - J1 is the 6 x 6 state transition matrix between t0 and t1.
 *                  - COORD_J_RCM is the jacobian matrix dz/ds, where z is the state in
 *                    coord_type coordinates and s is the state in RCM coordinates,
 *                    inside EML2 center-unstable manifold. COORD_J_RCM is a 6 x 5 matrix.
 *
 *         Inputs:  - Phi0, the 6 x 5 gsl_matrix to update.
 *                  - J1, the the 6 x 6 state transition matrix between t0 and t1.
 *                  - orbit_EM, the Orbit object that contains the center-unstable
 *                    manifold at EML2, and therefore can be used to evaluate COORD_J_RCM.
 *                  - t0_EM, the initial time t0, in EM units.
 *                  - coord_type, the final desired coordinate system.
 *
 *         Note: no testing is made on the inputs (dimensions, manifold type...). These
 *         tests should be made by the user higher in the code.
 **/
int ftc_compute_phi0(gsl_matrix* Phi0, gsl_matrix* J1, Orbit& orbit_EM, double t0_EM, int coord_type);

/**
 *  \brief Selection of a differential corrector (msft3d, msvt3d, etc).
 **/
diffcorrptr ftc_select_diffcorr(RefSt& refSt);

/**
 *  \brief Selection of a differential predictor (ufvarft3d, ufvarvt3d, etc).
 *         Must be coherent with ftc_select_diffcorr, just above.
 **/
predictorptr ftc_select_predictor(RefSt& refSt);

/**
 *  \brief Yields the number of free variables necessary to compute the refinment procedure.
 **/
int nfreevariables(RefSt refSt, int man_grid_size);


#endif // OOLENCONREF_H_INCLUDED
