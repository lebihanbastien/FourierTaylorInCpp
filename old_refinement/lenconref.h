#ifndef LENCONREF_H_INCLUDED
#define LENCONREF_H_INCLUDED

#include "lencon_io.h"
#include "ephemerides.h"
#include "Orbit.h"
#include "ftc_errno.h"

#define SI_NORM_EM_MAX  56.0
#define SI_NORM_SEM_MAX 1.0

//Parameters for the type of refinements
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

//Precision
#define PREC_GSM    5e-12

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

//========================================================================================
//
//          SUBROUTINES:
//
//========================================================================================
/**
 *  \brief Solve a definite-positive linear system:
 *         - Decompose a ncs x ncs matrix M that is definite-positive, using GSL routines.
 *         - Then inverse the system Fv = M*K3.
 **/
int ftc_inv_dfls(gsl_matrix* M, gsl_vector*Fv, gsl_vector* K3, int ncs);


/**
 *  \brief Computes the correction vector associated to the minimum norm solution.
 *         Given:
 *              - an ncs x 1   error vector Fv
 *              - an nfv x ncs Jacobian DF,
 *         This routine computes the correction vector associated to
 *         the minimum norm solution:
 *
 *              DQv = DF^T x (DF x DF^T) Fv.
 **/
int ftc_corrvec_mn(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int nfv, int ncs);


//========================================================================================
//
//          DIFFCORR BASED ON GOMEZ ET AL. 1998
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Based on a recursive scheme in Gomez et al., "Quasihalo orbits associated with libration points", 1998, JAS.
 *        The Multiple shooting is denoted GMS for Gomez Multiple Shooting in the commentaries.
 **/
int multiple_shooting_gomez(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                            int N, int man_grid_size,
                            int isPlotted, gnuplot_ctrl *h1);
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct(double **ymd, double *tmd, double **ymdn, double *tmdn,
                             int N, int man_grid_size, int coord_type,
                             int isPlotted, gnuplot_ctrl *h1);
/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        Note that, contrary to the Gomez et al. article, the times are free to vary, leading to additionnal free variables and constraints.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd, double **ymdn, double *tmdn,
                                           int N, int man_grid_size, int coord_type, double prec,
                                           int isPlotted, gnuplot_ctrl *h1);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        An additionnal free variable epsilon is added to the state, for continuation.
 **/
int multiple_shooting_direct_deps(double **ymd, double *tmd,
                                  double **ymdn, double *tmdn,
                                  double *nullvector, int isFirst,
                                  int N, int mgs, int coord_type,
                                  int isPlotted, gnuplot_ctrl *h1);

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CM
//
//========================================================================================
/**
 *  \brief Yields the number of free variables necessary to compute the refinment procedure.
 **/
int nfreevariables(RefSt refSt, int man_grid_size);

//===========================================================
// FIXED TIMES
//===========================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN is free to vary.
 *        - The times t0,..., tN are fixed.
 **/
int msd_CM_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
               int N, int man_grid_size, int coord_type,
               int isPlotted, gnuplot_ctrl *h1,
               matrix<Ofsc>  &Mcoc_EM,
               matrix<Oftsc> &DCM_EM_TFC,
               gsl_matrix_complex *CCM_R_RCM,
               SingleOrbit &orbit);

//===========================================================
// FREE TIMES
//===========================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN is free to vary.
 *        - The times t0,..., tN are free to vary.
 **/
int msdvt_CM_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                 int N,int man_grid_size, int coord_type,
                 int isPlotted, gnuplot_ctrl *h1,
                 matrix<Ofsc>  &Mcoc_EM, matrix<Ofsc>  &Mcoc_SEM, vector<Ofsc>  &Vcoc_SEM,
                 matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                 vector<Oftsc> &CM_SEM_TFC, gsl_matrix_complex *CCM_R_RCM,
                 SingleOrbit &orbit, SingleOrbit &orbit_SEM,
                 int isDebug);

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: UPDATE FREE VARIABLES
//
//========================================================================================
int ufvarft3d(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
              SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
              int man_grid_size, int coord_type);

int ufvarvt3d(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
              SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
              int man_grid_size, int coord_type);

int ufvarftplan(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
                SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                int man_grid_size, int coord_type);

int ufvarvtplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                int man_grid_size, int coord_type);


//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 **/
int msdvt_CMS_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn,
                  int number_of_variables, int man_grid_size, int coord_type, double precision,
                  matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                  gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                  SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                  gnuplot_ctrl *h1, int isPlotted, int isDebug);
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 **/
int msvt3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                       int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                       matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                       gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                       SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                       gnuplot_ctrl *h1, int isPlotted, int isDebug);
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msft3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                           int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                           matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                           gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                           SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                           gnuplot_ctrl *h1, int isPlotted, int isUserDefined, int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 **/
int msvtplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                              int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                              matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                              gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                              SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                              gnuplot_ctrl *h1, int isPlotted, int isDebug);
/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                                  int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                                  matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                                  gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                                  SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                                  gnuplot_ctrl *h1, int isPlotted, int isUserDefined, int isDebug);

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS with PSEUDO-ARCLENGTH CONSTRAINT
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 *        - The pseudo-arclength constraint is added to the constraints. Therefore, the system is SQUARED.
 *
 *        Note: does not yield satisfactory results for the moment (02/09/2016).
 **/
int msdvt_CMS_RCM_deps_planar_pac_ATF(double **ymd, double *tmd, double **ymdn, double *tmdn,
                                      double *nullvector, double *conv_free_var, double ds,
                                      int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                                      matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                                      gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                                      SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                                      gnuplot_ctrl *h1, int isPlotted, int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 *        - The pseudo-arclength constraint is added to the constraints. Therefore, the system is SQUARED.
 *
 *        Note: does not yield satisfactory results for the moment (02/09/2016).
 **/
int msdvt_CMS_RCM_deps_planar_pac(double **ymd, double *tmd, double **ymdn, double *tmdn,
                                  double *nullvector, double *conv_free_var, double ds,
                                  int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                                  matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                                  gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                                  SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                                  gnuplot_ctrl *h1, int isPlotted, int isDebug);

//========================================================================================
//
//          DIFFCORR BASED ON LEVEL II Differential Corrector (Howell & Barden)
//
//========================================================================================
/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_level_I(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                    int N, int man_grid_size,
                                    int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);
#endif // LENCONREF_H_INCLUDED
