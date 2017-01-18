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

#define REF_SINGLE     2

//CAREFUL: because of the routine isCont(), we cannot create a continuation procedure
//that has a number above 39. We're good for now (only 30, 31, and 32 are taken)
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
    int dim;              //planar or 3d
    int grid;             //variable or fixed grid
    int time;             //variable or fixed time
    int type;             //single solution or continuation procedure
    int cont_step_max;    //maximum number of steps in cont procedure, if necessary
    int cont_step_max_vt; //maximum number of steps in cont procedure with variable time
    int isDirUD;          //the direction of refinement is fixed by the user if true
    int Dir;              //the direction of refinement if isDirUD = false
    int isFlagOn;         //are the "press enter to go on" active
    int isLimUD;          //the domain of search for first guess fixed by the user if true
    int isPlotted;        //some additionnal plots are made if true
    int isSaved;          //some additionnal storage is made if true
    int isJPL;            //refinement to JPL ephemerides
    int isDebug;          //debugging flag
    int isFromServer;     //data are from server
    int coord_type;       //desired type of coordinates for certain applications
    int gridSize;         //desired grid size, for certain applications

    int termination;      //type of termination for time varying continuation

    //Sampling frequencies
    int sf_eml2;
    int sf_man;
    int sf_seml2;

    int djplcoord; //Default JPL coordinates

    //Limits for domain of research of the first guess
    double s1_CMU_EM_MIN;
    double s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN;
    double s3_CMU_EM_MAX;

    //Poincaré section x = cst
    double xps;

    //Desired time
    double t0_des;

    //Limits for integration time
    double tspan_EM;
    double tspan_SEM;

    //Check if the type is of continuation type
    // We check that the type is between 30 (REF_CONT) and 39 (maximum allowed)
    bool isCont(){return (type/10*10 == REF_CONT);}
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
int nfreevariables(RefSt refst, int man_grid_size);

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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
 *        - The final state zN vary in the center-stable manifold of SEML2.
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
