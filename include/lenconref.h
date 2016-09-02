#ifndef LENCONREF_H_INCLUDED
#define LENCONREF_H_INCLUDED

#include "lencon_io.h"

#define SI_NORM_EM_MAX  56.0
#define SI_NORM_SEM_MAX 1.0

//=======================================================================================================================================
//
//          DIFFCORR BASED ON GOMEZ ET AL. 1998
//
//=======================================================================================================================================
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
int multiple_shooting_direct(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                             int N, int man_grid_size, int coord_type,
                             int isPlotted, gnuplot_ctrl *h1);
/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        Note that, contrary to the Gomez et al. article, the times are free to vary, leading to additionnal free variables and constraints.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
        int N, int man_grid_size, int coord_type,
        int isPlotted, gnuplot_ctrl *h1);
//=======================================================================================================================================
//
//          DIFFCORR CUSTOM: CMU to CM
//
//=======================================================================================================================================

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

//=======================================================================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS
//
//=======================================================================================================================================
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
int msdvt_CMS_RCM_deps(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
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
int msdvt_CMS_RCM_deps_ATF(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
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
int msdvt_CMS_RCM_deps_planar(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
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
int msdvt_CMS_RCM_deps_planar_ATF(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                                  int number_of_variables, int man_grid_size, int coord_type, double precision, int isFirst,
                                  matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                                  gsl_matrix_complex *CCM_R_RCM_EM, gsl_matrix_complex *CCM_R_RCM_SEM,
                                  SingleOrbit &orbit_EM, SingleOrbit &orbit_SEM,
                                  gnuplot_ctrl *h1, int isPlotted, int isUserDefined, int isDebug);

//=======================================================================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS with PSEUDO-ARCLENGTH CONSTRAINT
//
//=======================================================================================================================================
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

//=======================================================================================================================================
//
//          DIFFCORR BASED ON LEVEL II Differential Corrector (Howell & Barden)
//
//=======================================================================================================================================
/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_level_I(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                    int N, int man_grid_size,
                                    int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);
#endif // LENCONREF_H_INCLUDED
