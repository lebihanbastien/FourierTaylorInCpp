#ifndef LENCONREF_H_INCLUDED
#define LENCONREF_H_INCLUDED

#include "lencon_io.h"

#define SI_NORM_EM_MAX  56.0
#define SI_NORM_SEM_MAX 1.0

/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_level_I(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                    int N, int man_grid_size,
                                    int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);

/**
 * \brief Multiple shooting scheme with no boundary conditions
 **/
int multiple_shooting_gomez(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                            int N, int man_grid_size,
                            int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                             int N, int man_grid_size, int coord_type,
                             int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);

/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
        int N, int man_grid_size, int coord_type,
        int isPlotted, int isTimeFixed, gnuplot_ctrl *h1);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msd_CM_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                int N, int man_grid_size, int coord_type,
                                int isPlotted, int isTimeFixed, gnuplot_ctrl *h1,
                                matrix<Ofsc>  &Mcoc_EM,
                                matrix<Oftsc> &DCM_EM_TFC,
                                gsl_matrix_complex *CCM_R_RCM,
                                SingleOrbit &orbit);
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CM_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                 int N, int man_grid_size, int coord_type,
                 int isPlotted, int isTimeFixed, gnuplot_ctrl *h1,
                 matrix<Ofsc>  &Mcoc_EM, matrix<Ofsc>  &Mcoc_SEM,
                 vector<Ofsc>  &Vcoc_SEM,
                 matrix<Oftsc> &DCM_EM_TFC, matrix<Oftsc> &DCM_SEM_TFC,
                 vector<Oftsc> &CM_SEM_TFC,
                 gsl_matrix_complex *CCM_R_RCM,
                 SingleOrbit &orbit,
                 SingleOrbit &orbit_SEM,
                 int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM(double **ymd,
                  double *tmd,
                  double **ymdn,
                  double *tmdn,
                  int number_of_variables,
                  int man_grid_size,
                  int coord_type,
                  matrix<Oftsc> &DCM_EM_TFC,
                  matrix<Oftsc> &DCM_SEM_TFC,
                  gsl_matrix_complex *CCM_R_RCM_EM,
                  gsl_matrix_complex *CCM_R_RCM_SEM,
                  SingleOrbit &orbit,
                  SingleOrbit &orbit_SEM,
                  gnuplot_ctrl *h1,
                  int isPlotted,
                  int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_ATF(double **ymd,
                           double *tmd,
                           double **ymdn,
                           double *tmdn,
                           double *nullvector,
                           int number_of_variables,
                           int man_grid_size,
                           int coord_type,
                           double precision,
                           int isFirst,
                           matrix<Oftsc> &DCM_EM_TFC,
                           matrix<Oftsc> &DCM_SEM_TFC,
                           gsl_matrix_complex *CCM_R_RCM_EM,
                           gsl_matrix_complex *CCM_R_RCM_SEM,
                           SingleOrbit &orbit_EM,
                           SingleOrbit &orbit_SEM,
                           gnuplot_ctrl *h1,
                           int isPlotted,
                           int isUserDefined,
                           int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions and fixed last time.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_ATF(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isUserDefined,
    int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions + pseuod-arclength constraint
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_pac(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    double *conv_free_var,
    double ds,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug);

/**
 * \brief Multiple shooting scheme with no boundary conditions and fixed last time + pseuod-arclength constraint.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_pac_ATF(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    double *conv_free_var,
    double ds,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug);

#endif // LENCONREF_H_INCLUDED
