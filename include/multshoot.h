#include "single_orbit.h"

//Precision
#define PREC_GSM    5e-12

using namespace std;

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
 *              DQv = DF^T x (DF x DF^T)^{-1} Fv.
 **/
int ftc_corrvec_mn(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int nfv, int ncs);

/**
 *  \brief Computes the correction vector for a square system.
 *         Given:
 *              - an ncs x 1   error vector Fv
 *              - an ncs x ncs Jacobian DF,
 *         This routine computes the correction vector given by:
 *
 *              DQv = DF^{-1} x Fv.
 **/
int ftc_corrvec_square(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int ncs);


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
                             int N, int man_grid_size, int coord_type, double prec,
                             int isPlotted, gnuplot_ctrl *h1, int strConv);
/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        Note that, contrary to the Gomez et al. article, the times are free to vary, leading to additionnal free variables and constraints.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd, double **ymdn, double *tmdn,
                                           int N, int man_grid_size, int coord_type, double prec,
                                           int isPlotted, gnuplot_ctrl *h1, int strConv);

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

/**
 * \brief Multiple shooting scheme with boundary conditions at the origin (makes the system square)
 *        Contrary to multiple_shooting_gomez,
 *        no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct_square(double **ymd, double *tmd,
                                    double **ymdn, double *tmdn,
                                    int N, int mgs, int coord_type, double prec,
                                    int isPlotted, gnuplot_ctrl *h1, int strConv);
