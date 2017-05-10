#ifndef OOLENCONREF_H_INCLUDED
#define OOLENCONREF_H_INCLUDED

#include "lencon_io.h"
#include "ephemerides.h"
#include "Orbit.h"
#include "ftc_errno.h"



//----------------------------------------------------------------------------------------
// Pointer structures
//----------------------------------------------------------------------------------------
/**
 *  \brief  Pointer to a given differential corrector. Examples: msft3d, msvt3d, etc.
 **/
typedef int (*diffcorrptr)(double**, double*, double**, double*, double*,
                           int, int, int, int, Orbit&, Orbit&,
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
int msft3d(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           int isFirst, Orbit& orbit_EM, Orbit& orbit_SEM,
           gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msftmixed(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
              int nov, int mgs, int coord_type,  int isFirst,
              Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1,
              RefSt& refSt, int* niter);


int msvt3d(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           int isFirst,
           Orbit& orbit_EM, Orbit& orbit_SEM,
           gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msvtmixed(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
              int nov, int mgs, int coord_type,
              int isFirst,
              Orbit& orbit_EM, Orbit& orbit_SEM,
              gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msftplan(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type,  int isFirst,
             Orbit& orbit_EM, Orbit& orbit_SEM,
             gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msvtplan(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type,  int isFirst,
             Orbit& orbit_EM, Orbit& orbit_SEM,
             gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msvltplan(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
              int number_of_variables, int man_grid_size,
              int coord_type,  int isFirst,
              Orbit& orbit_EM, Orbit& orbit_SEM,
              gnuplot_ctrl* h1, RefSt& refSt, int* niter);

int msvftplan(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
              int nov, int mgs, int coord_type,  int isFirst,
              Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1,
              RefSt& refSt, int* niter);


int msftplan_pa(double** ymd, double* tmd, double** ymdn, double* tmdn,
                double* nullvector, int nov, int mgs, int coord_type,  int isFirst,
                Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt,
                int* niter);

int msvftplan_dH(double** ymd, double* tmd, double** ymdn, double* tmdn,
                 double* nullvector,int nov, int mgs, int coord_type,  int isFirst,
                 Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt,
                 int* niter);

int msvftplan_dte(double** ymd, double* tmd, double** ymdn, double* tmdn,
                  double* nullvector, int nov, int mgs, int coord_type,  int isFirst,
                  Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt,
                  int* niter);

;//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: JACOBIAN MATRICES
//
//========================================================================================
int jacvftplan(int k, gsl_matrix* DF, int dim, int mgs, gsl_matrix** Ji, gsl_matrix* Phi0, gsl_matrix* PhiN, gsl_vector* K4, gsl_matrix* Id);
int jacftplan(int k, gsl_matrix* DF, int mgs, gsl_matrix** Ji, gsl_matrix* Phi0, gsl_matrix* PhiN, gsl_matrix* Id);
int jacftplan_dt0(int k, gsl_matrix* DF, int dim, int mgs, gsl_matrix** Ji, gsl_matrix* Phi0, gsl_matrix* PhiN, gsl_vector* K4, gsl_matrix* Id);
int nullvectorjac(double* nullvector, gsl_matrix* DF, int ncs, int nfv);
int t0jac(gsl_vector* Kout, double* tmdn, double** ymdn, Orbit& orbit_EM, Orbit& orbit_SEM,
          OdeParams& odeParams, vfptr vf, gsl_matrix* Phi0, gsl_vector* Ktemp);
//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: UPDATE FREE VARIABLES with NEW IMPLEMENTATION
//
//========================================================================================

int ufvarft3d(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
              double* nullvector,
              Orbit& orbit_EM, Orbit& orbit_SEM,
              int mgs, int coord_type,  RefSt& refSt);

int ufvarvt3d(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
              double* nullvector,
              Orbit& orbit_EM, Orbit& orbit_SEM,
              int mgs, int coord_type,  RefSt& refSt);

int ufvarftmixed(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                 double* nullvector,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int mgs, int coord_type,  RefSt& refSt);

int ufvarvtmixed(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                 double* nullvector,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int mgs, int coord_type,  RefSt& refSt);


int ufvarftplan(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                double* nullvector,
                Orbit& orbit_EM, Orbit& orbit_SEM,
                int mgs, int coord_type,  RefSt& refSt);

int ufvarvtplan(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                double* nullvector,
                Orbit& orbit_EM, Orbit& orbit_SEM,
                int man_grid_size, int coord_type,  RefSt& refSt);

int ufvarvltplan(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                 double* nullvector,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int man_grid_size, int coord_type,  RefSt& refSt);

int ufvarvftplan(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                 double* nullvector,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int mgs, int coord_type,  RefSt& refSt);

int ufvarvftplan_dH(double** y_traj_n, double* t_traj_n, double* ds, double ds0,
                    double* nullvector,
                    Orbit& orbit_EM, Orbit& orbit_SEM,
                    int mgs, int coord_type,  RefSt& refSt);

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
