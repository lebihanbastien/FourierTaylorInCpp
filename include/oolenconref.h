#ifndef OOLENCONREF_H_INCLUDED
#define OOLENCONREF_H_INCLUDED

#include "lencon_io.h"
#include "lenconref.h"
#include "ephemerides.h"
#include "Orbit.h"
#include "ftc_errno.h"


//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS with new implementation
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEML2.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msft3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           double precision, int isFirst,
           Orbit &orbit_EM, Orbit &orbit_SEM,
           gnuplot_ctrl *h1, RefSt &refst, int *niter);

int msvt3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
           int number_of_variables, int man_grid_size, int coord_type,
           double precision, int isFirst,
           Orbit &orbit_EM, Orbit &orbit_SEM,
           gnuplot_ctrl *h1, RefSt &refst, int *niter);

int msftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type, double precision, int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refst, int *niter);

int msvtplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int number_of_variables, int man_grid_size,
             int coord_type, double precision, int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refst, int *niter);

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: UPDATE FREE VARIABLES with NEW IMPLEMENTATION
//
//========================================================================================

int ufvarft3d(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
              Orbit &orbit_EM, Orbit &orbit_SEM,
              int man_grid_size, int coord_type);

int ufvarvt3d(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
              Orbit &orbit_EM, Orbit &orbit_SEM,
              int man_grid_size, int coord_type);

int ufvarftplan(double **y_traj_n, double *t_traj_n, double ds, double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int man_grid_size, int coord_type);

int ufvarvtplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int man_grid_size, int coord_type);


#endif // OOLENCONREF_H_INCLUDED
