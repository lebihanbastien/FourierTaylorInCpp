#ifndef OOLENCON_H_INCLUDED
#define OOLENCON_H_INCLUDED


#include <vector>

#include "ode.h"
#include "env.h"
#include "Oftsc.h"
#include "Orbit.h"
#include "matrix.h"
#include "lencon_io.h"
#include "lenconref.h"

extern "C" {
    #include "nrutil.h"
}

//========================================================================================
//
//          Computation of the CMU about EML2
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a 5-dimensional
 *         box: one dimension for the starting time, four dimensions for the
 *         parameterization of the Center Manifold (s1 to s4 coordinates). The RCM
 *         coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates on the unstable direction s5.
 *  \param tlim_CMU_EM:    the min/max starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param si_LIM_CMU_RCM: the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                         in the IC box.
 *  \param si_grid_size:   the number of points on the  s1, s2, s3, s4 values  grids
 *                         in the IC box.
 *  \param invman:         the center-unstable manifold, if the type of manifold provided
 *                         is not center-unstable, a warning message is displayed and
 *                         nothing is done.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/EM/L2/cu_3d_order_16.bin"
 **/
int oo_compute_grid_CMU_EM_3D(double dist_to_cm, double *tlim_CMU_EM,
                              int t_grid_size, double si_LIM_CMU_RCM[4][2],
                              int *si_grid_size, Invman &invman, bool isPar);

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         three-dimensional box: one dimension for the starting time, two dimensions for
 *         the parameterization of the Center Manifold (s1 and s3 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates on the unstable direction s5.
 *  \param tmin_CMU_EM:    the minimum starting time (in EM units) in the IC box.
 *  \param tmax_CMU_EM:    the maximum starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param s1_MIN_CMU_RCM: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param s1_MAX_CMU_RCM: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param s3_MIN_CMU_RCM: the minimum s3 value (in RCM coordinates) in the IC box.
 *  \param s3_MAX_CMU_RCM: the maximum s3 value (in RCM coordinates) in the IC box.
 *  \param s1_grid_size:   the number of points on the s1 grid in the IC box.
 *  \param s3_grid_size:   the number of points on the s3 grid in the IC box.
 *  \param invman:         the center-unstable manifold, if the type of manifold provided
 *                         is not center-unstable, a warning message is displayed and
 *                         nothing is done.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int oo_compute_grid_CMU_EM(double dist_to_cm, double tmin_CMU_EM, double tmax_CMU_EM,
                           int t_grid_size, double s1_MIN_CMU_RCM, double s1_MAX_CMU_RCM,
                           double s3_MIN_CMU_RCM, double s3_MAX_CMU_RCM,
                           int s1_grid_size, int s3_grid_size, Invman &invman, bool isPar);

//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEML2
//
//========================================================================================
int oo_int_proj_CMU_EM_on_CM_SEM_3D(double tmax_on_manifold_EM,
                                    int man_grid_size, int nod,
                                    int isPar, double ynormMax,
                                    double snormMax);

int oo_int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM, int t_grid_size_x,
                                 int s1_grid_size_x, int s3_grid_size_x,
                                 int man_grid_size, int NsortMin, int nod, int isPar,
                                 double ynormMax, double snormMax);

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - general routines
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2.
 *         The final conditions vary in the paramerization of the CMS of SEML2.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 **/
int oorefeml2seml(int man_grid_size, int coord_type, Invman &invman, RefSt &refst);

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - subroutines
//
//========================================================================================
/**
 *  \brief Continuation of a single of EML2-to-SEML2 connection,
 *         between orbit_EM and orbit_SEM.
 *         The Jacobian of the parameterization of the manifolds are
 *         contained in  DCM_EM_TFC and DCMS_SEM_TFC.
**/
void oosrefeml2seml(Orbit &orbit_EM, Orbit &orbit_SEM,
                  int dcs, int coord_type, int man_grid_size_t,
                  RefSt &refst, gnuplot_ctrl *h2);

//========================================================================================
//
//         Refinement of solutions: Complete trajectory
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the WHOLE trajectory
 *         (EML2 orbit + manifold leg + SEML2 orbit).
 **/
int oocomprefft3d(int man_grid_size_t, int coord_type,
                  Orbit &orbit_EM, Orbit &orbit_SEM,
                  RefSt refst);

/**
 *  \brief Refine a given output of oocomprefft3d into JPL ephemerides.
 **/
int oojplrefft3d(int coord_type);

//----------------------------------------------------------------------------------------
// Text format, read
//----------------------------------------------------------------------------------------
void toCelestiaFormat(string filename);

#endif // OOLENCON_H_INCLUDED
