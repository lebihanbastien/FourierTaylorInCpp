#ifndef LENCON_H_INCLUDED
#define LENCON_H_INCLUDED

#include "ode.h"
#include "lencon_io.h"
#include "lenconref.h"



//========================================================================================
//
//          Main routines (1): connections in 3D
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about EML2, in the QBCP model.
 *         The initial conditions (IC) are computed in a 5-dimensional box: one dimension for the starting time,
 *         four dimensions for the parameterization of the Center Manifold (s1 to s4 coordinates). The RCM coordinate s5 along the unstable
 *         direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates applied on the unstable direction s5.
 *  \param tlim_CMU_EM:    the min/max starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param si_LIM_CMU_RCM: the min/max of s1, s2, s3, s4 values (in RCM coordinates) in the IC box.
 *  \param si_grid_size:   the number of points on the  s1, s2, s3, s4 values  grids in the IC box.
 *  \param CM_TFC:         the Fourier-Taylor representation of the Center-Unstable Manifold about EML2, in TFC coordinates
 *  \param Mcoc:           the Matrix that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param MIcoc:          the invarse of Mcoc.
 *  \param Vcoc:           the vector that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int compute_grid_CMU_EM_3D(double dist_to_cm,
                           double *tlim_CMU_EM,
                           int t_grid_size,
                           double si_LIM_CMU_RCM[4][2],
                           int *si_grid_size,
                           vector<Oftsc> &CM_TFC,
                           matrix<Ofsc>  &Mcoc,
                           vector<Ofsc>  &Vcoc,
                           bool isPar);

//========================================================================================
//
//          Main routines (1): connections in the plane
//
//========================================================================================
/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file (obtained from gridOrbit).
 **/
int compute_grid_CMU_EM(double epsilon,
                        double tmin_CMU_EM,
                        double tmax_CMU_EM,
                        int t_grid_size,
                        double s1_MIN_CMU_RCM,
                        double s1_MAX_CMU_RCM,
                        double s3_MIN_CMU_RCM,
                        double s3_MAX_CMU_RCM,
                        int s1_grid_size,
                        int s3_grid_size,
                        vector<Oftsc> &CM_TFC,
                        matrix<Ofsc>  &Mcoc,
                        vector<Ofsc>  &Vcoc,
                        bool isPar);

/**
 *  \brief Computes the manifold branches from a discrete set of unstable directions obtained from compute_grid_CMU_EM
 **/
int int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM,
                              int t_grid_size_x,
                              int s1_grid_size_x,
                              int s3_grid_size_x,
                              int man_grid_size,
                              int NsortMin,
                              int nod,
                              int isPar,
                              double ynormMax,
                              double snormMax);

int int_proj_CMU_EM_on_CM_SEM_3D(double tmax_on_manifold_EM,
                                 int man_grid_size, int nod,
                                 int isPar, double ynormMax, double snormMax);

/**
 *  \brief Computes the manifold branches from a discrete set of unstable directions obtained from compute_grid_CMU_EM
 **/
int int_proj_CMU_EM_on_CMS_SEM(double tmax_on_manifold_EM,
                               int ofts_order,
                               int t_grid_size_x,
                               int s1_grid_size_x,
                               int s3_grid_size_x,
                               int man_grid_size,
                               int NsortMin,
                               int isPar,
                               double ynormMax,
                               double snormMax);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int int_sorted_sol_CMU_EM_to_CM_SEM(int ofts_order,
                                    int man_grid_size,
                                    int isPar,
                                    vector<Oftsc> &CM_EM_NC,
                                    vector<Oftsc> &CM_EM_TFC,
                                    matrix<Ofsc>  &Mcoc_EM,
                                    matrix<Ofsc>  &MIcoc_EM,
                                    vector<Ofsc>  &Vcoc_EM);
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD(int ofts_order,
                             int man_grid_size,
                             int isPar);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP(int ofts_order,
                                  int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP(int man_grid_size,
                                  int coord_type,
                                  SingleOrbit &orbit_EM,
                                  SingleOrbit &orbit_SEM);
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID(int coord_type,
                                                SingleOrbit &orbit_EM,
                                                SingleOrbit &orbit_SEM);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID_TEST(int coord_type,
                                                     SingleOrbit &orbit_SEM);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int comprefft3d(int man_grid_size,
                int coord_type,
                SingleOrbit &orbit_EM,
                SingleOrbit &orbit_SEM,
                RefSt refst);
/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_seml_dimjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_SEM,
                     RefSt refst);
/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_seml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_SEM,
                     RefSt refst);

/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_eml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_EM,
                     RefSt refst);

/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_eml2seml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_EM,
                     SingleOrbit &orbit_SEM,
                     RefSt refst);

/**
 *  \brief Computes the whole trajectory in INSEM coordinates
 **/
int comprefft3d_test_eml2seml_insem(int man_grid_size_t,
                                    int coord_type,
                                    SingleOrbit& orbit_EM,
                                    SingleOrbit& orbit_SEM,
                                    RefSt refst);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_DEP(int ofts_order,
                                 int man_grid_size,
                                 int coord_type,
                                 vector<Oftsc> &CM_EM_NC,
                                 vector<Oftsc> &CM_EM_TFC,
                                 matrix<Oftsc> &DCM_EM_TFC,
                                 matrix<Ofsc>  &Mcoc_EM,
                                 matrix<Ofsc>  &MIcoc_EM,
                                 vector<Ofsc>  &Vcoc_EM,
                                 int isDebug);


/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_PART(int ofts_order,
                                  int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM);


void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(SingleOrbit &orbit_EM,
                                    SingleOrbit &orbit_SEM,
                                    matrix<Oftsc> &DCM_EM_TFC,
                                    matrix<Oftsc> &DCMS_SEM_TFC,
                                    int dcs,
                                    int coord_type,
                                    int man_grid_size,
                                    gnuplot_ctrl *h2);


/**
 *  \brief Continuation of a single of EML2-to-SEML2 connection, between orbit_EM and orbit_SEM.
 *         The Jacobian of the parameterization of the manifolds are contained in  DCM_EM_TFC and DCMS_SEM_TFC.
**/
void srefeml2seml(SingleOrbit &orbit_EM,
                  SingleOrbit &orbit_SEM,
                  matrix<Oftsc> &DCM_EM_TFC,
                  matrix<Oftsc> &DCMS_SEM_TFC,
                  int dcs,
                  int coord_type,
                  int man_grid_size,
                  RefSt refst,
                  gnuplot_ctrl *h2);

void crefvtplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2);
/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D, in the planar case.
 **/
void crefftplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        int cont_steps_MAX,
        int isSaved,
        int isUserDefined,
        gnuplot_ctrl *h2);

/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D, in the planar case.
 **/
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT_PAC_AFT(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        int cont_steps_MAX,
        gnuplot_ctrl *h2);

void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT_PAC(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2);


void srefvtplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2);

/**
 *  \brief Computes ONE trajectory from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
void crefft3d(SingleOrbit &orbit_EM,
                                             SingleOrbit &orbit_SEM,
                                             matrix<Oftsc> &DCM_EM_TFC,
                                             matrix<Oftsc> &DCMS_SEM_TFC,
                                             int dcs,
                                             int coord_type,
                                             int man_grid_size,
                                             int cont_steps_MAX,
                                             int isUserDefined,
                                             gnuplot_ctrl *h2);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE(int cont_grid_size,
                                                int coord_type,
                                                vector<Oftsc> &CM_EM_NC,
                                                vector<Oftsc> &CM_EM_TFC,
                                                matrix<Oftsc> &DCM_EM_TFC,
                                                matrix<Ofsc>  &Mcoc_EM,
                                                matrix<Ofsc>  &MIcoc_EM,
                                                vector<Ofsc>  &Vcoc_EM,
                                                int isPlanar);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM,
                                  int isPlanar,
                                  int isSaved,
                                  int isUserDefined);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int refeml2seml(int man_grid_size,
                int coord_type,
                vector<Oftsc> &CM_EM_NC,
                vector<Oftsc> &CM_EM_TFC,
                matrix<Oftsc> &DCM_EM_TFC,
                matrix<Ofsc>  &Mcoc_EM,
                matrix<Ofsc>  &MIcoc_EM,
                vector<Ofsc>  &Vcoc_EM,
                RefSt refst);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_T(int proj_grid_size,
        int man_grid_size,
        int coord_type,
        double st_EM[],
        double t0_EM,
        double tmax_on_manifold_EM,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        double ynormMax,
        double snormMax,
        int isPlanar);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_OPT(int proj_grid_size,
        int man_grid_size,
        int coord_type,
        double st_EM[],
        double t0_EM,
        double tmax_on_manifold_EM,
        double s3_MIN_CMU_RCM,
        double s3_MAX_CMU_RCM,
        int s3_grid_size,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        double ynormMax,
        double snormMax,
        int isPlanar,
        int isPar);
//========================================================================================
//========================================================================================

/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file (obtained from gridOrbit).
 *         Interpolation is used to build a function s3 = f(s1). WARNING: the data is hardcoded for now: work in progress!
 **/
int compute_grid_CMU_EM_Interpolated(int ofts_order,
                                     double epsilon,
                                     double tmin_CMU_EM,
                                     double tmax_CMU_EM,
                                     int t_grid_size,
                                     double s1_MIN_CMU_RCM,
                                     double s1_MAX_CMU_RCM,
                                     double s3_MIN_CMU_RCM,
                                     double s3_MAX_CMU_RCM,
                                     int s1_grid_size,
                                     int s3_grid_size,
                                     vector<Oftsc> &CM_TFC,
                                     matrix<Ofsc>  &Mcoc,
                                     matrix<Ofsc>  &MIcoc,
                                     vector<Ofsc>  &Vcoc,
                                     bool isPar);


#endif // LENCON_H_INCLUDED
