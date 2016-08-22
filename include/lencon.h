#ifndef LENCON_H_INCLUDED
#define LENCON_H_INCLUDED

#include "lencon_io.h"
#include "lenconref.h"



//=======================================================================================================================================
//
//          Main routines (1): connections in the plane
//
//=======================================================================================================================================
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
                        matrix<Ofsc>  &MIcoc,
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
                                    matrix<Ofsc>  &Pcoc_EM,
                                    matrix<Ofsc>  &MIcoc_EM,
                                    matrix<Ofsc>  &PIcoc_EM,
                                    vector<Ofsc>  &Vcoc_EM);
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD(int ofts_order,
                             int man_grid_size,
                             int isPar,
                             vector<Oftsc> &CM_EM_NC,
                             vector<Oftsc> &CM_EM_TFC,
                             matrix<Ofsc>  &Mcoc_EM,
                             matrix<Ofsc>  &Pcoc_EM,
                             matrix<Ofsc>  &MIcoc_EM,
                             matrix<Ofsc>  &PIcoc_EM,
                             vector<Ofsc>  &Vcoc_EM);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP(int ofts_order,
                                  int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &Pcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  matrix<Ofsc>  &PIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM);
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
                                 matrix<Ofsc>  &Pcoc_EM,
                                 matrix<Ofsc>  &MIcoc_EM,
                                 matrix<Ofsc>  &PIcoc_EM,
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
                                  matrix<Ofsc>  &Pcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  matrix<Ofsc>  &PIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM);


void ref1_CMU_EM_to_CMS_SEM_MSD_RCM(SingleOrbit &orbit_EM,
                                    SingleOrbit &orbit_SEM,
                                    matrix<Oftsc> &DCM_EM_TFC,
                                    matrix<Oftsc> &DCMS_SEM_TFC,
                                    int dcs,
                                    int coord_type,
                                    int man_grid_size,
                                    gnuplot_ctrl *h2);


void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2);
/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM, in the planar case.
 **/
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT_AFT(SingleOrbit &orbit_EM,
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
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM, in the planar case.
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


void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR(SingleOrbit &orbit_EM,
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
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_CONT_AFT(SingleOrbit &orbit_EM,
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
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &Pcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  matrix<Ofsc>  &PIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM,
                                  int isPlanar,
                                  int isSaved,
                                  int isUserDefined);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM(int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &Pcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  matrix<Ofsc>  &PIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM,
                                  int isPlanar);

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ(int proj_grid_size,
        int man_grid_size,
        int coord_type,
        double st_EM[],
        double t0_EM,
        double tmax_on_manifold_EM,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &Pcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        matrix<Ofsc>  &PIcoc_EM,
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
        double tmin_on_manifold_EM,
        double tmax_on_manifold_EM,
        double s3_MIN_CMU_RCM,
        double s3_MAX_CMU_RCM,
        int s3_grid_size,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &Pcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        matrix<Ofsc>  &PIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        double ynormMax,
        double snormMax,
        int isPlanar,
        int isPar);
//=======================================================================================================================================
//=======================================================================================================================================

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
