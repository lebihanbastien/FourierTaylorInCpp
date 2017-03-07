#ifndef OOLENCON_H_INCLUDED
#define OOLENCON_H_INCLUDED

#include "oolenconref.h"

#include "oolencon.h"


/**
 *  \struct ProjSt
 *  \brief  Define a given projection structure, with a set of parameters
 **/
typedef struct ProjSt ProjSt;
struct ProjSt
{
    double TM, TMIN, TMAX, TLIM[2], GLIM_SI[4][2];
    int    TSIZE, GSIZE_SI[4], MSIZE, NSMIN, NOD, ISPAR;
    double YNMAX, SNMAX;
};

#include "oolencon.h"


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
 *  \param dist_to_cm:      the value in RCM coordinates on the unstable direction s5.
 *  \param projSt.TLIM:     the min/max starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:    the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI:  the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                          in the IC box.
 *  \param projSt.GSIZE_SI: the number of points on the  s1, s2, s3, s4 values  grids
 *                          in the IC box.
 *  \param projSt.ISPAR;    if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/EM/L2/cu_3d_order_16.bin"
 **/
int compute_grid_CMU_EM_3D(double dist_to_cm, ProjSt &projSt);

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         three-dimensional box: one dimension for the starting time, two dimensions for
 *         the parameterization of the Center Manifold (s1 and s3 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:           the value in RCM coordinates on the unst. direction s5.
 *  \param projSt.TLIM[0]:       the minimum starting time (in EM units) in the IC box.
 *  \param projSt.TLIM[1]:       the maximum starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:         the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI[0][0]: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[0][1]: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[2][0]: the minimum s3 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[2][1]: the maximum s3 value (in RCM coordinates) in the IC box.
 *  \param projSt.GSIZE_SI[0]:   the number of points on the s1 grid in the IC box.
 *  \param projSt.GSIZE_SI[2]:   the number of points on the s3 grid in the IC box.
 *  \param projSt.ISPAR:         if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int compute_grid_CMU_EM(double dist_to_cm, ProjSt &projSt);


//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEMLi
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM. Then, each point on the
 *         integration grid is projected on the Center Manifold CM_SEM_NC about SEMLi.
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM: the maximum integration time on each leg, in EM units.
 *
 *  \param projSt.MSIZE:       the number of points on each manifold leg.
 *
 *  \param projSt.NOD:         the number of dimensions on which the distance of
 *                             projection is computed (usually either 3 (the physical
 *                             distance) or 6 (the whole phase space)).
 *
 *  \param projSt.ISPAR:       if TRUE, the computation is parallelized.
 *
 *  \param projSt.YNMAX:       the maximum norm in NCSEM coordinates for which a given
 *                             state on the integration grid is projected on CM_SEM_NC
 *                             More precisely: for a given state y along the manifold leg,
 *                             if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                             Otherwise, it is considered too far away from SEMLi to be
 *                             a good candidate for projection.
 *
 *  \param projSt.SNMAX:       the maximum norm in RCM SEM coordinates for which a given
 *                             projection state on the CM of SEMLi (CM_SEM_NC) is
 *                             computed back in NCSEM coordinates. More precisely, for a
 *                             given state y in NCSEM coordinates, the result of the
 *                             projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                             coordinates.
 *                             if norm(sproj, 4) < projSt.SNMAX, the computation
 *                             yproj = CM_SEM_NC(sproj, t) is performed. Otherwise, the
 *                             state sproj is considered too far away from the RCM origin
 *                             to be a good candidate - it is out of the domain of
 *                             practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form:
 *          "plot/QBCP/EM/L2/projcu_3d_order_16.bin".
 **/
int int_proj_CMU_EM_on_CM_SEM_3D(ProjSt &projSt);

/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM.
 *         Then, each point on the integration grid is projected on the center Manifold
 *         about SEMLi (denoted here CM_SEM_NC).
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM:......the maximum integration time on each leg, in EM units.
 *  \param projSt.MSIZE:...the number of points on each manifold leg.
 *  \param projSt.NSMIN:...the number of best solutions that are kept
 *  \param projSt.NOD:.....the number of dimensions on which the distance of
 *                         projection is computed - usually either 3
 *                         (the physical distance) or 6 (the whole phase space).
 *  \param projSt.ISPAR:...if TRUE, the computation is parallelized.
 *  \param projSt.YNMAX:...the maximum norm in NCSEM coordinates for which a given
 *                         state on the integration grid is projected on CM_SEM_NC.
 *                         More precisely: for a given state y along the manifold leg,
 *                         if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                         Otherwise, it is considered too far away from SEMLi to
 *                         be a good candidate for projection.
 *  \param projSt.SNMAX:...the maximum norm in RCM SEM coordinates for which a given
 *                         projection state on the CM of SEMLi (CM_SEM_NC) is
 *                         computed back in NCSEM coordinates. More precisely, for a
 *                         given state y in NCSEM coordinates, the result of the
 *                         projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                         coordinates. If norm(sproj, 4) < projSt.SNMAX, the computation
 *                         yproj = CM_SEM_NC(sproj, t) is performed. Otherwise,
 *                         the state sproj is considered too far away from the RCM
 *                         origin to be a good candidate - it is out of the domain of
 *                         practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form
 * "plot/QBCP/EM/L2/projcu_order_16.bin", and
 * "plot/QBCP/EM/L2/sortprojcu_order_16.bin" for the projSt.NSMIN best solutions.
 **/
int int_proj_CMU_EM_on_CM_SEM(ProjSt &projSt);

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - general routine
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2.
 *         The final conditions vary in the paramerization of the CMS of SEMLi.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,
 *         the refSt.coord_type must be NCSEM. However, the user can put
 *         other coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines
 *         are not used.
 **/
int refemlisemli(RefSt& refSt);

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - subroutines
//
//========================================================================================
/**
 *  \brief Continuation of a single of EML2-to-SEMLi connection, between orbit_EM and
 *         orbit_SEM.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,the
 *         coord_type must be NCSEM. However, the user can put other
 *         coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines are
 *         not used.
 *
 *         In general, we recommend to use NCSEM at least for precision
 *         reasons. In this way, a warning is issued when coord_type is different
 *         from NCSEM.
 **/
int subrefemlisemli(Orbit& orbit_EM, Orbit& orbit_SEM, double** y_traj, double* t_traj,
                   int dcs, int coord_type, int* man_grid_size_t,
                   RefSt& refSt, gnuplot_ctrl* h2);

//----------------------------------------------------------------------------------------
//         Brick A: select the good IC for EML2-to-SEMli connections in data files
//----------------------------------------------------------------------------------------
/**
 *  \brief Selects good initial conditions for EML2-to-SEMLi connections, searching
 *         through data files produced by int_proj_CMU_EM_on_CM_SEM_3D or
 *         int_proj_CMU_EM_on_CM_SEM.
 **/
int selectemlisemli(RefSt& refSt, double st_EM[5], double st_SEM[5], double t_EM[2],
                  double* t0_SEM, double* pmin_dist_SEM_out);

//----------------------------------------------------------------------------------------
//         Brick B: generate a first guess (either a single unstable manifold leg or
//              a complete trajectory EML2 orbit + man leg + SEMLi orbit).
//----------------------------------------------------------------------------------------
/**
 *  \brief Computing the first guess for the connection leg between the orbit orbit_EM and
 *         and the orbit orbit_SEM. Only the manifold leg is computed.
 *         The grid size is returned.
 **/
int icmanemlisemli(double** y_traj, double* t_traj,
                    Orbit& orbit_EM, Orbit& orbit_SEM,
                    int dcs, int coord_type, int man_grid_size,
                    RefSt& refSt);

/**
 *  \brief Computing the first guess for the connection leg between the orbit orbit_EM and
 *         and the orbit orbit_SEM. The complete trajectory (EML2 + man leg + SEMLi)
 *         is computed.
 *         The grid size is returned.
 **/
int iccompemlisemli(double** y_traj, double* t_traj,
                 double** y_traj_comp, double* t_traj_comp,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int dcs, int coord_type, int grid_points_des[3],
                 int grid_points_eff[3], int max_grid,
                 RefSt& refSt, gnuplot_ctrl* h2, gnuplot_ctrl* h3);

//----------------------------------------------------------------------------------------
//         Brick C: Find the intersection of a EML2-SEMLi connection
//                     with a certain Pk section x = cst
//----------------------------------------------------------------------------------------
/**
 *  \brief Find the intersection of a EML2-SEMLi connection contained in
 *         y_traj_n/t_traj_n with a certain Pk section x = cst defined by refSt.
 **/
int xpkemlisemli(double ye[6], double* te, double* t_traj_n, double** y_traj_n,
                  int man_index, RefSt& refSt);

/**
 *  \brief Get the complementary coordinates associated to the coordinates coord_type.
 **/
int comp_coord_typ(int coord_type);

//========================================================================================
//
//         Refining trajectories from continuation procedures
//
//========================================================================================
/**
 *  \brief Refine complete EMLi-SEMLi connections from semi-analyticalcontinuation results
 *         The orbits at both ends are included in the refinement process. Each trajectory
 *         is then refined in a higher-fidelity model (JPL DE430).
 **/
int reffromcontemlisemli(RefSt& refSt);

//========================================================================================
//
//         Refinement of solutions: Complete trajectory
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the WHOLE trajectory
 *         (EML2 orbit + manifold leg + SEMLi orbit).
 **/
int comprefemlisemli3d(int grid_freq_days[3], int coord_type,
                  Orbit& orbit_EM, Orbit& orbit_SEM,
                  RefSt& refSt, int label, int isFirst);

//========================================================================================
//
//         Refinement of solutions: to JPL
//
//========================================================================================
/**
 *  \brief Refine a given output of comprefemlisemli3d into JPL ephemerides.
 **/
int jplref3d(int coord_type, RefSt& refSt, int label, int isFirst);

/**
 *  \brief Refine a given output of comprefemlisemli3d into Inertial Coordinates, then into
 *         JPL coordinates.
 **/
int comptojplref3d(int coord_type, RefSt& refSt);

//----------------------------------------------------------------------------------------
//         Refinement of solutions: to JPL - Subroutines
//----------------------------------------------------------------------------------------
/**
 * \brief Computes a first guess in ephemerides coordinates, switching between the
 *        Earth-Moon and Sun-Earth plane of motion when the discrepancy between the two
 *        coordinates system is minimal.
 **/
int jplfg3d_switch(double** y_traj_n, double* t_traj_n,
                     double** y_traj_jpl, double* t_traj_jpl,
                     int final_index, int coord_type,
                     int comp_type, int coord_int,
                     double et0, double tsys0, double tsys0_comp);

/**
 * \brief Same as jplfg3d_switch, with a refinement of the position of the minimum.
 **/
int jplfg3d_super_switch(double** y_traj_n, double* t_traj_n,
                           double** y_traj_jpl, double* t_traj_jpl,
                           int final_index, int coord_type,
                           int comp_type, int coord_int,
                           int mRef,
                           double et0, double tsys0, double tsys0_comp);

/**
 * \brief Same as jplfg3d_switch, with an interpolation before and after the switching
 *        point.
 **/
int jplfg3d_interpolation(double** y_traj_n, double* t_traj_n,
                            double** * y_traj_jpl, double** t_traj_jpl,
                            int final_index, int coord_type,
                            int comp_type, int coord_int,
                            int mRef,
                            double et0, double tsys0, double tsys0_comp);

//----------------------------------------------------------------------------------------
//         Refinement of solutions: to JPL - Tests
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes only a SEMLi orbit and test a JPL refinement.
 **/
int compref3d_test_seml_synjpl(int man_grid_size_t,
                                   int coord_type,
                                   Orbit& orbit_SEM,
                                   RefSt refSt);

/**
 *  \brief Computes only an EML2 orbit and test a JPL refinement.
 **/
int compref3d_test_eml_synjpl(int man_grid_size_t,
                                  int coord_type,
                                  Orbit& orbit_EM,
                                  RefSt refSt);


/**
 *  \brief Computes only a EML2-SEMLi connection and test a JPL refinement, in synodical coordinates
 **/
int compref3d_test_eml2seml_synjpl(int coord_type);

//========================================================================================
//
//         PLOTTING TRAJECTORIES
//
//========================================================================================
/**
 *  \brief Plotting the notable points in the SEM system, in coord_type coordinates
 **/
int notablePoints_sem(gnuplot_ctrl* h2, int coord_type);

/**
 *  \brief Plotting the notable points in the SEM & EM systems in coord_type coordinates,
 *         as well as complementary coordinates.
 **/
int notablePoints(gnuplot_ctrl* h2, gnuplot_ctrl* h3, int coord_type);

/**
 *  \brief Plot a trajectory, segment by segment, in to complementary coordinate systems.
 **/
int plottrajsegbyseg(double** y_traj, double* t_traj,
                     int final_index, int mPlot, int coord_int,
                     double et0,     double tsys0,
                     int coordsys1,  gnuplot_ctrl* h2,
                     int coordsys2,  gnuplot_ctrl* h3,
                     int color, string title);

/**
 *  \brief Save a trajectory, segment by segment, in to complementary coordinate systems.
 **/
int savetrajsegbyseg(double** y_traj, double* t_traj,
                     int final_index, int mPlot,
                     int coord_int,
                     double et0,      double tsys0,
                     int coordsys1,   int coordsys2,
                     string filename, int label, bool isFirst);

//========================================================================================
//
//          I/O (to set in oolencon_io.cpp)
//
//========================================================================================
/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void writeCONT_txt(Orbit& orbit_EM, Orbit& orbit_SEM,
                   double* te_NCSEM, double* ye_NCSEM,  int isFirst);

/**
 *  \brief Get the length the results of the continuation procedure, in txt file.
 **/
int getLengthCONT_txt(double t0xT);

/**
 *  \brief Reads the results of the continuation procedure, in txt file.
 **/
int readCONT_txt(double*  t0_CMU_EM, double*   tf_CMU_EM,
                 double** si_CMU_EM, double** si_CMS_SEM,
                 double** z0_CMU_NCEM, double** z0_CMS_NCSEM,
                 double* tethae, double** ye_NCSEM,
                 double* H0_NCEM, double* He_NCEM,
                 double* H0_NCSEM, double* He_NCSEM,
                 double tr0, int fsize);

/**
 *  \brief Save a given solution as a complete trajectory
 **/
int writeCONT_bin(RefSt& refSt, string filename_traj, int dcs, int coord_type,
                   double** y_traj_n, double* t_traj_n, int man_index, int mPlot,
                   Orbit &orbit_EM, Orbit &orbit_SEM, int label,
                   bool isFirst, int comp_orb_eml, int comp_orb_seml);
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void toCelestiaFormat(string filename);

#endif // OOLENCON_H_INCLUDED
