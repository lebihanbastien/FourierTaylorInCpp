#ifndef OOLENCON_H_INCLUDED
#define OOLENCON_H_INCLUDED

#include "oolenconref.h"
#include "oolencon.h"
#include "tinyfiledialogs.h"


//========================================================================================
//
//          Structure for projection parameters
//
//========================================================================================
/**
 *  \struct ProjSt
 *  \brief  Define a given projection structure, with a set of parameters
 **/
typedef struct ProjSt ProjSt;
struct ProjSt
{
    //------------------------------------------------------------------------------------
    // ... filename handling is private
    //------------------------------------------------------------------------------------
    private:

    //filename for the data file output of routines such as compute_grid_CMU_EM_3D/int_proj_CMU_EM_on_CM_SEM_3D
    string filename_output;


    //------------------------------------------------------------------------------------
    // Most of the structure is public to ease use...
    //------------------------------------------------------------------------------------
    public:

    int    OFTS_ORDER, LI_EM, LI_SEM, LI_START, LI_TARGET;
    int    IO_HANDLING, ISPAR;
    double TM, RMIN, RMAX, TMIN, TMAX, TLIM[2], GLIM_SI[4][2];
    int    TSIZE, GSIZE_SI[4], MSIZE, NSMIN, NOD, PRIMARY;
    double YNMAX, SNMAX, dHd, dt;
    double hyp_epsilon_eml2, hyp_epsilon_seml2;
    string plot_folder;
    string FILE_CU, FILE_PCU;
    double CENTER[3], RPS_NCEM;

    /**
     *  \brief Constructor for ProjSt
     **/
     ProjSt(int OFTS_ORDER_, int LI_EM_, int LI_SEM_,
            int LI_START_, int LI_TARGET_,
            int IO_HANDLING_, int ISPAR_,
            double hyp_epsilon_eml2_, double hyp_epsilon_seml2_,
            double RPS, CSYS *cs):
            OFTS_ORDER(OFTS_ORDER_), LI_EM(LI_EM_), LI_SEM(LI_SEM_),
            LI_START(LI_START_), LI_TARGET(LI_TARGET_),
            IO_HANDLING(IO_HANDLING_), ISPAR(ISPAR_),
            dHd(-1.0), dt(0.001),
            hyp_epsilon_eml2(hyp_epsilon_eml2_),
            hyp_epsilon_seml2(hyp_epsilon_eml2_),
            plot_folder(cs->F_PLOT)
    {
        // filename_output is initialized empty
        filename_output = "";

        //Center
        CENTER[0] = (cs->li == 1)? 1:-1;
        CENTER[1] = 0; CENTER[2] = 0;

        //Pk section
        RPS_NCEM = RPS;
        if(RPS_NCEM < 0) RPS_NCEM = cs->r3BSOI;
    }

    /**
     *  \brief Update the filename and return it.
     **/
    string get_and_update_filename(string filename_bash, int type, int mode)
    {
        // If TSIZE == 0, there is only one time in the time vector, and we can add the
        // initial time to the default filename

        double t0 = (TSIZE == 0)? RMIN:-1;
        this->filename_output = get_filenameCUM(this->IO_HANDLING, this->plot_folder,
                                                filename_bash, this->OFTS_ORDER, type,
                                                this->LI_TARGET, t0, dHd, mode);
        return this->filename_output;
    }

    string get_filename(int type, int mode)
    {
        //Updating the filename if it does not exist.
        if(this->filename_output.empty())
        {
            this->get_and_update_filename("", type, mode);
        }

        //Return its value
        return this->filename_output;
    }

};



//========================================================================================
//
//          Computation of the CMU about SEMLi
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about SEMLi,
 *         in the QBCP model. The initial conditions (IC) are computed in a 5-dimensional
 *         box: one dimension for the starting time, four dimensions for the
 *         parameterization of the Center Manifold (s1 to s4 coordinates). The RCM
 *         coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:      the value in RCM coordinates on the unstable direction s5.
 *  \param projSt.TLIM:     the min/max starting time (in SEM units) in the IC box.
 *  \param projSt.TSIZE:    the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI:  the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                          in the IC box.
 *  \param projSt.GSIZE_SI: the number of points on the  s1, s2, s3, s4 values  grids
 *                          in the IC box.
 *  \param projSt.ISPAR;    if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/SEM/L2/cu_3d_order_16.bin"
 **/
int compute_grid_CMU_SEM_3D(double dist_to_cm, ProjSt& projSt);


//========================================================================================
//
//         Projection on the CM/CMS/CMU of EMLi
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM. Then, each point on the
 *         integration grid is projected on the Center Manifold CM_EM_NC about EMLi.
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM:          the maximum integration time on each leg, in SEM units.
 *
 *  \param projSt.MSIZE:       the number of points on each manifold leg.
 *
 *  \param projSt.NOD:         the number of dimensions on which the distance of
 *                             projection is computed (usually either 3 (the physical
 *                             distance) or 6 (the whole phase space)).
 *
 *  \param projSt.ISPAR:       if TRUE, the computation is parallelized.
 *
 *  \param projSt.YNMAX:       the maximum norm in NCEM coordinates for which a given
 *                             state on the integration grid is projected on CM_EM_NC
 *                             More precisely: for a given state y along the manifold leg,
 *                             if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                             Otherwise, it is considered too far away from EMLi to be
 *                             a good candidate for projection.
 *
 *  \param projSt.SNMAX:       the maximum norm in RCM EM coordinates for which a given
 *                             projection state on the CM of EMLi (CM_SEM_NC) is
 *                             computed back in NCSEM coordinates. More precisely, for a
 *                             given state y in NCSEM coordinates, the result of the
 *                             projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                             coordinates.
 *                             if norm(sproj, 4) < projSt.SNMAX, the computation
 *                             yproj = CM_EM_NC(sproj, t) is performed. Otherwise, the
 *                             state sproj is considered too far away from the RCM origin
 *                             to be a good candidate - it is out of the domain of
 *                             practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form:
 *          "plot/QBCP/SEM/L2/projcu_3d_order_16.bin".
 **/
int int_proj_CMU_SEM_on_CM_EM_3D(ProjSt& projSt);


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

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         two-dimensional box: one dimension for the starting time (t0), one dimension
 *         for the parameterization of the Center Manifold (s1 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *          Two possibilities exist:
 *          (i) if projSt.dHd > 0, then a fixed value of energy
 *          at departure is desired by the user. Thereore, the value s3 = f(s1, t0)
 *          is computed and the IC are of the form (s1, 0, f(s1, t0), 0, dist_to_cm)
 *
 *          (ii) if projSt.dHd <= 0, the IC are of the form (s1, 0, 0 0, dist_to_cm)
 *
 *
 *  \param dist_to_cm:           the value in RCM coordinates on the unst. direction s5.
 *  \param projSt.TLIM[0]:       the minimum starting time (in EM units) in the IC box.
 *  \param projSt.TLIM[1]:       the maximum starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:         the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI[0][0]: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[0][1]: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GSIZE_SI[0]:   the number of points on the s1 grid in the IC box.
 *  \param projSt.ISPAR:         if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16_dH_0.005.bin"
 **/
int compute_grid_CMU_EM_dH(double dist_to_cm, ProjSt& projSt);

//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEMLi
//
//========================================================================================
/**
 *  \brief Projection subroutine, used in all the subsequent routines in this section.
 *         All outputs in projResSt are updated
 *         (see the declaration of this structure for details).
 *
 *         This routine performs the following steps:
 *
 *          - It integrates the initial state at emli, stored in projResSt
 *          - It compute the minimum and the argminimum distance of projection along the
 *            corresponding trajectory, targeting the center manifold invman_target
 *          - The information relative to this min/argmin are stored in projResSt.
 *
 *         The coordinates within which the results are saved are encoded in ProjResSt.
 **/
int proj_subroutine(ProjResSt& projResSt, Invman& invman_target, ProjSt& projSt);

/**
 *  \brief Projection subroutine, used in all the subsequent routines in this section.
 *         All outputs in projResSt are updated
 *         (see the declaration of this structure for details).
 *
 *         This routine performs the following steps:
 *
 *          - It integrates the initial state at emli, stored in projResSt
 *          - It compute the minimum and the argminimum distance of projection along the
 *            corresponding trajectory, targeting the center manifold invman_SEM
 *          - The information relative to this min/argmin are stored in projResSt.
 *
 *
 *         Note: this routine is the old version from the EML to SEML case.
 **/
int proj_subroutine_old_from_EML(ProjResSt& projResSt, Invman& invman_SEM, ProjSt& projSt);

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
 * "plot/QBCP/EM/L2/projcu_order_16.bin"
 **/
int int_proj_CMU_EM_on_CM_SEM(ProjSt &projSt);

/**
 *  \brief Same as int_proj_CMU_EM_on_CM_SEM but for the outputs from the routine
 *         compute_grid_CMU_EM_dH.
 *
 *         Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM_dH.
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
 *      "plot/QBCP/EM/L2/projcu_order_16_dH_0.005.bin"
 **/
int int_proj_CMU_EM_on_CM_SEM_dH(ProjSt &projSt);

/**
 *  \brief Detection of the possible connections between a set of orbits about EMLi
 *         and the center manifold at SEMLj.
 *         The routine works with a system of seefs in RCM coordinates.
 *
 *         For all projSt.GSIZE_SI[0] values of s1 in
 *         [projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1]], we define the following initial
 *         conditions:
 *                       s0 = (s1, projSt.GLIM_SI[1][0], -s1, projSt.GLIM_SI[3][0])
 *
 *         Note that s3 = -s1 guarantees that the IC are on the y = 0 plane in NC coord.
 *         (at least in the planar case).
 *
 *
 *
 *         Then, for all projSt.TSIZE value of initial time t0 in
 *         [projSt.TLIM[0], projSt.TLIM[1] ], we cen compute z0 = W(s0, t0). We can then
 *         project z0 on the center manifold, in order to get
 *                  si = (s1, s3, s3, s4) = W^{-1}(z0, t0),
 *
 *         Roughly speaking, si provides IC for a similar orbit but a different starting
 *         phase for the Sun-Earth-Moon system. These IC define a given orbit
 *         whose period is estimated via cmu_orbit_estimate_period. We can then integrate
 *         this orbit on one of its period and look for connections on a fine grid along
 *         the resulting trajectory
 *
 *  The output data are saved in a binary file of the form
 *      "plot/QBCP/EM/L2/projcu_order_16_Orbit.bin"
 **/
int int_proj_ORBIT_EM_on_CM_SEM(ProjSt& projSt, int N);

/**
 *  \brief Detection of the possible connections between a set of orbits about EMLi
 *         and the center manifold at SEMLj.
 *         We define the following initial
 *         conditions:
 *         s0 = (projSt.GLIM_SI[0][0], projSt.GLIM_SI[1][0], projSt.GLIM_SI[2][0], projSt.GLIM_SI[3][0])
 *
 *         Then, for all projSt.TSIZE value of initial time t0 in
 *         [projSt.TLIM[0], projSt.TLIM[1] ], we cen compute z0 = W(s0, t0). We can then
 *         project z0 on the center manifold, in order to get
 *                  si = (s1, s3, s3, s4) = W^{-1}(z0, t0),
 *
 *         Roughly speaking, si provides IC for a similar orbit but a different starting
 *         phase for the Sun-Earth-Moon system. These IC define a given orbit
 *         whose period is estimated via cmu_orbit_estimate_period. We can then integrate
 *         this orbit on one of its period and look for connections on a fine grid along
 *         the resulting trajectory
 *
 *  The output data are saved in a binary file of the form
 *      "plot/QBCP/EM/L2/projcu_order_16_Orbit.bin"
 **/
int int_proj_SINGLE_ORBIT_EM_on_CM_SEM(ProjSt& projSt, int Nperiods);

//========================================================================================
//
//         Initial conditions for projection of single orbits
//
//========================================================================================
/**
 *  \brief Computes the stroboscopic map on N+1 points along a trajectory in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] = t0 + k*SEML.us->T, for all k in [0, N],
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_strob(double* tNCE, double** yNCE, double** sRCM, double st0[], double t0, int N, int isPar, double hyp_epsilon);

/**
 *  \brief Computes N+1 points along a trajectory in the QBCP, on NPeriods periods T
 *         of the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] is the time, in [0 NPeriods*T]
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_orbit(double* tNCE, double** yNCE, double** sRCM, double st0[], double t0, int N, int NPeriods,  int isPar, double hyp_epsilon);

/**
 *  \brief Estimates the period T (in adimensionalized units) of a given orbit Orbit
 *         in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *         The period is estimated as follows:
 *         The initial condition are z0v = (x0, y0, z0, px0, py0, pz0) = W(s0, t0).
 *         The polar argument of (x0, y0) is computed.
 *         And the integration is stopped when a similar argument is obtained.
 *         The time of the event is the desired time T.
 *
 *         Moreover, the number of points N is computed as the desired number of points to
 *         achieve a frequency of dt in [0 T].
 *         Hence, N = floor(T/dt);
 **/
int cmu_orbit_estimate_period(const double st0[], double t0, double* T, int* N, double dt, Orbit& orbit);

/**
 *  \brief Computes N+1 points along a trajectory in the QBCP, on one period T
 *         of the orbit Orbit, in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] is the time, in [0 T]
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_orbit_on_one_period(Orbit& orbit, double* tNCE, double** yNCE, double** sRCM, const double st0[], double t0, double T, int N, int isPar, double hyp_epsilon);

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

;/**
 *  \brief Computes the best trajectories from int_proj_ORBIT_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of i.
 *         The final conditions vary in the paramerization of the CMS of SEMLi.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,
 *         the refSt.coord_type must be NCSEM. However, the user can put
 *         other coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines
 *         are not used.
 *
 *         so stands for single orbit.
 **/
int sorefemlisemli(RefSt& refSt);

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

/**
 *  \brief Selects good initial conditions for EML2-to-SEMLi connections, searching
 *         through data files produced by int_proj_ORBIT_EM_on_CM_SEM
 **/
int soselectemlisemli(RefSt& refSt, ProjResClass& subSt);

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
 *  \brief Find the intersection of a EML2-SEMLi connection contained in
 *         y_traj_n/t_traj_n with a certain Pk section x = cst defined by refSt.
 *         Once the intersection is found, it is incorporated in the sequence of patch points,
 *         in place of a given point, at position newpos
 **/
int xpk_emlis2emli(double ye[6], double* te, double* t_traj_n, double** y_traj_n,
                int *newpos, int man_index, RefSt& refSt);

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
int jplref3d(int coord_type, RefSt& refSt, int label, int isFirst, string filename_in);

/**
 *  \brief Refine a given output of comprefemlisemli3d into Inertial Coordinates, then into
 *         JPL coordinates.
 **/
int comptojplref3d(int coord_type, RefSt& refSt, string filename_in);

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
int compref3d_test_eml2seml_synjpl(int coord_type, string filename_in);

//========================================================================================
//
//         PLOTTING TRAJECTORIES
//
//========================================================================================
/**
 *  \brief Plotting the notable points in the SEM system, in coord_type coordinates
 **/
int notablePoints_sem(gnuplot_ctrl* h2, int coord_type, int isPlot);

/**
 *  \brief Plotting the notable points in the SEM & EM systems in coord_type coordinates,
 *         as well as complementary coordinates.
 **/
int notablePoints(gnuplot_ctrl* h2, gnuplot_ctrl* h3, int coord_type, int isPlot);

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


/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void toCelestiaFormat(string filename);

#endif // OOLENCON_H_INCLUDED
