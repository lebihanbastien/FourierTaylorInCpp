//========================================================================================
// Includes
//========================================================================================
// C++
#include <omp.h>
#include <iostream>
// Custom
#include "FTA.h"
#include "Ofsc.h"
#include "Oftsc.h"
#include "pmcoc.h"
#include "ode.h"
#include "vf.h"
#include "env.h"
#include "timec.h"
#include "single_orbit.h"
#include "oolencon.h"
#include "ephemerides.h"
#include "Invman.h"
#include "Orbit.h"
#include "Config.h"

#include "tinyfiledialogs.h"

//C
extern "C" {
#include "gnuplot_i.h"
#include "SpiceUsr.h"
}

//Standard namespaces
using namespace std;

//========================================================================================
// Constants
//========================================================================================
// TYPE OF COMPUTATIONS
#define COMP_CM_EML2_TO_CM_SEML           0    //EMLi Center Manifold to SEMLi Center Manifold
#define COMP_CM_EML2_TO_CMS_SEML          1    //EMLi Center Manifold to SEMLi Center-Stable Manifold
#define COMP_SINGLE_ORBIT                 2    //just some example of solutions
#define COMP_CM_EML2_TO_CMS_SEML_READ     3    //Read
#define COMP_VF_TEST                      4    //Test of the QBCP vector field. Should probably be put in OOFTDA
#define COMP_CM_EML2_TO_CM_SEML_REFINE    5    //EML2 Center Manifold to SEMLi Center Manifold
#define COMP_EPHEMERIDES_TEST             6    //Ephemerides test
#define COMP_CM_EML2_TO_CM_SEML_3D        7    //EML2 Center Manifold to SEMLi Center Manifold in 3D
#define COMP_VOFTS_TO_VOFTS               8    //Store CS/CU into one-dimensionnal series to gain memory
#define COMP_test_INVMAN                  9    //Test of the new invariant manifold implementation
#define COMP_REF_JPL                      10   //Refine to JPL ephemerides
#define COMP_CM_EML2_TO_CM_SEML_H         11   //Planar EMLi Center Manifold to SEMLi Center Manifold, at a given energy
#define COMP_ORBIT_EML2_TO_CM_SEML        12   //EMLi Center Manifold to SEMLi Center Manifold, on a given set of orbits
#define COMP_SINGLE_ORBIT_EML2_TO_CM_SEML 13   //EMLi Center Manifold to SEMLi Center Manifold, on a single orbit

#define COMP_CMU_SEML_TO_CM_EML           14   //SEMLi Center-Unstable Manifold to EMLj Center Manifold
#define COMP_CMU_SEML_TO_CMS_EML          15   //SEMLi Center-Unstable Manifold to EMLj Center-Stable Manifold

/************* NOTES *********************************************************************
 Notes from Reunion with Josep (15/12/2015)

+ Book I for an example of continuation.

+ @TODO: see if we can fasten the JPL and QBCP vf codes! Started to be sloooow. Can we also fasten ode78?
+ @TODO: in reffromcontemlisemli, try to move around the et0. Is there a connection with the size of the final refined EML2 orbit?
+ @TODO: in reffromcontemlisemli, what happens if we search for et0 in more than 3 month span? we can probably find a month for which the size is roughly equivalent

*****************************************************************************************/
int main(int argc, char** argv)
{
    //====================================================================================
    //Tests
    //====================================================================================


    //====================================================================================
    //Declare configuration parameters
    //====================================================================================
    // General parameters (orders, etc)
    int COMP_TYPE, NUM_THREADS, MODEL_TYPE, LI_EM, LI_SEM, ISPAR, IO_HANDLING;
    double HYP_EPSILON_EML2, HYP_EPSILON_SEML2;

    //------------------------------------------------------------------------------------
    // The variable index contains the index of the current argument of
    // the main routine. First index is 1 because index 0 is the application path
    //------------------------------------------------------------------------------------
    int index = 1;

    //====================================================================================
    //Update the general parameters
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Check if arguments have been passed
    //------------------------------------------------------------------------------------
    if(argc == 1) //no argument were passed
    {
        //--------------------------------------------------------------------------------
        // Type of computation
        //--------------------------------------------------------------------------------
        COMP_TYPE   = COMP_REF_JPL;

        //--------------------------------------------------------------------------------
        // Model and libration points
        //--------------------------------------------------------------------------------
        MODEL_TYPE  = M_QBCP;
        LI_EM       = 2;
        LI_SEM      = 2;

        //--------------------------------------------------------------------------------
        //Orders
        //--------------------------------------------------------------------------------
        OFTS_ORDER  = 16;
        if(MODEL_TYPE == M_RTBP) OFS_ORDER  = 0;
        else OFS_ORDER  = 30;

        //--------------------------------------------------------------------------------
        // Parallel computation parameters
        //--------------------------------------------------------------------------------
        ISPAR       = 1;
        NUM_THREADS = 4;

        //--------------------------------------------------------------------------------
        // Hyperbolic component (default values)
        //--------------------------------------------------------------------------------
        HYP_EPSILON_EML2  = HYP_EPSILON_EML2_DEFAULT;
        HYP_EPSILON_SEML2 = HYP_EPSILON_SEML2_DEFAULT;

        //--------------------------------------------------------------------------------
        // I/O handling
        //--------------------------------------------------------------------------------
        IO_HANDLING = IO_DEFAULT;
    }
    else  //arguments were passed
    {
        //--------------------------------------------------------------------------------
        //Orders
        //--------------------------------------------------------------------------------
        OFTS_ORDER  = atoi(argv[index++]);
        OFS_ORDER   = atoi(argv[index++]);

        //--------------------------------------------------------------------------------
        // Type of computation
        //--------------------------------------------------------------------------------
        COMP_TYPE   = atoi(argv[index++]);

        //--------------------------------------------------------------------------------
        // Model and libration points
        //--------------------------------------------------------------------------------
        MODEL_TYPE  = atoi(argv[index++]);
        LI_EM       = atoi(argv[index++]);
        LI_SEM      = atoi(argv[index++]);

        //--------------------------------------------------------------------------------
        // Parallel computation parameters
        //--------------------------------------------------------------------------------
        ISPAR       = atoi(argv[index++]);
        NUM_THREADS = atoi(argv[index++]);

        //--------------------------------------------------------------------------------
        // Hyperbolic component
        //--------------------------------------------------------------------------------
        HYP_EPSILON_EML2  = atof(argv[index++]);
        HYP_EPSILON_SEML2 = atof(argv[index++]);

        //--------------------------------------------------------------------------------
        // I/O handling
        //--------------------------------------------------------------------------------
        IO_HANDLING = atoi(argv[index++]);
    }

    //====================================================================================
    // General settings
    //====================================================================================
    cout << setiosflags(ios::scientific) << setprecision(15);
    omp_set_num_threads(NUM_THREADS);

    //Target is LI_SEM by default
    int LI_TARGET = LI_SEM;
    //Start is LI_EM by default
    int LI_START  = LI_EM;


    //====================================================================================
    // Parameters that depend on the type of computation
    //====================================================================================
    int model = 0, fwrk = 0, isNorm = 0;
    int pms_EM = 0, pms_SEM = 0, mType_EM = 0, mType_SEM = 0, reduced_nv = 0;

    switch(COMP_TYPE)
    {
        //================================================================================
        // Test and additionnal computations
        //================================================================================
    case COMP_test_INVMAN:
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER_S;
        break;

    case COMP_EPHEMERIDES_TEST:
    case COMP_VOFTS_TO_VOFTS:
        model     = MODEL_TYPE;
        fwrk      = F_SEM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER_S;
        break;

        //================================================================================
        // Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML:
    case COMP_CM_EML2_TO_CM_SEML_3D:
    case COMP_CM_EML2_TO_CM_SEML_REFINE:
    case COMP_CM_EML2_TO_CM_SEML_H:
    case COMP_ORBIT_EML2_TO_CM_SEML:
    case COMP_SINGLE_ORBIT_EML2_TO_CM_SEML:
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER;
        break;

        //================================================================================
        // Projection CMU SEMLi to CM EMLj
        //================================================================================
    case COMP_CMU_SEML_TO_CM_EML:
        model     = MODEL_TYPE;
        fwrk      = F_SEM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER_U;
        // Target & start are inversed
        LI_TARGET = LI_EM;
        LI_START  = LI_SEM;
        break;

        //================================================================================
        // Refinement CMU EML2 to CMS SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CMS_SEML:
    case COMP_CM_EML2_TO_CMS_SEML_READ:
    case COMP_REF_JPL:
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER_S;
        break;

        //================================================================================
        // Refinement CMU SEMLi to CMS EMLj
        //================================================================================
    case COMP_CMU_SEML_TO_CMS_EML:
        model     = MODEL_TYPE;
        fwrk      = F_SEM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_S;
        mType_SEM = MAN_CENTER_U;
        break;

        //================================================================================
        // Just some examples of solutions
        //================================================================================
    case COMP_SINGLE_ORBIT:
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER;
        break;

        //================================================================================
        // Test of the vector fields
        //================================================================================
    case COMP_VF_TEST:
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER;
        break;

    default:
        cout << "main. Warning: unknown COMP_TYPE!" << endl;
    }


    //====================================================================================
    // Initialization of the environnement
    // Mandatory to perform any computation
    //====================================================================================
    initialize_environment(LI_EM, LI_SEM, isNorm, model, fwrk,
                           pms_EM, pms_SEM, mType_EM, mType_SEM);


    //====================================================================================
    //Update the parameters specific to the type of computation (COMP_TYPE)
    //====================================================================================

    //------------------------------------------------------------------------------------
    // Structures
    //------------------------------------------------------------------------------------
    // Projection parameters (in structure)
    ProjSt projSt(OFTS_ORDER, LI_EM, LI_SEM, LI_START, LI_TARGET, IO_HANDLING, ISPAR,
                  HYP_EPSILON_EML2_DEFAULT, HYP_EPSILON_SEML2_DEFAULT, SEML.cs->F_PLOT);
    // Refinement parameters (in structure)
    RefSt refSt(OFTS_ORDER, LI_EM, LI_SEM, LI_START, LI_TARGET, IO_HANDLING, SEML.cs->F_PLOT);


    //------------------------------------------------------------------------------------
    //Check if arguments have been passed
    //------------------------------------------------------------------------------------
    if(argc == 1) //no argument were passed
    {
        //================================================================================
        // Projection parameters
        //================================================================================
        //--------------------------------------------------------------------------------
        // Hyperbolic component (default values)
        //--------------------------------------------------------------------------------
        projSt.hyp_epsilon_eml2  = HYP_EPSILON_EML2;
        projSt.hyp_epsilon_seml2 = HYP_EPSILON_SEML2;

        //--------------------------------------------------------------------------------
        //Time grid: min, max and number of points on the grid
        //--------------------------------------------------------------------------------
        projSt.RMIN  = 0.0;
        projSt.RMAX  = 0.5;
        projSt.TMIN  = projSt.RMIN*SEML.us->T;
        projSt.TMAX  = projSt.RMAX*SEML.us->T;

        projSt.TSIZE    = 5;
        projSt.TLIM[0]  = projSt.TMIN;
        projSt.TLIM[1]  = projSt.TMAX;

        //--------------------------------------------------------------------------------
        // Configuration (s1, s2, s3, s4) grid
        //--------------------------------------------------------------------------------
        projSt.GLIM_SI[0][0] = 10;
        projSt.GLIM_SI[0][1] = +30;

        projSt.GLIM_SI[1][0] = +0;
        projSt.GLIM_SI[1][1] = +0;

        projSt.GLIM_SI[2][0] = -30;
        projSt.GLIM_SI[2][1] = +30;

        projSt.GLIM_SI[3][0] = +0;
        projSt.GLIM_SI[3][1] = +0;

        projSt.GSIZE_SI[0]   = 0;
        projSt.GSIZE_SI[1]   = 0;
        projSt.GSIZE_SI[2]   = 50;
        projSt.GSIZE_SI[3]   = 0;

        //--------------------------------------------------------------------------------
        // Primary family
        //--------------------------------------------------------------------------------
        projSt.PRIMARY = 0;

        //--------------------------------------------------------------------------------
        // Energy
        //--------------------------------------------------------------------------------
        projSt.dHd = -1.0;

        //--------------------------------------------------------------------------------
        // Time frequency, in %T
        //--------------------------------------------------------------------------------
        projSt.dt = 0.005;

        //--------------------------------------------------------------------------------
        //Stable parameters (are not supposed to change)
        //--------------------------------------------------------------------------------
        projSt.TM    = 12.0*SEML.us->T; // Maximum integration time
        projSt.MSIZE = 500;             // Number of points on each trajectory
        projSt.NSMIN = 20;              // Number of sorted solutions
        projSt.YNMAX = 0.6;             // The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
        projSt.SNMAX = 0.6;             // The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
        projSt.NOD   = 6;               // Number of dimensions on which we compute the norm of the projection

        //--------------------------------------------------------------------------------
        // Filenames (used only if IO_HANDLING==$IO_BASH)
        //--------------------------------------------------------------------------------
        projSt.FILE_CU  = "cu.bin";
        projSt.FILE_PCU = "projcu.bin";


        //================================================================================
        // Refinement parameters
        //================================================================================

        //--------------------------------------------------------------------------------
        // Parameters that change often
        //--------------------------------------------------------------------------------
        //rk: set REF_CONT_D_HARD_CASE for difficult cases
        //with REF_CONT_D (ex: EML2-SEMLi via SEML1...)
        refSt.type          = REF_CONT;                       // Type of refinement
        refSt.dim           = REF_PLANAR;                      // Type of dimensions planar or 3d?
        refSt.t0xT_des      = 0.99;                            // Initial time (xT)
        refSt.t0_des        = refSt.t0xT_des*SEML.us->T;     // Initial time

        // Direction of the continuation procedure
        refSt.isDirUD       = 0;                  // is it user defined?
        refSt.Dir           = +1;                 // if not, +1 or -1

        // Domain of search for the first guess
        refSt.si_CMU_EM_MIN[0] = -40;
        refSt.si_CMU_EM_MAX[0] = +40;

        refSt.si_CMU_EM_MIN[1] = +0;
        refSt.si_CMU_EM_MAX[1] = +0;

        refSt.si_CMU_EM_MIN[2] = -40;
        refSt.si_CMU_EM_MAX[2] = +40;

        refSt.si_CMU_EM_MIN[3] = +4;
        refSt.si_CMU_EM_MAX[3] = +4;

        // Or, if we want the user to define such domain:
        refSt.isLimUD       =  0;

        // Domain of search for the seed of the first guess
        refSt.si_SEED_EM_MIN[0] = -40;
        refSt.si_SEED_EM_MAX[0] = +40;

        refSt.si_SEED_EM_MIN[1] = +0;
        refSt.si_SEED_EM_MAX[1] = +0;

        refSt.si_SEED_EM_MIN[2] = -40;
        refSt.si_SEED_EM_MAX[2] = +40;

        refSt.si_SEED_EM_MIN[3] = 0;
        refSt.si_SEED_EM_MAX[3] = 0;

        //Limits for the time of flight during transfers - not used if -1
        refSt.tof_MIN       = -1;
        refSt.tof_MAX       = -1;

        // Values for crossings
        refSt.crossings     = -1;

        // Maximum projection distance allowed during subselection
        refSt.pmax_dist_SEM = 1e-0;

        // Number of steps in the continuation procedure
        refSt.cont_step_max    = +450;            // with fixed times
        refSt.cont_step_max_vt = +2;              // with variable times

        // Initial step in the continuation procedure
        refSt.ds0    = 1e-3;                      //with fixed time
        refSt.ds0_vt = (LI_EM ==1)?  3e-2:1e-2;   //with variable time

        // Desired number of iterations in Newton's method in the continuation procedure
        refSt.nu0 = 2;          //with fixed time
        refSt.nu0_vt = 3;       //with variable time

        //User parameters
        refSt.isFlagOn      = 1;                  // do we have steps in the procedure - asking the user to press enter to go on?
        refSt.isPlotted     = 1;                  // do we plot the results during the computation?
        refSt.isSaved       = 1;                  // do we save the results in data files?
        refSt.isFromServer  = 0;                  // does the raw data comes from server files?
        refSt.isPar         = 0;                  //is parallel computation allowed?

        //Maximum angle around SEMLi if REF_COND_T is used (in degrees)
        refSt.thetaMax      = 240;                //should be a multiple of 90°

        //--------------------------------------------------------------------------------
        // Parameters that are stable
        //--------------------------------------------------------------------------------
        refSt.isDebug       = 0;                        // if yes, additionnal tests are performed
        refSt.gridSize      = 20;                       // number of points on the refinement grid. 20 is taken by heuristics.
        refSt.mplot         = 200;                      // number of points per plot between to pach points (e.g. total plot points is gridSize*mplot)

        refSt.time          = REF_VAR_TN;               // type of constraints on the times in REF_CONT
        refSt.grid          = REF_FIXED_GRID;           // type of grid
        refSt.termination   = REF_COND_T;               // termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
        refSt.coord_type    = NCSEM;                    // coordinates system in the refinement procedure

        // Maximum/Minimum step in the continuation procedure
        refSt.dsmin         = 1e-4;                     //with fixed time
        refSt.dsmin_vt      = 1e-4;                     //with variable time
        refSt.dsmax         = 1e-1;                       //with fixed time
        refSt.dsmax_vt      = 1e-1;                       //with variable time

        refSt.xps           = (LI_SEM == 1)? +0.7:-0.7; // position of the poincaré section in NCSEM coordinates
        refSt.isJPL         = 1;                        // is the JPL refinement performed when possible?
        refSt.djplcoord     = NJ2000;                   // coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with NJ2000
        refSt.sidim         = 0;                        // 0 or 2 - component of s0 that stays constant when t0 is free

        // Sampling frequencies in REF_COMP (complete trajectory) in days
        refSt.sf_eml2  = 2;                             // orbit at EML2
        refSt.sf_man   = 5;                             // transfer leg
        refSt.sf_seml2 = 10;                            // orbit at SEML2

        // Integration window for each orbit
        refSt.tspan_EM      = +10*SEML.us_em.T;
        refSt.tspan_SEM     = +10*SEML.us_sem.T;

        // Storing the orbits at each step?
        refSt.isSaved_EM    = 0;      //0: don't save, 1: save using projection method
        refSt.isSaved_SEM   = 0;      //0: don't save, 1: save using projection method,
        //2: save using integration in reduced coordinates

        //Energy
        refSt.dH     = 0.0;
        refSt.pkpos  = 0;

        //Type of time selection
        refSt.typeOfTimeSelection = TIME_SELECTION_RATIO;

        //--------------------------------------------------------------------------------
        // Filenames (used only if IO_HANDLING==$IO_BASH)
        //--------------------------------------------------------------------------------
        refSt.FILE_PCU       ="projcu_order_16_dest_L2_t0_0.995.bin";
        refSt.FILE_CONT      ="cont_atf.txt";
        refSt.FILE_CONT_TRAJ ="cont_atf_traj.bin";
        refSt.FILE_JPL_BIN   ="cont_jpl.bin";
        refSt.FILE_JPL_TXT   ="cont_jpl.txt";

    }
    else  //arguments were passed
    {
        //================================================================================
        // Parameters that depends on the type of computation
        //================================================================================
        switch(COMP_TYPE)
        {

            //============================================================================
            // 3D Projection CMU EML2 to CM SEMLi
            //============================================================================
        case COMP_CM_EML2_TO_CM_SEML_3D:
        case COMP_CM_EML2_TO_CM_SEML:
        case COMP_CMU_SEML_TO_CM_EML:
        case COMP_CM_EML2_TO_CM_SEML_H:
        case COMP_ORBIT_EML2_TO_CM_SEML:
        case COMP_SINGLE_ORBIT_EML2_TO_CM_SEML:
        {
            //----------------------------------------------------------------------------
            //Time grid: min, max and number of points on the grid
            //----------------------------------------------------------------------------
            projSt.RMIN  = atof(argv[index++]);
            projSt.RMAX  = atof(argv[index++]);
            projSt.TMIN  = projSt.RMIN*SEML.us->T;
            projSt.TMAX  = projSt.RMAX*SEML.us->T;
            projSt.TM    = atof(argv[index++])*SEML.us->T;

            projSt.TSIZE   = atoi(argv[index++]);
            projSt.TLIM[0] = projSt.TMIN;
            projSt.TLIM[1] = projSt.TMAX;

            //----------------------------------------------------------------------------
            // Configuration (s1, s2, s3, s4) grid
            //----------------------------------------------------------------------------
            projSt.GLIM_SI[0][0] = atof(argv[index++]);
            projSt.GLIM_SI[0][1] = atof(argv[index++]);

            projSt.GLIM_SI[1][0] = atof(argv[index++]);
            projSt.GLIM_SI[1][1] = atof(argv[index++]);

            projSt.GLIM_SI[2][0] = atof(argv[index++]);
            projSt.GLIM_SI[2][1] = atof(argv[index++]);

            projSt.GLIM_SI[3][0] = atof(argv[index++]);
            projSt.GLIM_SI[3][1] = atof(argv[index++]);

            projSt.GSIZE_SI[0]   = atoi(argv[index++]);
            projSt.GSIZE_SI[1]   = atoi(argv[index++]);
            projSt.GSIZE_SI[2]   = atoi(argv[index++]);
            projSt.GSIZE_SI[3]   = atoi(argv[index++]);

            //----------------------------------------------------------------------------
            // Primary family
            //----------------------------------------------------------------------------
            projSt.PRIMARY = atoi(argv[index++]);

            //----------------------------------------------------------------------------
            //Stable parameters (are not supposed to change)
            //----------------------------------------------------------------------------
            projSt.MSIZE = atoi(argv[index++]);
            projSt.NSMIN = atoi(argv[index++]);
            projSt.YNMAX = atof(argv[index++]);
            projSt.SNMAX = atof(argv[index++]);
            projSt.NOD   = atoi(argv[index++]);

            //----------------------------------------------------------------------------
            // Energy
            //----------------------------------------------------------------------------
            projSt.dHd  = atof(argv[index++]);

            //----------------------------------------------------------------------------
            // Time frequency, in %T
            //----------------------------------------------------------------------------
            projSt.dt   = atof(argv[index++]);

            //----------------------------------------------------------------------------
            // Filenames (used only if IO_HANDLING==$IO_BASH)
            //----------------------------------------------------------------------------
            projSt.FILE_CU  = argv[index++];
            projSt.FILE_PCU = argv[index++];

            break;
        }


        //================================================================================
        // Refinement CMU EML2 to CMS SEMLi
        //================================================================================
        case COMP_CM_EML2_TO_CMS_SEML:
        case COMP_CMU_SEML_TO_CMS_EML:
        case COMP_REF_JPL:
        {
            //----------------------------------------------------------------------------
            // Parameters that change often
            //----------------------------------------------------------------------------
            //rk: set REF_CONT_D_HARD_CASE for difficult cases
            //with REF_CONT_D (ex: EML2-SEMLi via SEML1...)
            refSt.type          = atoi(argv[index++]);            // Type of refinement
            refSt.dim           = atoi(argv[index++]);            // Type of dimensions planar or 3d?
            refSt.t0xT_des      = atof(argv[index++]);            // Initial time (xT)
            refSt.t0_des        = refSt.t0xT_des*SEML.us->T;      // Initial time

            // Direction of the continuation procedure
            refSt.isDirUD       = atoi(argv[index++]);    // is it user defined?
            refSt.Dir           = atoi(argv[index++]);    // if not, +1 or -1

            // Domain of search for the first guess
            refSt.si_CMU_EM_MIN[0] = atof(argv[index++]);
            refSt.si_CMU_EM_MAX[0] = atof(argv[index++]);

            refSt.si_CMU_EM_MIN[1] = atof(argv[index++]);
            refSt.si_CMU_EM_MAX[1] = atof(argv[index++]);

            refSt.si_CMU_EM_MIN[2] = atof(argv[index++]);
            refSt.si_CMU_EM_MAX[2] = atof(argv[index++]);

            refSt.si_CMU_EM_MIN[3] = atof(argv[index++]);
            refSt.si_CMU_EM_MAX[3] = atof(argv[index++]);

            // Or, if we want the user to define such domain:
            refSt.isLimUD       = atoi(argv[index++]);

            // Domain of search for the seed of the first guess
            refSt.si_SEED_EM_MIN[0] = atof(argv[index++]);
            refSt.si_SEED_EM_MAX[0] = atof(argv[index++]);

            refSt.si_SEED_EM_MIN[1] = atof(argv[index++]);
            refSt.si_SEED_EM_MAX[1] = atof(argv[index++]);

            refSt.si_SEED_EM_MIN[2] = atof(argv[index++]);
            refSt.si_SEED_EM_MAX[2] = atof(argv[index++]);

            refSt.si_SEED_EM_MIN[3] = atof(argv[index++]);
            refSt.si_SEED_EM_MAX[3] = atof(argv[index++]);

            //Limits for the time of flight during transfers - not used if negative
            refSt.tof_MIN       = atof(argv[index++])*SEML.us->T;
            refSt.tof_MAX       = atof(argv[index++])*SEML.us->T;

            // Values for crossings
            refSt.crossings     = atof(argv[index++]);

            // Maximum projection distance allowed during subselection
            refSt.pmax_dist_SEM = atof(argv[index++]);

            // Number of steps in the continuation procedure
            refSt.cont_step_max    = atoi(argv[index++]);  // with fixed times
            refSt.cont_step_max_vt = atoi(argv[index++]);  // with variable times

            // Initial step in the continuation procedure
            refSt.ds0    = atof(argv[index++]);   //with fixed time
            refSt.ds0_vt = atof(argv[index++]);   //with variable time

            // Desired number of iterations in Newton's method in the continuation procedure
            refSt.nu0    = atoi(argv[index++]);   //with fixed time
            refSt.nu0_vt = atoi(argv[index++]);   //with variable time

            //User parameters
            refSt.isFlagOn      = atoi(argv[index++]);     // do we have steps in the procedure - asking the user to press enter to go on?
            refSt.isPlotted     = atoi(argv[index++]);     // do we plot the results during the computation?
            refSt.isSaved       = atoi(argv[index++]);     // do we save the results in data files?
            refSt.isFromServer  = atoi(argv[index++]);     // does the raw data comes from server files?
            refSt.isPar         = ISPAR;                   //is parallel computation allowed?

            //Maximum angle around SEMLi if REF_COND_T is used (in degrees)
            refSt.thetaMax      = atof(argv[index++]);     //should be a multiple of 90°

            //----------------------------------------------------------------------------
            // Parameters that are stable
            //----------------------------------------------------------------------------
            refSt.isDebug       = atoi(argv[index++]);  // if yes, additionnal tests are performed
            refSt.gridSize      = atoi(argv[index++]);  // number of points on the refinement grid. 20 is taken by heuristics.
            refSt.mplot         = atoi(argv[index++]);  // number of points per plot between to pach points (e.g. total plot points is gridSize*mplot)

            refSt.time          = atoi(argv[index++]);  // type of constraints on the times in REF_CONT
            refSt.grid          = atoi(argv[index++]);  // type of grid
            refSt.termination   = atoi(argv[index++]);  // termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
            refSt.coord_type    = atoi(argv[index++]);  // coordinates system in the refinement procedure

            // Maximum/Minimum step in the continuation procedure (for now, not changed by the user)
            refSt.dsmin         = 1e-6;                 //with fixed time
            refSt.dsmin_vt      = 1e-6;                 //with variable time
            refSt.dsmax         = 10;                   //with fixed time
            refSt.dsmax_vt      = 10;                   //with variable time

            refSt.xps           = atof(argv[index++]);  // position of the poincaré section in NCSEM coordinates
            refSt.xps *= (LI_SEM == 1)? +1:-1;
            refSt.isJPL         = atoi(argv[index++]);  // is the JPL refinement performed when possible?
            refSt.djplcoord     = atoi(argv[index++]);  // coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with NJ2000
            refSt.sidim         = atoi(argv[index++]);  // 0 or 2 - component of s0 that stays constant when t0 is free

            // Sampling frequencies in REF_COMP (complete trajectory) in days
            refSt.sf_eml2  = atof(argv[index++]);       // orbit at EML2
            refSt.sf_man   = atof(argv[index++]);       // transfer leg
            refSt.sf_seml2 = atof(argv[index++]);       // orbit at SEML2

            // Integration window for each orbit
            refSt.tspan_EM      = atof(argv[index++])*SEML.us_em.T;
            refSt.tspan_SEM     = atof(argv[index++])*SEML.us_sem.T;

            // Storing the orbits at each step?
            refSt.isSaved_EM    = atoi(argv[index++]);  //0: don't save, 1: save using projection method
            refSt.isSaved_SEM   = atoi(argv[index++]);  //0: don't save, 1: save using projection method,
            //2: save using integration in reduced coordinates

            //Type of time selection
            refSt.typeOfTimeSelection = atoi(argv[index++]);

            //----------------------------------------------------------------------------
            // Filenames (used only if IO_HANDLING==$IO_BASH)
            //----------------------------------------------------------------------------
            refSt.FILE_PCU       = argv[index++];
            refSt.FILE_CONT      = argv[index++];
            refSt.FILE_CONT_TRAJ = argv[index++];
            refSt.FILE_JPL_BIN   = argv[index++];
            refSt.FILE_JPL_TXT   = "cont_jpl_temp.txt";
            break;
        }


        //================================================================================
        // Refinement of the whole trajectory CMU EML2 to CM SEMLi
        //================================================================================
        case COMP_CM_EML2_TO_CM_SEML_REFINE:
        case COMP_EPHEMERIDES_TEST:
        case COMP_CM_EML2_TO_CMS_SEML_READ:
        case COMP_SINGLE_ORBIT:
        case COMP_VF_TEST:
        case COMP_test_INVMAN:
        case COMP_VOFTS_TO_VOFTS:
        {
            // Nothing for now, no additionnal parameters are needed!
            break;
        }

        }
    }


    //====================================================================================
    // Switch
    //====================================================================================
    switch(COMP_TYPE)
    {
        //================================================================================
        // 3D Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML_3D:
    {
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //--------------------------------------------------------------------------------
        tic();
        compute_grid_CMU_EM_3D(HYP_EPSILON_EML2, projSt);
        cout << "End of in compute_grid_CMU_EM_3D in " << toc() << endl;


        //--------------------------------------------------------------------------------
        // Integrate those initial conditions and project them on the CM SEMLi
        //--------------------------------------------------------------------------------
        tic();
        int_proj_CMU_EM_on_CM_SEM_3D(projSt);
        cout << "End of in int_proj_CMU_EM_on_CM_SEM_3D in " << toc() << endl;
        break;
    }

        //================================================================================
        // 3D Projection CMU SEMLi to CM EMLj
        //================================================================================
    case COMP_CMU_SEML_TO_CM_EML:
    {
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU SEMLi on a given grid
        //--------------------------------------------------------------------------------
        tic();
        compute_grid_CMU_SEM_3D(HYP_EPSILON_SEML2, projSt);
        cout << "End of in compute_grid_CMU_SEM_3D in " << toc() << endl;


        //--------------------------------------------------------------------------------
        // Integrate those initial conditions and project them on the CM EMLj
        //--------------------------------------------------------------------------------
        tic();
        int_proj_CMU_SEM_on_CM_EM_3D(projSt);
        cout << "End of in int_proj_CMU_SEM_on_CM_EM_3D in " << toc() << endl;
        break;
    }

        //================================================================================
        // Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML:
    {
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //--------------------------------------------------------------------------------
        tic();
        compute_grid_CMU_EM(HYP_EPSILON_EML2, projSt);
        cout << "End of in compute_grid_CMU_EM in " << toc() << endl;

        //--------------------------------------------------------------------------------
        // Integrate those initial conditions and project them on the CM SEMLi
        //--------------------------------------------------------------------------------
        tic();
        int_proj_CMU_EM_on_CM_SEM(projSt);
        cout << "End of in int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;

        break;
    }

        //================================================================================
        // Projection CMU EML2 to CM SEMLi, on a given set of orbits
        //================================================================================
    case COMP_ORBIT_EML2_TO_CM_SEML:
    {
        tic();
        int_proj_ORBIT_EM_on_CM_SEM(projSt, 1);
        cout << "End of in int_proj_ORBIT_EM_on_CM_SEM in " << toc() << endl;
        break;
    }

        //================================================================================
        // Projection CMU EML2 to CM SEMLi, on a given orbit
        //================================================================================
    case COMP_SINGLE_ORBIT_EML2_TO_CM_SEML:
    {
        tic();
        int_proj_SINGLE_ORBIT_EM_on_CM_SEM(projSt, 1);
        cout << "End of in int_proj_ORBIT_EM_on_CM_SEM in " << toc() << endl;
        break;
    }

        //================================================================================
        // Projection CMU EML2 to CM SEMLi, at a fixed energy
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML_H:
    {
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //--------------------------------------------------------------------------------
        tic();
        compute_grid_CMU_EM_dH(HYP_EPSILON_EML2, projSt);
        cout << "End of in compute_grid_CMU_EM in " << toc() << endl;

        //--------------------------------------------------------------------------------
        // Integrate those initial conditions and project them on the CM SEMLi
        //--------------------------------------------------------------------------------
        tic();
        int_proj_CMU_EM_on_CM_SEM_dH(projSt);
        cout << "End of in int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;

        break;
    }
        //================================================================================
        // Refinement CMU EML2 to CMS SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CMS_SEML:
    {
        //--------------------------------------------------------------------------------
        // Complete routine: new version
        //--------------------------------------------------------------------------------
        switch(refSt.type)
        {
        case REF_ORBIT:
        case REF_CONT_ORBIT:
            sorefemlisemli(refSt);
            break;

        default:
            refemlisemli(refSt);
            break;
        }

        break;

        //--------------------------------------------------------------------------------
        // Celestia format (for movies)
        //--------------------------------------------------------------------------------
        //toCelestiaFormat("jpltraj.xyz");
        break;
    }

        //================================================================================
        // Refinement CMU SEMLi to CMS EMLj
        //================================================================================
    case COMP_CMU_SEML_TO_CMS_EML:
    {
        //--------------------------------------------------------------------------------
        // Complete routine: new version
        //--------------------------------------------------------------------------------

        break;
    }

        //================================================================================
        // Refinement of the whole trajectory to JPL ephemerides
        //================================================================================
    case COMP_REF_JPL:
    {
        reffromcontemlisemli(refSt);

        //jplref3d(refSt.coord_type, refSt);
        //comptojplref3d(refSt.coord_type, refSt);
        //compref3d_test_eml2seml_synjpl(refSt.coord_type);
        break;
    }

        //================================================================================
        // Test on ephemerides
        //================================================================================
    case COMP_EPHEMERIDES_TEST:
    {
        //----------------------
        //load kernels
        //----------------------
        furnsh_c("spice/kernels/metakernel.furnsh");

        //----------------------
        //TEST 1: cocs
        //----------------------
        //double B[3];
        //double C[3][3], k, kCprim[3][3], Bprim[3];
        //init_ecl2synpos("2000 JAN 31 01:00", B, C, &k);
        //init_ecl2synstate("1949 DEC 14 00:00:00.000", B, C, &k, Bprim, kCprim);

        //----------------------
        //TEST 2: Info
        //----------------------
        //displayKernelFeaturesOneBody();

        //----------------------
        //TEST 3: Search for the best fit
        //----------------------
        double et;
        qbcp2jpl_disp(0.2, &et, VSEM);

        //----------------------
        //TEST 4: Compare the numerical constants
        //----------------------
        //comp_num_const();

        //----------------------
        //TEST 5: vector field
        //----------------------
        //test_asteroid();
        //test_coc(VSEM);

        break;
    }
        //================================================================================
        // READ CMU EML2 to CMS SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CMS_SEML_READ:
    {
        //Filename;
        string filename = get_filenameCUM(IO_DEFAULT, SEML.cs->F_PLOT, "", OFTS_ORDER, TYPE_MAN_PROJ, SEML.li_SEM, -1, -1, ios::in);

        //Structure to store data;
        ProjResClass sortSt;

        //Read data file
        sortSt.readProjRes_t0(filename, -1, TIME_SELECTION_ABSOLUTE);


        //Display
        coutmp();
        cout << "--------------------------------------" << endl;
        cout << "The first entry for this time is:" << endl;
        sortSt.displayFirstEntry();

        cout << "--------------------------------------" << endl;
        cout << "The last entry for this time is:" << endl;
        sortSt.displayLastEntry();

        break;
    }

        //================================================================================
        // Just some examples of solutions
        //================================================================================
    case COMP_SINGLE_ORBIT:
    {
        //Reduced number of variables in the default invariant manifolds
        reduced_nv = Invman::compRNV(SEML.cs_em);
        //--------------------------------------------------------------------------------
        // Initial conditions
        //--------------------------------------------------------------------------------
        double st0[reduced_nv];
        switch(MODEL_TYPE)
        {
        case M_QBCP:
        {
            switch(fwrk)
            {
            case F_EM:

                switch(LI_EM)
                {
                case 1:
                    st0[0] = 0.0;//-0.145454545454546;
                    st0[1] = 0.5;//0.194601266014795;
                    st0[2] = 0.0;//-1.309090909090909;
                    st0[3] = 0.0;//0.194601266014795;
                    if(reduced_nv == 5) st0[4] = 0.0;
                    break;
                case 2:
                    st0[0] = 28;//-20 ;
                    st0[1] = 1.74347452709299;//-6.316522019280152e-03;
                    st0[2] = 36;//-36;
                    st0[3] = 1.74347452709299;//-1.681003648125090e-03;
                    if(reduced_nv == 5) st0[4] = 0.0;

                    break;
                }
                break;

            case F_SEM:
                //Values for Earth capture via Moon slingshot
                st0[0] =  0.1;
                st0[1] =  0.2;
                st0[2] = -0.1;
                st0[3] =  0.0;
                if(reduced_nv == 5) st0[4] =  0.0;
                break;
            }
        }
        break;
        case M_RTBP:
        {
            st0[0] =  0.17;
            st0[1] =  0.17;
            st0[2] =  0.17;
            st0[3] =  0.17;
            if(reduced_nv == 5) st0[4] =  0.0;
            break;
        }
        }

        //--------------------------------------------------------------------------------
        // Computation
        //--------------------------------------------------------------------------------
        int isFlagOn = 1;
        int isPlot   = 1;
        double t0 = 0.96*SEML.us->T;
        double  N = 50;

        //gridOrbit_si(st0, t0, t0 + N*SEML.us->T, 1e-2*SEML.us->T, isFlagOn, isPlot);
        gridOrbit_strob(st0, t0, N, isFlagOn, isPlot);

        //--------------------------------------------------------------------------------
        // Mean orbit of the Moon, in NCSEM coordinates
        //--------------------------------------------------------------------------------
        N = 1000;
        double tk = 0;
        double pm[3];
        double qpm, mean_qpm = 0.0;
        for(int k = 0; k < N; k++)
        {
            tk = (1.0*k)/N*SEML.us->T;
            evaluateCoef(pm, tk, SEML.us->n, SEML.nf, SEML.cs->pm, 3);
            qpm = sqrt((-1.0 - pm[0]) * (-1.0 - pm[0]) + (0.0 - pm[1]) * (0.0 - pm[1]) + (0.0 - pm[2]) * (0.0 - pm[2]));
            mean_qpm += qpm/N;
        }
        cout << "mean_qpm = " << mean_qpm << endl;

        break;
    }
        //================================================================================
        // Test of the vector fields. Better in OOFTDA??
        // @todo set this routine (qbtbp_test) in OOFTDA
        //================================================================================
    case COMP_VF_TEST:
    {
        qbtbp_test(SEML.us->T, SEML);
        break;
    }

        //================================================================================
        // Test of the new invariant manifold representation
        //================================================================================
    case COMP_test_INVMAN:
    {
        test_evalCCMtoTFC();
        test_evalRCMtoNC();
        test_evalDCCMtoTFC();
        test_evalDRCMtoTFC();
    }
    break;

        //================================================================================
        // Store CS/CU into one-dimensionnal series to gain memory.
        // DEPRECATED: NOW DONE IN OOFTDA
        //================================================================================
    case COMP_VOFTS_TO_VOFTS:
    {
        cout << "-------------------------------------------------------" << endl;
        cout << "Store CS/CU into one-dimensionnal series to gain memory" << endl;
        cout << "-------------------------------------------------------" << endl;

        //--------------------------------------------------------------------------------
        //Storing
        //--------------------------------------------------------------------------------
        Oftsc W  = Oftsc(reduced_nv, OFTS_ORDER, OFS_NV, OFS_ORDER);
        Oftsc W1 = Oftsc(1,          OFTS_ORDER, OFS_NV, OFS_ORDER);
        fromVOFTStoVOFTS_bin(W, W1, SEML.cs->F_PMS+"W/Wh", SEML.cs->F_PMS+"W/Wh1");

        //--------------------------------------------------------------------------------
        //Test, only the 3 first orders
        //--------------------------------------------------------------------------------
        vector<Oftsc> W1v;
        W1v.reserve(2);
        for(int i = 0; i < 2; i++) W1v.push_back(Oftsc(1, min(3, OFTS_ORDER), OFS_NV, OFS_ORDER));
        readVOFTS_bin(W1v, SEML.cs->F_PMS+"W/Wh1");

        cout << "Test:" << endl;
        for(int i = 0; i < 2; i++)
        {
            cout << "Wh1[" << i << "] = " << endl;
            cout << W1v[i] << endl;
            cout << "-------------------------------------------------------" << endl;
        }
        break;
    }
    }

    return FTC_SUCCESS;
}
