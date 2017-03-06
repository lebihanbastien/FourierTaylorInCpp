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
#include "lencon.h"
#include "oolencon.h"
#include "ephemerides.h"
#include "Invman.h"
#include "Orbit.h"
#include "Config.h"
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
#define COMP_CM_EML2_TO_CM_SEML        0    //EML2 Center Manifold to SEMLi Center Manifold
#define COMP_CM_EML2_TO_CMS_SEML       1    //EML2 Center Manifold to SEMLi Center-Stable Manifold
#define COMP_SINGLE_ORBIT              2    //just some example of solutions
#define COMP_CM_EML2_TO_CMS_SEML_READ  3    //Read
#define COMP_VF_TEST                   4    //Test of the QBCP vector field. Should probably be put in OOFTDA
#define COMP_CM_EML2_TO_CM_SEML_REFINE 5    //EML2 Center Manifold to SEMLi Center Manifold
#define COMP_EPHEMERIDES_TEST          6    //Ephemerides test
#define COMP_CM_EML2_TO_CM_SEML_3D     7    //EML2 Center Manifold to SEMLi Center Manifold in 3D
#define COMP_VOFTS_TO_VOFTS            8    //Store CS/CU into one-dimensionnal series to gain memory
#define COMP_test_INVMAN               9    //Test of the new invariant manifold implementation
#define COMP_REF_JPL                   10   //Refine to JPL ephemerides


/************* NOTES *********************************************************************
 Notes from Reunion with Josep (15/12/2015)

 Step 3: Fix a given Pk section and fix the phase @ this section
+ what about refinment in the JPL eph?
+ Book I for an example of continuation.

+ @TODO: see if we can fasten the JPL and QBCP vf codes! Started to be sloooow. Can we also fasten ode78?

+ @TODO: in ooconteml2seml, try to move around the et0. Is there a connection with the size of the final refined EML2 orbit?
+ @TODO: in ooconteml2seml, what happens if we search for et0 in more than 3 month span? we can probably find a month for which the size is roughly equivalent

+ @TODO: do NOT use function that select pointers to function (does not work in C/C++)
but pass the desired function as argument. Should replace ftc_select_vf that does
NOT work!!! Ou alors dans le return?? Yes, by creating a typedef

+ @TODO: look for the other kinds of switch... typically in oosrefeml2seml (call to mft3d...)


*****************************************************************************************/
int main(int argc, char** argv)
{
    //====================================================================================
    //Declare configuration parameters
    //====================================================================================
    // General parameters (orders, etc)
    int COMP_TYPE, NUM_THREADS, MODEL_TYPE, LI_EM, LI_SEM, ISPAR;
    // Projection parameters (in structure)
    ProjSt projSt;
    // Refinement parameters (in structure)
    RefSt refSt;

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
        COMP_TYPE   = COMP_CM_EML2_TO_CMS_SEML;

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
    }

    //====================================================================================
    // General settings
    //====================================================================================
    cout << setiosflags(ios::scientific) << setprecision(15);
    omp_set_num_threads(NUM_THREADS);


    //====================================================================================
    // Some other parameters depend on the type of computation
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
        model     = MODEL_TYPE;
        fwrk      = F_EM;
        isNorm    = 1;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER;
        break;

        //================================================================================
        // Projection CMU EML2 to CMS SEMLi
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
    //Check if arguments have been passed
    //------------------------------------------------------------------------------------
    if(argc == 1) //no argument were passed
    {
        //================================================================================
        // Projection parameters
        //================================================================================
        //--------------------------------------------------------------------------------
        //Time grid: min, max and number of points on the grid
        //--------------------------------------------------------------------------------
        projSt.TMIN  = 0.88*SEML.us.T;
        projSt.TMAX  = 1.00*SEML.us.T;

        projSt.TSIZE    = 0;
        projSt.TLIM[0]  = projSt.TMIN;
        projSt.TLIM[1]  = projSt.TMAX;

        //--------------------------------------------------------------------------------
        // Configuration (s1, s2, s3, s4) grid
        //--------------------------------------------------------------------------------
        projSt.GLIM_SI[0][0] = -30;
        projSt.GLIM_SI[0][1] = +30;

        projSt.GLIM_SI[1][0] = +0;
        projSt.GLIM_SI[1][1] = +10;

        projSt.GLIM_SI[2][0] = -30;
        projSt.GLIM_SI[2][1] = +30;

        projSt.GLIM_SI[3][0] = +0;
        projSt.GLIM_SI[3][1] = +10;

        projSt.GSIZE_SI[0]   = 50;
        projSt.GSIZE_SI[1]   = 5;
        projSt.GSIZE_SI[2]   = 50;
        projSt.GSIZE_SI[3]   = 5;

        //--------------------------------------------------------------------------------
        //Stable parameters (are not supposed to change)
        //--------------------------------------------------------------------------------
        projSt.TM    = 12.0*SEML.us.T; // Maximum integration time
        projSt.MSIZE = 500;            // Number of points on each trajectory
        projSt.NSMIN = 20;             // Number of sorted solutions
        projSt.YNMAX = 0.6;            // The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
        projSt.SNMAX = 0.6;            // The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
        projSt.NOD   = 6;              // Number of dimensions on which we compute the norm of the projection
        projSt.ISPAR = ISPAR;          // Boolean for parallel computation

        //================================================================================
        // Refinement parameters
        //================================================================================

        //--------------------------------------------------------------------------------
        // Parameters that change often
        //--------------------------------------------------------------------------------
        //rk: set REF_CONT_D_HARD_CASE for difficult cases
        //with REF_CONT_D (ex: EML2-SEMLi via SEML1...)
        refSt.type          = REF_CONT_D;         // Type of refinement
        refSt.dim           = REF_MIXED;         // Type of dimensions planar or 3d?
        refSt.t0_des        = 0.995*SEML.us_em.T;  // Initial time

        // Direction of the continuation procedure
        refSt.isDirUD       = 0;                  // is it user defined?
        refSt.Dir           = +1;                 // if not, +1 or -1

        // Domain of search for the first guess
        refSt.s1_CMU_EM_MIN = -35;
        refSt.s1_CMU_EM_MAX = +35;

        refSt.s2_CMU_EM_MIN = +0;
        refSt.s2_CMU_EM_MAX = +0;

        refSt.s3_CMU_EM_MIN = -35;
        refSt.s3_CMU_EM_MAX = +35;

        refSt.s4_CMU_EM_MIN = +4;
        refSt.s4_CMU_EM_MAX = +4;


        // Or, if we want the user to define such domain:
        refSt.isLimUD       =  0;

        //Limits for the time of flight during transfers - not used if -1
        refSt.tof_MIN       = -1;
        refSt.tof_MAX       = -1;

        // Number of steps in the continuation procedure
        refSt.cont_step_max    = +450;            // with fixed times
        refSt.cont_step_max_vt = +150;            // with variable times

        // Initial step in the continuation procedure
        refSt.ds0    = 8e-2;                      //with fixed time
        refSt.ds0_vt = (LI_EM ==1)?  3e-2:1e-2;   //with variable time

        // Desired number of iterations in Newton's method in the continuation procedure
        refSt.nu0 = 2;          //with fixed time
        refSt.nu0_vt = 3;       //with variable time

        //User parameters
        refSt.isFlagOn      = 1;                  // do we have steps in the procedure - asking the user to press enter to go on?
        refSt.isPlotted     = 1;                  // do we plot the results during the computation?
        refSt.isSaved       = 1;                  // do we save the results in data files?
        refSt.isFromServer  = 1;                  // does the raw data comes from server files?

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
        refSt.dsmin         = 1e-6;                     //with fixed time
        refSt.dsmin_vt      = 1e-6;                     //with variable time
        refSt.dsmax         = 10;                 //with fixed time
        refSt.dsmax_vt      = 10;                 //with variable time

        refSt.xps           = (LI_SEM == 1)? +0.6:-0.6; // position of the poincaré section in NCSEM coordinates
        refSt.isJPL         = 1;                        // is the JPL refinement performed when possible?
        refSt.djplcoord     = -1;                       // coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with NJ2000
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
        {
            //----------------------------------------------------------------------------
            //Time grid: min, max and number of points on the grid
            //----------------------------------------------------------------------------
            projSt.TMIN  = atof(argv[index++])*SEML.us.T;
            projSt.TMAX  = atof(argv[index++])*SEML.us.T;
            projSt.TM    = atof(argv[index++])*SEML.us.T;

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
            //Stable parameters (are not supposed to change)
            //----------------------------------------------------------------------------
            projSt.MSIZE = atoi(argv[index++]);
            projSt.NSMIN = atoi(argv[index++]);
            projSt.YNMAX = atof(argv[index++]);
            projSt.SNMAX = atof(argv[index++]);
            projSt.NOD   = atoi(argv[index++]);
            projSt.ISPAR = ISPAR;

            break;
        }


        //================================================================================
        // Refinement CMU EML2 to CMS SEMLi
        //================================================================================
        case COMP_CM_EML2_TO_CMS_SEML:
        case COMP_REF_JPL:
        {
            //------------------------------------------------------------------------
            // Parameters that change often
            //------------------------------------------------------------------------
            //rk: set REF_CONT_D_HARD_CASE for difficult cases
            //with REF_CONT_D (ex: EML2-SEMLi via SEML1...)
            refSt.type          = atoi(argv[index++]);            // Type of refinement
            refSt.dim           = atoi(argv[index++]);            // Type of dimensions planar or 3d?
            refSt.t0_des        = atof(argv[index++])*SEML.us.T;  // Initial time

            // Direction of the continuation procedure
            refSt.isDirUD       = atoi(argv[index++]);    // is it user defined?
            refSt.Dir           = atoi(argv[index++]);    // if not, +1 or -1

            // Domain of search for the first guess
            refSt.s1_CMU_EM_MIN = atof(argv[index++]);
            refSt.s1_CMU_EM_MAX = atof(argv[index++]);

            refSt.s2_CMU_EM_MIN = atof(argv[index++]);
            refSt.s2_CMU_EM_MAX = atof(argv[index++]);

            refSt.s3_CMU_EM_MIN = atof(argv[index++]);
            refSt.s3_CMU_EM_MAX = atof(argv[index++]);

            refSt.s4_CMU_EM_MIN = atof(argv[index++]);
            refSt.s4_CMU_EM_MAX = atof(argv[index++]);


            // Or, if we want the user to define such domain:
            refSt.isLimUD       = atoi(argv[index++]);

            //Limits for the time of flight during transfers - not used if negative
            refSt.tof_MIN       = atof(argv[index++])*SEML.us.T;
            refSt.tof_MAX       = atof(argv[index++])*SEML.us.T;

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

            //Maximum angle around SEMLi if REF_COND_T is used (in degrees)
            refSt.thetaMax      = atof(argv[index++]);     //should be a multiple of 90°

            //------------------------------------------------------------------------
            // Parameters that are stable
            //------------------------------------------------------------------------
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
            refSt.dsmax         = 2e-1;                 //with fixed time
            refSt.dsmax_vt      = 2e-1;                 //with variable time

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
        // New version
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid
        tic();
        oo_compute_grid_CMU_EM_3D(PROJ_EPSILON, projSt);
        cout << "End of in oo_compute_grid_CMU_EM_3D in " << toc() << endl;



        // Integrate those initial conditions and project them on the CM SEMLi
        tic();
        oo_int_proj_CMU_EM_on_CM_SEM_3D(projSt);
        cout << "End of in oo_int_proj_CMU_EM_on_CM_SEM_3D in " << toc() << endl;
        break;

        //--------------------------------------------------------------------------------
        // Old version
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid

        tic();
        //compute_grid_CMU_EM_3D(PROJ_EPSILON, TLIM, TSIZE, GLIM_SI, GSIZE_SI, CM_TFC, Mcoc, Vcoc, ISPAR);
        cout << "End of in compute_grid_CMU_EM_3D in " << toc() << endl;



        // Integrate those initial conditions and project them on the CM SEMLi
        tic();
        //int_proj_CMU_EM_on_CM_SEM_3D(TM, MSIZE , NOD, ISPAR, YNMAX, SNMAX);
        cout << "End of in int_proj_CMU_EM_on_CM_SEM_3D in " << toc() << endl;

        break;
    }

        //================================================================================
        // Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML:
    {
        //--------------------------------------------------------------------------------
        // New version
        //--------------------------------------------------------------------------------
        // Compute initial conditions in CMU EML2 on a given grid
        tic();
        oo_compute_grid_CMU_EM(PROJ_EPSILON, projSt);
        cout << "End of in oo_compute_grid_CMU_EM in " << toc() << endl;


        // Integrate those initial conditions and project them on the CM SEMLi
        tic();
        oo_int_proj_CMU_EM_on_CM_SEM(projSt);
        cout << "End of in oo_int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;
        break;

        //--------------------------------------------------------------------------------
        // Old version
        //--------------------------------------------------------------------------------
        tic();
        //compute_grid_CMU_EM(PROJ_EPSILON, TMIN, TMAX, TSIZE,
        //                    GMIN_S1, GMAX_S1, GMIN_S3, GMAX_S3, GSIZE_S1, GSIZE_S3,
        //                    CM_TFC, Mcoc, Vcoc, ISPAR);
        cout << "End of in compute_grid_CMU_EM in " << toc() << endl;




        // Integrate those initial conditions and project them on the CM SEMLi

        tic();
        //int_proj_CMU_EM_on_CM_SEM(TM, TSIZE, GSIZE_S1, GSIZE_S3, MSIZE, NSMIN, NOD,
        //                          ISPAR, YNMAX, SNMAX);
        cout << "End of in int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;



        //----------------------
        // Compute the best solutions from previous computation
        // No equivalent in new implementation for now
        //----------------------
        //int_sorted_sol_CMU_EM_to_CM_SEM(OFTS_ORDER, MSIZE, ISPAR, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);
        //cout << "End of int_sorted_sol_CMU_EM_to_CM_SEM in " << toc() << endl;
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
        oorefeml2seml(refSt);
        //ooconteml2seml(refSt);
        break;

        //--------------------------------------------------------------------------------
        // Complete routine: old version
        //--------------------------------------------------------------------------------
        //refeml2seml(20, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, MIcoc, Vcoc, refSt);
        break;

        //--------------------------------------------------------------------------------
        // Complete routine: original version
        //--------------------------------------------------------------------------------
        //int ISPLANAR = 1;
        //int ISSAVED = 1;
        //int ISUSERDEFINED = 1;
        // Continuation: CMU to CMS, with variable times
        //ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(20, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, MIcoc, Vcoc, ISPLANAR, ISSAVED, ISUSERDEFINED);
        break;

        // Multiple shooting: CMU to CMS, with variable times. DEPRECATED?
        //ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE(20, 60, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, ISPLANAR, ISSAVED, ISUSERDEFINED);
        //break;

        //--------------------------------------------------------------------------------
        // Celestia format (for movies)
        //--------------------------------------------------------------------------------
        //toCelestiaFormat("jpltraj.xyz");
        break;
    }

        //================================================================================
        // Refinement of the whole trajectory to JPL ephemerides
        //================================================================================
    case COMP_REF_JPL:
    {
        //oojplrefft3d(refSt.coord_type, refSt);
        oointojplrefft3d(refSt.coord_type, refSt);
        //oocomprefft3d_test_eml2seml_synjpl(refSt.coord_type);
        break;
    }
        //================================================================================
        // Refinement of the whole trajectory CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML_REFINE:
    {

        //----------------------
        // Multiple shooting on the whole trajectory (base for JPL refinment)
        //ref_CMU_EM_to_CM_SEM_MSD_COMP(OFTS_ORDER, 60, NCSEM, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);

        //----------------------
        // Multiple shooting on the manifold leg only
        //ref_CMU_EM_to_CM_SEM_MSD(OFTS_ORDER, 10, 1, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to free boundary (SEMLi orbit included)
        //ref_CMU_EM_to_CM_SEM_MSD_PART(OFTS_ORDER, 60, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to CM (DEPRECATED)
        //ref_CMU_EM_to_CM_SEM_MSD_DEP(OFTS_ORDER, 50, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, false);

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
        string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ, SEML.li_SEM);

        //Indix
        int ind = 0;
        vector<size_t> sortId;

        //Vector to store the values of t0v_EM;
        vector<double> t0_CMU_EM;
        vector<double> tf_man_EM;
        vector<double> s1_CMU_EM;
        vector<double> s2_CMU_EM;
        vector<double> s3_CMU_EM;
        vector<double> s4_CMU_EM;
        vector<double> s5_CMU_EM;
        vector<double> pmin_dist_SEM;
        vector<double> s1_CM_SEM;
        vector<double> s2_CM_SEM;
        vector<double> s3_CM_SEM;
        vector<double> s4_CM_SEM;


        //Read data file
        readIntProjCU_bin(filename, t0_CMU_EM, tf_man_EM, s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM, s5_CMU_EM,
                          pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM, sortId);

        //Display
        coutmp();
        ind = sortId[0];
        cout << "--------------------------------------" << endl;
        cout << "The first entry for this time is:" << endl;
        cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
        cout << "tf_EM    = " << tf_man_EM[ind] << endl;
        cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
        cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
        cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;

        ind = sortId[t0_CMU_EM.size()-1];
        cout << "--------------------------------------" << endl;
        cout << "The last entry for this time is:" << endl;
        cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
        cout << "tf_EM    = " << tf_man_EM[ind] << endl;
        cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
        cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
        cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;
        break;
    }

        //================================================================================
        // Just some examples of solutions
        //================================================================================
    case COMP_SINGLE_ORBIT:
    {
        //Reduced number of variables in the default invariant manifolds
        reduced_nv = Invman::compRNV(SEML.cs);
        //----------------------------------------------------------------------------
        // Initial conditions
        //----------------------------------------------------------------------------
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
                    st0[0] =  5.901025202644599e+00;//-20 ;
                    st0[1] = -6.316522019280152e-03;
                    st0[2] = -3.240398486261306e+01;//-36;
                    st0[3] = -1.681003648125090e-03;
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

        //----------------------------------------------------------------------------
        // Computation
        //----------------------------------------------------------------------------
        //oo_gridOrbit(st0, +6.016997770510415e+00, -2.793897158902917e+01, 1e-2);
        oo_gridOrbit(st0, 0.0, 5*SEML.us.T, 1e-3*SEML.us.T);
        //gridOrbit(st0, +6.016997770510415e+00, -2.793897158902917e+01, 1e-2, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);

        break;
    }
        //================================================================================
        // Test of the vector fields. Better in OOFTDA??
        // @todo set this routine (qbtbp_test) in OOFTDA
        //================================================================================
    case COMP_VF_TEST:
    {
        qbtbp_test(SEML.us.T, SEML);
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
        fromVOFTStoVOFTS_bin(W, W1, SEML.cs.F_PMS+"W/Wh", SEML.cs.F_PMS+"W/Wh1");

        //--------------------------------------------------------------------------------
        //Test, only the 3 first orders
        //--------------------------------------------------------------------------------
        vector<Oftsc> W1v;
        W1v.reserve(2);
        for(int i = 0; i < 2; i++) W1v.push_back(Oftsc(1, min(3, OFTS_ORDER), OFS_NV, OFS_ORDER));
        readVOFTS_bin(W1v, SEML.cs.F_PMS+"W/Wh1");

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
