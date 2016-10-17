#include <omp.h>
#include <iostream>

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
extern "C" {
    #include "gnuplot_i.h"
    #include "SpiceUsr.h"
}

// TYPE OF COMPUTATIONS
#define COMP_CM_EML2_TO_CM_SEML        0    //EML2 Center Manifold to SEML2 Center Manifold
#define COMP_CM_EML2_TO_CMS_SEML       1    //EML2 Center Manifold to SEML2 Center-Stable Manifold
#define COMP_SINGLE_ORBIT              2    //just some example of solutions
#define COMP_CM_EML2_TO_CMS_SEML_READ  3    //Read
#define COMP_VF_TEST                   4    //Test of the QBCP vector field. Should probably be put in OOFTDA
#define COMP_CM_EML2_TO_CM_SEML_REFINE 5    //EML2 Center Manifold to SEML2 Center Manifold
#define COMP_EPHEMERIDES_TEST          6    //Ephemerides test
#define COMP_CM_EML2_TO_CM_SEML_3D     7    //EML2 Center Manifold to SEML2 Center Manifold in 3D
#define COMP_VOFTS_TO_VOFTS            8    //Store CS/CU into one-dimensionnal series to gain memory
#define COMP_test_INVMAN               9    //Test of the new invariant manifold implementation


using namespace std;

/************* NOTES ********************
 Notes from Reunion with Josep (15/12/2015)

 Step 1: Find some solution "by hand"
 Step 2: Refine to get a final heteroclinic connection
         Look @ collocation methods (Doedel) or Simo methods
 Step 3: Fix a given Pk section and fix the phase @ this section

    + what about refinment in the JPL eph?
    + Book I for an example of continuation.

Look at the domain of validity of the (un)stable manifolds of EML1,2
Can we use it at its limit?
Is it close to the Moon?
****************************************/


int main()
{
    //====================================================================================
    // Type of computation
    //====================================================================================
    int COMPTYPE = COMP_CM_EML2_TO_CMS_SEML;//COMP_CM_EML2_TO_CM_SEML_3D;//COMP_CM_EML2_TO_CMS_SEML;

    //====================================================================================
    // cout settings
    //====================================================================================
    cout << setiosflags(ios::scientific) << setprecision(15);

    //====================================================================================
    // openMP settings
    //====================================================================================
    int num_threads = 4;
    omp_set_num_threads(num_threads);
    cout << "num_threads  = "  << num_threads << endl;

    //====================================================================================
    // Define the global parameters
    //====================================================================================
    MODEL_TYPE = M_QBCP;
    OFTS_ORDER = 16;
    if(MODEL_TYPE == M_RTBP) OFS_ORDER  = 0;
    else OFS_ORDER  = 30;

    //====================================================================================
    // Configuration parameters
    //====================================================================================
    int model = 0, coordsys = 0, isNorm = 0, li_EM = 0, li_SEM = 0;
    int pms_EM = 0, pms_SEM = 0, mType_EM = 0, mType_SEM = 0, reduced_nv = 0;
    switch(COMPTYPE)
    {
        //================================================================================
        // Test and additionnal computations
        //================================================================================
     case COMP_EPHEMERIDES_TEST:
     case COMP_VOFTS_TO_VOFTS:
     case COMP_test_INVMAN:
        model     = MODEL_TYPE;
        coordsys  = F_SEM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
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
        coordsys  = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
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
        model     = MODEL_TYPE;
        coordsys  = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
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
        coordsys  = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
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
        coordsys  = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER;
        break;

    default:
        cout << "main. Warning: unknown COMPTYPE!" << endl;
    }

    //====================================================================================
    //Define the number of reduced variables
    //depending on the type of manifold
    //====================================================================================
    if(coordsys == F_EM)
    {
        switch(mType_EM)
        {
        case MAN_CENTER:
            reduced_nv = 4;
            break;
        case MAN_CENTER_S:
        case MAN_CENTER_U:
            reduced_nv = 5;
            break;
        case MAN_CENTER_US:
            reduced_nv = 6;
            break;
        default:
            cout << "Warning: unknown type of manifold. reduced_nv = 4 by default." << endl;
            reduced_nv = 4;
            break;
        }
    }
    else
    {
        switch(mType_SEM)
        {
        case MAN_CENTER:
            reduced_nv = 4;
            break;
        case MAN_CENTER_S:
        case MAN_CENTER_U:
            reduced_nv = 5;
            break;
        case MAN_CENTER_US:
            reduced_nv = 6;
            break;
        default:
            cout << "Warning: unknown type of manifold. REDUCED_NV = 4 by default." << endl;
            reduced_nv = 4;
            break;
        }
    }

    //====================================================================================
    // Initialization of the environnement
    // Mandatory to perform any computation
    //====================================================================================
    initialize_environment(li_EM, li_SEM, isNorm, model, coordsys, pms_EM, pms_SEM, mType_EM, mType_SEM);

    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc> CM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_TFC;     //center manifold in TFC coordinates
    //Init routine
    initOFTS(CM_NC, CM_TFC, reduced_nv, OFTS_ORDER, OFS_ORDER, SEML.cs);

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc(6);      //COC vector

    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //--------------------------------------
    //DCM_TFC : Jacobian of CM_TFC
    //--------------------------------------
    matrix<Oftsc> DCM_TFC(6, reduced_nv, reduced_nv, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_TFC, SEML.cs.F_PMS+"DWf/DWhc");

    //====================================================================================
    // Second init, with new implementation!!
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs);

    //====================================================================================
    // Initial conditions
    //====================================================================================
    double st0[reduced_nv];
    switch(coordsys)
    {
    case F_EM:

        switch(li_EM)
        {
        case 1:
            st0[0] = -0.145454545454546;
            st0[1] =  0.194601266014795;
            st0[2] = -1.309090909090909;
            st0[3] =  0.194601266014795;
            if(reduced_nv == 5) st0[4] = 0.0;
            break;
        case 2:
            st0[0] = 5.901025202644599e+00;//-20 ;
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

    //====================================================================================
    // Computation parameters: Projection CMU EML2 to CM SEMLi
    //====================================================================================
    //TIME
    double TM      = 12*SEML.us.T;
    double TMIN    = 0.99*SEML.us.T;
    double TMAX    = 1.0*SEML.us.T;
    int    TSIZE   = 0;
    double TLIM[2] = {TMIN, TMAX};

    //UNSTABLE MANIFOLD AT EML2
    double GMIN_S1  = -30;
    double GMAX_S1  = +30;
    int    GSIZE_S1 = 80;

    double GMIN_S3  = -30;
    double GMAX_S3  = +30;
    int    GSIZE_S3 = 80;

    //In the form of a matrix & a vector
    double  GMIN_SI[4][2] = {{-20, 20} , {10, 10} , {-20, 20} , {10, 10} };
    int     GSIZE_SI[4]   = {30, 0, 30, 0};


    //TRAJECTORTY REFINMENT
    int MSIZE = 500;    //number of points on each trajectory
    int NSMIN = 20;     //number of sorted solutions

    //LIMITS OF DETECTION
    double YNMAX = 0.6;  //The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
    double SNMAX = 0.6;  //The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi

    //PARALLEL CPU
    int ISPAR = 1;

    //PLANAR CASE?
    int ISPLANAR = 1;

    //Store?
    int ISSAVED = 1;

    //User-defined direction?
    int ISUSERDEFINED = 1;

    //Number of dimensions on which we compute the norm of the projection error
    int NOD = 6;

    //Refinement structure
    RefSt refst;
    refst.type          = REF_CONT_D;
    refst.dim           = REF_PLANAR;
    refst.time          = REF_VAR_TIME;
    refst.grid          = REF_FIXED_GRID;
    refst.isDirUD       = 1;
    refst.isLimUD       = 1;
    refst.isPlotted     = 1;
    refst.isSaved       = 1;
    refst.isJPL         = 1;
    refst.isDebug       = 0;
    refst.isFromServer  = 0;
    refst.s1_CMU_EM_MIN = +24;
    refst.s1_CMU_EM_MAX = +30;
    refst.s3_CMU_EM_MIN = -18;
    refst.s3_CMU_EM_MAX = +6;
    refst.cont_step_max = +150;
    refst.tspan_EM      = +5*SEML.us_em.T;
    refst.tspan_SEM     = +5*SEML.us_sem.T;

    //====================================================================================
    // For COMP_CM_EML2_TO_CMS_SEML_READ
    //====================================================================================
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

    //Filename;
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ);

    //Indix
    int ind = 0;
    vector<size_t> sortId;

    //====================================================================================
    // Switch
    //====================================================================================
    switch(COMPTYPE)
    {
        //================================================================================
        // 3D Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML_3D:
        //----------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //----------------------
//                tic();
//                compute_grid_CMU_EM_3D(PROJ_EPSILON, TLIM, TSIZE, GMIN_SI, GSIZE_SI, CM_TFC, Mcoc, Vcoc, ISPAR);
//                cout << "End of in compute_grid_CMU_EM_3D in " << toc() << endl;
//        tic();
        oo_compute_grid_CMU_EM_3D(PROJ_EPSILON, TLIM, TSIZE, GMIN_SI, GSIZE_SI, invman, ISPAR);
        cout << "End of in oo_compute_grid_CMU_EM_3D in " << toc() << endl;


        //----------------------
        // Integrate those initial conditions and project them on the CM SEML2
        //----------------------
//                tic();
//                int_proj_CMU_EM_on_CM_SEM_3D(TM, MSIZE , NOD, ISPAR, YNMAX, SNMAX);
//                cout << "End of in int_proj_CMU_EM_on_CM_SEM_3D in " << toc() << endl;
//        tic();
        oo_int_proj_CMU_EM_on_CM_SEM_3D(TM, MSIZE, NOD, ISPAR, YNMAX, SNMAX);
        cout << "End of in oo_int_proj_CMU_EM_on_CM_SEM_3D in " << toc() << endl;
        break;

        //================================================================================
        // Projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML:
        tic();
        //----------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //----------------------
        //        tic();
        //        compute_grid_CMU_EM(PROJ_EPSILON, TMIN, TMAX, TSIZE,
        //                            GMIN_S1, GMAX_S1, GMIN_S3, GMAX_S3, GSIZE_S1, GSIZE_S3,
        //                            CM_TFC, Mcoc, Vcoc, ISPAR);
        //        cout << "End of in compute_grid_CMU_EM in " << toc() << endl;

        tic();
        oo_compute_grid_CMU_EM(PROJ_EPSILON, TMIN, TMAX, TSIZE,
                            GMIN_S1, GMAX_S1, GMIN_S3, GMAX_S3, GSIZE_S1, GSIZE_S3,
                            invman, ISPAR);
        cout << "End of in oo_compute_grid_CMU_EM in " << toc() << endl;

        //----------------------
        // Integrate those initial conditions and project them on the CM SEML2
        //----------------------
        //        tic();
        //        int_proj_CMU_EM_on_CM_SEM(TM, TSIZE, GSIZE_S1, GSIZE_S3, MSIZE, NSMIN, NOD,
        //                                  ISPAR, YNMAX, SNMAX);
        //        cout << "End of in int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;

        tic();
        oo_int_proj_CMU_EM_on_CM_SEM(TM, TSIZE, GSIZE_S1, GSIZE_S3, MSIZE, NSMIN, NOD,
                                  ISPAR, YNMAX, SNMAX);
        cout << "End of in oo_int_proj_CMU_EM_on_CM_SEM in " << toc() << endl;

        //----------------------
        // Compute the best solutions from previous computation
        // No equivalent in new implementation for now
        //----------------------
        //int_sorted_sol_CMU_EM_to_CM_SEM(OFTS_ORDER, MSIZE, ISPAR, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);
        //cout << "End of int_sorted_sol_CMU_EM_to_CM_SEM in " << toc() << endl;
        break;

        //================================================================================
        // Projection CMU EML2 to CMS SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CMS_SEML:
        //----------------------
        // Test: complete routine
        //----------------------
        //refeml2seml(20, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, MIcoc, Vcoc, refst);
        //toCelestiaFormat("jpltraj.xyz");
        refeml2seml(20, NCSEM, invman, refst);
        break;

        //----------------------
        // Multiple shooting: CMU to CMS, with variable times
        //----------------------
        //ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE(20, 60, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, ISPLANAR, ISSAVED, ISUSERDEFINED);
        //break;

        //----------------------
        // Continuation: CMU to CMS, with variable times
        //----------------------
        ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(20, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, MIcoc, Vcoc, ISPLANAR, ISSAVED, ISUSERDEFINED);
        break;

        //================================================================================
        // Refinement of the projection CMU EML2 to CM SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CM_SEML_REFINE:

        //----------------------
        // Multiple shooting on the whole trajectory (base for JPL refinment)
        ref_CMU_EM_to_CM_SEM_MSD_COMP(OFTS_ORDER, 60, NCSEM, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);

        //----------------------
        // Multiple shooting on the manifold leg only
        //ref_CMU_EM_to_CM_SEM_MSD(OFTS_ORDER, 10, 1, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to free boundary (SEML2 orbit included)
        //ref_CMU_EM_to_CM_SEM_MSD_PART(OFTS_ORDER, 60, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to CM (DEPRECATED)
        //ref_CMU_EM_to_CM_SEM_MSD_DEP(OFTS_ORDER, 50, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, false);

        break;

        //================================================================================
        // Test on ephemerides
        //================================================================================
    case COMP_EPHEMERIDES_TEST:

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

        //================================================================================
        // READ CMU EML2 to CMS SEMLi
        //================================================================================
    case COMP_CM_EML2_TO_CMS_SEML_READ:
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

        //================================================================================
        // Just some examples of solutions
        //================================================================================
    case COMP_SINGLE_ORBIT:
        //gridOrbit(st0, +6.016997770510415e+00, -2.793897158902917e+01, 1e-2, CM_NC, CM_TFC, Mcoc, MIcoc, Vcoc);
        oo_gridOrbit(st0, +6.016997770510415e+00, -2.793897158902917e+01, 1e-2, invman);
        break;

        //================================================================================
        // Test of the vector fields. Better in OOFTDA?? €€TODO
        //================================================================================
    case COMP_VF_TEST:
        qbtbp_test(0.25*SEML.us.T, SEML);
        break;

        //================================================================================
        // Test of the new invariant manifold representation
        //================================================================================
    case COMP_test_INVMAN:

        //Test in
        //test_evalCCMtoTFC();
        //test_evalRCMtoNC();
        //test_evalDCCMtoTFC();
        test_evalDRCMtoTFC();

        break;

        //================================================================================
        // Store CS/CU into one-dimensionnal series to gain memory.
        // DEPRECATED: NOW DONE IN OOFTDA
        //================================================================================
    case COMP_VOFTS_TO_VOFTS:

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


    return 0;
}
