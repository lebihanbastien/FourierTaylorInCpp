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

extern "C" {
#include "gnuplot_i.h"
}

// TYPE OF COMPUTATIONS
#define COMP_CM_EML2_TO_CM_SEML  0
#define COMP_CM_EML2_TO_CMS_SEML 1
#define COMP_SINGLE_ORBIT        2
#define COMP_VF_TEST             3
#define COMP_CM_EML2_TO_CMS_SEML_READ  4

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
    //=====================================================================
    // Type of computation
    //=====================================================================
    int COMPTYPE = COMP_CM_EML2_TO_CMS_SEML;//COMP_CM_EML2_TO_CMS_SEML;//COMP_SINGLE_ORBIT;//

    //=====================================================================
    // cout settings
    //=====================================================================
    cout << setiosflags(ios::scientific) << setprecision(15);

    //=====================================================================
    // openMP settings
    //=====================================================================
    int num_threads = 4;
    omp_set_num_threads(num_threads);
    cout << "num_threads  = "  << num_threads << endl;

    //=====================================================================
    // Define the global parameters
    //=====================================================================
    //MODEL
    MODEL_TYPE = M_QBCP;
    //OFTS_ORDER
    OFTS_ORDER = 16; //IF 16 is set, the results came from the server (OFTS_ORDER = 20...) NOT the best practice!
    //OFS_ORDER
    if(MODEL_TYPE == M_RTBP) OFS_ORDER  = 0;
    else OFS_ORDER  = 30;


    //=====================================================================
    // Configuration parameters
    //=====================================================================
    int model = 0, dcs = 0, isNorm = 0, li_EM = 0, li_SEM = 0, pms_EM = 0, pms_SEM = 0, mType_EM = 0, mType_SEM = 0, reduced_nv = 0;
    switch(COMPTYPE)
    {
        //================================
        // Projection CMU EML2 to CM SEMLi
        //================================
    case COMP_CM_EML2_TO_CM_SEML:
        model     = MODEL_TYPE;
        dcs       = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
        pms_EM    = PMS_MIXED;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER;
        break;

        //================================
        // Projection CMU EML2 to CMS SEMLi
        //================================
    case COMP_CM_EML2_TO_CMS_SEML:
    case COMP_CM_EML2_TO_CMS_SEML_READ:
        model     = MODEL_TYPE;
        dcs       = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
        pms_EM    = PMS_MIXED;
        pms_SEM   = PMS_MIXED;
        mType_EM  = MAN_CENTER_U;
        mType_SEM = MAN_CENTER_S;
        break;

        //================================
        // Just some examples of solutions
        //================================
    case COMP_SINGLE_ORBIT:
        model     = MODEL_TYPE;
        dcs       = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER;
        break;

        //================================
        // Test of the vector fields
        //================================
    case COMP_VF_TEST:
        model     = MODEL_TYPE;
        dcs       = F_EM;
        isNorm    = 1;
        li_EM     = 2;
        li_SEM    = 2;
        pms_EM    = PMS_GRAPH;
        pms_SEM   = PMS_GRAPH;
        mType_EM  = MAN_CENTER;
        mType_SEM = MAN_CENTER;
        break;
    }

    //=====================================================================
    //Define the number of reduced variables
    //depending on the type of manifold
    //=====================================================================
    if(dcs == F_EM)
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

    //=====================================================================
    // Initialization of the environnement
    // Mandatory to perform any computation
    //=====================================================================
    initialize_environment(li_EM, li_SEM, isNorm, model, dcs, pms_EM, pms_SEM, mType_EM, mType_SEM);

    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc> CM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_TFC;     //center manifold in TFC coordinates
    CM_NC.reserve(6);
    for(int i = 0; i < 6; i++) CM_NC.push_back(Oftsc(reduced_nv, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_TFC.reserve(6);
    for(int i = 0; i < 6; i++) CM_TFC.push_back(Oftsc(reduced_nv, OFTS_ORDER, OFS_NV, OFS_ORDER));

    //Read from file
    readVOFTS_bin(CM_NC,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CM_TFC, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc(6,6);    ///COC matrix
    matrix<Ofsc>  Pcoc(6,6);    ///COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc(6,6);   ///COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc(6,6);   ///COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc(6);      ///COC vector

    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //--------------------------------------
    //DCM_TFC : Jacobian of CM_TFC
    //--------------------------------------
    matrix<Oftsc> DCM_TFC(6, 4, reduced_nv, OFTS_ORDER, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_TFC, SEML.cs.F_PMS+"DWf/DWhc", OFS_ORDER);

    //=====================================================================
    // Initial conditions
    //=====================================================================
    double st0[reduced_nv];
    switch(dcs)
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


    //=====================================================================
    // Computation parameters: Projection CMU EML2 to CM SEMLi
    //=====================================================================
    //TIME
    double TM    = 12*SEML.us.T;
    double TMIN  = 0.9*SEML.us.T;
    double TMAX  = 1.0*SEML.us.T;
    int    TSIZE = 0;

    //UNSTABLE MANIFOLD AT EML2
    double GMIN_S1  = -35;
    double GMAX_S1  = +35;
    int    GSIZE_S1 = 100;

    double GMIN_S3  = -35;
    double GMAX_S3  = +35;
    int    GSIZE_S3 = 100;

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
    int ISUSERDEFINED = 0;

    //Number of dimensions on which we compute the norm of the projection error
    int NOD = 6;

    //=====================================================================
    // For COMP_CM_EML2_TO_CMS_SEML_READ
    //=====================================================================
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

    //=====================================================================
    // Switch
    //=====================================================================
    switch(COMPTYPE)
    {
        //================================
        // Projection CMU EML2 to CM SEMLi
        //================================
    case COMP_CM_EML2_TO_CM_SEML:

        tic();

        //----------------------
        // Compute initial conditions in CMU EML2 on a given grid
        //----------------------
        compute_grid_CMU_EM(PROJ_EPSILON, TMIN, TMAX, TSIZE,
                            GMIN_S1, GMAX_S1, GMIN_S3, GMAX_S3, GSIZE_S1, GSIZE_S3,
                            CM_TFC, Mcoc, MIcoc, Vcoc, ISPAR);
        cout << "End of in compute_grid_CMU_EM" <<  endl;

        //----------------------
        // Integrate those initial conditions and project them on the CM SEML2
        //----------------------
        int_proj_CMU_EM_on_CM_SEM(TM, TSIZE, GSIZE_S1, GSIZE_S3, MSIZE, NSMIN, NOD,
                                  ISPAR, YNMAX, SNMAX);
        cout << "int_proj_CMU_EM_on_CM_SEM" <<  endl;

        //----------------------
        // Compute the best solutions from previous computation
        //----------------------
        //int_sorted_sol_CMU_EM_to_CM_SEM(OFTS_ORDER, MSIZE, ISPAR, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
        cout << "End of int_sorted_sol_CMU_EM_to_CM_SEM in " << toc() << endl;
        break;

        //================================
        // Projection CMU EML2 to CMS SEMLi
        //================================
    case COMP_CM_EML2_TO_CMS_SEML:


        //----------------------
        // Multiple shooting on the manifold leg only
        //ref_CMU_EM_to_CM_SEM_MSD(OFTS_ORDER, 10, 1, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting on the whole trajectory (base for JPL refinment)
        //ref_CMU_EM_to_CM_SEM_MSD_COMP(OFTS_ORDER, 100, NCSEM, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to CM (DEPRECATED)
        //ref_CMU_EM_to_CM_SEM_MSD_DEP(OFTS_ORDER, 50, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, false);

        //----------------------
        // Multiple shooting: CMU to free boundary (SEML2 orbit included)
        //ref_CMU_EM_to_CM_SEM_MSD_PART(OFTS_ORDER, 60, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

        //----------------------
        // Multiple shooting: CMU to CMS, with variable times
        ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(20, NCSEM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, ISPLANAR, ISSAVED, ISUSERDEFINED);

        //----------------------
        // Search for minimum of projection + Multiple shooting: CMU to CMS, with variable times
        //        double st_EM[5];
        //        st_EM[0] = 35;
        //        st_EM[1] = 0.0;
        //        st_EM[2] = 10.0;
        //        st_EM[3] = 0.0;
        //        st_EM[4] = PROJ_EPSILON;

        // 1. the IC in CMS of EML2 are fixed
        //ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ(MSIZE, 20, NCSEM, st_EM, T0EM, TM, CM_NC, CM_TFC, DCM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, YNMAX, SNMAX, ISPLANAR);

        // 2. Only s1_EM is fixed
        //ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_OPT(MSIZE, 20, NCSEM, st_EM, 0.9*SEML.us.T, TM, TM, GMIN_S3, GMAX_S3,
        //                                                      GSIZE_S3, CM_NC, CM_TFC, DCM_TFC,
        //                                                      Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, YNMAX, SNMAX, ISPLANAR, ISPAR);

        break;

    //================================
    // READ CMU EML2 to CMS SEMLi
    //================================
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

        //================================
        // Just some examples of solutions
        //================================
    case COMP_SINGLE_ORBIT:
        gridOrbit(st0, +6.016997770510415e+00, -2.793897158902917e+01, 1e-2, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
        break;
        //================================
        // Test of the vector fields
        //================================
    case COMP_VF_TEST:
        qbtbp_test(0.25*SEML.us.T, SEML);

        break;

    }


    return 0;
}
