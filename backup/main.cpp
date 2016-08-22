#include <iostream>
#include <omp.h>
#include "FTA.h"
#include "Ofsc.h"
#include "Oftsc.h"
#include "env.h"
#include "pmcoc.h"
#include "ode.h"
#include "vf.h"
#include "single_orbit.h"
extern "C" {
#include "gnuplot_i.h"
}

using namespace std;

/************* NOTES ********************
 Notes from Reunion with Josep (15/12/2015)

 Step 1: Find some solution "by hand"
 Step 2: Refine to get a final heteroclinic connection
         Look @ colloquial methods (Doedel) or Simo methods
 Step 3: Fix a given Pk section and fix the phase @ this section

    + what about refinment in the JPL eph?
    + Book I for an example of continuation.

Look at the domain of validity of the (un)stable manifolds of EML1,2
Can we use it at its limit?
Is it close to the Moon?

Try a "fenetre glissante" of one period to check the distance to SEMLi

Project only if dist to Li is < threshold
****************************************/


int main()
{
    cout << setiosflags(ios::scientific) << setprecision(15);

    //---------------------------------------------------------------------
    // openMP settings
    //---------------------------------------------------------------------
    int num_threads = 50;
    omp_set_num_threads(num_threads);
    cout << "num_threads       = "  << num_threads << endl;

    //---------------------------------------------------------------------
    // Define the global parameters
    //---------------------------------------------------------------------
    MODEL_TYPE = M_QBCP;
    //MODEL_TYPE = M_RTBP;
    if(MODEL_TYPE == M_RTBP) OFS_ORDER  = 0;
    else OFS_ORDER  = 30;
    OFTS_ORDER = 20;

    //---------------------------------------------------------------------
    // Retrieve the configuration parameters
    //---------------------------------------------------------------------
    int model     = MODEL_TYPE;
    int dcs       = F_EM;
    int isNorm    = 1;
    int li_EM     = 2;
    int li_SEM    = 2;
    int pms_EM    = PMS_MIXED;
    int pms_SEM   = PMS_GRAPH;
    int mType_EM  = MAN_CENTER_U;
    int mType_SEM = MAN_CENTER;


    //---------------------------------------------------------------------
    //Define the number of reduced variables
    //depending on the type of manifold
    //---------------------------------------------------------------------
    if(dcs == F_EM)
    {
        switch(mType_EM)
        {
        case MAN_CENTER:
            REDUCED_NV = 4;
            break;
        case MAN_CENTER_S:
        case MAN_CENTER_U:
            REDUCED_NV = 5;
            break;
        case MAN_CENTER_US:
            REDUCED_NV = 6;
            break;
        default:
            cout << "Warning: unknown type of manifold. REDUCED_NV = 4 by default." << endl;
            REDUCED_NV = 4;
            break;
        }
    }
    else
    {
        switch(mType_SEM)
        {
        case MAN_CENTER:
            REDUCED_NV = 4;
            break;
        case MAN_CENTER_S:
        case MAN_CENTER_U:
            REDUCED_NV = 5;
            break;
        case MAN_CENTER_US:
            REDUCED_NV = 6;
            break;
        default:
            cout << "Warning: unknown type of manifold. REDUCED_NV = 4 by default." << endl;
            REDUCED_NV = 4;
            break;
        }
    }

    //---------------------------------------------------------------------
    // Initialization of the environnement
    // Mandatory to perform any computation except qbtbp(int)
    //---------------------------------------------------------------------
    init_env(li_EM, li_SEM, isNorm, model, dcs, pms_EM, pms_SEM, mType_EM, mType_SEM);

    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc>  CM(6);     ///center manifold in NC coordinates
    vector<Oftsc> CMh(6);     ///center manifold in TFC coordinates
    //Read from file
    cout << SEML.cs.F_PMS << endl;
    readVOFTS_bin(CM,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CMh, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);

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

    //---------------------------------------------------------------------
    // Initial conditions
    //---------------------------------------------------------------------
    double st0[REDUCED_NV];
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
            if(REDUCED_NV == 5) st0[4] = 0.0;
            break;

        case 2:
            st0[0] = -2;//-20 ;
            st0[1] = 0.0;//2.47827452767491;
            st0[2] = -2;//-36;
            st0[3] = 0.0;//2.47827452767491;
            if(REDUCED_NV == 5) st0[4] = 0.0;

            break;
        }
        break;

    case F_SEM:
        //Values for Earth capture via Moon slingshot
        st0[0] = -0.0471846673120537;
        st0[1] =  0.0;
        st0[2] =  0.461268865028523;
        st0[3] =  0.0;
        if(REDUCED_NV == 5) st0[4] =  0.0;
        break;
    }

    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------
    double tMin = -1.0/2*SEML.us.T, tMax = 1.0/2*SEML.us.T;
    double gMin = -30, gMax = +30;
    int tSize = 9, gSize = 20;
    double tm = 5*SEML.us.T;
    int Nman = 2000;
    cout << "Maximum integration time: = " << tm*SEML.cs.cr3bp.T/(2*M_PI*86400*24) << " months" << endl;

    switch(dcs)
    {
    case F_EM:

        switch(mType_EM)
        {
        case MAN_CENTER_U:

            //----------------------
            // Finding connections with SEMLi (NEW)
            //----------------------
            cusMan(OFTS_ORDER, 1e-6, tMin, tMax, tSize, gMin, gMax, gSize, CMh, Mcoc, MIcoc, Vcoc, 1);
            intMan(tm, OFTS_ORDER, -1, -1, Nman, 1, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            projMan(OFTS_ORDER, OFTS_ORDER, -1, -1, -1, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, 1);

            //----------------------
            // Finding connections with SEMLi (NEW)
            //----------------------
            //            int Size = floor(fabs(st0[0]));
            //            int Nt = 20;
            //            int Nman = 1000;

            //            //1. Building an EML2 orbit on a grid
            //            gridOrbit(st0, 0, 2*SEML.us.T, +1e-2, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //            //2. Computing the unstable directions along this orbit
            //            cusOrbit(OFTS_ORDER, Size, +1e-5, CMh, Mcoc, MIcoc, Vcoc, 1);
            //            //3. Integrating forward these unstable directions on a grid
            //            manOrbit(5*SEML.us.T, OFTS_ORDER, Size, -1, Nman, 1, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //            //4. Sort the manifold branches wrt to a certain criterion
            //            manOrbitPostProcess(OFTS_ORDER, Size, -1, 1);
            //            //5. Project the manifold branches on the center manifold of SEML1,2
            //            manOrbitProj(OFTS_ORDER, Size, TYPE_MAN_SORT_DR, OFTS_ORDER, -1, Nman, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, 1);
            //            //6. Diplay the solution, or...
            //            //manOrbitPlot(OFTS_ORDER, Size, Nt, Nman, TYPE_MAN_SORT_DR);
            //            //6. Refine the solutions to actually target the projected state on the center manifold
            //            manOrbitLambert(OFTS_ORDER, Size, TYPE_MAN_PROJ, OFTS_ORDER, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);

            //----------------------
            // Just some examples of solutions
            //----------------------
            //orbit_cus(st0, 0, 5*SEML.us.T, 8*SEML.us.T, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //----------------------
            // Finding connections with SEMLi (OLD)
            //----------------------
            //eml2Tosemli2(st0, 6.481139365237600e+00-1e-1, 6.481139365237600e+00+1e-1, 9*SEML.us.T, 1e-3, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //eml2Tosemli2(st0, 0, SEML.us.T, 9*SEML.us.T, 1e-3, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //eml2Tosemli_torb_vs_tman(st0, -0.1*SEML.us.T, 0.9*SEML.us.T, 8*SEML.us.T, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);


            break;
        case MAN_CENTER_S:
            //----------------------
            // Just some examples of solutions
            //----------------------
            //orbit_cus(st0, 0, 5*SEML.us.T, -2*SEML.us.T, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //----------------------
            // Finding connections with SEMLi
            //----------------------
            eml2Tosemli2(st0, 0, 1*SEML.us.T, -5*SEML.us.T, 5e-2, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            break;
        case MAN_CENTER:
            cout << "Nothing is done here" << endl;
            break;
        }
        break;


    case F_SEM:

        //Values of time on the manifold
        if(SEML.model == M_QBCP) tm = 7*SEML.us.T;
        if(SEML.model == M_RTBP) tm = 12*SEML.us.T;
        switch(mType_SEM)
        {
        case MAN_CENTER_U:
            //----------------------
            // Just some examples of solutions
            //----------------------
            //orbit_cu_sem(st0, 0, 10*SEML.us.T, tm, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            //----------------------
            // Finding connections with EMLi
            //----------------------
            //semliToseml2(st0, 0, 10*SEML.us.T, 12*SEML.us.T, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            break;
        case MAN_CENTER_S:
            //----------------------
            // Just some examples of solutions
            //----------------------
            orbit_cu_sem(st0, 0, 10*SEML.us.T, -tm, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            break;
        case MAN_CENTER:
            //----------------------
            // Just some examples of solutions
            //----------------------
            gridOrbit(st0, 0, 4*SEML.us.T, 1e-2, CM, CMh, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
            break;
        }
        break;
    }

    return 0;
}
