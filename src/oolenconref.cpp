#include "oolenconref.h"




//========================================================================================
//
//          SUBROUTINES:
//
//========================================================================================
/**
 *  \brief Computation of the matrix Phi0. This matrix gives at first order the linear
 *         relationship between the parameterization at EML2 (4 dimensions) at t = t0 and
 *         the final state at the first patch point at t = t1, in coord_type coordinates
 *         (6 dimensions).
 *         It is computed thanks to the following equation:
 *
 *              Phi0 =  J1 * COORD_J_RCM
 *
 *         Where:   - J1 is the 6 x 6 state transition matrix between t0 and t1.
 *                  - COORD_J_RCM is the jacobian matrix dz/ds, where z is the state in
 *                    coord_type coordinates and s is the state in RCM coordinates,
 *                    inside EML2 center-unstable manifold. COORD_J_RCM is a 6 x 5 matrix.
 *
 *         Inputs:  - Phi0, the 6 x 5 gsl_matrix to update.
 *                  - J1, the the 6 x 6 state transition matrix between t0 and t1.
 *                  - orbit_EM, the Orbit object that contains the center-unstable
 *                    manifold at EML2, and therefore can be used to evaluate COORD_J_RCM.
 *                  - t0_EM, the initial time t0, in EM units.
 *                  - coord_type, the final desired coordinate system.
 *
 *         Note: no testing is made on the inputs (dimensions, manifold type...). These
 *         tests should be made by the user higher in the code.
 **/
int ftc_compute_phi0(gsl_matrix* Phi0, gsl_matrix* J1, Orbit& orbit_EM, double t0_EM, int coord_type)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    gsl_matrix* COORD_J_RCM  = gsl_matrix_calloc(6,5);

    //====================================================================================
    // 2. Computation
    //====================================================================================
    //Evaluate COORD_J_RCM(orbit_EM.si, t0), in EM units, in R(6,5)
    orbit_EM.getInvman()->evalDRCMtoCOORD(orbit_EM.getSi(), t0_EM, COORD_J_RCM, OFTS_ORDER, OFS_ORDER, coord_type);

    //Phi0 = J1*COORD_J_RCM, in R(6,6)*R(6,5) = R(6,5)
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, J1, COORD_J_RCM, 0.0, Phi0);


    //====================================================================================
    // 3. Free
    //====================================================================================
    gsl_matrix_free(COORD_J_RCM);

    return(GSL_SUCCESS);
}

/**
 *  \brief Selection of a differential corrector (msft3d, msvt3d, etc).
 **/
diffcorrptr ftc_select_diffcorr(RefSt& refSt)
{
    switch(refSt.dim)
    {
    case REF_3D:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return msft3d;

        case REF_VAR_TIME:
        case REF_VAR_TN:
            return msvt3d;

        default:
            cerr << "ftc_select_diffcorr" << ". Unknown dcs. refSt.time = " << ". msft3d is returned by default." << endl;
            return msft3d;
        }
        break;
    case REF_MIXED:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return msftmixed;

        case REF_VAR_TIME:
        case REF_VAR_TN:
            return msvtmixed;

        default:
            cerr << "ftc_select_diffcorr" << ". Unknown dcs. refSt.time = " << ". msftmixed is returned by default." << endl;
            return msftmixed;
        }
        break;
    case REF_PLANAR:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return msftplan;
            //return msvftplan;
            //return msftplan_pa;
            //return msvftplan_dH;
            //return msvftplan_dte;

        case REF_VAR_TIME:
            return msvtplan;

        case REF_VAR_TN:
            return msvltplan;

        default:
            cerr << "ftc_select_diffcorr" << ". Unknown dcs. refSt.time = " << ". msftplan is returned by default." << endl;
            return msftplan;
        }
        break;

        default:
            cerr << "ftc_select_diffcorr" << ". Unknown dcs. refSt.dim = " << ". msftplan is returned by default." << endl;
            return msftplan;
    }
}


/**
 *  \brief Selection of a differential predictor (ufvarft3d, ufvarvt3d, etc).
 *         Must be coherent with ftc_select_diffcorr, just above.
 **/
predictorptr ftc_select_predictor(RefSt& refSt)

{
    switch(refSt.dim)
    {
    case REF_3D:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return ufvarft3d;

        case REF_VAR_TIME:
        case REF_VAR_TN:
            return ufvarvt3d;

        cerr << "ftc_select_predictor" << ". Unknown dcs. refSt.time = " << ". ufvarft3d is returned by default." << endl;
            return ufvarft3d;
        }
        break;

    case REF_MIXED:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return ufvarftmixed;

        case REF_VAR_TIME:
        case REF_VAR_TN:
            return ufvarvtmixed;

        cerr << "ftc_select_predictor" << ". Unknown dcs. refSt.time = " << ". ufvarftmixed is returned by default." << endl;
            return ufvarftmixed;
        }
        break;

    case REF_PLANAR:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
            return ufvarftplan;
            //return ufvarvftplan;
            //return ufvarvftplan_dH;

        case REF_VAR_TIME:
            return ufvarvtplan;

        case REF_VAR_TN:
            return ufvarvltplan;

        cerr << "ftc_select_predictor" << ". Unknown dcs. refSt.time = " << ". ufvarftplan is returned by default." << endl;
            return ufvarftplan;
        }
        break;

        default:
            cerr << "ftc_select_predictor" << ". Unknown dcs. refSt.dim = " << ". ufvarftplan is returned by default." << endl;
            return ufvarftplan;
    }

    return ufvarftplan;
}

/**
 *  \brief Yields the number of free variables necessary to compute the refinment procedure.
 **/
int nfreevariables(RefSt refSt, int mgs)
{
    int nfv = 0;
    switch(refSt.dim)
    {
    case REF_3D:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
                nfv = 6*mgs+3;
            break;

        case REF_VAR_TIME:
        case REF_VAR_TN:
                nfv = 7*mgs+3;
            break;
        default:
            perror("nfreevariables. Unknown refSt.time.");
        break;
        }
        break;
    case REF_MIXED:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
                nfv = 6*mgs+1;
            break;

        case REF_VAR_TIME:
        case REF_VAR_TN:
                nfv = 7*mgs+1;
            break;
        default:
            perror("nfreevariables. Unknown refSt.time.");
        break;
        }
        break;
    case REF_PLANAR:
        switch(refSt.time)
        {
        case REF_FIXED_TIME:
               nfv = 4*mgs+1;
            break;

        case REF_VAR_TN:
               nfv = 4*mgs+2;
            break;

        case REF_VAR_TIME:
               nfv = 5*mgs+1;
            break;

        default:
            perror("nfreevariables. Unknown refSt.time.");
        break;
        }
         break;
    default:
            perror("nfreevariables. Unknown refSt.dim.");
        break;
    }

    return nfv;
}

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS with new implementation
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msft3d(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
           int nov, int mgs, int coord_type,  int isFirst,
           Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msft3d";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE)
    {
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE)
    {
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double** ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double* tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*mgs+3;  //free variables
    int ncs = 6*mgs;    //constraints

    // Correction vector at patch points
    gsl_vector* DQv = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector* Fv  = gsl_vector_calloc(ncs);
    // Error isolated at final point
    gsl_vector* Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix** Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix* DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix* Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix* Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix* PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {

                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == mgs - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //------------------------------------------------------------------------
                //DF/DS
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i + 6*k, j + 6*k-2,       gsl_matrix_get(Ji[k], i, j));
                    for(int j = 0; j < 5; j++) gsl_matrix_set(DF, i + 6*k, j + 6*(k+1)-2,  -gsl_matrix_get(PhiN, i, j));
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 6*k-2,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 6*(k+1)-2, -gsl_matrix_get(Id, i, j));
                    }
                }

            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        // Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }


        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }


        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_EM.addSi(gsl_vector_get(DQv, i), i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 6*k-2);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.addSi(gsl_vector_get(DQv, i+ 6*mgs-2), i);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        //cout << "--------------------" << endl;
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }

        //--------------------------------------------------------------------------------
        // Update number of iterations
        //--------------------------------------------------------------------------------
        iter++;
    }


    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1,  ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    //------------------------------------------------------------------------------------
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(refSt.isDirUD && refSt.isCont())
        {
            do
            {
                cout << "-------------------------------------------------------" << endl;
                cout << "After refinement: s1_CMU_EM = " << orbit_EM.getSi()[0]   << endl;
                cout << "Please choose a direction for the cont. procedure:"      << endl;
                cout << "+1: s1_CMU_EM is increasing"                             << endl;
                cout << "-1: s1_CMU_EM is decreasing"                             << endl;
                cout << " to select a specific starting time"                     << endl;
                cin >> dti;
                ti = (int) dti;
            }
            while(ti != 1 && ti != -1);
        }
        else ti = refSt.Dir;

        //Here, we want to make s_EM[1] "grow"
        sign = nullvector_temp[1] > 0? ti:-ti;
    }
    else
    {
        //OR
        //Here, we want to go "in the same direction" for some components of Q
        dotNV += nullvector_temp[0]*nullvector[0];                  //CMU of  EML2
        //dotNV += gsl_matrix_get(Q, 1, nfv-1)*nullvector[1];       //CMU of  EML2
        dotNV += nullvector_temp[2]*nullvector[2];                  //CMU of  EML2
        //dotNV += gsl_matrix_get(Q, 3, nfv-1)*nullvector[3];       //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEMLi
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);

    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(Fvn);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 *
 *        Difference with msft3d: only s1 and s3 are corrected at EML2
 **/
int msftmixed(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
               int nov, int mgs, int coord_type,  int isFirst,
               Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msftmixed";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE)
    {
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE)
    {
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double** ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double* tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*mgs+1;  //free variables
    int ncs = 6*mgs;    //constraints

    // Correction vector at patch points
    gsl_vector* DQv = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector* Fv  = gsl_vector_calloc(ncs);
    // Error isolated at final point
    gsl_vector* Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix** Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix* DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix* Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix* Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix* PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {

                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == mgs - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //------------------------------------------------------------------------
                //DF/DS
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                    gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));

                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i + 6*k, j + 6*k-4,       gsl_matrix_get(Ji[k], i, j));
                    for(int j = 0; j < 5; j++) gsl_matrix_set(DF, i + 6*k, j + 6*(k+1)-4,  -gsl_matrix_get(PhiN, i, j));
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 6*k-4,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 6*(k+1)-4, -gsl_matrix_get(Id, i, j));
                    }
                }

            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        // Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }


        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }


        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 6*k-4);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.addSi(gsl_vector_get(DQv, i+ 6*mgs-4), i);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        //cout << "--------------------" << endl;
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }

        //--------------------------------------------------------------------------------
        // Update number of iterations
        //--------------------------------------------------------------------------------
        iter++;
    }


    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1,  ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(refSt.isDirUD && refSt.isCont())
        {
            do
            {
                cout << "-------------------------------------------------------" << endl;
                cout << "After refinement: s1_CMU_EM = " << orbit_EM.getSi()[0]   << endl;
                cout << "Please choose a direction for the cont. procedure:"      << endl;
                cout << "+1: s1_CMU_EM is increasing"                             << endl;
                cout << "-1: s1_CMU_EM is decreasing"                             << endl;
                cout << " to select a specific starting time"                     << endl;
                cin >> dti;
                ti = (int) dti;
            }
            while(ti != 1 && ti != -1);
        }
        else ti = refSt.Dir;

        //Here, we want to make s_EM[1] "grow"
        sign = nullvector_temp[1] > 0? ti:-ti;
    }
    else
    {
        //OR
        //Here, we want to go "in the same direction for some components of Q
        for(int i = 0; i < 2; i++) dotNV += nullvector_temp[i]*nullvector[i];    //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //        dotNV += gsl_matrix_get(Q, nfv-4, nfv-1)*nullvector[nfv-4];                //CMS of  SEMLi
        //        dotNV += gsl_matrix_get(Q, nfv-2, nfv-1)*nullvector[nfv-2];                //CMS of  SEMLi
        //        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEMLi
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];


    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);

    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(Fvn);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}


/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 **/
int msvt3d(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
           int nov, int mgs, int coord_type,
            int isFirst,
           Orbit &orbit_EM, Orbit &orbit_SEM,
           gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //Name of the routine
    string fname = "msvt3d";

    //====================================================================================
    // 0. Test on coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cerr << fname << ". Wrong coord_type. Must be NCSEM." << endl;
        return FTC_EDOM;
    }


    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov], f[6], te;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 7*mgs+3;  //free variables
    int ncs = 6*mgs;    //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //For time derivatives
    double z1[6];

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];
            te = tm[1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == mgs - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            //----------------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &odeParams);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &odeParams);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == mgs-1)
            {
                // z1 = dCM_SEM_NC/dt[k+1]
                orbit_SEM.getInvman()->evaldotRCMtoNC(orbit_SEM.getSi(), tmdn[mgs], z1, OFTS_ORDER, OFS_ORDER);

                // Then f = f - z1
                for(int i = 0; i < 6; i++) f[i] -= z1[i];
            }


            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //------------------------------------------------------------------------
                //DF/DS
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i + 6*k, j + 7*k-3,       gsl_matrix_get(Ji[k], i, j));
                    for(int j = 0; j < 5; j++) gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-3,  -gsl_matrix_get(PhiN, i, j));
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 7*k-3,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-3, -gsl_matrix_get(Id, i, j));
                    }
                }


                //------------------------------------------------------------------------
                //DF/DT
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    //--------------------------------------------------------------------
                    //dF[0]/dt[0]
                    //--------------------------------------------------------------------
                    //NOTHING IS DONE FOR NOW

                    //--------------------------------------------------------------------
                    //dF[0]/dt[1]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }
                else if(k == mgs-1)
                {
                    //--------------------------------------------------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*mgs-4, gsl_vector_get(K4, i));

                    //--------------------------------------------------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*mgs+2, f[i]);
                }
                else
                {
                    //--------------------------------------------------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+1)-4, gsl_vector_get(K4, i));

                    //--------------------------------------------------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }

            }


        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }


        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_EM.addSi(gsl_vector_get(DQv, i), i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < nov; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 7*k-3);
            tmdn[k] = max(tmdn[k] + gsl_vector_get(DQv, 7*(k+1)-4), tmdn[k-1]);
            //tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-4);
        }
        //Last time:
        tmdn[mgs] = max(tmdn[mgs]+ gsl_vector_get(DQv, 7*mgs+2), tmdn[mgs-1]);
        //tmdn[mgs] += gsl_vector_get(DQv, 7*mgs+2);


        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.addSi(gsl_vector_get(DQv, i + 7*mgs-3), i);

        //Updating CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }


    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;


    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);


    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2]  > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }
    else
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //------------------------------------------------------------------------------------
    // 4. Free
    //------------------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);

    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(Fvn);
    gsl_vector_free(Kf);
    gsl_vector_free(K4);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);


    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        Difference with msvt3d: only s1 and s3 are corrected at EML2
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 **/
int msvtmixed(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
               int nov, int mgs, int coord_type,
                int isFirst,
               Orbit &orbit_EM, Orbit &orbit_SEM,
               gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //Name of the routine
    string fname = "msvtmixed";

    //====================================================================================
    // 0. Test on coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cerr << fname << ". Wrong coord_type. Must be NCSEM." << endl;
        return FTC_EDOM;
    }


    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);


    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov], f[6], te;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 7*mgs+1;  //free variables
    int ncs = 6*mgs;    //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //For time derivatives
    double z1[6];

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];
            te = tm[1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == mgs - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            //----------------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &odeParams);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &odeParams);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == mgs-1)
            {
                // z1 = dCM_SEM_NC/dt[k+1]
                orbit_SEM.getInvman()->evaldotRCMtoNC(orbit_SEM.getSi(), tmdn[mgs], z1, OFTS_ORDER, OFS_ORDER);

                // Then f = f - z1
                for(int i = 0; i < 6; i++) f[i] -= z1[i];
            }


            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //------------------------------------------------------------------------
                //DF/DS
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    //for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                    gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));

                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i + 6*k, j + 7*k-5,       gsl_matrix_get(Ji[k], i, j));
                    for(int j = 0; j < 5; j++) gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-5,  -gsl_matrix_get(PhiN, i, j));
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 7*k-5,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-5, -gsl_matrix_get(Id, i, j));
                    }
                }


                //------------------------------------------------------------------------
                //DF/DT
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    //--------------------------------------------------------------------
                    //dF[0]/dt[0]
                    //--------------------------------------------------------------------
                    //NOTHING IS DONE FOR NOW

                    //--------------------------------------------------------------------
                    //dF[0]/dt[1]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-6, f[i]);
                }
                else if(k == mgs-1)
                {
                    //--------------------------------------------------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*mgs-6, gsl_vector_get(K4, i));

                    //--------------------------------------------------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*mgs, f[i]);
                }
                else
                {
                    //--------------------------------------------------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+1)-6, gsl_vector_get(K4, i));

                    //--------------------------------------------------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------------------------------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-6, f[i]);
                }

            }


        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }


        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < nov; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 7*k-5);
            tmdn[k] = max(tmdn[k] + gsl_vector_get(DQv, 7*(k+1)-6), tmdn[k-1]);
            //tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-6);
        }
        //Last time:
        tmdn[mgs] = max(tmdn[mgs]+ gsl_vector_get(DQv, 7*mgs), tmdn[mgs-1]);
        //tmdn[mgs] += gsl_vector_get(DQv, 7*mgs);


        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.addSi(gsl_vector_get(DQv, i + 7*mgs-5), i);

        //Updating CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }


    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;


    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);


    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    //double dotNV = 0.0;
    if(isFirst)
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }
    else
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //------------------------------------------------------------------------------------
    // 4. Free
    //------------------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);

    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(Fvn);
    gsl_vector_free(Kf);
    gsl_vector_free(K4);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}


/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 **/
int msftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int nov, int mgs,
             int coord_type,  int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM,
             gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msftplan";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];


    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs + 1;    //free variables
    int ncs = 4*mgs;        //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            jacftplan(k, DF, mgs, Ji, Phi0, PhiN, Id);
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        //Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //------------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //------------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates. Again, we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-2), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }


        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1,  ymdn[0],       ymdn[1],       ymdn[2],     mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs], &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0],   &ymdn[1][0],   &ymdn[2][0],   1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    //3. Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(refSt.isDirUD && refSt.isCont())
        {
            do
            {
                cout << "-------------------------------------------------------" << endl;
                cout << "After refinement: s1_CMU_EM = " << orbit_EM.getSi()[0]   << endl;
                cout << "Please choose a direction for the cont. procedure:"      << endl;
                cout << "+1: s1_CMU_EM is increasing"                             << endl;
                cout << "-1: s1_CMU_EM is decreasing"                             << endl;
                cout << " to select a specific starting time"                     << endl;
                cin  >> dti;
                ti = (int) dti;
            }
            while(ti != 1 && ti != -1);
        }
        else ti = refSt.Dir;

        //Here, we want to make s_EM[0] "grow"
        sign = nullvector_temp[0] > 0? ti:-ti;
    }
    else
    {
        //Here, we want to go "in the same direction for some components of Q
        for(int i = 0; i < 2; i++) dotNV += nullvector_temp[i]*nullvector[i];    //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEMLi
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. FINALLY, we update orbit_EM.tf, even if it is not necessary
    //====================================================================================
    orbit_EM.setTf(tmdn[mgs]/SEML.us_em.ns);

    //====================================================================================
    // 5. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t1,..., tN are free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 **/
int msvtplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
             int nov, int mgs, int coord_type,  int isFirst,
             Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //Name of the routine
    string fname = "msvtplan";

    //====================================================================================
    // 0. Test on coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cerr << fname << ". Wrong coord_type. Must be NCSEM." << endl;
        return FTC_EDOM;
    }

    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov], f[6], te;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 5*mgs+1;  //number of free variables
    int ncs = 4*mgs;    //number of constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //For time derivatives
    double z1[6];

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];
            te = tm[1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            //----------------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &odeParams);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &odeParams);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == mgs-1)
            {
                // z1 = dCM_SEM_NC/dt[k+1]
                orbit_SEM.getInvman()->evaldotRCMtoNC(orbit_SEM.getSi(), tmdn[mgs], z1, OFTS_ORDER, OFS_ORDER);


                // Then f = f - z1
                for(int i = 0; i < 6; i++) f[i] -= z1[i];
            }


            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                //---------------------------
                //DF/DS
                //---------------------------
                if(k == 0)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                        gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i+1, 0));
                        gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i+1, 2));
                    }

                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-8,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-7,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-6,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-5,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*mgs-3,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-2,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-1,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-8,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-7,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-6,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-5,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*mgs-3,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-2,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*mgs-1,  -gsl_matrix_get(PhiN, i+1, 4));
                    }
                }
                else
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*k-3,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*k-2,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*k-1,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*k-0,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*k+2,  -gsl_matrix_get(Id, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*k+3,  -gsl_matrix_get(Id, i, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*k+4,  -gsl_matrix_get(Id, i, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*k+5,  -gsl_matrix_get(Id, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*k-3,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*k-2,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*k-1,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*k-0,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*k+2,  -gsl_matrix_get(Id, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*k+3,  -gsl_matrix_get(Id, i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*k+4,  -gsl_matrix_get(Id, i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*k+5,  -gsl_matrix_get(Id, i+1, 4));
                    }

                }


                //---------------------------
                //DF/DT
                //---------------------------
                if(k == 0)
                {
                    //--------------------------
                    //dF[0]/dt[0]
                    //--------------------------
                    //NOTHING IS DONE FOR NOW

                    //--------------------------
                    //dF[0]/dt[1]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i, 6, f[i]);
                    else  gsl_matrix_set(DF, i, 6, f[i+1]);
                }
                else if(k == mgs-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*mgs-4, gsl_vector_get(K4, i));
                    else gsl_matrix_set(DF, i + 4*k, 5*mgs-4, gsl_vector_get(K4, i+1));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*mgs, f[i]);
                    else gsl_matrix_set(DF, i + 4*k, 5*mgs, f[i+1]);
                }
                else
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*k+1, gsl_vector_get(K4, i));
                    else gsl_matrix_set(DF, i + 4*k, 5*k+1, gsl_vector_get(K4, i+1));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*k+6, f[i]);
                    else gsl_matrix_set(DF, i + 4*k, 5*k+6, f[i+1]);
                }

            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 0 + 5*k-3);
            ymdn[1][k] += gsl_vector_get(DQv, 1 + 5*k-3);
            ymdn[3][k] += gsl_vector_get(DQv, 2 + 5*k-3);
            ymdn[4][k] += gsl_vector_get(DQv, 3 + 5*k-3);
            tmdn[k] = max(tmdn[k] + gsl_vector_get(DQv, 5*k+1), tmdn[k-1]);
            //tmdn[k] =tmdn[k] + gsl_vector_get(DQv, 5*k+1);
        }
        //Last time:
        tmdn[mgs] = max(tmdn[mgs]+ gsl_vector_get(DQv, 5*mgs), tmdn[mgs-1]);
        //tmdn[mgs] = tmdn[mgs] + gsl_vector_get(DQv, 5*mgs);


        //--------------------------------------------------------------------------------
        //Last 3 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 5*mgs-3), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 5*mgs-2), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 5*mgs-1), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];


        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }

        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    // Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }
    else
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 5*mgs
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);


    return FTC_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The time tN alone is free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 **/
int msvltplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
              int nov, int mgs, int coord_type,  int isFirst,
              Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //Name of the routine
    string fname = "msvltplan";

    //====================================================================================
    // 0. Test on coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cerr << fname << ". Wrong coord_type. Must be NCSEM." << endl;
        return FTC_EDOM;
    }

    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov], f[6], te;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs+2;  //number of free variables
    int ncs = 4*mgs;    //number of constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //For time derivatives
    double z1[6];

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];
            te = tm[1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            // Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //------------------------------------------------------------------------
                //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
                //------------------------------------------------------------------------
                //Computing f[Q[k+1], t[k+1])
                vf(te, ye, f, &odeParams);

                //------------------------------------------------------------------------
                // Then f = f - z1, with z1 = dCM_SEM_NC/dt[k+1]
                //------------------------------------------------------------------------
                orbit_SEM.getInvman()->evaldotRCMtoNC(orbit_SEM.getSi(), tmdn[mgs], z1, OFTS_ORDER, OFS_ORDER);
                for(int i = 0; i < 6; i++) f[i] -= z1[i];
            }


            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                //------------------------------------------------------------------------
                //DF/DS
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                        gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i+1, 0));
                        gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i+1, 2));
                    }

                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
                }
                else if(k == mgs-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i+1, 4));
                    }
                }
                else
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i+1, 4));
                    }
                }


                //------------------------------------------------------------------------
                //DF/DT
                //------------------------------------------------------------------------
                if(k == 0)
                {
                    //NOTHING IS DONE FOR NOW
                }
                else if(k == mgs-1)
                {
                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 4*mgs+1, f[i]);
                    else gsl_matrix_set(DF, i + 4*k, 4*mgs+1, f[i+1]);
                }
                else
                {
                    //NOTHING IS DONE FOR NOW
                }

            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        // Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }
        //Last time:
        tmdn[mgs] = max(tmdn[mgs] + gsl_vector_get(DQv, 4*mgs+1), tmdn[mgs-1]);

        //--------------------------------------------------------------------------------
        //Last 3 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-2), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];


        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }

        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_temp[nfv];
    nullvectorjac(nullvector_temp, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = 4*mgs+1
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }
    else
    {
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            //----------------------------------------------------------------------------
            // First type of condition: stop when we are close enough to the
            // center manifold (unstable component is small enough)
            // So, here, we want to make s_SEM[4] "decrease"
            //----------------------------------------------------------------------------
            if(orbit_SEM.getSi()[4] > 0) sign = nullvector_temp[nfv-2] < 0? 1:-1;
            else sign = nullvector_temp[nfv-2] > 0? 1:-1;
            break;
        }
        case REF_COND_T:
        {
            //----------------------------------------------------------------------------
            // Another possible condition: enough turns around SEMLi. So, here,
            // we just want to increase the last time, at position nfv-1 = nfv-1 = 4*mgs+1
            //----------------------------------------------------------------------------
            sign = nullvector_temp[nfv-1] > 0? 1:-1;
            break;
        }
        }
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_temp[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);


    return FTC_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The time t0 alone is free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 **/
int msvftplan(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
              int nov, int mgs, int coord_type,  int isFirst,
              Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msftplan";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);


    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];


    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs + 1;    //free variables
    int ncs = 4*mgs;        //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *DF2  = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            // Special case of the first point
            //----------------------------------------------------------------------------
            if(k == 0) t0jac(K4, tmdn, ymdn, orbit_EM, orbit_SEM, odeParams, vf, Ji[0], Kf);

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            jacftplan(k, DF, mgs, Ji, Phi0, PhiN, Id);

            //----------------------------------------------------------------------------
            // Update DF2
            //----------------------------------------------------------------------------
            jacvftplan(k, DF2, refSt.sidim, mgs, Ji, Phi0, PhiN, K4, Id);
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        //Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //------------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //------------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates. Again, we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //------------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //------------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }

        //------------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //------------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-2), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }


        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    // 3. Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_new[nfv];
    nullvectorjac(nullvector_new, DF2, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 1
        //--------------------------------------------------------------------------------
        sign = nullvector_new[1] > 0? -1:+1;
    }
    else
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 1
        //--------------------------------------------------------------------------------
        sign = nullvector_new[1] > 0? -1:+1;
    }

    //Finally, we update nullvector with the right direction
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_new[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. FINALLY, we update orbit_EM.tf, even if it is not necessary
    //====================================================================================
    orbit_EM.setTf(tmdn[mgs]/SEML.us_em.ns);

    //====================================================================================
    // 5. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The times t0,..., tN are fixed.
 *        - The null vector associated to the solution is computed.
 *
 *   Contrary to msftplan, the pseudo-arclentgh constraint is added to the vector of constraints when !isFirst.
 **/
int msftplan_pa(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                int nov, int mgs, int coord_type,  int isFirst,
                Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msftplan_pa";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];


    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs + 1;    //free variables
    int ncs; //constraints
    if(isFirst)  ncs = 4*mgs;
    else  ncs = 4*mgs + 1;

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    double vecfv[nfv];
    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Update the vector of free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 2 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        vecfv[0] = orbit_EM.getSi(0);
        vecfv[1] = orbit_EM.getSi(2);

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            vecfv[4*k-2] = ymdn[0][k];
            vecfv[4*k-1] = ymdn[1][k];
            vecfv[4*k-0] = ymdn[3][k];
            vecfv[4*k+1] = ymdn[4][k];
        }

        //--------------------------------------------------------------------------------
        //Last 3 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        vecfv[4*mgs-2] = orbit_SEM.getSi(0);
        vecfv[4*mgs-1] = orbit_SEM.getSi(2);
        vecfv[4*mgs-0] = orbit_SEM.getSi(4);


        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            jacftplan(k, DF, mgs, Ji, Phi0, PhiN, Id);
        }

        //--------------------------------------------------------------------------------
        // Last component is the pseudo-arclength constraint, if !isFirst
        //--------------------------------------------------------------------------------
        if(!isFirst)
        {
            double res = 0.0;
            for(int i = 0; i < nfv; i++) res += (vecfv[i] - refSt.goodvector[i])*nullvector[i];
            res -=  refSt.dsc;
            gsl_vector_set(Fv, ncs-1, res);

            cout << "goodvector = " << endl;
            vector_printf_prec(refSt.goodvector, 2);

            cout << "nullvector = " << endl;
            vector_printf_prec(nullvector, 2);

            cout << "refSt.dsc = " << refSt.dsc << endl;

            cout << "gsl_vector_get(Fv, ncs-1) = " << gsl_vector_get(Fv, ncs-1)  << endl;

            pressEnter(true);

            // The last column is the transpose of the null vector
            for(int i = 0; i < nfv; i++) gsl_matrix_set(DF, ncs-1, i, nullvector[i]);
        }


        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        //Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        if(isFirst) status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        else ftc_corrvec_square(DQv, Fv, DF, ncs);

        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates. Again, we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-2), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }


        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    // Update the good vector of free variables
    //====================================================================================
    for(int i = 0; i < nfv; i++) refSt.goodvector[i] = vecfv[i];

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    // 3. Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    //Select a 4*mgs x 4*mgs+1 view in DF
    gsl_matrix_view DFview = gsl_matrix_submatrix (DF, 0, 0, 4*mgs, 4*mgs+1);

    //QR elements
    gsl_vector *work  = gsl_vector_calloc(4*mgs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,4*mgs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,4*mgs);

    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, &DFview.matrix);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(refSt.isDirUD && refSt.isCont())
        {
            do
            {
                cout << "-------------------------------------------------------" << endl;
                cout << "After refinement: s1_CMU_EM = " << orbit_EM.getSi()[0]   << endl;
                cout << "Please choose a direction for the cont. procedure:"      << endl;
                cout << "+1: s1_CMU_EM is increasing"                             << endl;
                cout << "-1: s1_CMU_EM is decreasing"                             << endl;
                cout << " to select a specific starting time"                     << endl;
                cin  >> dti;
                ti = (int) dti;
            }
            while(ti != 1 && ti != -1);
        }
        else ti = refSt.Dir;

        //Here, we want to make s_EM[0] "grow"
        sign = gsl_matrix_get(Q, 0, nfv-1) > 0? ti:-ti;
    }
    else
    {
        //Here, we want to go "in the same direction for some components of Q
        for(int i = 0; i < 2; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];    //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEMLi
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    gsl_vector_view nvview = gsl_matrix_subcolumn(Q, nfv-1, 0, nfv);
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_vector_get(&nvview.vector, i)/gsl_blas_dnrm2(&nvview.vector);

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. FINALLY, we update orbit_EM.tf, even if it is not necessary
    //====================================================================================
    orbit_EM.setTf(tmdn[mgs]/SEML.us_em.ns);

    //====================================================================================
    // 5. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);
    gsl_vector_free(work);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_matrix_free(DFT);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The time t0 alone is free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 *
 *        A constraint H(0) = cst is added.
 **/
int msvftplan_dH(double **ymd, double *tmd, double **ymdn, double *tmdn, double *nullvector,
                 int nov, int mgs, int coord_type,  int isFirst,
                 Orbit &orbit_EM, Orbit &orbit_SEM, gnuplot_ctrl *h1, RefSt &refSt, int *niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msvftplan_dH";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE){
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);


    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];


    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs + 2;    //free variables
    int ncs = 4*mgs + 1;    //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //Final column: energy at the origin
            if(refSt.dH > 0.0)
            {
                double z[6], H, H0;
                for(int i = 0; i < 6; i++) z[i] = 0.0;
                H0 = qbcp_H_complete(tmdn[0], z, coord_type, coord_type);
                for(int i = 0; i < 6; i++) z[i] = ymdn[i][0];
                H = qbcp_H_complete(tmdn[0], z, coord_type, coord_type);
                gsl_vector_set(Fv, ncs-1, H - H0 - refSt.dH);
            }
            else gsl_vector_set(Fv, ncs-1, 0.0); //constraint automatically satisfied

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            // Special case of the first point
            //----------------------------------------------------------------------------
            if(k == 0) t0jac(K4, tmdn, ymdn, orbit_EM, orbit_SEM, odeParams, vf, Ji[0], Kf);

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            jacftplan_dt0(k, DF, refSt.sidim, mgs, Ji, Phi0, PhiN, K4, Id);

            //Final column: energy at the origin
            if(k == 0)
            {
                gsl_matrix* COORD_J_RCM  = gsl_matrix_calloc(6,5);
                gsl_vector* K5  = gsl_vector_calloc(5);
                gsl_vector* dH  = gsl_vector_calloc(6);

                //Evaluate COORD_J_RCM(orbit_EM.si, t0), in EM units, in R(6,5)
                orbit_EM.getInvman()->evalDRCMtoCOORD(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns, COORD_J_RCM, OFTS_ORDER, OFS_ORDER, coord_type);

                //K5 = COORD_J_RCM' * dH/dz   in R(5,6) * R(6) = R(5)
                gsl_blas_dgemv (CblasTrans, 1.0, COORD_J_RCM, dH, 0.0, K5);


                //Hamiltonian derivatives:
                double z[6];
                for(int i = 0; i < 6; i++) z[i] = 0.0;
                double Hdot0 = qbcp_Hn_SEM_dot(tmdn[0], z, &SEML);

                for(int i = 0; i < 6; i++) z[i] = ymdn[i][0];
                double Hdot = qbcp_Hn_SEM_dot(tmdn[0], z, &SEML);

                //Updating the last row
                gsl_matrix_set(DF, ncs-1, 0, gsl_vector_get(K5, 0));
                gsl_matrix_set(DF, ncs-1, 1, gsl_vector_get(K5, 2));
                gsl_matrix_set(DF, ncs-1, 2, Hdot - Hdot0);


                //Free
                gsl_matrix_free(COORD_J_RCM);
                gsl_vector_free(K5);
                gsl_vector_free(dH);

            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        //Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Updating the initial time
        tmdn[0] += gsl_vector_get(DQv, 2);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates. Again, we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k+1);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+2);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs+1), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }


        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    // Energy at the origin
    //====================================================================================
    if(isFirst)
    {
        double z[6], H, H0;
        for(int i = 0; i < 6; i++) z[i] = 0.0;
        H0 = qbcp_H_complete(tmdn[0], z, coord_type, coord_type);
        for(int i = 0; i < 6; i++) z[i] = ymdn[i][0];
        H = qbcp_H_complete(tmdn[0], z, coord_type, coord_type);
        refSt.dH = H - H0;
    }

    //====================================================================================
    // 3. Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_new[nfv];
    nullvectorjac(nullvector_new, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 1
        //--------------------------------------------------------------------------------
        sign = nullvector_new[1] > 0? -1:+1;
    }
    else
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 1
        //--------------------------------------------------------------------------------
        sign = nullvector_new[1] > 0? -1:+1;
    }

    //Finally, we update nullvector with the right direction
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_new[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. FINALLY, we update orbit_EM.tf, even if it is not necessary
    //====================================================================================
    orbit_EM.setTf(tmdn[mgs]/SEML.us_em.ns);

    //====================================================================================
    // 5. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. PLANAR CASE.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        - The initial conditions z0 vary in the center-unstable manifold of EML2.
 *        - The final state zN vary in the center-stable manifold of SEMLi.
 *        - The time t0 alone is free to vary.
 *        - The null vector associated to the solution is computed.
 *
 *        For now, the computation is limited to coord_type == NCSEM. The general structure
 *        Of the code leaves room for an extension to other types of coordinates. To do so,
 *        One should adapt the part that computes dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k]).
 *
 *        A constraint x = xe at t = te is added
 **/
int msvftplan_dte(double** ymd, double* tmd, double** ymdn, double* tmdn, double* nullvector,
                  int nov, int mgs, int coord_type,  int isFirst,
                  Orbit& orbit_EM, Orbit& orbit_SEM, gnuplot_ctrl* h1, RefSt& refSt, int* niter)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status = 0;
    //Name of the routine
    string fname = "msvftplan_dH";

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE)
    {
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    int fwrk = default_framework(coord_type);
    if(fwrk == FTC_FAILURE)
    {
        cerr << fname << ". The selection of fwrk failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Selection of the vector field vf
    //------------------------------------------------------------------------------------
    vfptr vf  = ftc_select_vf(dcs, 6);

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // Other initialization
    //------------------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC = 0.0;
    //Current state along the trajectory
    double** ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double* tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[nov], ye[nov];


    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 4*mgs + 2;    //free variables
    int ncs = 4*mgs + 1;    //constraints

    // Correction vector at patch points
    gsl_vector* DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector* Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix** Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix* DF   = gsl_matrix_calloc(ncs, nfv);

    //Identity matrix eye(6)
    gsl_matrix* Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector* Kf  = gsl_vector_calloc(6);
    gsl_vector* K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix* Phi0 = gsl_matrix_calloc(6,5);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix* PhiN = gsl_matrix_calloc(6,5);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < nov; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    int ode78coll = 0, tempcoll = 0;
    while(iter <  itermax)
    {
        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        ode78coll = 0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            tempcoll = 0;
            for(int i = 0; i < nov; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &tempcoll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Collisionner. If a collision occured, we save it in ode78coll
            //----------------------------------------------------------------------------
            if(tempcoll && !ode78coll) ode78coll = tempcoll;

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < nov; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update Phi0
            //----------------------------------------------------------------------------
            if(k == 0)
            {
                //Phi0 = Ji[0] x COORD_J_RCM(orbit_EM.si, t0)
                ftc_compute_phi0(Phi0, Ji[0], orbit_EM, tmdn[0]/SEML.us_em.ns, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }

            //----------------------------------------------------------------------------
            // Update PhiN
            //----------------------------------------------------------------------------
            if(k == mgs-1)
            {
                //PhiN = COORD_J_RCM(orbit_SEM.si, tf), in SEM units, in R(6,5)
                orbit_SEM.getInvman()->evalDRCMtoCOORD(orbit_SEM.getSi(), tmdn[mgs], PhiN, OFTS_ORDER, OFS_ORDER, coord_type);

                if(refSt.isDebug)
                {
                    cout << fname << ". PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //Final column: x = xe at k = refSt.pkpos, if pkpos has been computed (i.e. > 0)
            if(refSt.pkpos > 0)
            {
                if(k == refSt.pkpos)
                {
                    //gsl_vector_set(Fv, ncs-1, ye[0] - refSt.xps);
                    //cout << "gsl_vector_set(Fv, ncs-1, ye[0] - refSt.xps) = " << gsl_vector_get(Fv, ncs-1) << endl;
                    gsl_vector_set(Fv, ncs-1, ymdn[0][k] - refSt.xps);
                    //cout << "gsl_vector_set(Fv, ncs-1, ymdn[0][k] - refSt.xps) = " << gsl_vector_get(Fv, ncs-1) << endl;
                }
            }
            else gsl_vector_set(Fv, ncs-1, 0.0); //constraint automatically satisfied

            //----------------------------------------------------------------------------
            // Update the derivatives wrt to time
            // Special case of the first point
            //----------------------------------------------------------------------------
            if(k == 0) t0jac(K4, tmdn, ymdn, orbit_EM, orbit_SEM, odeParams, vf, Ji[0], Kf);

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            jacftplan_dt0(k, DF, refSt.sidim, mgs, Ji, Phi0, PhiN, K4, Id);

            //Final column: x = xe at k = refSt.pkpos
            if(refSt.pkpos > 0)
            {
                //Identity: x[k] = xe
                if(k == refSt.pkpos)
                {
                    gsl_matrix_set(DF, ncs-1, 4*k-1,  gsl_matrix_get(Id, 0, 0));
                }

                //STM: Phi_x(z[k]) = xe. Not necessary, because we ensure the continuity otherwise
                if(k == refSt.pkpos-1)
                {
                    //                    gsl_matrix_set(DF, ncs-1, 4*k-1,  gsl_matrix_get(Ji[k], 0, 0));
                    //                    gsl_matrix_set(DF, ncs-1, 4*k-0,  gsl_matrix_get(Ji[k], 0, 1));
                    //                    gsl_matrix_set(DF, ncs-1, 4*k+1,  gsl_matrix_get(Ji[k], 0, 3));
                    //                    gsl_matrix_set(DF, ncs-1, 4*k+2,  gsl_matrix_get(Ji[k], 0, 4));
                }
            }
            else
            {
                //No constraint, some ones on the row to avoid non-inversible matrix
                gsl_matrix_set(DF, ncs-1, 0,  1.0);
                gsl_matrix_set(DF, ncs-1, 1,  1.0);
            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        //Norm
        normC  = gsl_blas_dnrm2(Fv);

        //Display current status
        cout << fname << ". nerror = " << normC << endl;

        //Check that all points are under a given threshold
        if(normC < refSt.inner_prec)
        {
            cout << fname << ". Desired precision was reached. " << endl; cout << "nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //================================================================================
        // Update the free variables
        //================================================================================
        //--------------------------------------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //--------------------------------------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.addSi(gsl_vector_get(DQv, 0), 0);
        orbit_EM.addSi(gsl_vector_get(DQv, 1), 2);

        //Updating the initial time
        tmdn[0] += gsl_vector_get(DQv, 2);

        cout << "gsl_vector_get(DQv, 1) = " << gsl_vector_get(DQv, 1) << endl;
        cout << "gsl_vector_get(DQv, 2) = " << gsl_vector_get(DQv, 2) << endl;
        cout << "gsl_vector_get(DQv, 3) = " << gsl_vector_get(DQv, 3) << endl;

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), tmdn[0]/SEML.us_em.ns);

        //To CM_EM_NCSEM coordinates. Again, we need to normalize the time
        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //--------------------------------------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //--------------------------------------------------------------------------------
        for(int k = 1; k < mgs; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k+1);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+2);
        }

        //--------------------------------------------------------------------------------
        //Last 4 correction variables is orbit.si
        //--------------------------------------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-1), 0);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs-0), 2);
        orbit_SEM.addSi(gsl_vector_get(DQv, 4*mgs+1), 4);

        //Updating in CM_SEM_NCSEM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), tmdn[mgs]);

        //Updating in CM_SEM_NCSEM coordinates
        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
        qbcp_coc(tmdn[mgs], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][mgs] = ye[i];

        //================================================================================
        // Norm check: we need to know if we are in the DPC of the semi-analytical tools
        //================================================================================
        // First check at EML2
        si_norm_EM  = ENorm(orbit_EM.getSi(), 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cerr << fname << ". si_norm_EM has reached its limits: " << endl;
            cout << " si_norm_EM = "       << si_norm_EM << " >";
            cout << " SI_NORM_EM_MAX = "  << SI_NORM_EM_MAX << endl;
            return REF_EOUTOFDPC;
        }

        // Second check at SEMLi
        si_norm_SEM = ENorm(orbit_SEM.getSi(), 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cerr << fname << ". si_norm_SEM has reached its limits: " << endl;
            cout << " si_norm_SEM = "      << si_norm_SEM << " >";
            cout << " SI_NORM_SEM_MAX = "  << SI_NORM_SEM_MAX << endl;
            return REF_EOUTOFDPC;
        }


        //--------------------------------------------------------------------------------
        // Norm display
        //--------------------------------------------------------------------------------
        if(refSt.isDebug)
        {
            cout << fname << ". si_norm_EM = "   << si_norm_EM << endl;
            cout << fname << ". si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }

    //====================================================================================
    //Collision check: just a warning (for now). In the long run we need to give it to
    // the upper level!
    //====================================================================================
    if(ode78coll) cout << fname << ". A collision has occured with " << ode78coll << endl;

    //------------------------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_xyz(refSt.isPlotted, h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "2", "2", 4);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][mgs], &ymdn[1][mgs],  &ymdn[2][mgs], 1, (char*)"", "points", "2", "2", 0);
    gnuplot_plot_xyz(refSt.isPlotted, h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);

    //====================================================================================
    // 3. Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    double nullvector_new[nfv];
    nullvectorjac(nullvector_new, DF, ncs, nfv);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 2
        //--------------------------------------------------------------------------------
        sign = nullvector_new[2] < 0? -1:+1;
    }
    else
    {
        //--------------------------------------------------------------------------------
        // we just want to increase the first time, at position 2
        //--------------------------------------------------------------------------------
        sign = nullvector_new[2] < 0? -1:+1;
    }

    //Finally, we update nullvector with the right direction
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector_new[i];

    //------------------------------------------------------------------------------------
    //Number of iterations
    //------------------------------------------------------------------------------------
    *niter = iter;

    //------------------------------------------------------------------------------------
    //Last error
    //------------------------------------------------------------------------------------
    refSt.last_error = normC;

    //====================================================================================
    // 4. FINALLY, we update orbit_EM.tf, even if it is not necessary
    //====================================================================================
    orbit_EM.setTf(tmdn[mgs]/SEML.us_em.ns);

    //====================================================================================
    // 5. Free
    //====================================================================================
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_matrix_free(DF);
    gsl_matrix_free(Id);
    gsl_matrix_free(Phi0);
    gsl_matrix_free(PhiN);

    return GSL_SUCCESS;
}

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: JACOBIAN MATRICES
//
//========================================================================================
/**
 *  \brief Update of the jacobian matrix at step k for the Newton method in the routine msvftplan
 **/
int jacvftplan(int k, gsl_matrix* DF, int dim, int mgs, gsl_matrix **Ji, gsl_matrix *Phi0, gsl_matrix *PhiN, gsl_vector *K4, gsl_matrix *Id)
{
    //------------------------------------------------------------------------------------
    // Update DF
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        //--------------------------------------------------------------------------------
        //DF/DS
        //--------------------------------------------------------------------------------
        if(k == 0)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i, 0, gsl_matrix_get(Phi0, i, dim));
            }
            else
            {
                gsl_matrix_set(DF, i, 0, gsl_matrix_get(Phi0, i+1, dim));
            }

            for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
        }
        else if(k == mgs-1)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i+1, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i+1, 4));
            }
        }
        else
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i+1, 4));
            }
        }


        //--------------------------------------------------------------------------------
        //DF/DT
        //--------------------------------------------------------------------------------
        if(k == 0)
        {
            //--------------------------
            //dF[0]/dt[0]
            //--------------------------
            if(i < 2) gsl_matrix_set(DF, i, 1, gsl_vector_get(K4, i));
            else gsl_matrix_set(DF, i, 1, gsl_vector_get(K4, i+1));
        }
        else if(k == mgs-1)
        {
            //NOTHING IS DONE FOR NOW
        }
        else
        {
            //NOTHING IS DONE FOR NOW
        }

    }

    return GSL_SUCCESS;
}

/**
 *  \brief Update of the jacobian matrix at step k for the Newton method in the routine msftplan
 **/
int jacftplan(int k, gsl_matrix* DF, int mgs, gsl_matrix** Ji, gsl_matrix* Phi0, gsl_matrix* PhiN, gsl_matrix *Id)
{
    //------------------------------------------------------------------------------------
    // Update DF
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        //--------------------------------------------------------------------------------
        //DF/DS
        //--------------------------------------------------------------------------------
        if(k == 0)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));
            }
            else
            {
                gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i+1, 0));
                gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i+1, 2));
            }

            for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j+2, -gsl_matrix_get(Id, i, j));
        }
        else if(k == mgs-1)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-6,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  -gsl_matrix_get(PhiN, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i+1, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i+1, 4));
            }
        }
        else
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-2,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+2,  -gsl_matrix_get(Id, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i+1, 4));
            }

        }
    }

    return(GSL_SUCCESS);
}

/**
 *  \brief Update of the jacobian matrix at step k for the Newton method in the routine vith variable t0
 **/
int jacftplan_dt0(int k, gsl_matrix* DF, int dim, int mgs, gsl_matrix **Ji, gsl_matrix *Phi0, gsl_matrix *PhiN, gsl_vector *K4, gsl_matrix *Id)
{
    //------------------------------------------------------------------------------------
    // Update DF
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        //--------------------------------------------------------------------------------
        //DF/DS
        //--------------------------------------------------------------------------------
        if(k == 0)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i, 0));
                gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i, 2));
            }
            else
            {
                gsl_matrix_set(DF, i, 0,    gsl_matrix_get(Phi0, i+1, 0));
                gsl_matrix_set(DF, i, 1,    gsl_matrix_get(Phi0, i+1, 2));
            }

            //THERE IS A HIATUS HERE: it is the column of t0 (see below)

            for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j+3, -gsl_matrix_get(Id, i, j));
        }
        else if(k == mgs-1)
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs+1,  -gsl_matrix_get(PhiN, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*mgs-5,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-4,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-3,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-2,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*mgs-1,  -gsl_matrix_get(PhiN, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*mgs-0,  -gsl_matrix_get(PhiN, i+1, 2));
                gsl_matrix_set(DF, i + 4*k, 4*mgs+1,  -gsl_matrix_get(PhiN, i+1, 4));
            }
        }
        else
        {
            if(i < 2)
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+2,  gsl_matrix_get(Ji[k], i, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+6,  -gsl_matrix_get(Id, i, 4));
            }
            else
            {
                gsl_matrix_set(DF, i + 4*k, 4*k-1,  gsl_matrix_get(Ji[k], i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k-0,  gsl_matrix_get(Ji[k], i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+1,  gsl_matrix_get(Ji[k], i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+2,  gsl_matrix_get(Ji[k], i+1, 4));

                gsl_matrix_set(DF, i + 4*k, 4*k+3,  -gsl_matrix_get(Id, i+1, 0));
                gsl_matrix_set(DF, i + 4*k, 4*k+4,  -gsl_matrix_get(Id, i+1, 1));
                gsl_matrix_set(DF, i + 4*k, 4*k+5,  -gsl_matrix_get(Id, i+1, 3));
                gsl_matrix_set(DF, i + 4*k, 4*k+6,  -gsl_matrix_get(Id, i+1, 4));
            }

        }

        //--------------------------------------------------------------------------------
        //DF/DT
        //--------------------------------------------------------------------------------
        if(k == 0)
        {
            //----------------------------------------------------------------------------
            //dF[0]/dt[0]
            //----------------------------------------------------------------------------
            if(i < 2) gsl_matrix_set(DF, i, 2, gsl_vector_get(K4, i));
            else gsl_matrix_set(DF, i, 2, gsl_vector_get(K4, i+1));
        }
        else if(k == mgs-1)
        {
            //NOTHING IS DONE FOR NOW
        }
        else
        {
            //NOTHING IS DONE FOR NOW
        }

    }

    return(GSL_SUCCESS);
}

/**
 *      \brief Computes the nullvector of a given Jacobian DF, of size ncs x nfv
 **/
int nullvectorjac(double *nullvector, gsl_matrix* DF, int ncs, int nfv)
{
    //------------------------------------------------------------------------------------
    // Init the GSL objects
    //------------------------------------------------------------------------------------
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);

    //------------------------------------------------------------------------------------
    // Compute the nullvector
    //------------------------------------------------------------------------------------
    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);

    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //Null vector is the last column of Q, normalized
    gsl_vector_view nvview = gsl_matrix_subcolumn(Q, nfv-1, 0, nfv);
    for(int i = 0; i < nfv; i++) nullvector[i] = gsl_vector_get(&nvview.vector, i)/gsl_blas_dnrm2(&nvview.vector);

    //------------------------------------------------------------------------------------
    // Free
    //------------------------------------------------------------------------------------
    gsl_vector_free(work);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_matrix_free(DFT);

    return(GSL_SUCCESS);
}

/**
 *  \brief Dependence of the final point zf with respect to the initial time t0,
 *         in NCSEM coordinates, when the initial conditions are taken
 *         in a FT series in RCM (EM) coordinates.
 **/
int t0jac(gsl_vector *Kout, double* tmdn, double **ymdn, Orbit &orbit_EM, Orbit &orbit_SEM,
          OdeParams &odeParams, vfptr vf, gsl_matrix* Phi0, gsl_vector *Ktemp)
{
    //====================================================================================
    // dzf/dt0 = + Phi(z0, t0, tf)*(-f[z0, t0) + Qt0)
    // where Qt0 = dCMU_EM_NC/dt0
    //====================================================================================
    double yv[6], f[6];

    //------------------------------------------------------------------------------------
    // Computing Qt0 = dCMU_EM_NC/dt0; i.e. the dependency of the semi-numerical
    // approximation of the manifold in orbit_EM BUT in NCSEM coordinates.
    // Hence, a big change of coordinates RMC (EM) -> NCEM -> NCSEM has to be done
    // in evaldotRCMEMtoNCSEM. Note that the initial time is set in EM units right below
    //------------------------------------------------------------------------------------
    gsl_vector* dzNCSEM = gsl_vector_alloc(6);
    //dCMU_EM_NC/dt[k] in NCSEM coordinates
    orbit_EM.getInvman()->evaldotRCMEMtoNCSEM(orbit_EM.getSi(),
    tmdn[0]/SEML.us_em.ns, dzNCSEM, OFTS_ORDER, OFS_ORDER, *orbit_SEM.getInvman());

    //------------------------------------------------------------------------------------
    //dzf/dt0 = + Phi(z0, t0, tf)*(-f[z0, t0) + Qt0)
    //------------------------------------------------------------------------------------
    //Computing f[z0, t0)
    for(int i = 0; i < 6; i++) yv[i] = ymdn[i][0];
    vf(tmdn[0], yv, f, &odeParams);

    //Ktemp = -f[z0, t0) + Qt0
    for(int i = 0; i < 6; i++)
    {
        gsl_vector_set(Ktemp, i, -f[i] + gsl_vector_get(dzNCSEM, i));
        //if(f[i] != 0.0) cout << "f[i]/dzNCSEM[i] = " << f[i]/gsl_vector_get(dzNCSEM, i) << endl;
    }

    //Kout = dF[k]/dt[k] = +Phi0*Ktemp
    gsl_blas_dgemv(CblasNoTrans, 1.0, Phi0, Ktemp, 0.0, Kout);

    return GSL_SUCCESS;
}

//========================================================================================
//
//          DIFFCORR CUSTOM: CMU to CMS: UPDATE FREE VARIABLES with NEW IMPLEMENTATION
//
//========================================================================================
int ufvarft3d(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //------------------------------------------------------------------------------------
    //Here, we simply update *ds = ds0, no additionnal constraint.
    //------------------------------------------------------------------------------------
    *ds = ds0;

    //------------------------------------------------------------------------------------
    //Updating CM_EM
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    for(int i = 0; i < 4; i++) orbit_EM.addSi(*ds*nullvector[i], i);
    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);
    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

    //------------------------------------------------------------------------------------
    //The middle (patch points) is classical cartesian coordinates at patch points
    //------------------------------------------------------------------------------------
    for(int k = 1; k < mgs; k++)
    {
        for(int i = 0; i < 6; i++) y_traj_n[i][k] += *ds*nullvector[i+6*k-2];
    }

    //------------------------------------------------------------------------------------
    //Updating CM_SEM
    //------------------------------------------------------------------------------------
    //Updating CM_SEM_RCM coordinates
    for(int i = 0; i < 5; i++) orbit_SEM.addSi(*ds*nullvector[i + 6*mgs-2], i);

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];


    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}

int ufvarftmixed(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //------------------------------------------------------------------------------------
    //Here, we simply update *ds = ds0, no additionnal constraint.
    //------------------------------------------------------------------------------------
    *ds = ds0;

    //------------------------------------------------------------------------------------
    //Updating CM_EM
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);
    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

    //------------------------------------------------------------------------------------
    //The middle (patch points) is classical cartesian coordinates at patch points
    //------------------------------------------------------------------------------------
    for(int k = 1; k < mgs; k++)
    {
        for(int i = 0; i < 6; i++) y_traj_n[i][k] += *ds*nullvector[i+6*k-4];
    }

    //------------------------------------------------------------------------------------
    //Updating CM_SEM
    //------------------------------------------------------------------------------------
    //Updating CM_SEM_RCM coordinates
    for(int i = 0; i < 5; i++) orbit_SEM.addSi(*ds*nullvector[i + 6*mgs-4], i);

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];


    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}


int ufvarvt3d(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //------------------------------------------------------------------------------------
    //Here, we simply update *ds = ds0, no additionnal constraint.
    //------------------------------------------------------------------------------------
    // The value of ds depends on the type of termination condition for the continuation
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        // The objective of this continuation is to bring orbit_SEM.si[4] to 0.0.
        // Prior to updating, we check that the continuation is not going "to far"
        // (orbit_SEM.si[4] may change sign).
        //--------------------------------------------------------------------------------
        double dkn = orbit_SEM.getSi()[4] + ds0*nullvector[7*mgs+1];
        if(dkn * orbit_SEM.getSi()[4] < 0) //if there is a change of sign, we reduce the stepsize
        {
            *ds = -orbit_SEM.getSi()[4]/nullvector[7*mgs+1];
        }
        else *ds = ds0;
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        *ds = ds0;
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //Updating CM_EM
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    for(int i = 0; i < 4; i++) orbit_EM.addSi(*ds*nullvector[i], i);
    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);
    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

    //------------------------------------------------------------------------------------
    //The middle (patch points) is classical cartesian coordinates at patch points
    //------------------------------------------------------------------------------------
    for(int k = 1; k < mgs; k++)
    {
        for(int i = 0; i < 6; i++) y_traj_n[i][k] += *ds*nullvector[i + 7*k-3];
        t_traj_n[k] += *ds*nullvector[7*(k+1)-4];
    }
    //Last time:
    t_traj_n[mgs] += *ds*nullvector[7*mgs+2];


    //------------------------------------------------------------------------------------
    //Last 4 correction variables is orbit.si
    //------------------------------------------------------------------------------------
    //Updating CM_SEM_RCM coordinates
    for(int i = 0; i < 4; i++) orbit_SEM.addSi(*ds*nullvector[i + 7*mgs-3], i);
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(max(0.0, orbit_SEM.getSi()[4] + *ds*nullvector[7*mgs+1]), 4);
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(orbit_SEM.getSi()[4] + *ds*nullvector[7*mgs+1], 4);
        break;
    }
    }

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];


    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}

int ufvarvtmixed(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //------------------------------------------------------------------------------------
    //Here, we simply update *ds = ds0, no additionnal constraint.
    //------------------------------------------------------------------------------------
    // The value of ds depends on the type of termination condition for the continuation
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        // The objective of this continuation is to bring orbit_SEM.si[4] to 0.0.
        // Prior to updating, we check that the continuation is not going "to far"
        // (orbit_SEM.si[4] may change sign).
        //--------------------------------------------------------------------------------
        double dkn = orbit_SEM.getSi()[4] + ds0*nullvector[7*mgs-1];
        if(dkn * orbit_SEM.getSi()[4] < 0) //if there is a change of sign, we reduce the stepsize
        {
            *ds = -orbit_SEM.getSi()[4]/nullvector[7*mgs-1];
        }
        else *ds = ds0;
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        *ds = ds0;
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //Updating CM_EM
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);
    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);
    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

    //------------------------------------------------------------------------------------
    //The middle (patch points) is classical cartesian coordinates at patch points
    //------------------------------------------------------------------------------------
    for(int k = 1; k < mgs; k++)
    {
        for(int i = 0; i < 6; i++) y_traj_n[i][k] += *ds*nullvector[i + 7*k-5];
        t_traj_n[k] += *ds*nullvector[7*(k+1)-6];
    }
    //Last time:
    t_traj_n[mgs] += *ds*nullvector[7*mgs];


    //------------------------------------------------------------------------------------
    //Last 4 correction variables is orbit.si
    //------------------------------------------------------------------------------------
    //Updating CM_SEM_RCM coordinates
    for(int i = 0; i < 4; i++) orbit_SEM.addSi(*ds*nullvector[i + 7*mgs-5], i);
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(max(0.0, orbit_SEM.getSi()[4] + *ds*nullvector[7*mgs-1]), 4);
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(orbit_SEM.getSi()[4] + *ds*nullvector[7*mgs-1], 4);
        break;
    }
    }

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];


    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}


int ufvarftplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);


    //------------------------------------------------------------------------------------
    //Here, we simply update *ds = ds0, no additionnal constraint.
    //------------------------------------------------------------------------------------
    *ds = ds0;


    //------------------------------------------------------------------------------------
    // Uncomment for pure "blind" t0 = t0 + dt0 continuation
    //------------------------------------------------------------------------------------
    //        t_traj_n[0] -= 1e-4**ds;
    //
    //        //Updating CM_EM_NCEM coordinates
    //        orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);
    //
    //        //To CM_EM_NCSEM coordinates
    //        for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    //        qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    //        for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];
    //
    //        //Updating CM_SEM_NCSEM coordinates
    //        orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);
    //
    //        //Updating in CM_SEM_NCSEM coordinates
    //        for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    //        qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    //        for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];
    //
    //        return GSL_SUCCESS;


    //====================================================================================
    //Updating the free variables
    //====================================================================================
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);

    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


    //The middle (patch points) is classical cartesian coordinates at patch points
    for(int k = 1; k < mgs; k++)
    {
        y_traj_n[0][k] += *ds*nullvector[0 + 4*k-2];
        y_traj_n[1][k] += *ds*nullvector[1 + 4*k-2];
        y_traj_n[3][k] += *ds*nullvector[2 + 4*k-2];
        y_traj_n[4][k] += *ds*nullvector[3 + 4*k-2];
    }

    //Last 3 correction variables is orbit.si
    //Updating CM_SEM_RCM coordinates
    orbit_SEM.addSi(*ds*nullvector[4*mgs-2], 0);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-1], 2);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-0], 4);

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];


    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}


int ufvarvtplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //====================================================================================
    //Updating the free variables
    //====================================================================================
    // The value of ds depends on the type of termination condition for the continuation
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        // The objective of this continuation is to bring orbit_SEM.si[4] to 0.0.
        // Prior to updating, we check that the continuation is not going "to far"
        // (orbit_SEM.si[4] may change sign).
        //--------------------------------------------------------------------------------
        double dkn = orbit_SEM.getSi()[4] + ds0*nullvector[5*mgs-1];
        if(dkn * orbit_SEM.getSi()[4] < 0) //if there is a change of sign, we reduce the stepsize
        {
            *ds = -orbit_SEM.getSi()[4]/nullvector[5*mgs-1];
        }
        else *ds = ds0;
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        *ds = ds0;
        break;
    }
    }

    //------------------------------------------------------------------------------------
    // Then we can go on
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);

    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


    //The middle (patch points) is classical cartesian coordinates at patch points
    for(int k = 1; k < mgs; k++)
    {
        y_traj_n[0][k] += *ds*nullvector[0 + 5*k-3];
        y_traj_n[1][k] += *ds*nullvector[1 + 5*k-3];
        y_traj_n[3][k] += *ds*nullvector[2 + 5*k-3];
        y_traj_n[4][k] += *ds*nullvector[3 + 5*k-3];
        t_traj_n[k]    += *ds*nullvector[5*k+1];
    }

    //Last time:
    t_traj_n[mgs] += *ds*nullvector[5*mgs];

    //Last 4 correction variables is orbit.si
    //Updating CM_SEM_RCM coordinates
    orbit_SEM.addSi(*ds*nullvector[5*mgs-3], 0);
    orbit_SEM.addSi(*ds*nullvector[5*mgs-2], 2);
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(max(0.0, orbit_SEM.getSi()[4] + *ds*nullvector[5*mgs-1]), 4);
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(orbit_SEM.getSi()[4] + *ds*nullvector[5*mgs-1], 4);
        break;
    }
    }

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];

    //====================================================================================
    //Free variables
    //====================================================================================
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}

int ufvarvltplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //====================================================================================
    //Updating the free variables
    //====================================================================================
    // The value of ds depends on the type of termination condition for the continuation
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        // The objective of this continuation is to bring orbit_SEM.si[4] to 0.0.
        // Prior to updating, we check that the continuation is not going "to far"
        // (orbit_SEM.si[4] may change sign).
        //--------------------------------------------------------------------------------
        double dkn = orbit_SEM.getSi()[4] + ds0*nullvector[4*mgs];
        if(dkn * orbit_SEM.getSi()[4] < 0) //if there is a change of sign, we reduce the stepsize
        {
            *ds = -orbit_SEM.getSi()[4]/nullvector[4*mgs];
        }
        else *ds = ds0;
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        *ds = ds0;
        break;
    }
    }

    //------------------------------------------------------------------------------------
    // Then we can go on
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);

    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


    //The middle (patch points) is classical cartesian coordinates at patch points
    for(int k = 1; k < mgs; k++)
    {
        y_traj_n[0][k] += *ds*nullvector[4*k-2];
        y_traj_n[1][k] += *ds*nullvector[4*k-1];
        y_traj_n[3][k] += *ds*nullvector[4*k-0];
        y_traj_n[4][k] += *ds*nullvector[4*k+1];
    }

    //Last time:
    t_traj_n[mgs] += *ds*nullvector[4*mgs+1];

    //Last 4 correction variables is orbit.si
    //Updating CM_SEM_RCM coordinates
    orbit_SEM.addSi(*ds*nullvector[4*mgs-2], 0);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-1], 2);
    switch(refSt.termination)
    {
    case REF_COND_S5:
    {
        //--------------------------------------------------------------------------------
        // First type of condition: stop when we are close enough to the
        // center manifold (unstable component is small enough)
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(max(0.0, orbit_SEM.getSi()[4] + *ds*nullvector[4*mgs]), 4);
        break;
    }
    case REF_COND_T:
    {
        //--------------------------------------------------------------------------------
        // Another possible condition: enough turns around SEMLi. So, here,
        // we just want to increase the last time, at position nfv-1
        //--------------------------------------------------------------------------------
        orbit_SEM.addSi(orbit_SEM.getSi()[4] + *ds*nullvector[4*mgs], 4);
        break;
    }
    }

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];

    //====================================================================================
    //Free variables
    //====================================================================================
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}


int ufvarvftplan(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //====================================================================================
    //Updating the free variables
    //====================================================================================
    *ds = ds0;

    //------------------------------------------------------------------------------------
    // Then we can go on
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], refSt.sidim);
    //First time:
    t_traj_n[0] += *ds*nullvector[1];

    cout << "nullvector[1] = " << nullvector[1] << endl;

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);

    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


    //The middle (patch points) is classical cartesian coordinates at patch points
    for(int k = 1; k < mgs; k++)
    {
        y_traj_n[0][k] += *ds*nullvector[4*k-2];
        y_traj_n[1][k] += *ds*nullvector[4*k-1];
        y_traj_n[3][k] += *ds*nullvector[4*k-0];
        y_traj_n[4][k] += *ds*nullvector[4*k+1];
    }

    //Last 4 correction variables is orbit.si
    //Updating CM_SEM_RCM coordinates
    orbit_SEM.addSi(*ds*nullvector[4*mgs-2], 0);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-1], 2);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-0], 4);

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];

    //====================================================================================
    //Free variables
    //====================================================================================
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}

int ufvarvftplan_dH(double **y_traj_n, double *t_traj_n, double *ds, double ds0,
                double *nullvector,
                Orbit &orbit_EM, Orbit &orbit_SEM,
                int mgs, int coord_type,  RefSt &refSt)
{
    //------------------------------------------------------------------------------------
    //Temp variables
    //------------------------------------------------------------------------------------
    double *yv = dvector(0, 5);
    double *ye = dvector(0, 5);

    //====================================================================================
    //Updating the free variables
    //====================================================================================
    *ds = ds0;

    //------------------------------------------------------------------------------------
    // Then we can go on
    //------------------------------------------------------------------------------------
    //Updating CM_EM_RCM coordinates
    orbit_EM.addSi(*ds*nullvector[0], 0);
    orbit_EM.addSi(*ds*nullvector[1], 2);
    //First time:
    t_traj_n[0] += *ds*nullvector[2];

    //Updating CM_EM_NCEM coordinates
    orbit_EM.update_ic(orbit_EM.getSi(), t_traj_n[0]/SEML.us_em.ns);

    //To CM_EM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_EM.getZ0()[i];
    qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


    //The middle (patch points) is classical cartesian coordinates at patch points
    for(int k = 1; k < mgs; k++)
    {
        y_traj_n[0][k] += *ds*nullvector[4*k-1];
        y_traj_n[1][k] += *ds*nullvector[4*k-0];
        y_traj_n[3][k] += *ds*nullvector[4*k+1];
        y_traj_n[4][k] += *ds*nullvector[4*k+2];
    }

    //Last 4 correction variables is orbit.si
    //Updating CM_SEM_RCM coordinates
    orbit_SEM.addSi(*ds*nullvector[4*mgs-1], 0);
    orbit_SEM.addSi(*ds*nullvector[4*mgs-0], 2);
    orbit_SEM.addSi(*ds*nullvector[4*mgs+1], 4);

    //Updating CM_SEM_NCSEM coordinates
    orbit_SEM.update_ic(orbit_SEM.getSi(), t_traj_n[mgs]);

    //Updating in CM_SEM_NCSEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = orbit_SEM.getZ0()[i];
    qbcp_coc(t_traj_n[mgs], yv, ye, NCSEM, coord_type);
    for(int i = 0; i < 6; i++) y_traj_n[i][mgs] = ye[i];

    //====================================================================================
    //Free variables
    //====================================================================================
    free_dvector(yv, 0, 5);
    free_dvector(ye, 0, 5);
    return GSL_SUCCESS;
}
