#include "multshoot.h"

//========================================================================================
//
//          SUBROUTINES:
//
//========================================================================================
/**
 *  \brief Solve a definite-positive linear system:
 *         - Decompose a ncs x ncs matrix M that is definite-positive, using GSL routines.
 *         - Then inverse the system Fv = M*K3.
 **/
int ftc_inv_dfls(gsl_matrix* M, gsl_vector*Fv, gsl_vector* K3, int ncs)
{
    //Name of the routine
    string fname = "ftc_inv_dfls";

    //====================================================================================
    //We start by turning off the handler in the next chunk of code, so that we can test
    //if M is truly definite-positive.
    //We will restore it at the end of the decomposition + inversion process.
    //====================================================================================
    gsl_error_handler_t* error_handler = gsl_set_error_handler_off();

    //====================================================================================
    //Since M is in theory definite-positive,
    //a Cholesky decomposition can be used, instead of a LU decomposition.
    //====================================================================================
    int status = gsl_linalg_cholesky_decomp (M);

    //====================================================================================
    //If the decomposition went bad, status is different from GSL_SUCCESS. For example,
    //if M is not strictly definite-positive, the constant GSL_EDOM is returned by
    //gsl_linalg_cholesky_decomp. In any case, we try a LU decomposition.
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Cholesky decomposition failed. A LU decomposition is tested. ref_errno = " << gsl_strerror(status) << endl;
        int s;
        gsl_permutation* p = gsl_permutation_alloc (ncs);
        status = gsl_linalg_LU_decomp (M, p, &s);
        if(!status) status = gsl_linalg_LU_solve(M, p, Fv, K3);
        gsl_permutation_free (p);

        //--------------------------------------------------------------------------------
        //If the decomposition or solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". The LU decomposition failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return FTC_FAILURE;
        }

    }
    else  //if not, the decomposition went fine, and we can go on with the inversion.
    {
        status = gsl_linalg_cholesky_solve(M, Fv, K3);
        //--------------------------------------------------------------------------------
        //If the solving went bad, status is different from GSL_SUCCESS.
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". Cholesky solving failed. ref_errno = " << gsl_strerror(status) << endl;
            //Finally, we return an error to the user
            gsl_set_error_handler (error_handler);
            return FTC_FAILURE;
        }
    }

    //We check that the result is not undefined
    if(gsl_isnan(gsl_vector_max(K3)))
    {
        cerr << fname << ". The decomposition + solving failed, result contains NaN" << endl;
        //Finally, we return an error to the user
        gsl_set_error_handler (error_handler);
        return FTC_FAILURE;
    }


    //====================================================================================
    //We restore the error handler at the end
    //====================================================================================
    gsl_set_error_handler (error_handler);
    return FTC_SUCCESS;
}

/**
 *  \brief Computes the correction vector associated to the minimum norm solution.
 *         Given:
 *              - an ncs x 1   error vector Fv
 *              - an nfv x ncs Jacobian DF,
 *         This routine computes the correction vector associated to
 *         the minimum norm solution:
 *
 *              DQv = DF^T x (DF x DF^T) Fv.
 **/
int ftc_corrvec_mn(gsl_vector* DQv, gsl_vector *Fv, gsl_matrix* DF, int nfv, int ncs)
{
    //Name of the routine
    string fname = "ftc_corrvec_mn";

    //====================================================================================
    // Initialize the local GSL objects
    //====================================================================================
    gsl_matrix *M   = gsl_matrix_calloc(ncs, ncs);
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //====================================================================================
    // Update M = DF x DF^T
    //====================================================================================
    int status = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);
    if(status)
    {
        cerr << fname << ". The computation of M failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return FTC_FAILURE;
    }

    //====================================================================================
    // Update correction vector DQv = DF^T x (DF x DF^T) Fv
    //====================================================================================
    //Inverse Fv = M*K3 via Cholesky decomposition of M
    status = ftc_inv_dfls(M, Fv, K3, ncs);
    //We check that the inversing went well
    if(status)
    {
        cerr << fname << ". The decomposition + solving failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return FTC_FAILURE;
    }
    //Then, DQv = DF^T*K3
    status = gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);
    if(status)
    {
        cerr << fname << ". The computation of DQv failed. ref_errno = " << gsl_strerror(status) << endl;
        gsl_matrix_free(M);
        gsl_vector_free(K3);
        return FTC_FAILURE;
    }

    //====================================================================================
    // Kill the local GSL objects and return FTC_SUCCESS
    //====================================================================================
    gsl_matrix_free(M);
    gsl_vector_free(K3);
    return FTC_SUCCESS;
}


//========================================================================================
//
//          DIFFCORR BASED ON GOMEZ ET AL. 1998
//
//========================================================================================
/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Based on a recursive scheme in Gomez et al., "Quasihalo orbits associated with libration points", 1998, JAS.
 *        The Multiple shooting is denoted GMS for Gomez Multiple Shooting in the commentaries.
 **/
int multiple_shooting_gomez(double **ymd, double *tmd,
                            double **ymdn, double *tmdn, double **yma,
                            int N, int mgs,
                            int isPlotted, gnuplot_ctrl *h1)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, mgs);
    //Current time along the trajectory
    double *tm   = dvector(0, mgs);
    //Various temporary states and times
    double yv[N], ye[N];

    //------------------------------------------------------------------------------------
    //GSL matrices and vectors
    //------------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector **DQ  = gslc_vector_array_calloc(6, mgs+1);
    // Error vector at patch points
    gsl_vector **F   = gslc_vector_array_calloc(6, mgs);
    //Intermediate variables to compute DQ (see Gomez et al.)
    gsl_vector **X   = gslc_vector_array_calloc(6, mgs);
    gsl_vector **Y   = gslc_vector_array_calloc(6, mgs);
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    //Intermediate matrices to compute DQ (see Gomez et al.)
    gsl_matrix **Di  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix **Li  = gslc_matrix_array_calloc(6, 6, mgs);

    //Temporary variables
    gsl_matrix *K1 = gsl_matrix_calloc(6,6);
    gsl_matrix *K2 = gsl_matrix_calloc(6,6);
    gsl_matrix *K4 = gsl_matrix_calloc(6,6);
    gsl_vector_view K2v;
    gsl_vector_view K4v;

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn: ymdn = ymd
    //--------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //----------------------------------------------------------------------
    // Arrival state: yma(1, :) = ymd(1,:);
    // Note that the arrival state does not make much sense for now.
    // Its definition needs to be clarified, in the case of the GMS.
    //----------------------------------------------------------------------
    for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    int ode78coll;
    while(iter < 20)
    {
        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        normC = 0.0;
        for(int k = 0; k <= mgs-1; k++)
        {
            //------------------------------------------------------------------
            // Integration between tmdn[k] and tmdn[k+1]
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78(ym, tm, &ode78coll, tmdn[k], tmdn[k+1], yv, 42, 1, I_VNCSEM, VNCSEM, VNCSEM);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Arrival state (including STM) = final position is at the end of ym
            // Note that the arrival state does not make much sense for now.
            // Its definition needs to be clarified, in the case of the GMS.
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Di and Li
            //------------------------------------------------------------------
            if(k == 0)
            {
                //D[0] = Id
                gsl_matrix_set_identity (Di[k]);
                //D[0] = Id + Ji[0]*Ji[0]^T
                gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Ji[k] , Ji[k], 1.0, Di[k]);
            }
            else
            {
                //----------------------
                // Li[k] = -Ji[k]*inv(Di[k-1])
                // We first compute Li[k]^T = - inv(Di[k-1])*Ji[k]^T
                //----------------------
                // K1 = inv(Di[k-1]).
                // Since Di is symmetric positive,
                // a Cholesky decomposition can be used
                gsl_matrix_transpose_memcpy(K1, Di[k-1]);
                gsl_linalg_cholesky_decomp(K1);

                //K2 = -Ji[k]^T
                gsl_matrix_transpose_memcpy(K2, Ji[k]);
                gsl_matrix_scale(K2, -1.0);

                //K4 = Li[k]^T = - inv(Di[k-1])*Ji[k]^T
                for(int i = 0; i < 6; i++)
                {
                    K2v = gsl_matrix_column (K2, i);

                    if(gsl_blas_dnrm2(&K2v.vector) > 0)
                    {
                        K4v = gsl_matrix_column (K4, i);
                        gsl_linalg_cholesky_solve(K1, &K2v.vector, &K4v.vector);
                    }
                }

                //Li[k] = K4^T
                gsl_matrix_transpose_memcpy(Li[k], K4);

                //----------------------
                //Di[k] = Id + Ji[k]*Ji[k]^T + Ji[k]*Li[k]^T
                //----------------------
                //Di[k]  = Id
                gsl_matrix_set_identity (Di[k]);
                //Di[k]  = Id + Ji[k]*Ji[k]^T
                gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Ji[k] , Ji[k], 1.0, Di[k]);

                //Di[k] = Id + Ji[k]*Ji[k]^T + Ji[k]*Li[k]^T
                gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, Ji[k], Li[k], 1.0, Di[k]);
            }

            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(F[k], i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update the norm of the error
            //------------------------------------------------------------------
            normC  += gsl_blas_dnrm2(F[k]);
        }

        //----------------------------------------------------------------------
        // Build the intermediate variables X
        //----------------------------------------------------------------------
        for(int k = 0; k <= mgs-1; k++)
        {
            //X[k] = F[k]
            gsl_vector_memcpy(X[k], F[k]);

            if(k > 0)
            {
                //X[k] = F[k] - Li[k]*X[k-1]
                gsl_blas_dgemv(CblasNoTrans, -1.0, Li[k], X[k-1], 1.0, X[k] );
            }
        }

        //----------------------------------------------------------------------
        // Build the intermediate variables X
        //----------------------------------------------------------------------
        for(int k = mgs-1; k >= 0; k--)
        {
            //---------------------------------
            //Y[k] = inv(Di[k])*X[k]
            //---------------------------------
            //K1 = Di[k]
            gsl_matrix_memcpy(K1, Di[k]);
            //Y[k] = inv(Di[k])*X[k]
            gsl_linalg_cholesky_decomp(K1);
            gsl_linalg_cholesky_solve(K1, X[k], Y[k]);
            //Y[k] = inv(Di[k])*X[k] - L[k+1]^T*Y[k+1]
            if(k < mgs-1)
            {
                //Y[k] = inv(Di[k])*X[k] - L[k+1]^T*Y[k+1]
                gsl_blas_dgemv(CblasTrans, -1.0, Li[k+1], Y[k+1], 1.0, Y[k]);
            }

        }

        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_gomez. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        for(int k = 0; k <= mgs; k++)
        {
            if(k == 0)
            {
                //DQ[0] = -Ji[0]^T*Y[0]
                gsl_blas_dgemv(CblasTrans, -1.0, Ji[k], Y[k], 0.0, DQ[k]);
            }
            else if(k == mgs)
            {
                //DQ[mgs] = Y[mgs-1]
                gsl_vector_memcpy(DQ[k], Y[k-1]);
            }
            else
            {
                //DQ[k] = Y[k-1]- Ji[k]^T*Y[k]
                gsl_vector_memcpy(DQ[k], Y[k-1]);
                gsl_blas_dgemv(CblasTrans, -1.0, Ji[k], Y[k], 1.0, DQ[k]);
            }
        }


        //----------------------------------------------------------------------
        // Update the free variables
        //----------------------------------------------------------------------
        for(int k = 0; k <= mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQ[k], i);
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }


    //--------------------------------------------------------------------------
    // Free
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, mgs);
    free_dvector(tm, 0, mgs);
    gslc_vector_array_free(DQ, mgs+1);
    gslc_vector_array_free(F, mgs);
    gslc_vector_array_free(X, mgs);
    gslc_vector_array_free(Y, mgs);
    gslc_matrix_array_free(Ji, mgs);
    gslc_matrix_array_free(Di, mgs);
    gslc_matrix_array_free(Li, mgs);
    gsl_matrix_free(K1);
    gsl_matrix_free(K2);
    gsl_matrix_free(K4);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct(double **ymd, double *tmd,
                             double **ymdn, double *tmdn,
                             int N, int mgs, int coord_type,
                             int isPlotted, gnuplot_ctrl *h1, int isPar)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Name of the routine
    string fname = "multiple_shooting_direct";

    //Cumulated norm of the error
    double normC;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    if(dcs == FTC_FAILURE){
        cerr << fname << ". The selection of dcs failed." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Get the correct integration routine
    //------------------------------------------------------------------------------------
    //by default, to avoid gcc warning
    int (*ode78_int)(double**, double*, int*, double, double, const double*y, int, int, int, int, int) = ode78;

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*(mgs+1);  //free variables
    int ncs = 6*mgs;      //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    while(iter < itermax)
    {

        //--------------------------------------------------------------------------------
        //Trajectory plot
        //--------------------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_X(h1, ymdn, mgs+1, (char*)"", "points", "1", "4", 4);

        //--------------------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //--------------------------------------------------------------------------------
        #pragma omp parallel for if(isPar)
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Initialization for // computing
            //----------------------------------------------------------------------------
            int ode78coll;
            double yv[N], ye[N];
            //Current state along the trajectory
            double **ym  = dmatrix(0, 41, 0, mgs);
            //Current time along the trajectory
            double *tm   = dvector(0, mgs);

            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_int(ym, tm, &ode78coll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 6; j++)
                {
                    gsl_matrix_set(DF, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                    gsl_matrix_set(DF, i + 6*k, j + 6*(k+1), -gsl_matrix_get(Id, i, j));
                }
            }

            //----------------------------------------------------------------------------
            //Display
            //----------------------------------------------------------------------------
            //            #pragma omp critical
            //            {
            //                cout << fname << ".Step " << k << endl;
            //            }

            //----------------------------------------------------------------------------
            //Free
            //----------------------------------------------------------------------------
            free_dmatrix(ym, 0, 41, 0, mgs);
            free_dvector(tm, 0, mgs);
        }


        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        normC  = gsl_blas_dnrm2(Fv);

        //--------------------------------------------------------------------------------
        // Check that all points are under a given threshold
        //--------------------------------------------------------------------------------
        cout << fname << ". n° " << iter+1 << "/" << itermax << ". nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            cout << fname << ". Desired precision was reached. break. nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        int status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }

        //        if(iter == 0)
        //        {
        //            for(int k = 0; k < DQv->size; k++) cout << gsl_vector_get(DQv, k) << endl;
        //        }

        //================================================================================
        // Update the free variables
        //================================================================================
        double factor;
        if(dcs == I_ECLI) factor = 1.0;
        else factor = 1.0;

        for(int k = 0; k <= mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += factor*gsl_vector_get(DQv, 6*k+i);
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }

    //====================================================================================
    // Free
    //====================================================================================
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return FTC_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        Note that, contrary to the Gomez et al. article, the times are free to vary,
 *        leading to additionnal free variables and constraints.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd,
                                           double **ymdn, double *tmdn,
                                           int N, int mgs, int coord_type,
                                           double prec,
                                           int isPlotted, gnuplot_ctrl *h1, int isPar)
{
    //Name of the routine
    string fname = "multiple_shooting_direct_variable_time";

    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Cumulated norm of the error
    double normC;

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
    // Get the correct integration routine
    //------------------------------------------------------------------------------------
    //by default, to avoid gcc warning
    int (*ode78_int)(double**, double*, int*, double, double, const double*y, int, int, int, int, int) = ode78;

    //------------------------------------------------------------------------------------
    // Create a local OdeParams
    //------------------------------------------------------------------------------------
    OdeParams odeParams(&SEML, dcs);

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 7*(mgs+1);    //free variables
    int ncs = 6*mgs;        //constraints

    // Correction vector at patch points
    gsl_vector *DQv  = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector *Fv   = gsl_vector_calloc(ncs);
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(6*mgs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6, 6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    while(iter <  itermax)
    {
        //================================================================================
        //Trajectory plot
        //================================================================================
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], mgs+1, (char*)"", "points", "1", "4", 4);

        //================================================================================
        // Build the Jacobian and other useful matrices
        //================================================================================
        #pragma omp parallel for if(isPar)
        for(int k = 0; k <= mgs-1; k++)
        {
            //----------------------------------------------------------------------------
            // Initialization for // computing
            //----------------------------------------------------------------------------
            int ode78coll;
            double yv[N], ye[N], f[6];
            //Current state along the trajectory
            double **ym  = dmatrix(0, 41, 0, mgs);
            //Current time along the trajectory
            double *tm   = dvector(0, mgs);

            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_int(ym, tm, &ode78coll, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //----------------------------------------------------------------------------
            // Final position is at the end of ym
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update the error vector: case of the time variables
            //----------------------------------------------------------------------------
            //--------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //--------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &odeParams);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //--------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //--------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(tmdn[k+1], ye, f, &odeParams);


            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 6; j++)
                {
                    //--------------------------
                    //dF[k]/dQ[k]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, j + 7*k, gsl_matrix_get(Ji[k], i, j));

                    //--------------------------
                    //dF[k]/dQ[k+1]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, j + 7*(k+1), -gsl_matrix_get(Id, i, j));
                }

                //--------------------------
                //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                //--------------------------
                gsl_matrix_set(DF, i + 6*k, 7*(k+1)-1, gsl_vector_get(K4, i));

                //--------------------------
                //dF[k]/dt[k+1]
                //--------------------------
                gsl_matrix_set(DF, i + 6*k, 7*(k+2)-1, f[i]);
            }

            //----------------------------------------------------------------------------
            //Display
            //----------------------------------------------------------------------------
            //            #pragma omp critical
            //            {
            //                cout << fname << ".Step " << k << endl;
            //            }

            //----------------------------------------------------------------------------
            //Free
            //----------------------------------------------------------------------------
            free_dmatrix(ym, 0, 41, 0, mgs);
            free_dvector(tm, 0, mgs);
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
        if(normC < prec)
        {
            cout << fname << ". Desired precision was reached. break. nerror = " << normC << endl;
            break;
        }

        //================================================================================
        //Compute the correction vector
        //================================================================================
        int status = ftc_corrvec_mn(DQv, Fv, DF, nfv, ncs);
        if(status)
        {
            cerr << fname << ". The computation of the correction vector failed."  << endl;
            return FTC_FAILURE;
        }


        //================================================================================
        // Update the free variables
        //================================================================================
        for(int k = 0; k <= mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 7*k+i);
            tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-1);
        }

        //================================================================================
        // Update number of iterations
        //================================================================================
        iter++;
    }


    //====================================================================================
    // Free
    //====================================================================================
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 *        An additionnal free variable epsilon is added to the state, for continuation.
 **/
int multiple_shooting_direct_deps(double **ymd, double *tmd,
                                  double **ymdn, double *tmdn,
                                  double *nullvector, int isFirst,
                                  int N, int mgs, int coord_type,
                                  int isPlotted, gnuplot_ctrl *h1)
{
    //====================================================================================
    // 1. Initialization
    //====================================================================================
    //Status along the computation
    int status;
    //Name of the routine
    string fname = "multiple_shooting_direct_deps";

    //Cumulated norm of the error
    double normC;
    //Various temporary states and times
    double yv[N], ye[N];

    //------------------------------------------------------------------------------------
    //Driver
    //------------------------------------------------------------------------------------
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, ECISEM);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 48, qbcp_ecisem_cont_necijpl, &odeParams);

    //------------------------------------------------------------------------------------
    // GSL matrices and vectors
    //------------------------------------------------------------------------------------
    int nfv = 6*(mgs+1);  //free variables
    int ncs = 6*mgs;      //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, mgs);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    gsl_matrix *DF1  = gsl_matrix_calloc(ncs, nfv+1);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //------------------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //------------------------------------------------------------------------------------
    for(int k = 0; k <= mgs; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //====================================================================================
    // 2. Loop correction
    //====================================================================================
    //Maximum number of iterations is retrieved from config manager
    int itermax = Config::configManager().G_DC_ITERMAX();
    int iter = 0;
    double t = 0;
    while(iter < itermax)
    {
        //--------------------------------------------------------------------------------
        //Trajectory plot
        //--------------------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_X(h1, ymdn, mgs+1, (char*)"", "points", "1", "3", 2);
        //pressEnter(true);

        //--------------------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //--------------------------------------------------------------------------------
        for(int k = 0; k <= mgs-1; k++)
        {
            //cout << "k = " << k << endl;
            //----------------------------------------------------------------------------
            // Initialization of yv
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            // For safety, the identity matrix is set at the end of yv STM(t = 0)
            gslc_matrixToVector(yv, Id, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Integration
            //----------------------------------------------------------------------------
            t = tmdn[k];
            status = gsl_odeiv2_driver_apply (odestruct.d, &t , tmdn[k+1] , yv);

            //----------------------------------------------------------------------------
            // Final position is at the end of yv
            //----------------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = yv[i];

            //----------------------------------------------------------------------------
            // Update the Jacobian
            //----------------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //----------------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }

            //----------------------------------------------------------------------------
            // Update DF
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 6; j++)
                {
                    gsl_matrix_set(DF, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                    gsl_matrix_set(DF, i + 6*k, j + 6*(k+1), -gsl_matrix_get(Id, i, j));
                }

                //Last column contains the variations with respect to epsilon!
                //gsl_matrix_set(DF, i + 6*k, nfv-1, ye[42+i]);
            }

            //----------------------------------------------------------------------------
            // Update DF1
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 6; j++)
                {
                    gsl_matrix_set(DF1, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                    gsl_matrix_set(DF1, i + 6*k, j + 6*(k+1), -gsl_matrix_get(Id, i, j));
                }

                //Last column contains the variations with respect to epsilon!
                gsl_matrix_set(DF1, i + 6*k, nfv, ye[42+i]);
                //cout << "ye[42+i] = " << ye[42+i] << endl;
            }
        }

        //================================================================================
        //Termination condition: if the desired precision is met,
        //the process is terminated.
        //================================================================================
        normC  = gsl_blas_dnrm2(Fv);

        //--------------------------------------------------------------------------------
        // Check that all points are under a given threshold
        //--------------------------------------------------------------------------------
        cout << fname << ". n° " << iter+1 << "/" << itermax << ". nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            cout << fname << ". Desired precision was reached. break" << endl;//. nerror = " << normC << endl;
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
        //Update the state
        for(int k = 0; k <= mgs; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 6*k+i);
        }

        //Update epsilon at the last position? Useful?
        //SEML.epsilon += gsl_vector_get(DQv, nfv-1);
        //cout << fname << ". SEML.epsilon = " << SEML.epsilon << endl;

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }

    //====================================================================================
    //Last plot
    //====================================================================================
    if(isPlotted) gnuplot_plot_X(h1, ymdn, mgs+1, (char*)"", "points", "1", "3", 2);

    //====================================================================================
    //Compute the null vector: QR decomposition of DP^T
    //====================================================================================
    //Add one free variable: epsilon
    nfv = nfv+1;

    //QR elements
    gsl_vector* work  = gsl_vector_calloc(ncs);
    gsl_matrix* Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix* R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix* DFT   = gsl_matrix_calloc(nfv,ncs);

    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF1);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //------------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //------------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        //If it is the first use, we go in the direction so that depsilon > 0
        sign = gsl_matrix_get(Q, nfv-1, nfv-1) > 0? ti:-ti;
    }
    else
    {
        //If not, we follow the last direction of the null vector.
        for(int i = 0; i < nfv; i++)
        {
            nullvector[i] = gsl_matrix_get(Q, i, nfv-1);
            dotNV += nullvector[i]*nullvector[i];
        }
        sign = dotNV > 0? 1:-1;
        for(int i = 0; i < nfv; i++) nullvector[i] = sign*nullvector[i];
    }
    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);


    //cout << "Nullvector = " << endl;
    //vector_printf_prec(nullvector, nfv);


    //====================================================================================
    // Free
    //====================================================================================
    gslc_matrix_array_free(Ji , mgs);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(DF1);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);
    gsl_vector_free(work);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    gsl_matrix_free(DFT);

    return FTC_SUCCESS;
}


