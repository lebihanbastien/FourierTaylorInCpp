#include "lenconref.h"

//=======================================================================================================================================
//
//          DIFFCORR
//
//=======================================================================================================================================
/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_level_I(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                    int N, int man_grid_size,
                                    int isPlotted, int isTimeFixed, gnuplot_ctrl *h1)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N], te, yedot[6];
    //Error vector
    double bv[3];

    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_alloc(6,6);
    gsl_vector *P    = gsl_vector_alloc(3);
    gsl_matrix *Prod = gsl_matrix_calloc(3,3);
    gsl_vector *P1   = gsl_vector_calloc(3);
    gsl_vector *Pn;
    gsl_matrix *DP;

    if(isTimeFixed)
    {
        DP  = gsl_matrix_alloc(3,3);
        Pn = gsl_vector_calloc(3);
    }
    else
    {
        DP  = gsl_matrix_alloc(3,4);
        Pn = gsl_vector_calloc(4);
    }

    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (3);
    int s;



    //--------------------------------------------------------------------------
    // Vector that will contain all corrections to make at all patch points
    //--------------------------------------------------------------------------
    double **kv;
    if(isTimeFixed) kv = dmatrix(0, 2, 0, man_grid_size-1);
    else kv = dmatrix(0, 3, 0, man_grid_size-1);


    //--------------------------------------------------------------------------
    // Copy the departure state
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter < 50)
    {

        //----------------------------------------------------------------------
        // Norm check
        //----------------------------------------------------------------------
        normC = 0.0;

        //----------------------------------------------------------------------
        // Arrival state: yma(1, :) = ymd(1,:);
        //----------------------------------------------------------------------
        for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

        //----------------------------------------------------------------------
        // Level one
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++) //TBC: the LAST point
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, F_VNCSEM, VNCSEM, VNCSEM);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Arrival state (including STM)
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            // STM: STM = vectorToMatrix(ye, 6, 6, 6);
            gslc_vectorToMatrix(STM, ye, 6, 6, 6);

            // Derivatives at tf
            qbfbp_vfn_novar_xv(te, ye, yedot, &SEML);

            if(isTimeFixed)
                //DP = STM(1:3,4:6);
                for(int i = 0; i < 3; i++) for(int j = 0; j <3; j++) gsl_matrix_set(DP, i, j, gsl_matrix_get(STM, i, j+3));
            else
            {
                //DP = [STM(1:3,4:6) yedot(1:3)];
                for(int i = 0; i < 3; i++)
                {
                    for(int i = 0; i < 3; i++)
                    {
                        for(int j = 0; j <3; j++) gsl_matrix_set(DP, i, j, gsl_matrix_get(STM, i, j+3));
                        gsl_matrix_set(DP, i, 3, yedot[i]);
                    }
                }
            }

            //------------------------------------------------------------------
            // Update the error vector: P = [ye(0:2) - ymdn(k+1,0:2)']
            //------------------------------------------------------------------
            for(int i = 0; i < 3; i++)
            {
                gsl_vector_set(P, i, ye[i] - ymdn[i][k+1]);
                bv[i] = ye[i] - ymdn[i][k+1];
            }

            //------------------------------------------------------------------
            // Norm
            //------------------------------------------------------------------
            normC += ENorm(bv, 3);

            //------------------------------------------------------------------
            // Minimum norm solution
            //------------------------------------------------------------------
            if(isTimeFixed) //inverse a 3x3 system
            {
                //Inverse
                gsl_linalg_LU_decomp (DP, p , &s);
                //Pn = DP^{-1}*P
                gsl_linalg_LU_solve(DP, p, P, Pn);
                //Update kv
                for(int i = 0; i < 3; i++) kv[i][k] = gsl_vector_get(Pn, i);
            }
            else //minimum norm solution on a 3x4 system
            {
                //Compute Prod = DP*DP^T
                gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DP , DP, 0.0, Prod);
                //Inverse
                gsl_linalg_LU_decomp (Prod, p , &s);
                //P1 = Prod^{-1}*P
                gsl_linalg_LU_solve(Prod, p, P, P1);
                //Pn = DP^T*P1
                gsl_blas_dgemv (CblasTrans, 1.0, DP, P1, 0.0, Pn);
                //Update kv
                for(int i = 0; i < 4; i++) kv[i][k] = gsl_vector_get(Pn, i);
            }
        }

        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "normC = " << normC << endl;
        if(normC < 1e-9) break;


        //----------------------------------------------------------------------
        // Update the free variables
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //The state at position k
            for(int i = 0; i < 3; i++) ymdn[i+3][k] -= kv[i][k];
            //The final time, at position k+1
            if(!isTimeFixed) tmdn[k+1] -= kv[3][k];
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;

        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "lines", "1", "3", 5);
    }


    //--------------------------------------------------------------------------
    // Free
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    if(isTimeFixed) free_dmatrix(kv, 0, 2, 0, man_grid_size-1);
    else free_dmatrix(kv, 0, 3, 0, man_grid_size-1);

    gsl_matrix_free(STM);
    gsl_matrix_free(Prod);
    gsl_matrix_free(DP);

    gsl_vector_free(P);
    gsl_vector_free(P1);
    gsl_vector_free(Pn);

    gsl_permutation_free(p);

    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions.
 *        Based on a recursive scheme in Gomez et al., "Quasihalo orbits associated with libration points", 1998, JAS.
 *        The Multiple shooting is denoted GMS for Gomez Multiple Shooting in the commentaries.
 **/
int multiple_shooting_gomez(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                            int N, int man_grid_size,
                            int isPlotted, int isTimeFixed, gnuplot_ctrl *h1)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N];

    //-------------------------------------------------------------------------------
    //GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector **DQ  = gslc_vector_array_calloc(6, man_grid_size+1);
    // Error vector at patch points
    gsl_vector **F   = gslc_vector_array_calloc(6, man_grid_size);
    //Intermediate variables to compute DQ (see Gomez et al.)
    gsl_vector **X   = gslc_vector_array_calloc(6, man_grid_size);
    gsl_vector **Y   = gslc_vector_array_calloc(6, man_grid_size);
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    //Intermediate matrices to compute DQ (see Gomez et al.)
    gsl_matrix **Di  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix **Li  = gslc_matrix_array_calloc(6, 6, man_grid_size);

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
    for(int k = 0; k <= man_grid_size; k++)
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
    while(iter < 20)
    {
        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        normC = 0.0;
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration between tmdn[k] and tmdn[k+1]
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, F_VNCSEM, VNCSEM, VNCSEM);

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
        for(int k = 0; k <= man_grid_size-1; k++)
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
        for(int k = man_grid_size-1; k >= 0; k--)
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
            if(k < man_grid_size-1)
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
        for(int k = 0; k <= man_grid_size; k++)
        {
            if(k == 0)
            {
                //DQ[0] = -Ji[0]^T*Y[0]
                gsl_blas_dgemv(CblasTrans, -1.0, Ji[k], Y[k], 0.0, DQ[k]);
            }
            else if(k == man_grid_size)
            {
                //DQ[man_grid_size] = Y[man_grid_size-1]
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
        for(int k = 0; k <= man_grid_size; k++)
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
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_vector_array_free(DQ, man_grid_size+1);
    gslc_vector_array_free(F, man_grid_size);
    gslc_vector_array_free(X, man_grid_size);
    gslc_vector_array_free(Y, man_grid_size);
    gslc_matrix_array_free(Ji, man_grid_size);
    gslc_matrix_array_free(Di, man_grid_size);
    gslc_matrix_array_free(Li, man_grid_size);
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
int multiple_shooting_direct(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                             int N, int man_grid_size, int coord_type,
                             int isPlotted, int isTimeFixed, gnuplot_ctrl *h1)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(6*(man_grid_size+1));
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(6*(man_grid_size));
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(6*man_grid_size, 6*(man_grid_size+1));
    gsl_matrix *M    = gsl_matrix_calloc(6*(man_grid_size), 6*(man_grid_size));

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //----------------------------------------------------------------------
    // Arrival state: yma(1, :) = ymd(1,:);
    //----------------------------------------------------------------------
    for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter < 20)
    {

        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Arrival state (including STM) = final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                for(int j = 0; j < 6; j++)
                {
                    gsl_matrix_set(DF, i + 6*k, j + 6*k,      gsl_matrix_get(Ji[k], i, j));
                    gsl_matrix_set(DF, i + 6*k, j + 6*(k+1), -gsl_matrix_get(Id, i, j));
                }
            }
        }

        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);


        //----------------------------------------------------------------------
        // Update the free variables
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 6*k+i);
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }


    //--------------------------------------------------------------------------
    // Free
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msd_CM_RCM(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
                                int N, int man_grid_size, int coord_type,
                                int isPlotted, int isTimeFixed, gnuplot_ctrl *h1,
                                matrix<Ofsc>  &Mcoc_EM,
                                matrix<Oftsc> &DCM_EM_TFC,
                                gsl_matrix_complex *CCM_R_RCM,
                                SingleOrbit &orbit)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(6*(man_grid_size+1)-2);
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(6*(man_grid_size));
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(6*man_grid_size, 6*(man_grid_size+1)-2);
    gsl_matrix *M    = gsl_matrix_calloc(6*(man_grid_size), 6*(man_grid_size));

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));

    //Phi0
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);

    gsl_matrix *K1   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2   = gsl_matrix_calloc(6,4);
    matrix<Ofsc>  DWh_ofsc(6, 4);
    gsl_matrix_complex *RDWhc = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *PC  = gsl_matrix_complex_calloc(6,6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);


    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //----------------------------------------------------------------------
    // Arrival state: yma(1, :) = ymd(1,:);
    //----------------------------------------------------------------------
    for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter < 20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Arrival state (including STM) = final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc, RDWhc, true);

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, Mcoc_EM, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc, CCM_R_RCM, gslc_complex(0.0,0.0), K1c);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c, gslc_complex(0.0,0.0), K2c);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c, K2);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2, 0.0, K1);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1, 0.0, Phi0);

            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
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


        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);


        //----------------------------------------------------------------------
        // Update the free variables
        //----------------------------------------------------------------------
        //-------------------------------
        //First 4 correction variables is orbit.si
        //-------------------------------
        for(int i = 0; i < 4; i++) orbit.si[i] += gsl_vector_get(DQv, i);
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        orbit_update_ic(orbit, orbit.si, tmdn[0]/SEML.us_em.ns);
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < N; i++) yv[i] = orbit.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //The rest is classical cartesian coordinates at patch points
        for(int k = 1; k <= man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 6*k+i-2);
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }


    //--------------------------------------------------------------------------
    // Free
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msdvt_CM_RCM(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double **yma,
    int N,
    int man_grid_size,
    int coord_type,
    int isPlotted,
    int isTimeFixed,
    gnuplot_ctrl *h1,
    matrix<Ofsc>  &Mcoc_EM,
    matrix<Ofsc>  &Mcoc_SEM,
    vector<Ofsc>  &Vcoc_SEM,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    vector<Oftsc> &CM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM,
    SingleOrbit &orbit,
    SingleOrbit &orbit_SEM,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N], f[6];
    double te;

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(7*(man_grid_size+1)-5);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(6*(man_grid_size));
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(6*man_grid_size, 7*(man_grid_size+1)-5);
    gsl_matrix *M    = gsl_matrix_calloc(6*(man_grid_size), 6*(man_grid_size));

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,4);

    //Intermediate variables
    gsl_matrix *K1   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc(6, 4);
    gsl_matrix_complex *RDWhc = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }

    //----------------------------------------------------------------------
    // Arrival state: yma(1, :) = ymd(1,:);
    //----------------------------------------------------------------------
    for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  10)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Arrival state (including STM) = final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc, RDWhc, true);

                //cout << "n*t0 (2) = " << tmdn[0]/SEML.us_em.ns*SEML.us_em.n << endl;

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, Mcoc_EM, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc, CCM_R_RCM, gslc_complex(0.0,0.0), K1c);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c, gslc_complex(0.0,0.0), K2c);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c, K2);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2, 0.0, K1);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,4)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 4, DCM_SEM_TFC, DWh_ofsc, RDWhc, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, Mcoc_SEM, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc, CCM_R_RCM, gslc_complex(0.0,0.0), K1c);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c, gslc_complex(0.0,0.0), K2c);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c, K2);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == man_grid_size - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update the derivatives wrt to time
            //------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &SEML);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == man_grid_size-1)
            {
                //------------------------------------------
                // RCM to TFC
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zOut[i].zero();
                }

                RCMtoTFC(orbit_SEM.si, OFTS_ORDER, OFS_ORDER, 4, CM_SEM_TFC, zIn, 1);

                //------------------------------------------
                // TFC to NC to build zOut = CM_SEM_NC
                //------------------------------------------
                applyCOC(Mcoc_SEM, Vcoc_SEM, zIn, zOut);

                if(isDebug)
                {
                    cout << "zIn[1] = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

                //------------------------------------------
                // zIn = dot(zOut)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zIn[i].dot(zOut[i], SEML.us_sem.n);
                }


                //------------------------------------------
                // Then f = f - zIn(tf)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    f[i] -= creal(zIn[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                    yv[i] = creal(zOut[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                }

                if(isDebug)
                {
                    cout << "zOut[1] = " << endl;
                    cout << zOut[1] << endl;
                    cout << "--------------------" << endl;

                    cout << "dot(zOut[1]) = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

            }


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //---------------------------
                //DF/DS
                //---------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == man_grid_size-1)
                {
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i + 6*k, j + 7*k-3,       gsl_matrix_get(Ji[k], i, j));
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-3,  -gsl_matrix_get(PhiN, i, j));
                }
                else
                {
                    for(int j = 0; j < 6; j++)
                    {
                        gsl_matrix_set(DF, i + 6*k, j + 7*k-3,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 6*k, j + 7*(k+1)-3, -gsl_matrix_get(Id, i, j));
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
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }
                else if(k == man_grid_size-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size+1, f[i]);
                }
                else
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+1)-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }

            }


        }

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit.si[i] += gsl_vector_get(DQv, i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit, orbit.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < N; i++) yv[i] = orbit.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 7*k-3);
            tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-4);
        }
        //Last time:
        tmdn[man_grid_size] += gsl_vector_get(DQv, 7*man_grid_size+1);


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_SEM.si[i] += 1e-1*gsl_vector_get(DQv, i + 7*man_grid_size-3);

        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 4);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        for(int i = 0; i < N; i++) yv[i] = orbit_SEM.z0[i];
        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = ye[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size-1; k++)
    {

        for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
        ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);
        if(isPlotted) gnuplot_plot_xyz(h1, ym[0], ym[1],  ym[2], 2, (char*)"", "lines", "2", "4", 5);
    }
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "4", 5);



    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msdvt_CMS_RCM(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables], f[6];
    double te;

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(7*(man_grid_size+1)-4);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(6*(man_grid_size));
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(6*man_grid_size, 7*(man_grid_size+1)-4);
    gsl_matrix *M    = gsl_matrix_calloc(6*(man_grid_size), 6*(man_grid_size));

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //cout << "n*t0 (2) = " << tmdn[0]/SEML.us_em.ns*SEML.us_em.n << endl;

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == man_grid_size - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update the derivatives wrt to time
            //------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &SEML);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == man_grid_size-1)
            {
                //------------------------------------------
                // RCM to TFC
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zOut[i].zero();
                }

                RCMtoTFC(orbit_SEM.si, OFTS_ORDER, OFS_ORDER, 5, *orbit_SEM.Wh, zIn, 1);

                //------------------------------------------
                // TFC to NC to build zOut = CM_SEM_NC
                //------------------------------------------
                applyCOC(*orbit_SEM.PC, *orbit_SEM.V, zIn, zOut);

                if(isDebug)
                {
                    cout << "zIn[1] = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

                //------------------------------------------
                // zIn = dot(zOut)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zIn[i].dot(zOut[i], SEML.us_sem.n);
                }


                //------------------------------------------
                // Then f = f - zIn(tf)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    f[i] -= creal(zIn[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                }

                if(isDebug)
                {
                    cout << "zOut[1] = " << endl;
                    cout << zOut[1] << endl;
                    cout << "--------------------" << endl;

                    cout << "dot(zOut[1]) = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

            }


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //---------------------------
                //DF/DS
                //---------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == man_grid_size-1)
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
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }
                else if(k == man_grid_size-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size+2, f[i]);
                }
                else
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+1)-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }

            }


        }

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Same with LU decomposition:
        //        int s;
        //        gsl_permutation * p = gsl_permutation_alloc (6*(man_grid_size));
        //        gsl_linalg_LU_decomp (M, p, &s);
        //        gsl_linalg_LU_solve(M, p, Fv, K3);
        //        gsl_permutation_free (p);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_EM.si[i] += gsl_vector_get(DQv, i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 7*k-3);
            tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-4);
        }
        //Last time:
        tmdn[man_grid_size] += gsl_vector_get(DQv, 7*man_grid_size+2);


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.si[i] += gsl_vector_get(DQv, i + 7*man_grid_size-3);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;

        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);

    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "4", 4);



    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msdvt_CMS_RCM_deps(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables], f[6];
    double te;

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 7*(man_grid_size+1)-4;  //free variables
    int ncs = 6*man_grid_size;        //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //cout << "n*t0 (2) = " << tmdn[0]/SEML.us_em.ns*SEML.us_em.n << endl;

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == man_grid_size - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update the derivatives wrt to time
            //------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &SEML);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == man_grid_size-1)
            {
                //------------------------------------------
                // RCM to TFC
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zOut[i].zero();
                }

                RCMtoTFC(orbit_SEM.si, OFTS_ORDER, OFS_ORDER, 5, *orbit_SEM.Wh, zIn, 1);

                //------------------------------------------
                // TFC to NC to build zOut = CM_SEM_NC
                //------------------------------------------
                applyCOC(*orbit_SEM.PC, *orbit_SEM.V, zIn, zOut);

                if(isDebug)
                {
                    cout << "zIn[1] = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

                //------------------------------------------
                // zIn = dot(zOut)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zIn[i].dot(zOut[i], SEML.us_sem.n);
                }


                //------------------------------------------
                // Then f = f - zIn(tf)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    f[i] -= creal(zIn[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                }

                if(isDebug)
                {
                    cout << "zOut[1] = " << endl;
                    cout << zOut[1] << endl;
                    cout << "--------------------" << endl;

                    cout << "dot(zOut[1]) = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

            }


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //---------------------------
                //DF/DS
                //---------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == man_grid_size-1)
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
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }
                else if(k == man_grid_size-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*man_grid_size+2, f[i]);
                }
                else
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+1)-4, gsl_vector_get(K4, i));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]]
                    //--------------------------
                    gsl_matrix_set(DF, i + 6*k, 7*(k+2)-4, f[i]);
                }

            }


        }

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_EM.si[i] += gsl_vector_get(DQv, i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 7*k-3);
            tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-4);
        }
        //Last time:
        tmdn[man_grid_size] += gsl_vector_get(DQv, 7*man_grid_size+2);


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.si[i] += gsl_vector_get(DQv, i + 7*man_grid_size-3);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;

        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);

    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        //Here, we want to make s_EM[0] "grow"
        sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }
    else
    {
        //Here, we want to go "in the same direction for s_EM = Q[0:3]": no u_turn!
        for(int i = 0; i < 4; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];
        //sign = dotNV > 0? 1:-1;

        //OR

        //Here, we want to make s_SEM[4] "decrease"
        sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);



    //    gsl_vector *Knv  = gsl_vector_calloc(nfv);
    //    gsl_vector *Knv2  = gsl_vector_calloc(ncs);
    //
    //    cout << "HERE" << endl;
    //    for(int i = 0; i < nfv; i++) gsl_vector_set(Knv, i, nullvector[i]);
    //
    //    gsl_blas_dgemv (CblasNoTrans, 1.0, DF, Knv, 0.0, Knv2);
    //
    //    cout << "Knv2 = " << endl;
    //    for(int i = 0; i < ncs; i++) cout << gsl_vector_get(Knv2, i) << endl;


    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msdvt_CMS_RCM_deps_ATF(double **ymd,
                           double *tmd,
                           double **ymdn,
                           double *tmdn,
                           double *nullvector,
                           int number_of_variables,
                           int man_grid_size,
                           int coord_type,
                           double precision,
                           int isFirst,
                           matrix<Oftsc> &DCM_EM_TFC,
                           matrix<Oftsc> &DCM_SEM_TFC,
                           gsl_matrix_complex *CCM_R_RCM_EM,
                           gsl_matrix_complex *CCM_R_RCM_SEM,
                           SingleOrbit &orbit_EM,
                           SingleOrbit &orbit_SEM,
                           gnuplot_ctrl *h1,
                           int isPlotted,
                           int isUserDefined,
                           int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 6*man_grid_size+3;  //free variables
    int ncs = 6*man_grid_size;    //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);
    //Error isolated at final point
    gsl_vector *Fvn  = gsl_vector_calloc(6);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    int iterMax = 20;
    while(iter <  iterMax)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //cout << "n*t0 (2) = " << tmdn[0]/SEML.us_em.ns*SEML.us_em.n << endl;

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
                if(k == man_grid_size - 1) gsl_vector_set(Fvn, i, ye[i] - ymdn[i][k+1]);
            }

            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                //---------------------------
                //DF/DS
                //---------------------------
                if(k == 0)
                {
                    for(int j = 0; j < 4; j++) gsl_matrix_set(DF, i, j,    gsl_matrix_get(Phi0, i, j));
                    for(int j = 0; j < 6; j++) gsl_matrix_set(DF, i, j+4, -gsl_matrix_get(Id, i, j));
                }
                else if(k == man_grid_size-1)
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

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        for(int i = 0; i < 4; i++) orbit_EM.si[i] += gsl_vector_get(DQv, i);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, i + 6*k-2);
        }

        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.si[i] += gsl_vector_get(DQv, i+ 6*man_grid_size-2);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);
    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][man_grid_size], &ymdn[1][man_grid_size],  &ymdn[2][man_grid_size], 1, (char*)"", "points", "2", "2", 0);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(isUserDefined)
        {
            do
            {
                cout << "After refinment: s1_CMU_EM = " << orbit_EM.si[0] << endl;
                cout << "Please choose a direction for the continuation procedure:" << endl;
                cout << "+1: s1_CMU_EM is increasing" << endl;
                cout << "-1: s1_CMU_EM is decreasing" << endl;
                cout << " to select a specific starting time" << endl;
                cin >> dti;
                ti = (int) dti;
            }while(ti != 1 && ti != -1);
        }else ti = 1;

        //Here, we want to make s_EM[1] "grow"
        sign = gsl_matrix_get(Q, 1, nfv-1) > 0? ti:-ti;
    }
    else
    {
        //OR
        //Here, we want to go "in the same direction" for some components of Q
        dotNV += gsl_matrix_get(Q, 0, nfv-1)*nullvector[0];                //CMU of  EML2
        dotNV += gsl_matrix_get(Q, 1, nfv-1)*nullvector[1];                //CMU of  EML2
        dotNV += gsl_matrix_get(Q, 2, nfv-1)*nullvector[2];                //CMU of  EML2
        dotNV += gsl_matrix_get(Q, 3, nfv-1)*nullvector[3];                //CMU of  EML2
        sign = dotNV > 0? 1:-1;
        //        dotNV += gsl_matrix_get(Q, nfv-4, nfv-1)*nullvector[nfv-4];                //CMS of  SEML2
        //        dotNV += gsl_matrix_get(Q, nfv-2, nfv-1)*nullvector[nfv-2];                //CMS of  SEML2
        //        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEML2
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);

    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
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
 **/
int msdvt_CMS_RCM_deps_planar(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables], f[6];
    double te;

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 5*man_grid_size+1;  //free variables
    int ncs = 4*man_grid_size;    //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //------------------------------------------------------------------
            // Update the derivatives wrt to time
            //------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &SEML);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == man_grid_size-1)
            {
                //------------------------------------------
                // RCM to TFC
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zOut[i].zero();
                }

                RCMtoTFC(orbit_SEM.si, OFTS_ORDER, OFS_ORDER, 5, *orbit_SEM.Wh, zIn, 1);

                //------------------------------------------
                // TFC to NC to build zOut = CM_SEM_NC
                //------------------------------------------
                applyCOC(*orbit_SEM.PC, *orbit_SEM.V, zIn, zOut);

                if(isDebug)
                {
                    cout << "zIn[1] = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

                //------------------------------------------
                // zIn = dot(zOut)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zIn[i].dot(zOut[i], SEML.us_sem.n);
                }


                //------------------------------------------
                // Then f = f - zIn(tf)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    f[i] -= creal(zIn[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                }

                if(isDebug)
                {
                    cout << "zOut[1] = " << endl;
                    cout << zOut[1] << endl;
                    cout << "--------------------" << endl;

                    cout << "dot(zOut[1]) = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

            }


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
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
                else if(k == man_grid_size-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-8,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-7,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-6,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-5,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-3,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-2,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-1,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-8,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-7,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-6,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-5,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-3,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-2,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-1,  -gsl_matrix_get(PhiN, i+1, 4));
                    }
                }
                else
                {
                    for(int j = 0; j < 4; j++)
                    {
                        gsl_matrix_set(DF, i + 4*k, j + 5*k-3,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 4*k, j + 5*(k+1)-3, -gsl_matrix_get(Id, i, j));
                    }

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
                else if(k == man_grid_size-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-4, gsl_vector_get(K4, i));
                    else gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-4, gsl_vector_get(K4, i+1));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*man_grid_size, f[i]);
                    else gsl_matrix_set(DF, i + 4*k, 5*man_grid_size, f[i+1]);
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

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += gsl_vector_get(DQv, 0);
        orbit_EM.si[2] += gsl_vector_get(DQv, 1);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 0 + 5*k-3);
            ymdn[1][k] += gsl_vector_get(DQv, 1 + 5*k-3);
            ymdn[3][k] += gsl_vector_get(DQv, 2 + 5*k-3);
            ymdn[4][k] += gsl_vector_get(DQv, 3 + 5*k-3);
            //tmdn[k] = max(tmdn[k] + gsl_vector_get(DQv, 5*k+1), tmdn[k-1]);
            tmdn[k] =tmdn[k] + gsl_vector_get(DQv, 5*k+1);
        }
        //Last time:
        //tmdn[man_grid_size] = max(tmdn[man_grid_size]+ gsl_vector_get(DQv, 5*man_grid_size), tmdn[man_grid_size-1]);
        tmdn[man_grid_size] = tmdn[man_grid_size]+ gsl_vector_get(DQv, 5*man_grid_size);


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += gsl_vector_get(DQv, 5*man_grid_size-3);
        orbit_SEM.si[2] += gsl_vector_get(DQv, 5*man_grid_size-2);
        orbit_SEM.si[4] += gsl_vector_get(DQv, 5*man_grid_size-1);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][man_grid_size], &ymdn[1][man_grid_size],  &ymdn[2][man_grid_size], 1, (char*)"", "points", "2", "2", 0);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    if(isFirst)
    {
        //Here, we want to make s_SEM[4] "decrease"
        if(orbit_SEM.si[4] > 0) sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;
        else sign = gsl_matrix_get(Q, nfv-2, nfv-1) > 0? 1:-1;
    }
    else
    {
        //Here, we want to make s_SEM[4] "decrease"
        if(orbit_SEM.si[4] > 0) sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;
        else sign = gsl_matrix_get(Q, nfv-2, nfv-1) > 0? 1:-1;;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);



    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions and fixed last time.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_ATF(double **ymd,
                                  double *tmd,
                                  double **ymdn,
                                  double *tmdn,
                                  double *nullvector,
                                  int number_of_variables,
                                  int man_grid_size,
                                  int coord_type,
                                  double precision,
                                  int isFirst,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Oftsc> &DCM_SEM_TFC,
                                  gsl_matrix_complex *CCM_R_RCM_EM,
                                  gsl_matrix_complex *CCM_R_RCM_SEM,
                                  SingleOrbit &orbit_EM,
                                  SingleOrbit &orbit_SEM,
                                  gnuplot_ctrl *h1,
                                  int isPlotted,
                                  int isUserDefined,
                                  int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 4*(man_grid_size-1) + 5;    //free variables
    int ncs = 4*man_grid_size;            //constraints

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    int iterMax = 20;
    while(iter <  iterMax)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
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
                else if(k == man_grid_size-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-6,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-5,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-4,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-3,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-2,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-1,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-0,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-6,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-5,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-4,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-3,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-2,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-1,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-0,  -gsl_matrix_get(PhiN, i+1, 4));
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
        }

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }


        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += gsl_vector_get(DQv, 0);
        orbit_EM.si[2] += gsl_vector_get(DQv, 1);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            return GSL_FAILURE;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }

        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += gsl_vector_get(DQv, 4*man_grid_size-2);
        orbit_SEM.si[2] += gsl_vector_get(DQv, 4*man_grid_size-1);
        orbit_SEM.si[4] += gsl_vector_get(DQv, 4*man_grid_size-0);


        //-----------------------------------------------------
        //Display
        //-----------------------------------------------------
        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);
    }


    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][man_grid_size], &ymdn[1][man_grid_size],  &ymdn[2][man_grid_size], 1, (char*)"", "points", "2", "2", 0);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        if(isUserDefined)
        {
            do
            {
                cout << "After refinment: s1_CMU_EM = " << orbit_EM.si[0] << endl;
                cout << "Please choose a direction for the continuation procedure:" << endl;
                cout << "+1: s1_CMU_EM is increasing" << endl;
                cout << "-1: s1_CMU_EM is decreasing" << endl;
                cout << " to select a specific starting time" << endl;
                cin >> dti;
                ti = (int) dti;
            }while(ti != 1 && ti != -1);
        }else ti = 1;

        //Here, we want to make s_EM[0] "grow"
        sign = gsl_matrix_get(Q, 0, nfv-1) > 0? ti:-ti;
    }
    else
    {
        //OR
        //Here, we want to go "in the same direction for some components of Q
        for(int i = 0; i < 2; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];                //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEML2
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);


    //-----------------------------------------------------
    //FINALLY, we update orbit_EM.tf, even if it is not necessary
    //-----------------------------------------------------
    orbit_EM.tf = tmdn[man_grid_size]/SEML.us_em.ns;

    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);

    if(iter >= iterMax) return GSL_FAILURE;
    else return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions and fixed last time + pseuod-arclength constraint.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_pac_ATF(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    double *conv_free_var,
    double ds,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 4*man_grid_size+1;    //free variables
    int ncs = 4*man_grid_size+1;    //constraints

    //Free variable vector
    double *fvc = dvector(0, nfv-1);

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //For GSL inversion
    gsl_permutation * perm = gsl_permutation_alloc (nfv);
    int s;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }


            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }

            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
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
                else if(k == man_grid_size-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-6,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-5,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-4,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-3,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-2,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-1,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-0,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-6,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-5,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-4,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-3,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-2,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-1,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 4*man_grid_size-0,  -gsl_matrix_get(PhiN, i+1, 4));
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
        }

        //----------------------------------------------------------------------
        // Build the free variable vector
        //----------------------------------------------------------------------
        // CM of EML2
        fvc[0] = orbit_EM.si[0];
        fvc[1] = orbit_EM.si[2];

        //Patch point
        for(int p = 1; p < man_grid_size; p++)
        {
            fvc[0 + 4*p-2] = ymdn[0][p];
            fvc[1 + 4*p-2] = ymdn[1][p];
            fvc[2 + 4*p-2] = ymdn[3][p];
            fvc[3 + 4*p-2] = ymdn[4][p];
        }

        //CMS of SEML2
        fvc[4*man_grid_size-2] = orbit_SEM.si[0];
        fvc[4*man_grid_size-1] = orbit_SEM.si[2];
        fvc[4*man_grid_size-0] = orbit_SEM.si[4];

        //----------------------------------------------------------------------
        //Last row of the error vector is the pseudo-arclength constraint
        //----------------------------------------------------------------------
        double res = 0;
        for(int i = 0; i < nfv; i++)
        {
            fvc[i] -= conv_free_var[i];
            res += fvc[i]*nullvector[i];
        }
        gsl_vector_set(Fv, ncs-1, res - ds);

        //----------------------------------------------------------------------
        //Last row of the Jacobian is equal to the null vector of the previous step
        //----------------------------------------------------------------------
        for(int i = 0; i < nfv; i++) gsl_matrix_set(DF, ncs-1, i, nullvector[i]);


        //------------------------------------------------------------------
        // Display
        //------------------------------------------------------------------
        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }
        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        gsl_linalg_LU_decomp (DF, perm , &s);
        gsl_linalg_LU_solve(DF, perm, Fv, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += gsl_vector_get(DQv, 0);
        orbit_EM.si[2] += gsl_vector_get(DQv, 1);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 4*k-2);
            ymdn[1][k] += gsl_vector_get(DQv, 4*k-1);
            ymdn[3][k] += gsl_vector_get(DQv, 4*k-0);
            ymdn[4][k] += gsl_vector_get(DQv, 4*k+1);
        }

        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += gsl_vector_get(DQv, 4*man_grid_size-2);
        orbit_SEM.si[2] += gsl_vector_get(DQv, 4*man_grid_size-1);
        orbit_SEM.si[4] += gsl_vector_get(DQv, 4*man_grid_size-0);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);
    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][man_grid_size], &ymdn[1][man_grid_size],  &ymdn[2][man_grid_size], 1, (char*)"", "points", "2", "2", 0);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1, ti = 1;
    double dti = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        do
        {
            cout << "After refinment: s1_CMU_EM = " << orbit_EM.si[0] << endl;
            cout << "Please choose a direction for the continuation procedure:" << endl;
            cout << "+1: s1_CMU_EM is increasing" << endl;
            cout << "-1: s1_CMU_EM is decreasing" << endl;
            cout << " to select a specific starting time" << endl;
            cin >> dti;
            ti = (int) dti;
        }while(ti != 1 && ti != -1);

        //Here, we want to make s_EM[0] "grow"
        sign = gsl_matrix_get(Q, 0, nfv-1) > 0? ti:-ti;
    }
    else
    {
        //OR
        //Here, we want to go "in the same direction for some components of Q
        for(int i = 0; i < 2; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];                //CMU of  EML2
        sign = dotNV > 0? 1:-1;

        //OR
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        //for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR always decrease the stable component at SEML2
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;

        //OR always increase the s1 component at EML2
        //sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);



    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}


/**
 * \brief Multiple shooting scheme with no boundary conditions + pseuod-arclength constraint
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int msdvt_CMS_RCM_deps_planar_pac(
    double **ymd,
    double *tmd,
    double **ymdn,
    double *tmdn,
    double *nullvector,
    double *conv_free_var,
    double ds,
    int number_of_variables,
    int man_grid_size,
    int coord_type,
    double precision,
    int isFirst,
    matrix<Oftsc> &DCM_EM_TFC,
    matrix<Oftsc> &DCM_SEM_TFC,
    gsl_matrix_complex *CCM_R_RCM_EM,
    gsl_matrix_complex *CCM_R_RCM_SEM,
    SingleOrbit &orbit_EM,
    SingleOrbit &orbit_SEM,
    gnuplot_ctrl *h1,
    int isPlotted,
    int isDebug)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[number_of_variables], ye[number_of_variables], f[6];
    double te;

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    int nfv = 5*man_grid_size+1;    //free variables
    int ncs = 4*man_grid_size+1;    //constraints

    //Free variable vector
    double *fvc = dvector(0, nfv-1);

    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(nfv);

    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(ncs);

    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(ncs, nfv);
    gsl_matrix *M    = gsl_matrix_calloc(ncs, ncs);

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(ncs);
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //Phi0: Jacobian wrt to EM RCM variables
    gsl_matrix *Phi0 = gsl_matrix_calloc(6,4);
    //PhiN: Jacobian wrt to SEM RCM variables
    gsl_matrix *PhiN = gsl_matrix_calloc(6,5);

    //Intermediate variables
    gsl_matrix *K1_EM   = gsl_matrix_calloc(6,4);
    gsl_matrix *K2_EM   = gsl_matrix_calloc(6,4);
    matrix<Ofsc> DWh_ofsc_EM(6, 4);
    gsl_matrix_complex *RDWhc_EM = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K1c_EM   = gsl_matrix_complex_calloc(6,4);
    gsl_matrix_complex *K2c_EM   = gsl_matrix_complex_calloc(6,4);

    gsl_matrix *K2_SEM    = gsl_matrix_calloc(6,5);
    matrix<Ofsc> DWh_ofsc_SEM (6, 5);
    gsl_matrix_complex *RDWhc_SEM  = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K1c_SEM    = gsl_matrix_complex_calloc(6,5);
    gsl_matrix_complex *K2c_SEM    = gsl_matrix_complex_calloc(6,5);


    gsl_matrix_complex *PC    = gsl_matrix_complex_calloc(6,6);
    //For time derivatives
    vector<Ofsc> zIn(6), zOut(6);

    //NCEM to coord_type (e.g. NCSEM_R_NCEM)
    gsl_matrix *COORD_R_NCEM = gsl_matrix_calloc(6,6);
    //NCSEM to coord_type (e.g. NCSEM_R_NCSEM)
    gsl_matrix *COORD_R_NCSEM = gsl_matrix_calloc(6,6);

    //Norms
    double si_norm_EM, si_norm_SEM;

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < number_of_variables; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter <  20)
    {
        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < number_of_variables; i++) ye[i] = ym[i][1];
            te = tm[1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update Phi0
            //------------------------------------------------------------------
            if(k == 0)
            {
                //----------------------------
                //RDWhc = DCM_EM_TFC(orbit_EM.si, t0), in EM units, in C(6,4)
                //----------------------------
                //Here we suppose that the default framework is SEM, so we need to normalize the time
                RCMtoTFC_JAC(orbit_EM.si, tmdn[0]/SEML.us_em.ns, SEML.us_em.n, OFTS_ORDER, OFS_ORDER, 4, DCM_EM_TFC, DWh_ofsc_EM, RDWhc_EM, true);

                //----------------------------
                //PC = Mcoc_EM(t0), in EM units, in C(6,6)
                //----------------------------
                //Here we again suppose that the default framework is SEM, so we need to normalize the time
                evaluate(tmdn[0]/SEML.us_em.ns, SEML.us_em.n, *orbit_EM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,4)*C(4,4) = C(6,4)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_EM, CCM_R_RCM_EM, gslc_complex(0.0,0.0), K1c_EM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_EM, gslc_complex(0.0,0.0), K2c_EM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_EM, K2_EM);

                //----------------------------
                //COORD_R_NCEM
                //----------------------------
                rot_mat_coc(tmdn[0]/SEML.us_em.ns, COORD_R_NCEM, NCEM, coord_type);

                //----------------------------
                //K1 = COORD_R_NCEM*K2, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCEM, K2_EM, 0.0, K1_EM);

                //----------------------------
                //Phi0 = Ji[0]*K1, in R(6,6)*R(6,4) = R(6,4)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Ji[0], K1_EM, 0.0, Phi0);

                if(isDebug)
                {
                    cout << "Phi0 = " << endl;
                    gslc_matrix_printf(Phi0);
                }
            }


            //------------------------------------------------------------------
            // Update PhiN
            //------------------------------------------------------------------
            if(k == man_grid_size-1)
            {
                //----------------------------
                //RDWhc = DCM_SEM_TFC(orbit_SEM.si, tN), in SEM units, in C(6,5)
                //----------------------------
                RCMtoTFC_JAC(orbit_SEM.si, tmdn[man_grid_size], SEML.us_sem.n, OFTS_ORDER, OFS_ORDER, 5, DCM_SEM_TFC, DWh_ofsc_SEM, RDWhc_SEM, true);

                //----------------------------
                //PC = Mcoc_SEM(tN), in SEM units, in C(6,6)
                //----------------------------
                evaluate(tmdn[man_grid_size], SEML.us_sem.n, *orbit_SEM.PC, PC);

                //----------------------------
                //K2c = PC*RDWhc*CCM_R_RCM, in C(6,6)*C(6,5)*C(5,5) = C(6,5)
                //----------------------------
                //K1c = RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), RDWhc_SEM, CCM_R_RCM_SEM, gslc_complex(0.0,0.0), K1c_SEM);
                //K2c = PC*RDWhc*CCM_R_RCM
                gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c_SEM, gslc_complex(0.0,0.0), K2c_SEM);

                //----------------------------
                //K2 = real(K2c), in R(6,4)
                //----------------------------
                gslc_matrix_complex_to_matrix(K2c_SEM, K2_SEM);

                //----------------------------
                //COORD_R_NCSEM (useless for now, but may be needed in future release)
                //----------------------------
                rot_mat_coc(tmdn[man_grid_size], COORD_R_NCSEM, NCSEM, coord_type);

                //----------------------------
                //PhiN = COORD_R_NCSEM*K2, in R(6,6)*R(6,5) = R(6,5)
                //----------------------------
                gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NCSEM, K2_SEM, 0.0, PhiN);

                if(isDebug)
                {
                    cout << "PhiN = " << endl;
                    gslc_matrix_printf(PhiN);
                }
            }

            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            // EXCEPT THE LAST ROW
            //------------------------------------------------------------------
            for(int i = 0; i < 4; i++)
            {
                if(i < 2)  gsl_vector_set(Fv, 4*k+i, ye[i] - ymdn[i][k+1]);
                if(i >= 2) gsl_vector_set(Fv, 4*k+i, ye[i+1] - ymdn[i+1][k+1]);
            }


            //------------------------------------------------------------------
            // Update the derivatives wrt to time
            //------------------------------------------------------------------
            //------------------------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //------------------------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //K4 = dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //------------------------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //------------------------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(te, ye, f, &SEML);

            //Special case of the last point: dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]) - dCM_SEM_NC/dt[k+1]
            if(k == man_grid_size-1)
            {
                //------------------------------------------
                // RCM to TFC
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zOut[i].zero();
                }

                RCMtoTFC(orbit_SEM.si, OFTS_ORDER, OFS_ORDER, 5, *orbit_SEM.Wh, zIn, 1);

                //------------------------------------------
                // TFC to NC to build zOut = CM_SEM_NC
                //------------------------------------------
                applyCOC(*orbit_SEM.PC, *orbit_SEM.V, zIn, zOut);

                if(isDebug)
                {
                    cout << "zIn[1] = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

                //------------------------------------------
                // zIn = dot(zOut)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    zIn[i].zero();
                    zIn[i].dot(zOut[i], SEML.us_sem.n);
                }


                //------------------------------------------
                // Then f = f - zIn(tf)
                //------------------------------------------
                for(int i = 0; i < 6; i++)
                {
                    f[i] -= creal(zIn[i].evaluate(tmdn[man_grid_size]*SEML.us_sem.n));
                }

                if(isDebug)
                {
                    cout << "zOut[1] = " << endl;
                    cout << zOut[1] << endl;
                    cout << "--------------------" << endl;

                    cout << "dot(zOut[1]) = " << endl;
                    cout << zIn[1] << endl;
                    cout << "--------------------" << endl;
                }

            }


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
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
                else if(k == man_grid_size-1)
                {
                    if(i < 2)
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-8,  gsl_matrix_get(Ji[k], i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-7,  gsl_matrix_get(Ji[k], i, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-6,  gsl_matrix_get(Ji[k], i, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-5,  gsl_matrix_get(Ji[k], i, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-3,  -gsl_matrix_get(PhiN, i, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-2,  -gsl_matrix_get(PhiN, i, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-1,  -gsl_matrix_get(PhiN, i, 4));
                    }
                    else
                    {
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-8,  gsl_matrix_get(Ji[k], i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-7,  gsl_matrix_get(Ji[k], i+1, 1));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-6,  gsl_matrix_get(Ji[k], i+1, 3));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-5,  gsl_matrix_get(Ji[k], i+1, 4));

                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-3,  -gsl_matrix_get(PhiN, i+1, 0));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-2,  -gsl_matrix_get(PhiN, i+1, 2));
                        gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-1,  -gsl_matrix_get(PhiN, i+1, 4));
                    }
                }
                else
                {
                    for(int j = 0; j < 4; j++)
                    {
                        gsl_matrix_set(DF, i + 4*k, j + 5*k-3,      gsl_matrix_get(Ji[k], i, j));
                        gsl_matrix_set(DF, i + 4*k, j + 5*(k+1)-3, -gsl_matrix_get(Id, i, j));
                    }

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
                else if(k == man_grid_size-1)
                {
                    //--------------------------
                    //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-4, gsl_vector_get(K4, i));
                    else gsl_matrix_set(DF, i + 4*k, 5*man_grid_size-4, gsl_vector_get(K4, i+1));

                    //--------------------------
                    //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1]] - dCM_SEM_NC/dt[k+1]
                    //--------------------------
                    if(i < 2) gsl_matrix_set(DF, i + 4*k, 5*man_grid_size, f[i]);
                    else gsl_matrix_set(DF, i + 4*k, 5*man_grid_size, f[i+1]);
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

        //----------------------------------------------------------------------
        // Build the free variable vector
        //----------------------------------------------------------------------
        // CM of EML2
        fvc[0] = orbit_EM.si[0];
        fvc[1] = orbit_EM.si[2];

        //Patch point
        for(int p = 1; p < man_grid_size; p++)
        {
            fvc[0 + 5*p-3] = ymdn[0][p];
            fvc[1 + 5*p-3] = ymdn[1][p];
            fvc[2 + 5*p-3] = ymdn[3][p];
            fvc[3 + 5*p-3] = ymdn[4][p];
            fvc[5*p+1] = tmdn[p];
        }

        //Last time:
        fvc[5*man_grid_size] = tmdn[man_grid_size];

        //CMS of SEML2
        fvc[5*man_grid_size-3] = orbit_SEM.si[0];
        fvc[5*man_grid_size-2] = orbit_SEM.si[2];
        fvc[5*man_grid_size-1] = orbit_SEM.si[4];

        //----------------------------------------------------------------------
        //Last row of the error vector is the pseudo-arclength constraint
        //----------------------------------------------------------------------
        double res = 0;
        for(int i = 0; i < nfv; i++)
        {
            fvc[i] -= conv_free_var[i];
            res += fvc[i]*nullvector[i];
        }
        gsl_vector_set(Fv, ncs-1, res - ds);

        //----------------------------------------------------------------------
        //Last row of the Jacobian is equal to the null vector of the previous step
        //----------------------------------------------------------------------
        for(int i = 0; i < nfv; i++) gsl_matrix_set(DF, ncs-1, i, nullvector[i]);


        //------------------------------------------------------------------
        // Display
        //------------------------------------------------------------------
        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "F = " << endl;
            for(int j = 0; j < (int) Fv->size; j++) cout << gsl_vector_get(Fv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "---------------------------------------------------------------" << endl;
            cout << "DF = " << endl;
            gslc_matrix_printf_approx(DF);
            cout << "---------------------------------------------------------------" << endl;
        }

        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "M = " << endl;
            gslc_matrix_printf_approx(M);
            cout << "---------------------------------------------------------------" << endl;
        }


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < precision)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);

        //Then, DQv = DF*inv(M)*DF
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);

        if(isDebug && man_grid_size == 2)
        {
            cout << "---------------------------------------------------------------" << endl;
            cout << "DQv = " << endl;
            for(int j = 0; j < (int) DQv->size; j++) cout << gsl_vector_get(DQv, j) << endl;
            cout << "---------------------------------------------------------------" << endl;
        }


        //======================================================================
        // Update the free variables
        //======================================================================
        //-----------------------------------------------------
        //First 4 correction variables is orbit_EM.si
        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += gsl_vector_get(DQv, 0);
        orbit_EM.si[2] += gsl_vector_get(DQv, 1);

        //Here we suppose that the default framework is SEM, so we need to normalize the time
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, tmdn[0]/SEML.us_em.ns);

        // Norm check
        si_norm_EM  = ENorm(orbit_EM.si, 4);
        if(si_norm_EM > SI_NORM_EM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_EM has reached its limits. break"  << endl;
            cout << "si_norm_EM = "   << si_norm_EM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //To CM_EM_NCSEM coordinates
        //Here we suppose that the default framework is SEM, so we need to normalize the time
        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(tmdn[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            ymdn[0][k] += gsl_vector_get(DQv, 0 + 5*k-3);
            ymdn[1][k] += gsl_vector_get(DQv, 1 + 5*k-3);
            ymdn[3][k] += gsl_vector_get(DQv, 2 + 5*k-3);
            ymdn[4][k] += gsl_vector_get(DQv, 3 + 5*k-3);
            //tmdn[k] = max(tmdn[k] + gsl_vector_get(DQv, 5*k+1), tmdn[k-1]);
            tmdn[k] =tmdn[k] + gsl_vector_get(DQv, 5*k+1);
        }
        //Last time:
        //tmdn[man_grid_size] = max(tmdn[man_grid_size]+ gsl_vector_get(DQv, 5*man_grid_size), tmdn[man_grid_size-1]);
        tmdn[man_grid_size] = tmdn[man_grid_size]+ gsl_vector_get(DQv, 5*man_grid_size);


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += gsl_vector_get(DQv, 5*man_grid_size-3);
        orbit_SEM.si[2] += gsl_vector_get(DQv, 5*man_grid_size-2);
        orbit_SEM.si[4] += gsl_vector_get(DQv, 5*man_grid_size-1);

        if(isDebug)
        {
            cout << "si_SEM = " << endl;
            vector_printf(orbit_SEM.si, 5);
        }


        // Norm check
        si_norm_SEM = ENorm(orbit_SEM.si, 5);
        if(si_norm_SEM > SI_NORM_SEM_MAX)
        {
            cout << "#########################################" << endl;
            cout << "si_norm_SEM has reached its limits. break" << endl;
            cout << "si_norm_SEM = " << si_norm_SEM             << endl;
            cout << "#########################################" << endl;
            break;
        }

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tmdn[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        //may change in a future release.
        //        for(int i = 0; i < number_of_variables; i++) yv[i] = orbit_SEM.z0[i];
        //        qbcp_coc(tmdn[man_grid_size], yv, ye, NCSEM, coord_type);
        for(int i = 0; i < 6; i++) ymdn[i][man_grid_size] = orbit_SEM.z0[i];

        //----------------------------------------------------------------------
        // Norm display
        //----------------------------------------------------------------------
        //cout << "nerror[end]*gamma  = " << gsl_blas_dnrm2(Fvn)*SEML.cs_sem.gamma << endl;
        //cout << "--------------------" << endl;
        if(isDebug)
        {
            cout << "si_norm_EM = "   << si_norm_EM << endl;
            cout << "si_norm_SEM = "  << si_norm_SEM << endl;
        }


        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        //if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
        //        char ch;
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);
    }

    //------------------------------------------------------------------
    //Last plot
    //------------------------------------------------------------------
    if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "2", "2", 4);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][man_grid_size], &ymdn[1][man_grid_size],  &ymdn[2][man_grid_size], 1, (char*)"", "points", "2", "2", 0);
    if(isPlotted) gnuplot_plot_xyz(h1, &ymdn[0][0], &ymdn[1][0],  &ymdn[2][0], 1, (char*)"", "points", "2", "2", 0);


    //-------------------------------------------------------------------------------
    //Compute the null vector: QR decomposition of DP^T
    //-------------------------------------------------------------------------------
    //QR elements
    gsl_vector *work  = gsl_vector_calloc(ncs);
    gsl_matrix *Q     = gsl_matrix_calloc(nfv,nfv);
    gsl_matrix *R     = gsl_matrix_calloc(nfv,ncs);
    gsl_matrix *DFT   = gsl_matrix_calloc(nfv,ncs);


    //DPT = transpose(DP)
    gsl_matrix_transpose_memcpy(DFT, DF);
    //QR decomposition
    gsl_linalg_QR_decomp (DFT, work);
    gsl_linalg_QR_unpack (DFT, work, Q, R);

    //-------------------------------------------------------------------------------
    //Null vector is the last column of Q
    //-------------------------------------------------------------------------------
    //Sign of the null vector ?
    int sign = 1;
    double dotNV = 0.0;
    if(isFirst)
    {
        //Here, we want to make s_EM[0] "grow"
        sign = gsl_matrix_get(Q, 0, nfv-1) > 0? 1:-1;
    }
    else
    {
        //Here, we want to go "in the same direction for the whole Q vector": no u_turn!
        for(int i = 0; i < nfv-1; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];

        //OR
        //Here, we want to go "in the same direction for somr components of Q
        for(int i = 0; i < 2; i++) dotNV += gsl_matrix_get(Q, i, nfv-1)*nullvector[i];              //CMU of  EML2
        for(int i = 1; i < 4; i++) dotNV += gsl_matrix_get(Q, nfv-i-1, nfv-1)*nullvector[nfv-i-1];  //CMS of SEML2
        dotNV += gsl_matrix_get(Q, nfv-1, nfv-1)*nullvector[nfv-1];                                 //tN
        sign = dotNV > 0? 1:-1;
        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;
    }

    //Null vector is the last column of Q
    for(int i = 0; i < nfv; i++) nullvector[i] = sign*gsl_matrix_get(Q, i, nfv-1);


    cout << "-----------------------------" << endl;
    cout << "Q(end) = " << endl;
    cout << gsl_matrix_get(Q, 0, nfv-1) << endl;
    cout << gsl_matrix_get(Q, 1, nfv-1) << endl;

    cout << gsl_matrix_get(Q, nfv-4, nfv-1) << endl;
    cout << gsl_matrix_get(Q, nfv-3, nfv-1) << endl;
    cout << gsl_matrix_get(Q, nfv-2, nfv-1) << endl;
    cout << gsl_matrix_get(Q, nfv-1, nfv-1) << endl;

    //    //-------------------------------------------------------------------------------
    //    //Reduced Nullvector
    //    //-------------------------------------------------------------------------------
    //    int ncs_r = ncs;
    //    int nfv_r = nfv-1;
    //
    //    //QR elements
    //    gsl_vector *work_r  = gsl_vector_calloc(ncs_r);
    //    gsl_matrix *Q_r       = gsl_matrix_calloc(nfv_r,nfv_r);
    //    gsl_matrix *R_r       = gsl_matrix_calloc(nfv_r,ncs_r);
    //    gsl_matrix *DFT_r     = gsl_matrix_calloc(nfv_r,ncs_r);
    //    gsl_matrix *DF_r      = gsl_matrix_calloc(ncs_r,nfv_r);
    //
    //    //Copy DF into DF_r
    //    for(int i = 0; i < ncs_r; i++)
    //    {
    //        for(int j = 0; j < nfv_r; j++) gsl_matrix_set(DF_r, i, j, gsl_matrix_get(DF, i, j));
    //    }
    //
    //    //DPT = transpose(DP)
    //    gsl_matrix_transpose_memcpy(DFT_r, DF_r);
    //    //QR decomposition
    //    gsl_linalg_QR_decomp (DFT_r, work_r);
    //    gsl_linalg_QR_unpack (DFT_r, work, Q_r, R_r);
    //
    //    //-------------------------------------------------------------------------------
    //    //Display
    //    //-------------------------------------------------------------------------------
    //    cout << "-----------------------------" << endl;
    //    cout << "Q_r(end) = " << endl;
    //    cout << gsl_matrix_get(Q_r, 0, nfv_r-1) << endl;
    //    cout << gsl_matrix_get(Q_r, 1, nfv_r-1) << endl;
    //
    //    cout << gsl_matrix_get(Q_r, nfv_r-3, nfv_r-1) << endl;
    //    cout << gsl_matrix_get(Q_r, nfv_r-2, nfv_r-1) << endl;
    //    cout << gsl_matrix_get(Q_r, nfv_r-1, nfv_r-1) << endl;
    //
    //
    //    //-------------------------------------------------------------------------------
    //    //Null vector is the last column of Q
    //    //-------------------------------------------------------------------------------
    //    if(isFirst)
    //    {
    //        //Here, we want to make s_EM[0] "grow"
    //        sign = gsl_matrix_get(Q_r, 0, nfv_r-1) > 0? 1:-1;
    //    }
    //    else
    //    {
    //        //Here, we want to go "in the same direction for somr components of Q
    //        for(int i = 0; i < 2; i++) dotNV += gsl_matrix_get(Q_r, i, nfv_r-1)*nullvector[i];                //CMU of  EML2
    //        for(int i = 0; i < 3; i++) dotNV += gsl_matrix_get(Q_r, nfv_r-i-1, nfv_r-1)*nullvector[nfv_r-i-1];  //CMS of SEML2
    //        sign = dotNV > 0? 1:-1;
    //        //sign = gsl_matrix_get(Q, nfv-2, nfv-1) < 0? 1:-1;
    //    }
    //
    //    //Null vector is the last column of Q
    //    for(int i = 0; i < nfv_r; i++) nullvector[i] = sign*gsl_matrix_get(Q_r, i, nfv_r-1);
    //    nullvector[nfv - 1] = 0.0;

    //--------------------------------------------------------------------------
    // Free (TO UPDATE!)
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}

/**
 * \brief Multiple shooting scheme with no boundary conditions. The time vector is also variable.
 *        Contrary to multiple_shooting_gomez, no recursive scheme is used to compute the correction vector.
 **/
int multiple_shooting_direct_variable_time(double **ymd, double *tmd, double **ymdn, double *tmdn, double **yma,
        int N, int man_grid_size, int coord_type,
        int isPlotted, int isTimeFixed, gnuplot_ctrl *h1)
{
    //--------------------------------------------------------------------------
    // Init
    //--------------------------------------------------------------------------
    //Cumulated norm of the error
    double normC;
    //Current state along the trajectory
    double **ym  = dmatrix(0, 41, 0, man_grid_size);
    //Current time along the trajectory
    double *tm   = dvector(0, man_grid_size);
    //Various temporary states and times
    double yv[N], ye[N], f[6];

    //-------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //-------------------------------------------------------------------------------
    int dcs = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //--------------------------------------------------------------------------
    // Selection of the vector field
    //--------------------------------------------------------------------------
    int (*vf)(double, const double*, double*, void*) = qbfbp_vfn_novar; //by default, to avoid warning from gcc compiler
    switch(dcs)
    {
    case F_SEM:
    case F_EM:
        vf = qbfbp_vf; //vector field with a state (X, PX)
        break;
    case F_NCSEM:
    case F_NCEM:
        vf = qbfbp_vfn_novar; //vector field with a state (x, px) (default)
        break;
    case F_VNCSEM:
    case F_VNCEM:
        vf = qbfbp_vfn_novar_xv; //vector field with a state (x, vx)
        break;
    }


    //--------------------------------------------------------------------------
    // Check that the focus in SEML is in accordance with the dcs.
    //--------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);

    //-------------------------------------------------------------------------------
    // GSL matrices and vectors
    //-------------------------------------------------------------------------------
    // Correction vector at patch points
    gsl_vector *DQv = gsl_vector_calloc(7*(man_grid_size+1));
    // Error vector at patch points
    gsl_vector *Fv  = gsl_vector_calloc(6*(man_grid_size));
    //Jacobian at patch points
    gsl_matrix **Ji  = gslc_matrix_array_calloc(6, 6, man_grid_size);
    gsl_matrix *DF   = gsl_matrix_calloc(6*man_grid_size, 7*(man_grid_size+1));
    gsl_matrix *M    = gsl_matrix_calloc(6*(man_grid_size), 6*(man_grid_size));

    //Identity matrix eye(6)
    gsl_matrix *Id = gsl_matrix_calloc(6,6);
    gsl_matrix_set_identity (Id);

    //Intermediate variables
    gsl_vector *K3  = gsl_vector_calloc(6*(man_grid_size));
    gsl_vector *Kf  = gsl_vector_calloc(6);
    gsl_vector *K4  = gsl_vector_calloc(6);

    //--------------------------------------------------------------------------
    // Copy the departure state in ymdn
    //--------------------------------------------------------------------------
    for(int k = 0; k <= man_grid_size; k++)
    {
        for(int i = 0; i < N; i++) ymdn[i][k] = ymd[i][k];
        tmdn[k] = tmd[k];
    }


    //----------------------------------------------------------------------
    // Arrival state: yma(1, :) = ymd(1,:);
    //----------------------------------------------------------------------
    for(int i = 0; i < N; i++) yma[i][0] = ymd[i][0];

    //--------------------------------------------------------------------------
    // Loop correction
    //--------------------------------------------------------------------------
    int iter = 0;
    while(iter < 20)
    {

        //----------------------------------------------------------------------
        //Trajectory plot
        //----------------------------------------------------------------------
        if(isPlotted) gnuplot_plot_xyz(h1, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"", "points", "1", "4", 4);

        //----------------------------------------------------------------------
        // Build the Jacobian and other useful matrices
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size-1; k++)
        {
            //------------------------------------------------------------------
            // Integration
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yv[i] = ymdn[i][k];
            ode78_qbcp(ym, tm, tmdn[k], tmdn[k+1], yv, 42, 1, dcs, coord_type, coord_type);

            //------------------------------------------------------------------
            // Final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) ye[i] = ym[i][1];

            //------------------------------------------------------------------
            // Arrival state (including STM) = final position is at the end of ym
            //------------------------------------------------------------------
            for(int i = 0; i < N; i++) yma[i][k+1] = ym[i][1];

            //------------------------------------------------------------------
            // Update the Jacobian
            //------------------------------------------------------------------
            gslc_vectorToMatrix(Ji[k], ye, 6, 6, 6);

            //------------------------------------------------------------------
            // Update the error vector: F[k] = [ye[k] - ymdn[k+1]]
            //------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
            {
                gsl_vector_set(Fv, 6*k+i, ye[i] - ymdn[i][k+1]);
            }


            //--------------------------
            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            //--------------------------
            //Computing f[Q[k], t[k])
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            vf(tmdn[k], yv, f, &SEML);

            //Kf = -f[Q[k], t[k])
            for(int i = 0; i < 6; i++) gsl_vector_set(Kf, i, -f[i]);

            //dF[k]/dt[k] = - Ji[k]*f[Q[k], t[k])
            gsl_blas_dgemv(CblasNoTrans, 1.0, Ji[k], Kf, 0.0, K4);

            //--------------------------
            //dF[k]/dt[k+1] = +f[Q[k+1], t[k+1])
            //--------------------------
            //Computing f[Q[k+1], t[k+1])
            vf(tmdn[k+1], ye, f, &SEML);


            //------------------------------------------------------------------
            // Update DF
            //------------------------------------------------------------------
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
        }

        //------------------------------------------------------------------
        // Norm
        //------------------------------------------------------------------
        normC  = gsl_blas_dnrm2(Fv);

        //------------------------------------------------------------------
        // Update M
        //------------------------------------------------------------------
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DF , DF, 0.0, M);


        //----------------------------------------------------------------------
        // Check that all points are under a given threshold
        //----------------------------------------------------------------------
        cout << "multiple_shooting_direct. nerror = " << normC << endl;
        if(normC < PREC_GSM)
        {
            break;
        }

        //----------------------------------------------------------------------
        // Update correction vector
        //----------------------------------------------------------------------
        //Since M is definite-positive, a Cholesky decomposition can be used, instead of a LU decomposition.
        gsl_linalg_cholesky_decomp (M);
        gsl_linalg_cholesky_solve(M, Fv, K3);
        gsl_blas_dgemv(CblasTrans, -1.0, DF, K3, 0.0, DQv);


        //----------------------------------------------------------------------
        // Update the free variables
        //----------------------------------------------------------------------
        for(int k = 0; k <= man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) ymdn[i][k] += gsl_vector_get(DQv, 7*k+i);
            tmdn[k] += gsl_vector_get(DQv, 7*(k+1)-1);
        }

        //----------------------------------------------------------------------
        // Update number of iterations
        //----------------------------------------------------------------------
        iter++;
    }

    //--------------------------------------------------------------------------
    // Reset the focus in SEML, if necessary
    //--------------------------------------------------------------------------
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);


    //--------------------------------------------------------------------------
    // Free
    //--------------------------------------------------------------------------
    free_dmatrix(ym, 0, 41, 0, man_grid_size);
    free_dvector(tm, 0, man_grid_size);
    gslc_matrix_array_free(Ji , man_grid_size);
    gsl_vector_free(DQv);
    gsl_vector_free(Fv);
    gsl_vector_free(K3);
    gsl_matrix_free(DF);
    gsl_matrix_free(M);
    gsl_matrix_free(Id);


    return GSL_SUCCESS;
}
