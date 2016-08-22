#include "diffcorr.h"

/**
 * \brief Differential correction scheme
 **/
int differential_correction_mns(double *ystart, double *ydest, double t0, double *t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double tof;
    double dx;

    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *DP   = gsl_matrix_calloc(3,4);
    gsl_matrix *Prod = gsl_matrix_calloc(3,3);
    gsl_vector *P    = gsl_vector_calloc(3);
    gsl_vector *P1   = gsl_vector_calloc(3);
    gsl_vector *Pn   = gsl_vector_calloc(4);
    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (3);
    int s;

    //-------------------------------------------------------------------------------
    //Optionnal plotting
    //-------------------------------------------------------------------------------
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");
    double ydot[6];

    //-------------------------------------------------------------------------------
    //Correction loop
    //-------------------------------------------------------------------------------
    tof = *t1;
    do
    {
        //-------------------------------------------------------------------------------
        //Update the starting point (STM (0) = Id included)
        //-------------------------------------------------------------------------------
        for(i=0; i< N; i++) y[i] = ystart[i];

        //-------------------------------------------------------------------------------
        //Integration until t=t1 is reached
        //-------------------------------------------------------------------------------
        t = t0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, N, t0, tof, d, hc, 500, "DiffCorr", "lines", "1", "2", 8);
        //Integration
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , tof, y);

        cout << "Final state in diffcorr: " << endl;
        vector_printf(y, 6);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //-------------------------------------------------------------------------------
        // Update the Jacobian
        //-------------------------------------------------------------------------------
        //Update the matrix DP = DP(x^k), a 2x3 matrix
        gsl_matrix_set(DP, 0, 0, y[5+4]);  //Phi14
        gsl_matrix_set(DP, 0, 1, y[5+5]);  //Phi15
        gsl_matrix_set(DP, 0, 2, y[5+6]);  //Phi16

        gsl_matrix_set(DP, 1, 0, y[5+10]);  //Phi24
        gsl_matrix_set(DP, 1, 1, y[5+11]);  //Phi25
        gsl_matrix_set(DP, 1, 2, y[5+12]);  //Phi26

        gsl_matrix_set(DP, 2, 0, y[5+16]);  //Phi34
        gsl_matrix_set(DP, 2, 1, y[5+17]);  //Phi35
        gsl_matrix_set(DP, 2, 2, y[5+18]);  //Phi36

        //Derivative of y @ y= t12 in ydot
        qbfbp_vfn_novar(tof, y, ydot, d->sys->params);
        gsl_matrix_set(DP, 0, 3, ydot[0]); //xdot
        gsl_matrix_set(DP, 1, 3, ydot[1]); //ydot
        gsl_matrix_set(DP, 2, 3, ydot[2]); //zot


        //-------------------------------------------------------------------------------
        // Update the error vector = [x  - xd, y - yd, z - zd]
        //-------------------------------------------------------------------------------
        gsl_vector_set(P, 0, y[0] - ydest[0]);
        gsl_vector_set(P, 1, y[1] - ydest[1]);
        gsl_vector_set(P, 2, y[2] - ydest[2]);

        //-------------------------------------------------------------------------------
        //Minimum norm solution
        //-------------------------------------------------------------------------------
        //Update the matrix
        //Compute Prod = DP(x^k)*DP(x^k)^T, a 2x2 matrix
        gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, DP , DP, 0.0, Prod);
        //Inverse and product
        gsl_linalg_LU_decomp (Prod, p , &s);
        //P1(x^k) = Prod^{-1}*P(x^k), a 2x1 vector
        gsl_linalg_LU_solve(Prod, p, P, P1);
        //Pn(x^k) = DP(x^k)^T*P1(x^k), a 3x1 vector
        gsl_blas_dgemv (CblasTrans, 1.0, DP, P1, 0.0, Pn);

        //-------------------------------------------------------------------------------
        //Update the state
        //-------------------------------------------------------------------------------
        ystart[3] -= gsl_vector_get(Pn, 0);
        ystart[4] -= gsl_vector_get(Pn, 1);
        ystart[5] -= gsl_vector_get(Pn, 2);
        tof       -= gsl_vector_get(Pn, 3);

        //-------------------------------------------------------------------------------
        //Norm of the correction
        //-------------------------------------------------------------------------------
        dx = gsl_blas_dnrm2(P);

        cout << iter << "  dx = " <<  fabs(dx) << endl;
    }
    while((fabs(dx)> eps_diff) && (++iter) < itermax);


    //-------------------------------------------------------------------------------
    //End the computation if there was no convergence
    //-------------------------------------------------------------------------------
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dx));
        return GSL_FAILURE;
    }
    cout << "Time gap = " << *t1 - tof << endl;
    *t1 = tof;

    //-------------------------------------------------------------------------------
    //Free memory
    //-------------------------------------------------------------------------------
    gsl_matrix_free(DP);
    gsl_matrix_free(Prod);
    gsl_vector_free(P);
    gsl_vector_free(P1);
    gsl_vector_free(Pn);
    gnuplot_close(hc);


    return GSL_SUCCESS;
}

/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_ft(double *ystart, double *ydest, double t0, double *t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted)
{
    //-------------------------------------------------------------------------------
    //Init
    //-------------------------------------------------------------------------------
    double y[N];
    int i;
    int iter = 0;
    int status;
    int itermax = 70;

    double t;
    double tof;
    double dx;

    //-------------------------------------------------------------------------------
    //Gsl matrices and vectors
    //-------------------------------------------------------------------------------
    gsl_matrix *DP   = gsl_matrix_calloc(3,3);
    gsl_vector *P    = gsl_vector_calloc(3);
    gsl_vector *Pn   = gsl_vector_calloc(3);
    //For GSL inversion
    gsl_permutation * p = gsl_permutation_alloc (3);
    int s;

    //-------------------------------------------------------------------------------
    //Optionnal plotting
    //-------------------------------------------------------------------------------
    gnuplot_ctrl  *hc;
    hc = gnuplot_init();
    gnuplot_setstyle(hc, (char*)"lines");
    gnuplot_set_xlabel(hc, (char*)"x [-]");
    gnuplot_set_ylabel(hc, (char*)"y [-]");

    //-------------------------------------------------------------------------------
    //Correction loop
    //-------------------------------------------------------------------------------
    tof = *t1;
    do
    {
        //-------------------------------------------------------------------------------
        //Update the starting point (STM (0) = Id included)
        //-------------------------------------------------------------------------------
        for(i=0; i< N; i++) y[i] = ystart[i];

        //-------------------------------------------------------------------------------
        //Integration until t=t1 is reached
        //-------------------------------------------------------------------------------
        t = t0;
        //Optionnal plotting
        if(isPlotted) odePlotGen(y, N, t0, tof, d, hc, 500, "DiffCorr", "lines", "1", "2", 8);
        //Integration
        gsl_odeiv2_driver_reset(d);
        status = gsl_odeiv2_driver_apply (d, &t , tof, y);

        if(status != GSL_SUCCESS)
        {
            printf("WARNING: GSL driver failed to converge in differential_correction. Premature ending.\n");
            return GSL_FAILURE;
        }

        //-------------------------------------------------------------------------------
        // Update the Jacobian
        //-------------------------------------------------------------------------------
        //Update the matrix DP = DP(x^k), a 2x3 matrix
        gsl_matrix_set(DP, 0, 0, y[5+4]);  //Phi14
        gsl_matrix_set(DP, 0, 1, y[5+5]);  //Phi15
        gsl_matrix_set(DP, 0, 2, y[5+6]);  //Phi16

        gsl_matrix_set(DP, 1, 0, y[5+10]);  //Phi24
        gsl_matrix_set(DP, 1, 1, y[5+11]);  //Phi25
        gsl_matrix_set(DP, 1, 2, y[5+12]);  //Phi26

        gsl_matrix_set(DP, 2, 0, y[5+16]);  //Phi34
        gsl_matrix_set(DP, 2, 1, y[5+17]);  //Phi35
        gsl_matrix_set(DP, 2, 2, y[5+18]);  //Phi36


        //-------------------------------------------------------------------------------
        // Update the error vector = [x  - xd, y - yd, z - zd]
        //-------------------------------------------------------------------------------
        gsl_vector_set(P, 0, y[0] - ydest[0]);
        gsl_vector_set(P, 1, y[1] - ydest[1]);
        gsl_vector_set(P, 2, y[2] - ydest[2]);


        //-------------------------------------------------------------------------------
        //Minimum norm solution
        //-------------------------------------------------------------------------------
        //Inverse and product
        gsl_linalg_LU_decomp (DP, p , &s);
        //Pn(x^k) = DP^{-1}*P(x^k), a 2x1 vector
        gsl_linalg_LU_solve(DP, p, P, Pn);

        //-------------------------------------------------------------------------------
        //Update the state
        //-------------------------------------------------------------------------------
        ystart[3] -= gsl_vector_get(Pn, 0);
        ystart[4] -= gsl_vector_get(Pn, 1);
        ystart[5] -= gsl_vector_get(Pn, 2);

        //-------------------------------------------------------------------------------
        //Norm of the correction
        //-------------------------------------------------------------------------------
        dx = gsl_blas_dnrm2(P);
        cout << iter << "  dx = " <<  fabs(dx) << endl;
    }
    while((fabs(dx)> eps_diff) && (++iter) < itermax);


    //-------------------------------------------------------------------------------
    //End the computation if there was no convergence
    //-------------------------------------------------------------------------------
    if(iter>=itermax)
    {
        printf("WARNING: number of iter max exceeded during differential correction. Final precision is around %+5.2e. Premature ending\n", fabs(dx));
        return GSL_FAILURE;
    }
    cout << "Time gap = " << *t1 - tof << endl;
    *t1 = tof;

    //-------------------------------------------------------------------------------
    //Free memory
    //-------------------------------------------------------------------------------
    gsl_matrix_free(DP);
    gsl_vector_free(P);
    gsl_vector_free(Pn);
    gnuplot_close(hc);


    return GSL_SUCCESS;
}


/**
 *  \brief Same as odePlot, but with more flexibility on the parameters: can choose the title, the line style, the line type and the line color.
 **/
int odePlotGen(const double y[],
               int N,
               double t0,
               double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc)
{
    gsl_odeiv2_driver_reset(d);

    double xc[Npoints];
    double yc[Npoints];

    //Initial conditions
    double ys[N];
    for(int i=0; i<N; i++) ys[i] = y[i];
    double ti = 0;

    xc[0] = ys[0];
    yc[0] = ys[1];
    double t = t0;
    if(t1 <0) d->h = -d->h;
    for(int i =1; i<= Npoints; i++)
    {
        ti = t0 + i * (t1 - t0) / Npoints;
        gsl_odeiv2_driver_apply (d, &t, ti, ys);
        xc[i] = ys[0];
        yc[i] = ys[1];
    }
    if(t1 <0) d->h = -d->h;

    //------------------------------------------
    //Plot
    //------------------------------------------
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xy(h1, xc, yc, Npoints, title, ls, lt, lw, lc);

    return GSL_SUCCESS;
}


