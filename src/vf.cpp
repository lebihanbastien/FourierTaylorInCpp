#include "vf.h"

extern "C" {
#include "gnuplot_i.h"
}


//============================================================================================================================================
//
// Integration without STM
//
//============================================================================================================================================

//----------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------
//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc      = qbp->numberOfCoefs;
    int nf       = qbp->nf;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return 0;
}

//Update the normalized vector field of the state
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma)
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0]*y[3] + alpha[1]*y[0] + alpha[2]*y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0]*y[4] + alpha[1]*y[1] - alpha[2]*y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0]*y[5] + alpha[1]*y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1]*y[3] + alpha[2]*y[4] + alpha[14]*y[0]
           + alpha[12];
    if(me != 0) f[3] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);

    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1]*y[4] - alpha[2]*y[3] + alpha[14]*y[1]
           + alpha[13];
    if(me != 0) f[4] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);


    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1]*y[5] + alpha[14]*y[2];
    if(me != 0) f[5] += - alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] += - alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] += - alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);

    return 0;
}


//----------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------
//Update the normalized vector field of the state. The state is here (x, xdot), and not (x, px)
int vfn_state_xv(const double y[], double f[], double alpha[], double alphad[],
                 double ps[], double pe[], double pm[],
                 double qps2, double qpe2, double qpm2,
                 double ms, double me, double mm,
                 double gamma)
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = y[3];
    //f[1] = y'
    //------------------------------------
    f[1] = y[4];
    //f[2] = z'
    //------------------------------------
    f[2] = y[5];

    //-------------------------------------------------------------------------------
    //Intermediate coefficients
    //-------------------------------------------------------------------------------
    double alpha16 = alphad[1] - alphad[0]*alpha[1]/alpha[0] + alpha[1]*alpha[1] + alpha[2]*alpha[2]; //
    double alpha17 = alphad[2] - alphad[0]*alpha[2]/alpha[0]; //
    double alpha18 = alphad[0]/alpha[0]; //

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x'', y'', z''
    //-------------------------------------------------------------------------------
    //f[3] = x''
    //------------------------------------
    f[3] = alpha16*y[0] + alpha17*y[1] + alpha18*y[3] + 2*alpha[2]*y[4] + alpha[0]*alpha[12];
    if(me != 0) f[3] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);

    //f[4] = y''
    //------------------------------------
    f[4] = - alpha17*y[0] + alpha16*y[1] - 2*alpha[2]*y[3] + alpha18*y[4] + alpha[0]*alpha[13];
    if(me != 0) f[4] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);


    //f[5] = z''
    //------------------------------------
    f[5] = (alpha16 - alpha[2]*alpha[2])*y[2] + alpha18*y[5];
    if(me != 0) f[5] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] -=  alpha[0]*alpha[5]/pow(gamma,3.0)*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);

    return 0;
}

//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations). The state is here (x, xdot), and not (x, px)
int qbfbp_vfn_novar_xv(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc      = qbp->numberOfCoefs;
    int nf       = qbp->nf;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @ t, as well as the derivatives of the first 3 alphas
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, nf, qbp->cs.coeffs, noc);
    double alphad[3];
    evaluateCoefDerivatives(alphad, t, n, nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, nf, qbp->cs.pm, 3);


    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state_xv(y, f, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return 0;
}

//----------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbfbp_vf(double t, const double y[], double f[], void *params_void)
{

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vf_state(y, f, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    return GSL_SUCCESS;
}

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized)
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm )
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0]*y[3] + alpha[1]*y[0] + alpha[2]*y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0]*y[4] + alpha[1]*y[1] - alpha[2]*y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0]*y[5] + alpha[1]*y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1]*y[3] + alpha[2]*y[4] + alpha[14]*y[0]
           - alpha[3];
    if(me != 0) f[3] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[0] - pe[0]);
    if(mm != 0) f[3] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[0] - pm[0]);
    if(ms != 0) f[3] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[0] - ps[0]);


    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1]*y[4] - alpha[2]*y[3] + alpha[14]*y[1]
           - alpha[4];
    if(me != 0) f[4] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[1] - pe[1]);
    if(mm != 0) f[4] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[1] - pm[1]);
    if(ms != 0) f[4] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[1] - ps[1]);

    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1]*y[5] + alpha[14]*y[2];
    if(me != 0) f[5] += - alpha[5]*me/pow(qpe2, 3.0/2) * (y[2] - pe[2]);
    if(mm != 0) f[5] += - alpha[5]*mm/pow(qpm2, 3.0/2) * (y[2] - pm[2]);
    if(ms != 0) f[5] += - alpha[5]*ms/pow(qps2, 3.0/2) * (y[2] - ps[2]);


    return GSL_SUCCESS;
}



//============================================================================================================================================
//
// Integration with STM
//
//============================================================================================================================================
/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j)
{
    if(mc != 0.0)
    {
        if(i == j) return factor*(3*mc/pow(qpc2,5.0/2)*pow(y[i]-pc[i],2.0) - mc/pow(qpc2,3.0/2));
        else return factor*3*mc/pow(qpc2,5.0/2)*(y[i]-pc[i])*(y[j]-pc[j]);
    }
    else return 0.0;
}

//----------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha
 **/
void set_vareq_matrix(gsl_matrix *Q, double b[], double alpha[])
{
    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set(Q, 0, 0,  alpha[1]);
    gsl_matrix_set(Q, 0, 1,  alpha[2]);
    gsl_matrix_set(Q, 0, 3,  alpha[0]);
    gsl_matrix_set(Q, 1, 0, -alpha[2]);
    gsl_matrix_set(Q, 1, 1,  alpha[1]);
    gsl_matrix_set(Q, 1, 4,  alpha[0]);
    gsl_matrix_set(Q, 2, 2,  alpha[1]);
    gsl_matrix_set(Q, 2, 5,  alpha[0]);
    gsl_matrix_set(Q, 3, 0,  b[0]+alpha[14]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]+alpha[14]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]+alpha[14]);
    gsl_matrix_set(Q, 3, 3, -alpha[1]);
    gsl_matrix_set(Q, 3, 4,  alpha[2]);
    gsl_matrix_set(Q, 4, 3, -alpha[2]);
    gsl_matrix_set(Q, 4, 4, -alpha[1]);
    gsl_matrix_set(Q, 5, 5, -alpha[1]);
}

/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm(const double y[], gsl_matrix *Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[5]/pow(gamma,3.0);

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
             + Uij(y, pm, qpm2, mm, factor, 0, 0)
             + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *Q    = gsl_matrix_calloc(6,6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp  = (QBCP_L *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, qbp->nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, qbp->nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, qbp->nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);
    vfn_stm(y, Q, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha
 **/
void set_vareq_matrix_xv(gsl_matrix *Q, double b[], double alpha[], double alphad[])
{
    //-------------------------------------------------------------------------------
    //Intermediate coefficients
    //-------------------------------------------------------------------------------
    double alpha16 = alphad[1] - alphad[0]*alpha[1]/alpha[0] + alpha[1]*alpha[1] + alpha[2]*alpha[2]; //
    double alpha17 = alphad[2] - alphad[0]*alpha[2]/alpha[0]; //
    double alpha18 = alphad[0]/alpha[0]; //

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0]+alpha16);
    gsl_matrix_set(Q, 3, 1,  b[1]+alpha17);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]-alpha17);
    gsl_matrix_set(Q, 4, 1,  b[2]+alpha16);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]+alpha16-alpha[2]*alpha[2]);

    gsl_matrix_set(Q, 3, 3,  alpha18);
    gsl_matrix_set(Q, 3, 4,  2*alpha[2]);
    gsl_matrix_set(Q, 4, 3, -2*alpha[2]);
    gsl_matrix_set(Q, 4, 4,  alpha18);
    gsl_matrix_set(Q, 5, 5,  alpha18);
}



/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm_xv(const double y[], gsl_matrix *Q, double alpha[], double alphad[],
               double ps[], double pe[], double pm[],
               double qps2, double qpe2, double qpm2,
               double ms, double me, double mm,
               double gamma)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[0]*alpha[5]/pow(gamma,3.0);

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
             + Uij(y, pm, qpm2, mm, factor, 0, 0)
             + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix_xv(Q, b, alpha, alphad);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates. The state is here (x, xdot), and not (x, px)
 **/
int qbfbp_vfn_varnonlin_xv(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *Q    = gsl_matrix_calloc(6,6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp  = (QBCP_L *) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);
    double alphad[3];
    evaluateCoefDerivatives(alphad, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, qbp->nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, qbp->nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, qbp->nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0]-pe[0])*(y[0]-pe[0]) + (y[1]-pe[1])*(y[1]-pe[1]) + (y[2]-pe[2])*(y[2]-pe[2]);
    double qps2 = (y[0]-ps[0])*(y[0]-ps[0]) + (y[1]-ps[1])*(y[1]-ps[1]) + (y[2]-ps[2])*(y[2]-ps[2]);
    double qpm2 = (y[0]-pm[0])*(y[0]-pm[0]) + (y[1]-pm[1])*(y[1]-pm[1]) + (y[2]-pm[2])*(y[2]-pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state_xv(y, f, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);
    vfn_stm_xv(y, Q, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double y[], gsl_matrix *Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[5];

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
             + Uij(y, pm, qpm2, mm, factor, 0, 0)
             + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}


/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame. Variational equations included.
 **/
int qbfbp_vf_varnonlin(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix *STM  = gsl_matrix_calloc(6,6);
    gsl_matrix *STMd = gsl_matrix_calloc(6,6);
    gsl_matrix *B    = gsl_matrix_calloc(6,6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;
    int noc   = qbp->numberOfCoefs;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vf_state(y, f, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(B);
    vf_stm(y, B, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = B * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, STM, 0.0, STMd);

    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 ,6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(B);

    return GSL_SUCCESS;
}


//=======================================================================================================================================
//
//          Integration of the pm
//
//=======================================================================================================================================
/**
 *  \brief Computes the reduce vector field dot(s) = f(s,t)
 *  \param t the current time
 *  \param y the complex state s (CCM), with real and imag parts separated (CCM8 form)
 *  \param f the complex reduced vector field f(s,t), with real and imag parts separated (CCM8 form)
 *  \param params_void a (void) pointer towards a RVF structure that contains:
 *           (1) the Fourier-Taylor expansions of the rvf
 *           (2) the order of the evaluation
 *           (3) the pulsation n of the system
 *           (4) an OFS temp variable.
 *
 *     Note that the use of void *params_void instead of the direct use of an RVF structure is to comply with GSL ODE integration tools,
 *     that require a vector field of the exact form : int vf(double t, const double y[], double f[], void *params_void)
 *
 **/
int qbfbp_fh(double t, const double y[], double f[], void *params_void)
{
    //-------------------------------
    //Initialization
    //-------------------------------
    RVF *rvf = (RVF*) params_void;

    //-------------------------------
    //Evaluation of the reduced vector field at order rvf->order
    //-------------------------------
    CCM8toRVF8(y, t, rvf->n, rvf->order, rvf->ofs_order, rvf->reduced_nv, *rvf->fh, *rvf->ofs, f);

    return GSL_SUCCESS;
}




//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Hamiltonians
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian of the QBCP with SYS units and SYS coordinates
 **/
double qbfbp_H(double t, const double y[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //-------------------------------------------------------------------------------
    // Hamiltonian
    //-------------------------------------------------------------------------------
    double H = 0.5*alpha[0]*(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]) + alpha[1]*(y[3]*y[0] + y[4]*y[1] + y[5]*y[2])
               + alpha[2]*(y[3]*y[1] - y[4]*y[0])
               + alpha[3]*y[0] + alpha[4]*y[1]
               - 0.5*alpha[14]*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2])
               - alpha[5]*( me/pow(qPe2, 1.0/2)
                            + mm/pow(qPm2, 1.0/2)
                            + ms/pow(qPs2, 1.0/2) );

    return H;
}


/**
 *  \brief Hamiltonian of the QBCP with SEM units and SEM coordinates
 **/
double qbfbp_H_SEM(double t, const double y[], void *params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L *) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us_sem.ms;
    double me = qbp->us_sem.me;
    double mm = qbp->us_sem.mm;
    double n  = qbp->us_sem.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs_sem.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs_sem.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs_sem.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0]-Pe[0])*(y[0]-Pe[0]) + (y[1]-Pe[1])*(y[1]-Pe[1]) + (y[2]-Pe[2])*(y[2]-Pe[2]);
    double qPs2 = (y[0]-Ps[0])*(y[0]-Ps[0]) + (y[1]-Ps[1])*(y[1]-Ps[1]) + (y[2]-Ps[2])*(y[2]-Ps[2]);
    double qPm2 = (y[0]-Pm[0])*(y[0]-Pm[0]) + (y[1]-Pm[1])*(y[1]-Pm[1]) + (y[2]-Pm[2])*(y[2]-Pm[2]);

    //-------------------------------------------------------------------------------
    // Hamiltonian
    //-------------------------------------------------------------------------------
    double H = 0.5*alpha[0]*(y[3]*y[3] + y[4]*y[4] + y[5]*y[5]) + alpha[1]*(y[3]*y[0] + y[4]*y[1] + y[5]*y[2])
               + alpha[2]*(y[3]*y[1] - y[4]*y[0])
               + alpha[3]*y[0] + alpha[4]*y[1]
               - 0.5*alpha[14]*(y[0]*y[0] + y[1]*y[1] + y[2]*y[2])
               - alpha[5]*( me/pow(qPe2, 1.0/2)
                            + mm/pow(qPm2, 1.0/2)
                            + ms/pow(qPs2, 1.0/2) );

    return H;
}





//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Testing the QBTBP
//
//--------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Integrating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Derivatives of the QBTBP. To plug into GSL integrator.
 *
 * Note: the use of gsl livrary forces us to use double variables
 * As a consequence, in the case of z and Z, we need to use real and imaginary parts
 * as separate variables
 */
int qbtbp_derivatives(double t, const double y[], double f[], void *params)
{
    //Parameters for the qbtbp
    double mu = * (double * ) params;
    double ms = * ((double * ) (params)+1);

    //Reconstruction of z and Z
    cdouble z = y[0] + I*y[1];
    cdouble Z = y[2] + I*y[3];

    cdouble temp1 = Z-mu*z;
    temp1 = temp1/pow(cabs(temp1), 3.0);

    cdouble temp2 = Z+(1-mu)*z;
    temp2 = temp2/pow(cabs(temp2), 3.0);

    cdouble zdd = +0.0*I-z/pow(cabs(z), 3.0) + ms*(temp1-temp2);
    cdouble Zdd = -(1+ms)*(mu*temp2 + (1-mu)*temp1);

    //-------------------------------------------------------------------------------
    //Phase space derivatives
    //-------------------------------------------------------------------------------
    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];
    f[4] = creal(zdd);
    f[5] = cimag(zdd);
    f[6] = creal(Zdd);
    f[7] = cimag(Zdd);

    return GSL_SUCCESS;
}

/**
 * \brief Computes the semi-analytical position of the Earth, the Moon, and the Sun from the coefficients alpha and delta of the EM and SEM vector field.
 **/
void primary_analytical_position(int prim, double tem, double yEM_to_IN[], double ySEM_to_IN[], QBCP_L &qbcp_l)
{

    //====================================================================
    // 1. Init
    //====================================================================
    int noc     = qbcp_l.numberOfCoefs;
    int nf      = qbcp_l.nf;
    double tsem = tem*qbcp_l.us_em.ns;     //new tem in SEM units

    //====================================================================
    // 2. Evaluate the Fourier coefficients
    //====================================================================
    //-------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //-------------------------------------------------------------------------------
    double delta[noc];
    evaluateCoef(delta, tsem, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);
    double deltad[noc];
    evaluateCoefDerivatives(deltad, tsem, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, tem, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);
    double alphad[noc];
    evaluateCoefDerivatives(alphad, tem, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);


    //====================================================================
    // 3. Evaluate the semi-analytical position of the Primary
    //====================================================================
    double ySEM[6], yEM[6];
    switch(prim)
    {
    case EARTH:
        //Earth position & momenta in SEM ref and SEM units
        ySEM[0] = delta[8];
        ySEM[1] = delta[9];
        ySEM[2] = 0.0;
        ySEM[3] = deltad[8];
        ySEM[4] = deltad[9];
        ySEM[5] = 0.0;
        //Earth position & momenta in EM ref and EM units
        yEM[0] = qbcp_l.us_em.mu_EM;
        yEM[1] = 0.0;
        yEM[2] = 0.0;
        yEM[3] = 0.0;
        yEM[4] = 0.0;
        yEM[5] = 0.0;
        break;

    case MOON:
        //Moon position & velocity in SEM ref and SEM units
        ySEM[0] = delta[10];
        ySEM[1] = delta[11];
        ySEM[2] = 0.0;
        ySEM[3] = deltad[10];
        ySEM[4] = deltad[11];
        ySEM[5] = 0.0;
        //Moon position & velocity in EM ref and EM units
        yEM[0] = qbcp_l.us_em.mu_EM-1;
        yEM[1] = 0.0;
        yEM[2] = 0.0;
        yEM[3] = 0.0;
        yEM[4] = 0.0;
        yEM[5] = 0.0;
        break;

    case SUN:
        //Sun position & velocity in SEM ref and SEM units
        ySEM[0] = qbcp_l.us_sem.mu_SEM;
        ySEM[1] = 0.0;
        ySEM[2] = 0.0;
        ySEM[3] = 0.0;
        ySEM[4] = 0.0;
        ySEM[5] = 0.0;
        //Sun position & velocity in EM ref and EM units
        yEM[0] = alpha[6];
        yEM[1] = alpha[7];
        yEM[2] = 0.0;
        yEM[3] = alphad[6];
        yEM[4] = alphad[7];
        yEM[5] = 0.0;

        //            cout << "t*qbcp_l.us_em.n/(2*pi) = " << tem*qbcp_l.us_em.n/(2*M_PI) << endl;
        //            cout << "angle/(2*pi)  = " << atan2 (yEM[1],yEM[0])/(2*M_PI) << endl;
        //            cout << "Sun position = " << endl;
        //            vector_printf(yEM, 3);

        break;

    default:
        perror("Unknown primary");
    }



    //====================================================================
    // 4. Back in inertial coordinates,
    //    still with the use of the semi-analytical expressions
    //====================================================================
    double ytp[6];
    qbcp_coc(tsem, ySEM, ytp, VSEM, VEM);
    EMvtoIN(tem, ytp, ySEM_to_IN, &qbcp_l);
    EMvtoIN(tem, yEM, yEM_to_IN,  &qbcp_l);
}

/**
 * \brief Computes the "true" integrated position of the Earth, the Moon, and the Sun from the integrated bicircular movement contained in Z(t) and z(t).
 **/
void primary_integrated_position(int prim, double yIN[], cdouble z0, cdouble Z0, cdouble zdot0, cdouble Zdot0, QBCP_L &qbcp_l)
{
    //Initialization
    USYS us_em  = qbcp_l.us_em;

    //====================================================================
    // 1. Evaluate the integrated position of the Primary
    //====================================================================
    switch(prim)
    {
    case EARTH:
        yIN[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) + us_em.mu_EM*creal(z0);
        yIN[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) + us_em.mu_EM*cimag(z0);
        yIN[2] = 0.0;
        yIN[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) + us_em.mu_EM*creal(zdot0);
        yIN[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) + us_em.mu_EM*cimag(zdot0);
        yIN[5] = 0.0;
        break;

    case MOON:
        yIN[0] = -us_em.ms/(1.0+us_em.ms)*creal(Z0) - (1-us_em.mu_EM)*creal(z0);
        yIN[1] = -us_em.ms/(1.0+us_em.ms)*cimag(Z0) - (1-us_em.mu_EM)*cimag(z0);
        yIN[2] = 0.0;
        yIN[3] = -us_em.ms/(1.0+us_em.ms)*creal(Zdot0) - (1-us_em.mu_EM)*creal(zdot0);
        yIN[4] = -us_em.ms/(1.0+us_em.ms)*cimag(Zdot0) - (1-us_em.mu_EM)*cimag(zdot0);
        yIN[5] = 0.0;
        break;

    case SUN:
        yIN[0] = 1.0/(1.0+us_em.ms)*creal(Z0);
        yIN[1] = 1.0/(1.0+us_em.ms)*cimag(Z0);
        yIN[2] = 0.0;
        yIN[3] = 1.0/(1.0+us_em.ms)*creal(Zdot0);
        yIN[4] = 1.0/(1.0+us_em.ms)*cimag(Zdot0);
        yIN[5] = 0.0;
        break;

    default:
        perror("Unknown primary");
    }

}

/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, QBCP_L &qbcp_l)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                 Test of the                       " << endl;
    cout << "   Quasi-Bicircular Three-Body Problem (QBTBP)     " << endl;
    cout << "---------------------------------------------------" << endl;
    //Param
    double n   = qbcp_l.us_em.n;
    double ns  = qbcp_l.us_em.ns;
    double as  = qbcp_l.us_em.as;
    double ni  = qbcp_l.us_em.ni;
    double ai  = qbcp_l.us_em.ai;
    int noc    = qbcp_l.numberOfCoefs;
    int nf     = qbcp_l.nf;

    //====================================================================
    // 0. Gnuplot
    //====================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h1;
    h1 = gnuplot_init();
    gnuplot_cmd(h1, "set title \"Variations of the error in the Primaries' state vector: EM vs IN\" ");
    gnuplot_cmd(h1, "set grid");
    gnuplot_set_xlabel(h1, (char*)"t [EM units]");
    gnuplot_set_ylabel(h1, (char*)"Error [x Tsem]");
    gnuplot_cmd(h1, "set logscale y");
    gnuplot_cmd(h1, "set format y \"1e\%%L\"");


    gnuplot_ctrl *h2;
    h2 = gnuplot_init();
    gnuplot_cmd(h2, "set title \"Variations of the error in the Primaries' state vector: SEM vs IN\" ");
    gnuplot_cmd(h2, "set grid");
    gnuplot_set_xlabel(h2, (char*)"t [EM units]");
    gnuplot_set_ylabel(h2, (char*)"Error [x Tsem]");
    gnuplot_cmd(h2, "set logscale y");
    gnuplot_cmd(h2, "set format y \"1e\%%L\"");


    //====================================================================
    // 1. Init the integration tools
    //====================================================================
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    OdeStruct ode_s;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method
    //Ode solver parameters
    double param[2];
    param[0] = qbcp_l.us_em.mu_EM; //note that the qbtbp is computed in EM framework
    param[1] = qbcp_l.us_em.ms;    //note that the qbtbp is computed in EM framework
    //General structures
    init_ode_structure(&ode_s, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF,  8, PREC_HSTART,  qbtbp_derivatives, NULL, param);

    //====================================================================
    // 2. Initial conditions
    //====================================================================
    double t0  = 0.0;             //new t0 in EM units
    double t0c = t0*qbcp_l.us_em.ns;     //new t0 in SEM units

    //z(0) and Z(0)
    cdouble z0    = evz(qbcp_l.cs_em.zt, t0, n, ni, ai);
    cdouble Z0    = evz(qbcp_l.cs_em.Zt, t0, n, ns, as);
    cdouble zdot0 = evzdot(qbcp_l.cs_em.zt, qbcp_l.cs_em.ztdot, t0, n, ni, ai);
    cdouble Zdot0 = evzdot(qbcp_l.cs_em.Zt, qbcp_l.cs_em.Ztdot, t0, n, ns, as);

    //Put the IC in real form
    double yv[8];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << "Initial positions: " << endl;
    cout << "Internal motion: z(t=0.0) = " << creal(z0) << " + " << cimag(z0) << "i" <<  endl;
    cout << "External motion: Z(t=0.0) = " << creal(Z0) << " + " << cimag(Z0) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 2.2 Same for the Primaries
    //====================================================================
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    double yEM_to_IN[6], ySEM_to_IN[6], yIN[6];

    //-------------------------------------------------------------------------------
    //Initial Earth position
    //-------------------------------------------------------------------------------
    primary_analytical_position(EARTH, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(EARTH, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Earth position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    cout << "-------------------------------------------" << endl;

    //-------------------------------------------------------------------------------
    //Initial Moon position
    //-------------------------------------------------------------------------------
    primary_analytical_position(MOON, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(MOON, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Moon position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    cout << "-------------------------------------------" << endl;


    //-------------------------------------------------------------------------------
    //Initial Sun position
    //-------------------------------------------------------------------------------
    primary_analytical_position(SUN, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(SUN, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Sun position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 3. Initial conditions of the primaries
    //====================================================================
    //-------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //-------------------------------------------------------------------------------
    double delta[noc];
    evaluateCoef(delta, t0c, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);
    double deltad[noc];
    evaluateCoefDerivatives(deltad, t0c, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t0, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);
    double alphad[noc];
    evaluateCoefDerivatives(alphad, t0, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);


    //====================================================================
    //Loop
    //====================================================================
    int Nplot = 200;
    double t = t0, ti;
    double tvec[Nplot+1], yEM_vs_IN_E[Nplot+1], yEM_vs_IN_M[Nplot+1], yEM_vs_IN_S[Nplot+1];
    double ySEM_vs_IN_E[Nplot+1], ySEM_vs_IN_M[Nplot+1], ySEM_vs_IN_S[Nplot+1];

    for(int k = 0; k <= Nplot; k++)
    {
        //-------------------------------------------------------------------------------
        //Integration from t to ti
        //-------------------------------------------------------------------------------
        ti = t0 + (t1 - t0)*k/Nplot;
        gsl_odeiv2_driver_apply (ode_s.d, &t , ti, yv);
        tvec[k] = ti/qbcp_l.us_em.T;

        //-------------------------------------------------------------------------------
        //Current Earth position
        //-------------------------------------------------------------------------------
        primary_analytical_position(EARTH, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(EARTH, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);

        yEM_vs_IN_E[k]  = ENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_E[k] = ENorm(ySEM_to_IN, yIN, 6);

        //-------------------------------------------------------------------------------
        //Current Moon position
        //-------------------------------------------------------------------------------
        primary_analytical_position(MOON, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(MOON, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);

        yEM_vs_IN_M[k]  = ENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_M[k] = ENorm(ySEM_to_IN, yIN, 6);

        //-------------------------------------------------------------------------------
        //Current Sun position
        //-------------------------------------------------------------------------------
        primary_analytical_position(SUN, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(SUN, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);

        yEM_vs_IN_S[k]  = ENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_S[k] = ENorm(ySEM_to_IN, yIN, 6);
    }

    //-------------------------------------------------------------------------------
    //Plot: EM vs IN
    //-------------------------------------------------------------------------------
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_E, Nplot+1, (char*)"Earth", "lines", "3", "2", 3);
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_M, Nplot+1, (char*)"Moon", "lines", "3", "2", 4);
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_S, Nplot+1, (char*)"Sun", "lines", "3", "2", 5);

    gnuplot_fplot_xy(tvec, yEM_vs_IN_E, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_EM_vs_IN_E.txt").c_str());
    gnuplot_fplot_xy(tvec, yEM_vs_IN_M, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_EM_vs_IN_M.txt").c_str());
    gnuplot_fplot_xy(tvec, yEM_vs_IN_S, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_EM_vs_IN_S.txt").c_str());


    //-------------------------------------------------------------------------------
    //Plot: SEM vs IN
    //-------------------------------------------------------------------------------
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_E, Nplot+1, (char*)"Earth", "lines", "3", "2", 3);
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_M, Nplot+1, (char*)"Moon", "lines", "3", "2", 4);
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_S, Nplot+1, (char*)"Sun", "lines", "3", "2", 5);

    gnuplot_fplot_xy(tvec, ySEM_vs_IN_E, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_SEM_vs_IN_E.txt").c_str());
    gnuplot_fplot_xy(tvec, ySEM_vs_IN_M, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_SEM_vs_IN_M.txt").c_str());
    gnuplot_fplot_xy(tvec, ySEM_vs_IN_S, Nplot+1, (char*) (qbcp_l.cs.F_PLOT+"QBTBP_SEM_vs_IN_S.txt").c_str());

    //====================================================================
    //Final state
    //====================================================================
    cout << "-------------------------------------------" << endl;
    cout << "End of integration." << endl;
    cout << "Final t: " << t << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Results from integration" << endl;
    cout << "Internal motion: z(t=t1) = " << yv[0] << " + " << yv[1] << "i" << endl;
    cout << "External motion: Z(t=t1) = " << yv[2] << " + " << yv[3] << "i" << endl;
    cout << "-------------------------------------------" << endl;

    //Analytical final state
    cdouble zfinal = evz(qbcp_l.cs_em.zt, t1, n, ni, ai);
    cdouble Zfinal = evz(qbcp_l.cs_em.Zt, t1, n, ns, as);

    cdouble zdotfinal = evzdot(qbcp_l.cs_em.zt, qbcp_l.cs_em.ztdot, t1, n, ni, ai);
    cdouble Zdotfinal = evzdot(qbcp_l.cs_em.Zt, qbcp_l.cs_em.Ztdot, t1, n, ns, as);


    cout << "Analytical results" << endl;
    cout << "Internal motion: z(t=t1) = " << creal(zfinal) << " + " << cimag(zfinal) << "i" <<  endl;
    cout << "External motion: Z(t=t1) = " << creal(Zfinal) << " + " << cimag(Zfinal) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;


    cout << "Absolute delta between analytical and numerical results: " << endl;
    cout << "dz    = " << cabs(zfinal - yv[0]-I*yv[1])    << endl;
    cout << "dzdot = " << cabs(zdotfinal - yv[4]-I*yv[5]) << endl;
    cout << "dZ    = " << cabs(Zfinal - yv[2]-I*yv[3])    << endl;
    cout << "dZdot = " << cabs(Zdotfinal - yv[6]-I*yv[7]) << endl;
    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 2.2 Same for the Primaries
    //====================================================================
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);

    //-------------------------------------------------------------------------------
    //Current Earth position
    //-------------------------------------------------------------------------------
    primary_analytical_position(EARTH, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(EARTH, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);
    cout << "Final Earth position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    cout << "-------------------------------------------" << endl;

    //-------------------------------------------------------------------------------
    //Current Moon position
    //-------------------------------------------------------------------------------
    primary_analytical_position(MOON, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(MOON, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);
    cout << "Final Moon position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    cout << "-------------------------------------------" << endl;


    //-------------------------------------------------------------------------------
    //Current Sun position
    //-------------------------------------------------------------------------------
    primary_analytical_position(SUN, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(SUN, yIN, yv[0]+I*yv[1], yv[2]+I*yv[3], yv[4]+I*yv[5], yv[6]+I*yv[7], qbcp_l);
    cout << "Final Sun position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;
    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;
    char ch;
    gnuplot_cmd(h1, "set logscale y");
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);
    gnuplot_close(h2);
}


