#include "vf.h"
#include "eminsem.h"

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Integration without STM
//
//--------------------------------------------------------------------------------------------------------------------------------------------
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

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Integration with STM
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j)
{
    if(mc != 0.0)
    {
        if(i == j) return factor*(3*mc/pow(qpc2,5.0/2)*pow(y[i]-pc[i],2.0) - mc/pow(qpc2,3.0/2));
        else return factor*3*mc/pow(qpc2,5.0/2)*(y[i]-pc[i])*(y[j]-pc[j]);
    }else return 0.0;
}

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
    gsl_matrix *STM  = gsl_matrix_alloc(6,6);
    gsl_matrix *STMd = gsl_matrix_alloc(6,6);
    gsl_matrix *Q    = gsl_matrix_alloc(6,6);

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

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Change of coordinates: SEM <-> IN <-> EM
//
//--------------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Normalized-Centered coordinates to Earth-Moon coordinates
 **/
void NCtoEM(double t, const double yNC[], double yEM[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];

}

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}


//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------
/**
 *  \brief COC: Normalized-Centered coordinates to system coordinates. Use in priority instead of NCtoEM or NCtoSEM.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from system coordinates to Normalized-Centered coordinates
 **/
void SYStoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -ySEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -ySEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +ySEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -ySEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -ySEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +ySEM[5]/gamma;
}

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], QBCP_L *qbp)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //-------------------------------------------------------------------------------
    //CoC
    //-------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    ySEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    ySEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    ySEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    ySEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    ySEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    ySEM[5] = +gamma*yNC[5];

}

//-----------------------------------------------------------------------------
// COC: NCEM <--> NCSEM
//-----------------------------------------------------------------------------

/**
 *  \brief COC: from NC coordinates (SEM) to NC coordinates (EM)
 **/
void SNCtoENC(double t, const double ySNC[], double yENC[], QBCP_L *qbp)
{
    double yEMm[6];
    double ySEMm[6];

    //SNC to SEM
    NCtoSEM(t, ySNC, ySEMm, qbp);
    //SEM to EM
    SEMmtoEMm(t, ySEMm, yEMm, qbp);
    //EM to ENC
    EMtoNC(t, yEMm, yENC, qbp);
}


/**
 *  \brief COC: from NC coordinates (EM) to NC coordinates (SEM)
 **/
void ENCtoSNC(double t, const double yENC[], double ySNC[], QBCP_L *qbp)
{
    double yEMm[6];
    double ySEMm[6];

    //ENC to EM
    NCtoEM(t, yENC, yEMm, qbp);
    //EM to SEM
    EMmtoSEMm(t, yEMm, ySEMm, qbp);
    //SEM to SNC
    SEMtoNC(t, ySEMm, ySNC, qbp);
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
// Vector <--> matrix
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            gsl_matrix_set(m, i-1, j-1, y[shift+k-1]);
        }
}

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift)
{
    int i,j,k;
    for(i=1; i<=rows ; i++)
        for(j=1; j<=columns; j++)
        {
            k = rows*(i-1)+j;
            y[shift+k-1] = gsl_matrix_get(m, i-1, j-1);
        }
}

