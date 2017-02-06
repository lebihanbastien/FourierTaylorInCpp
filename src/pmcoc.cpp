#include "pmcoc.h"
#include <gsl/gsl_sf_trig.h>

/**
 * \file pmcoc.cpp
 * \brief Change of coordinates for the parameterization methods.
 * \author BLB.
 * \date 2016
 * \version 1.0
 *
 *  Two types of change of coordinates are implemented:
 *              1. Changes between different manifold coordinates (Real Center to Complex Center...).
 *              2. Evaluations of the parameterization (Real Center to Normalized-Centered and projection the other way around).
 *
 */


//----------------------------------------------------------------------------------------
//
//          Change of coordinates
//
//----------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void CCMtoRCM(const cdouble s1[], double si[], int nv)
{
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[2]*I));
    si[2] = creal(1.0/sqrt(2)*(s1[0]*I + s1[2]));
    si[1] = creal(1.0/sqrt(2)*(s1[1]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[1]*I + s1[3]));
    if(nv == 5)
    {
        si[4] = creal(s1[4]);
    }
    else if(nv == 6)
    {
        si[4] = creal(s1[4]);
        si[5] = creal(s1[5]);
    }
}

/**
 *  \brief from RCM to CCM coordinates
 **/
void RCMtoCCM(const double si[], cdouble s1[], int nv)
{
    //From real to complex TFC
    s1[0] = 1.0/sqrt(2)*(si[0] - si[2]*I);
    s1[2] = 1.0/sqrt(2)*(si[2] - si[0]*I);
    s1[1] = 1.0/sqrt(2)*(si[1] - si[3]*I);
    s1[3] = 1.0/sqrt(2)*(si[3] - si[1]*I);
    if(nv == 5)
    {
        s1[4] = si[4]+I*0.0;
    }
    else if(nv == 6)
    {
        s1[4] = si[4]+I*0.0;
        s1[5] = si[5]+I*0.0;
    }
}

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void RCMtoCCM8(const double si[], double s0d[], int nv)
{
    //From real to complex TFC
    cdouble s1[nv];
    RCMtoCCM(si, s1, nv);

    //Store real and imag part separately
    CCMtoCCM8(s1, s0d, nv);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void CCM8toRCM(const double s0d[], double si[], int nv)
{
    //CCM8 to CCM
    cdouble s1[nv];
    CCM8toCCM(s0d, s1, nv);
    //CCM to RCM
    CCMtoRCM(s1, si, nv);
}

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void CCM8toCCM(const double s0d[], cdouble s1[], int nv)
{
    int p2 = 0;
    for(int p = 0; p < nv; p++)
    {
        s1[p]  =   s0d[p2++];
        s1[p] += I*s0d[p2++];
    }
}

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void CCMtoCCM8(const cdouble s1[], double s0d[], int nv)
{
    //Store real and imag part separately
    int p2 = 0;
    for(int p = 0; p < nv; p++)
    {
        s0d[p2++] = creal(s1[p]);
        s0d[p2++] = cimag(s1[p]);
    }
}

/**
 *  \brief from TFC to TF coordinates
 **/
void TFCtoTF(const cdouble s1[6], double si[6])
{
    //First center
    si[0] = creal(1.0/sqrt(2)*(s1[0]   + s1[3]*I));
    si[3] = creal(1.0/sqrt(2)*(s1[0]*I + s1[3]));
    //Second center
    si[2] = creal(1.0/sqrt(2)*(s1[2]   + s1[5]*I));
    si[5] = creal(1.0/sqrt(2)*(s1[2]*I + s1[5]));
    //Hyperbolic dir
    si[1] = creal(s1[1]);
    si[4] = creal(s1[4]);
}


//----------------------------------------------------------------------------------------
//
//          Evaluation of the pm
//
//----------------------------------------------------------------------------------------
/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param W the Fourier-Taylor expansion that contains the parameterization W(s,t)
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute W.
 *          Namely: W[2] and W[5] are of order one, the rest of the components are complete polynomials.
 **/
void RCMtoNC(const double st0[],
             const double t,
             const double n,
             const int order,
             const int ofs_order,
             const int reduced_nv,
             vector<Oftsc> &W,
             Ofsc &ofs,
             double z1[],
             bool isGS)
{
    //------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------
    cdouble s0[reduced_nv];
    cdouble z0[6];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);
    //------------------------------------------
    // 2. Update z0
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            if(p == 2 || p == 5)
            {
                //order 1 is sufficient for p = 2,5
                ofs.zero();
                W[p].evaluate(s0, ofs, 1, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
            else
            {
                //For p = 0,1,3,4 normal computation
                ofs.zero();
                W[p].evaluate(s0, ofs, order, ofs_order);
                z0[p] = ofs.evaluate(n*t, ofs_order);
                z1[p] = creal(z0[p]);
            }
        }
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            W[p].evaluate(s0, ofs, order, ofs_order);
            z0[p] = ofs.evaluate(n*t, ofs_order);
            z1[p] = creal(z0[p]);
        }
    }

}

/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector:  z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoNCbyTFC(const double st0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  const int reduced_nv,
                  vector<Oftsc> &Wh,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // RCM to TFC
    //------------------------------------------
    RCMtoTFC(st0, order, ofs_order, reduced_nv, Wh, zIn, isGS);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}

/**
 *   \brief Evaluate the configuration z1 = W(g(st0), n*t) with the use of an intermediate TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param s0 an array of 4 complex double which gives the configuration to input in complex CM coordinates
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param PC COC matrix: z = PC*zh+V
 *   \param V COC vector: z = PC*zh+V
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoNCbyTFC(cdouble s0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  vector<Oftsc> &Wh,
                  matrix<Ofsc> &PC,
                  vector<Ofsc> &V,
                  double z1[],
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble z0[6];
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // CCM to TFC
    // If reduced_nv > 4, we cannot use the simplified
    // version of the parameterization, since the expansions are all full expansions.
    //------------------------------------------
    if(Wh[0].getNV() > 4) CCMtoTFC(s0, order, ofs_order, Wh, zIn, 0);
    else CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);


    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC(PC, V, zIn, zOut);
    for(int p = 0; p < 6; p++)
    {
        z0[p] = zOut[p].evaluate(n*t, ofs_order);
        z1[p] = creal(z0[p]);
    }
}


/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(st0), t)
 *   \param st0 an array of 4 double which gives the configuration to input in real CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoTFC(const double st0[],
              const int order,
              const int ofs_order,
              const int reduced_nv,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[reduced_nv];

    //------------------------------------------
    // 1. RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);

    //------------------------------------------
    // 2. Update zIn
    // If reduced_nv > 4, we cannot use the simplified
    // version of the parameterization, since the expansions are all full expansions.
    //------------------------------------------
    if(reduced_nv > 4) CCMtoTFC(s0, order, ofs_order, Wh, zIn, 0);
    else CCMtoTFC(s0, order, ofs_order, Wh, zIn, isGS);
}


/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(const matrix<Ofsc> &PC,
              const vector<Ofsc> &V,
              const vector<Ofsc> &zIn,
              vector<Ofsc> &zOut)
{
    //zeroing the target
    for(unsigned int i = 0; i < zOut.size(); i++) zOut[i].zero();
    //zOut = PC*zIn
    smvprod_ofs(PC, zIn, zOut);
    //zOut+=V(theta)
    for(int i = 0; i < (int) zOut.size(); i++) zOut[i] += V[i];
}


/**
 *   \brief Evaluate the TFC configuration zIn(t) = Wh(g(s0), t)
 *   \param s0 an array of 4 complex which gives the configuration to input in complex CM coordinates
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param Wh the Fourier-Taylor expansion that contains the parameterization Wh(s,t), in TFC coordinates
 *   \param ofs a side variable to compute OFS objects
 *   \param z1 the output array to update
 *   \param isGS if true, the special case of the QBCP is used to compute Wh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoTFC(cdouble s0[],
              const int order,
              const int ofs_order,
              vector<Oftsc> &Wh,
              vector<Ofsc> &zIn,
              bool isGS)
{
    //------------------------------------------
    // 1. Zeroing the results
    //------------------------------------------
    for(int i = 0; i < (int) zIn.size(); i++) zIn[i].zero();

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP
        //--------------------
        cdouble temp;
        // zIn[0]
        //---------------
        temp = Wh[0].getCA(1,0)->getCoef(0);
        zIn[0].setCoef(temp*s0[0],0);
        // zIn[1]
        //---------------
        Wh[1].evaluate(s0, zIn[1], order, ofs_order);
        // zIn[2]
        //---------------
        temp = Wh[2].getCA(1,1)->getCoef(0);
        zIn[2].setCoef(temp*s0[1],0);
        // zIn[3]
        //---------------
        temp = Wh[3].getCA(1,2)->getCoef(0);
        zIn[3].setCoef(temp*s0[2],0);
        // zIn[4]
        //---------------
        Wh[4].evaluate(s0, zIn[4], order, ofs_order);
        // zIn[5]
        //---------------
        temp = Wh[5].getCA(1,3)->getCoef(0);
        zIn[5].setCoef(temp*s0[3],0);
    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int p = 0; p < 6; p++)
        {
            Wh[p].evaluate(s0, zIn[p], order, ofs_order);
        }
    }
}

/**
 *  \brief Projection of the current NC state on the central manifold, via CCM coordinates
 **/
void NCprojCCM(const double z[], const double t, const double n, const int ofs_order,
               const matrix<Ofsc> &CQ, const vector<Ofsc> &V,
               double omega1, double omega3, cdouble sc[], int nv)
{
    //-------------------
    //Wh: TFC coordinates
    //-------------------
    cdouble zh[6];

    //-------------------
    //z - V
    //-------------------
    cdouble zd[6];
    for(int p = 0; p < 6; p++) zd[p] = z[p] - V[p].evaluate(n*t, ofs_order);

    //-------------------
    //Update Wh = CQ*(z - V)
    //-------------------
    for(int k = 0; k <6; k++)
    {
        zh[k] = 0.0+0.0*I;
        for(int p = 0; p <6; p++)
        {
            zh[k] += zd[p]* CQ.getCoef(k, p).evaluate(n*t, ofs_order);
        }
    }

    //-------------------
    //Projection on the center manifold
    //-------------------
    TFCprojCCM(zh, omega1, omega3, sc, nv);
}

/**
 *  \brief Projection of the current TFC state on the central manifold, via CCM coordinates
 **/
void TFCprojCCM(const cdouble zh[], double omega1, double omega3, cdouble sc[], int nv)
{
    //-------------------
    //Projection on the center manifold
    //-------------------
    //Wh1 = i*w1*s1 => s1 = -i/w1*Wh1
    sc[0] = -1.0*I/omega1*zh[0];

    //Wh3 = i*w3*s2 => s2 = -i/w3*Wh3
    sc[1] = -1.0*I/omega3*zh[2];

    //Wh4 = -i*w1*s3 => s3 = +i/w1*Wh4
    sc[2] = +1.0*I/omega1*zh[3];

    //Wh6 = -i*w3*s4 => s4 = +i/w3*Wh6
    sc[3] = +1.0*I/omega3*zh[5];

    if(nv > 4) sc[4] = 0.0;
    if(nv > 5) sc[5] = 0.0;
}

/**
 *   \brief Evaluate the Jacobian of the TFC configuration zIn(t) = Wh(g(s0), t), in mIn(t) = DWh(g(s0), t).
 *   \param isGS if true, the special case of the QBCP is used to compute DWh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoTFC_JAC(const double st0[],
                  const int order,
                  const int ofs_order,
                  const int reduced_nv,
                  matrix<Oftsc> &DWh,
                  matrix<Ofsc> &mIn,
                  bool isGS)
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[reduced_nv];

    //------------------------------------------
    // 1. RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);

    //------------------------------------------
    // 2. Update zIn
    //------------------------------------------
    if(reduced_nv > 4)  CCMtoTFC_JAC(s0, order, ofs_order, DWh, mIn, 0);
    else CCMtoTFC_JAC(s0, order, ofs_order, DWh, mIn, isGS);
}


/**
 *   \brief Evaluate the Jacobian of the TFC configuration zIn(t) = Wh(g(s0), t), in mIn(t) = DWh(g(s0), t).
 *   \param isGS if true, the special case of the QBCP is used to compute DWh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void RCMtoTFC_JAC(const double st0[],
                  const double t,
                  const double n,
                  const int order,
                  const int ofs_order,
                  const int reduced_nv,
                  matrix<Oftsc> &DWh,
                  matrix<Ofsc>  &mIn,
                  gsl_matrix_complex *m1,
                  bool isGS)
{
    //------------------------------------------
    // 1. RCM to TFC
    //------------------------------------------
    RCMtoTFC_JAC(st0, order, ofs_order, reduced_nv, DWh, mIn, isGS);

    //------------------------------------------
    // 3. Evaluate mIn in m2
    //------------------------------------------
    evaluate(n*t, mIn, m1);
}


/**
 *   \brief Evaluate the Jacobian of the TFC configuration zIn(t) = Wh(g(s0), t), in mIn(t) = DWh(g(s0), t).
 *   \param isGS if true, the special case of the QBCP is used to compute DWh.
 *          Namely: only Wh[1] and Wh[4] are complete polynomials, the rest of the components are of order one.
 **/
void CCMtoTFC_JAC(cdouble s0[],
                  const int order,
                  const int ofs_order,
                  matrix<Oftsc> &DWh,
                  matrix<Ofsc>  &mIn,
                  bool isGS)
{
    //------------------------------------------
    // 1. Check the size
    //------------------------------------------
    //    if(mIn.getSize(1) != DWh.getSize(1) || mIn.getSize(2) != DWh.getSize(2))
    //    {
    //        cout << "CCMtoTFC_JAC. Dimensions do not match:" << endl;
    //        cout << "mIn.getSize(1) = " << mIn.getSize(1) << ", DWh.getSize(1) = " << DWh.getSize(1) << endl;
    //        cout << "mIn.getSize(2) = " << mIn.getSize(2) << ", DWh.getSize(2) = " << DWh.getSize(2) << endl;
    //        cout << "return." << endl;
    //        return;
    //    }

    //--------------------
    // 2. Zeroing
    //--------------------
    for(int i = 0; i < mIn.getSize(1); i++)
    {
        for(int j = 0; j < mIn.getSize(2); j++)
        {
            mIn.getCA(i,j)->zero();
        }
    }

    //Temp variables
    Ofsc ofs_temp(mIn.getCA(0,0)->getOrder());
    cdouble temp;

    //------------------------------------------
    // 3. Update mIn
    //------------------------------------------
    if(isGS)
    {
        //--------------------
        // Using particular geometry of the QBCP:
        //--------------------
        // mIn[0][*] = [dWh[0]/ds[0] 0 0 ... 0]
        //---------------
        temp = DWh.getCoef(0,0).getCA(0,0)->getCoef(0);
        mIn.getCA(0,0)->setCoef(temp,0);


        // mIn[1][*] = full FT series
        //---------------
        for(int j = 0; j < mIn.getSize(2); j++)
        {
            DWh.getCoef(1,j).evaluate(s0, ofs_temp, order, ofs_order);
            mIn.getCA(1,j)->ccopy(ofs_temp);
        }

        // mIn[2][*] = [0 dWh[2]/ds[1] 0 ... 0]
        //---------------
        temp = DWh.getCoef(2,1).getCA(0,0)->getCoef(0);
        mIn.getCA(2,1)->setCoef(temp,0);

        // mIn[3][*] = [0 0 dWh[3]/ds[2] 0 ... 0]
        //---------------
        temp = DWh.getCoef(3,2).getCA(0,0)->getCoef(0);
        mIn.getCA(3,2)->setCoef(temp,0);

        // mIn[4][*] = full FT series
        //---------------
        for(int j = 0; j < mIn.getSize(2); j++)
        {
            DWh.getCoef(4,j).evaluate(s0, ofs_temp, order, ofs_order);
            mIn.getCA(4,j)->ccopy(ofs_temp);
        }

        // mIn[5][*] = [0 0 0 dWh[5]/ds[4] ... 0]
        //---------------
        temp = DWh.getCoef(5,3).getCA(0,0)->getCoef(0);
        mIn.getCA(5,3)->setCoef(temp,0);

    }
    else
    {
        //--------------------
        // General computation
        //--------------------
        for(int i = 0; i < mIn.getSize(1); i++)
        {
            for(int j = 0; j < mIn.getSize(2); j++)
            {
                DWh.getCoef(i,j).evaluate(s0, ofs_temp, order, ofs_order);
                mIn.getCA(i,j)->ccopy(ofs_temp);
            }
        }
    }
}


/**
 *   \brief Evaluate the reduced vector field (RVF) dot(s4) = f4(s4, t)
 *   \param s8 an array of 8 doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param f4 the RVF output array of 4 cdouble to update
 **/
void CCM8toRVF(const double s8[],
               const double t,
               const double n,
               const int order,
               const int ofs_order,
               const int reduced_nv,
               vector<Oftsc> &fh,
               Ofsc &ofs,
               cdouble f4[])
{
    //CCM8 to CCM
    //----------
    cdouble s[reduced_nv];
    CCM8toCCM(s8, s, reduced_nv);

    //Evaluation of fh
    //----------
    for(int p = 0; p < reduced_nv; p++)
    {
        fh[p].evaluate(s, ofs, order, ofs_order);
        f4[p] = ofs.evaluate(n*t, ofs_order);
    }
}


/**
 *   \brief Evaluate the reduced vector field (RVF) dot(s4) = f4(s4, t)
 *   \param s8 an array of 8 doubles which gives the current CCM state, with separated real and imag parts.
 *   \param t the current time
 *   \param n the pulsation of the system
 *   \param order the order of the Taylor expansions to use in the evaluation
 *   \param ofs_order the order of the Fourier expansions to use in the evaluation
 *   \param fh the Fourier-Taylor expansion that contains the RVF.
 *   \param ofs a side variable to compute OFS objects
 *   \param f8 the RVF output array of 8 double to update, with separated real and imag parts.
 **/
void CCM8toRVF8(const double s8[],
                const double t,
                const double n,
                const int order,
                const int ofs_order,
                const int reduced_nv,
                vector<Oftsc> &fh,
                Ofsc &ofs,
                double f8[])
{
    //CCM8 to RVF
    //----------
    cdouble sd[reduced_nv];
    CCM8toRVF(s8, t, n, order, ofs_order, reduced_nv, fh, ofs, sd);

    //Separation of real and imag parts
    //----------
    int p2 = 0;
    for(int p = 0; p < reduced_nv; p++)
    {
        f8[p2++] = creal(sd[p]);
        f8[p2++] = cimag(sd[p]);
    }
}


//----------------------------------------------------------------------------------------
//
//          Evaluation at time t
//
//----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i < order; i++)
    {
        //From trigo formulae (fastest)
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
        //From GSL
        //cR[i] =  gsl_sf_cos((i+1)*omega*t);
        //sR[i] =  gsl_sf_sin((i+1)*omega*t);
        //Native C
        //cR[i] =  cos((i+1)*omega*t);
        //sR[i] =  sin((i+1)*omega*t);
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13) alpha[l] = evaluateOdd(order, header, sR);    //Odd funtions (alpha_2,5,8,10,12,14)
        else  alpha[l] = evaluateEven(order, header, cR);                                                  //Even functions
    }
}

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number)
{
    double *header;
    int l;
    double cR[order];
    double sR[order];

    cR[0] = cos(omega*t);
    sR[0] = sin(omega*t);

    for(int i = 1; i< order; i++)
    {
        cR[i] =  cR[i-1]*cR[0] - sR[i-1]*sR[0];
        sR[i] =  sR[i-1]*cR[0] + cR[i-1]*sR[0];
    }

    for(l = 0, header = params; l < number ; l++, header+=(order+1))
    {
        alpha[l] = 0.0;
        if(l==1 || l== 4 || l==7 || l==9 || l==11 || l==13)
        {
            alpha[l] = evaluateOddDerivative(omega, order, header, cR); //Odd funtions (alpha_2,5,8,10,12,14)
        }
        else  alpha[l] = evaluateEvenDerivative(omega, order, header, sR);   //Even functions
    }
}


//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*cR[i-1];//even type
    result += coef[0];
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double omega, int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += -omega*i*coef[i]*sR[i-1];//even type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(int order, double *coef, double *sR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += coef[i]*sR[i-1]; //odd type
    return result;
}

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double omega, int order, double *coef, double *cR)
{
    double result = 0.0;
    for(int i= order; i>=1; i--) result += omega*i*coef[i]*cR[i-1];//odd type
    return result;
}

//-----------------------------------------------------------------------------
// Evaluating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(const Ofsc& zt, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*zt.evaluate(n*t);
}

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(const Ofsc& zt, const Ofsc& ztdot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*(ztdot.evaluate(n*t) + I*ni*zt.evaluate(n*t));
}

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(const Ofsc& zt, const Ofsc& ztdot, const Ofsc& ztddot, double t, double n, double ni, double ai)
{
    return ai*(cos(ni*t)+I*sin(ni*t))*( 2*I*ni*ztdot.evaluate(n*t) - ni*ni*zt.evaluate(n*t) + ztddot.evaluate(n*t));
}


//----------------------------------------------------------------------------------------
//
//          Matrices for COC
//
//----------------------------------------------------------------------------------------
/**
 *  \brief RCOC: from inputType to outputType. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void rot_mat_coc(double t, gsl_matrix* Rcoc, int inputType, int outputType)
{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    //Type of inputs
    if(inputType > INSEM)
        perror("rot_mat_coc. Unknown inputType");

    //Type of outputs
    if(outputType > INSEM)
        perror("rot_mat_coc. Unknown outputType");

    //====================================================================================
    // 2. Define the default framework wrt the inputType
    //====================================================================================
    int fwrk = default_framework(inputType);

    //====================================================================================
    // 2. Check that the focus in SEML is
    // in accordance with the inputType.
    //====================================================================================
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);


    //====================================================================================
    // 3. Updating output
    //====================================================================================

    //-----------------------------------
    // 3.1 Parameters in coordsys units
    //-----------------------------------
    double n     = SEML.us.n;
    double ns    = SEML.us.ns;
    double as    = SEML.us.as;
    double ni    = SEML.us.ni;
    double ai    = SEML.us.ai;
    double gamma = SEML.cs.gamma;

    //-----------------------------------
    // 3.2 Jacobi decomposition (r, R)
    //     of the QBTBP in coordsys units
    //-----------------------------------
    //r
    double r1 = creal(evz(SEML.cs.zt, t, n, ni, ai));
    double r2 = cimag(evz(SEML.cs.zt, t, n, ni, ai));
    double r  = sqrt(r1*r1 + r2*r2);

    //rdot
    double r1dot = creal(evzdot(SEML.cs.zt, SEML.cs.ztdot, t, n, ni, ai));
    double r2dot = cimag(evzdot(SEML.cs.zt, SEML.cs.ztdot, t, n, ni, ai));
    double rdot  = 1.0/r*(r1*r1dot + r2*r2dot);

    //R
    double R1 = creal(evz(SEML.cs.Zt, t, n, ns, as));
    double R2 = cimag(evz(SEML.cs.Zt, t, n, ns, as));
    double R  = sqrt(R1*R1 + R2*R2);

    //Rdot
    double R1dot = creal(evzdot(SEML.cs.Zt, SEML.cs.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(SEML.cs.Zt, SEML.cs.Ztdot, t, n, ns, as));
    double Rdot  = 1.0/R*(R1*R1dot + R2*R2dot);

    //Additionnal parameters
    double a1 = r1/(r*r);
    double a2 = r2/(r*r);
    double a6 = 1.0/r;
    double a1dot = +pow(r,-4.0)*(r1dot*r*r - 2*r1*(r1*r1dot + r2*r2dot)); //dot(r1/(r*r))
    double a2dot = +pow(r,-4.0)*(r2dot*r*r - 2*r2*(r1*r1dot + r2*r2dot)); //dot(r2/(r*r))
    double a6dot = -pow(r,-3.0)*(r1*r1dot + r2*r2dot);                    //dot(1/r)

    double A1 = R1/(R*R);
    double A2 = R2/(R*R);
    double A6 = 1.0/R;
    double A1dot = +pow(R,-4.0)*(R1dot*R*R - 2*R1*(R1*R1dot + R2*R2dot)); //dot(R1/(R*R))
    double A2dot = +pow(R,-4.0)*(R2dot*R*R - 2*R2*(R1*R1dot + R2*R2dot)); //dot(R2/(R*R))
    double A6dot = -pow(R,-3.0)*(R1*R1dot + R2*R2dot);                    //dot(1/R)


    //-----------------------------------
    // 3.3. Coefficients of the QBCP in coordsys units
    //-----------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, SEML.nf, SEML.cs.coeffs, 3);


    //-----------------------------------
    // 3.4. Temporary matrices
    //-----------------------------------
    gsl_matrix *RA = gsl_matrix_alloc(6,6);
    gsl_matrix *RB = gsl_matrix_alloc(6,6);
    gsl_matrix *RC = gsl_matrix_alloc(6,6);
    gsl_matrix *RD = gsl_matrix_alloc(6,6);
    gsl_matrix *RE = gsl_matrix_alloc(6,6);
    gsl_matrix *RF = gsl_matrix_alloc(6,6);


    //-----------------------------------
    // 3.5. Super switch (input/output)
    //-----------------------------------
    switch(inputType)
    {
        //-----------------------------------------------------------------
        // For inputType = VNCEM, NCEM, PEM, VEM, INEM, the focus of SEML is
        // on the EM SYSTEM
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // VEM
        //-----------------------------------------------------------------
    case VEM:
        switch(outputType)
        {

        case VEM:
            //-----------------------------------------------------
            // VEM -> VEM
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //VEM -> INEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0,  r1);
            gsl_matrix_set(Rcoc, 0, 1, -r2);
            gsl_matrix_set(Rcoc, 1, 0,  r2);
            gsl_matrix_set(Rcoc, 1, 1,  r1);

            gsl_matrix_set(Rcoc, 2, 2,  r);

            gsl_matrix_set(Rcoc, 3, 0,  r1dot);
            gsl_matrix_set(Rcoc, 3, 1, -r2dot);
            gsl_matrix_set(Rcoc, 4, 0,  r2dot);
            gsl_matrix_set(Rcoc, 4, 1,  r1dot);

            gsl_matrix_set(Rcoc, 3, 3,  r1);
            gsl_matrix_set(Rcoc, 3, 4, -r2);
            gsl_matrix_set(Rcoc, 4, 3,  r2);
            gsl_matrix_set(Rcoc, 4, 4,  r1);

            gsl_matrix_set(Rcoc, 5, 2,  rdot);
            gsl_matrix_set(Rcoc, 5, 5,  r);
            break;

        case VNCEM:
            //-----------------------------------------------------
            //VEM -> VNCEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +1.0/gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +1.0/gamma);
            break;

        case PEM:
            //-----------------------------------------------------
            // VEM -> PEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, -alpha[2]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 0, +alpha[2]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 0, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 1, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 2, -alpha[1]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 3, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, 1.0/alpha[0]);
            break;

        case NCEM:
            //-----------------------------------------------------
            // VEM -> NCEM
            //-----------------------------------------------------
            //RA = EM_R_VEM
            rot_mat_coc(t, RA, VEM, PEM);
            //RB =  NCEM_R_EM
            rot_mat_coc(t, RB, PEM, NCEM);
            //Rcoc = NCEM_R_EM*EM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB, RA, 0.0, Rcoc);
            break;

        case INSEM:
            //-----------------------------------------------------
            // VEM -> INSEM
            //-----------------------------------------------------
            //RA = INEM_R_VEM
            rot_mat_coc(t, RA, VEM, INEM);
            //RB =  INSEM_R_INEM
            rot_mat_coc(t, RB, INEM, INSEM);
            //Rcoc = INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);

        case VSEM:
            //-----------------------------------------------------
            // VEM -> VSEM
            //-----------------------------------------------------
            //RA = INEM_R_VEM
            rot_mat_coc(t, RA, VEM, INEM);
            //RB =  INSEM_R_INEM
            rot_mat_coc(t, RB, INEM, INSEM);

            //RC =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RC, INSEM, VSEM);

            //Rcoc = INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = VSEM_R_INSEM*RB = VSEM_R_INSEM*INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case PSEM:
            //-----------------------------------------------------
            // VEM -> PSEM
            //-----------------------------------------------------
            //RA = INEM_R_VEM
            rot_mat_coc(t, RA, VEM, INEM);
            //RB =  INSEM_R_INEM
            rot_mat_coc(t, RB, INEM, INSEM);

            //RC =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RC, INSEM, VSEM);
            //RD =  SEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RD, VSEM, PSEM);

            //Rcoc = INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= SEM_R_VSEM*VSEM_R_INSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case NCSEM:
            //-----------------------------------------------------
            // VEM -> NCSEM
            //-----------------------------------------------------
            //RA = PSEM_R_VEM
            rot_mat_coc(t, RA, VEM, PSEM);
            //RB =  NCSEM_R_PSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, PSEM, NCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCSEM:
            //-----------------------------------------------------
            // VEM -> VNCSEM
            //-----------------------------------------------------
            //RA = VSEM_R_VEM
            rot_mat_coc(t, RA, VEM, VSEM);
            //RB =  VNCSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, VSEM, VNCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // VNCEM
        //-----------------------------------------------------------------
    case VNCEM:
        switch(outputType)
        {
        case VNCEM:
            //-----------------------------------------------------
            //VNCEM -> VNCEM nothing to do
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VEM:
            //-----------------------------------------------------
            //VNCEM -> VEM
            //-----------------------------------------------------
            //Rcoc = VEM_R_VNCEM
            gsl_matrix_set(Rcoc, 0, 0,  -gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +gamma);
            break;

        case NCEM:
            //-----------------------------------------------------
            //VNCEM -> NCEM: equivalent to VEM -> PEM
            //-----------------------------------------------------
            //Rcoc = NCEM_R_VNCEM = PEM_R_VEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, -alpha[2]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 0, +alpha[2]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 0, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 1, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 2, -alpha[1]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 3, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, 1.0/alpha[0]);

            break;

        case INEM:
            //-----------------------------------------------------
            //VNCEM -> INEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VEM);
            //RB =  INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //Rcoc = INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case PEM:
            //-----------------------------------------------------
            //VNCEM -> PEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VEM);
            //RB =  PEM_R_VEM
            rot_mat_coc(t, RB, VEM, PEM);
            //Rcoc = PEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case INSEM:
            //-----------------------------------------------------
            //VNCEM -> INSEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VEM);
            //RB = INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);

            //Rcoc = RB*RA= INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INSEM_R_INEM*RB = INSEM_R_INEM*INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case VSEM:
            //-----------------------------------------------------
            //VNCEM -> VSEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VEM);
            //RB = INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);
            //RD =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RD, INSEM, VSEM);


            //Rcoc = INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VSEM_R_INSEM*INSEM_R_INEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VNCEM -> PSEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VEM);
            //RB = INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);
            //RD =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RD, INSEM, VSEM);
            //RE =  PSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RE, VSEM, PSEM);

            //Rcoc = INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VSEM_R_INSEM*INSEM_R_INEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = VSEM_R_INSEM*INSEM_R_INEM*INEM_R_VEM*VEM_R_VNCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;


        case NCSEM:
            //-----------------------------------------------------
            // VNCEM -> NCSEM
            //-----------------------------------------------------
            //RA = PSEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, PSEM);
            //RB =  NCSEM_R_PSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, PSEM, NCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCSEM:
            //-----------------------------------------------------
            // VNCEM -> VNCSEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCEM
            rot_mat_coc(t, RA, VNCEM, VSEM);
            //RB =  VNCSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, VSEM, VNCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;


        //-----------------------------------------------------------------
        // NCEM
        //-----------------------------------------------------------------
    case NCEM:
        switch(outputType)
        {

        case NCEM:
            //-----------------------------------------------------
            //NCEM -> NCEM nothing to do
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VNCEM:
            //-----------------------------------------------------
            //NCEM -> VNCEM
            //-----------------------------------------------------
            //Rcoc = VNCEM_R_NCEM = VEM_R_PEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, +alpha[2]);
            gsl_matrix_set(Rcoc, 4, 0, -alpha[2]);

            gsl_matrix_set(Rcoc, 3, 0, +alpha[1]);
            gsl_matrix_set(Rcoc, 4, 1, +alpha[1]);
            gsl_matrix_set(Rcoc, 5, 2, +alpha[1]);

            gsl_matrix_set(Rcoc, 3, 3, +alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, +alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, +alpha[0]);
            break;

        case PEM:
            //-----------------------------------------------------
            //NCEM -> PEM
            //-----------------------------------------------------
            //Rcoc = EM_R_NCEM
            gsl_matrix_set(Rcoc, 0, 0,  -gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +gamma);
            break;

        case VEM:
            //-----------------------------------------------------
            //NCEM -> VEM
            //-----------------------------------------------------
            //RA = VNCEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VNCEM);
            //RB =  VEM_R_VNCEM
            rot_mat_coc(t, RB, VNCEM, VEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //NCEM -> INEM
            //-----------------------------------------------------
            //RA = VNCEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VNCEM);
            //RB =  VEM_R_VNCEM
            rot_mat_coc(t, RB, VNCEM, VEM);
            //RC =  INEM_R_VEM
            rot_mat_coc(t, RC, VEM, INEM);

            //Rcoc = VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INEM_R_VEM*RB = INEM_R_VEM*VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case INSEM:
            //-----------------------------------------------------
            //NCEM -> INSEM
            //-----------------------------------------------------
            //RA = VNCEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VNCEM);
            //RB =  VEM_R_VNCEM
            rot_mat_coc(t, RB, VNCEM, VEM);
            //RC =  INEM_R_VEM
            rot_mat_coc(t, RC, VEM, INEM);
            //RD =  INSEM_R_INEM
            rot_mat_coc(t, RD, INEM, INSEM);

            //Rcoc = VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case VSEM:
            //-----------------------------------------------------
            //NCEM -> VSEM
            //-----------------------------------------------------
            //RA = VNCEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VNCEM);
            //RB =  VEM_R_VNCEM
            rot_mat_coc(t, RB, VNCEM, VEM);
            //RC =  INEM_R_VEM
            rot_mat_coc(t, RC, VEM, INEM);
            //RD =  INSEM_R_INEM
            rot_mat_coc(t, RD, INEM, INSEM);
            //RE =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RE, INSEM, VSEM);

            //Rcoc = VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = INSEM_R_INEM*INEM_R_VEM*VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;

        case PSEM:
            //-----------------------------------------------------
            //NCEM -> PSEM
            //-----------------------------------------------------
            //RA = VNCEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VNCEM);
            //RB =  VEM_R_VNCEM
            rot_mat_coc(t, RB, VNCEM, VEM);
            //RC =  INEM_R_VEM
            rot_mat_coc(t, RC, VEM, INEM);
            //RD =  INSEM_R_INEM
            rot_mat_coc(t, RD, INEM, INSEM);
            //RE =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RE, INSEM, VSEM);
            //RF =  PSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RF, VSEM, PSEM);


            //Rcoc = RB*RA = VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB = RC*RD = INSEM_R_INEM*INEM_R_VEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = INSEM_R_INEM*INEM_R_VEM*VEM_R_VNCEM*VNCEM_R_NCEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //RC = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, RC);

            //Rcoc = RF*RC = RE*RF*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RF , RC, 0.0, Rcoc);
            break;


        case NCSEM:
            //-----------------------------------------------------
            // NCEM -> NCSEM
            //-----------------------------------------------------
            //RA = PSEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, PSEM);
            //RB =  NCSEM_R_PSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, PSEM, NCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCSEM:
            //-----------------------------------------------------
            // NCEM -> VNCSEM
            //-----------------------------------------------------
            //RA = VSEM_R_NCEM
            rot_mat_coc(t, RA, NCEM, VSEM);
            //RB =  VNCSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, VSEM, VNCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // PEM
        //-----------------------------------------------------------------
    case PEM:
        switch(outputType)
        {
        case PEM:
            //-----------------------------------------------------
            //PEM -> PEM
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VEM:
            //-----------------------------------------------------
            // PEM -> VEM
            //-----------------------------------------------------
            //Rcoc = VEM_R_PEM = VNCEM_R_NCEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, +alpha[2]);
            gsl_matrix_set(Rcoc, 4, 0, -alpha[2]);

            gsl_matrix_set(Rcoc, 3, 0, +alpha[1]);
            gsl_matrix_set(Rcoc, 4, 1, +alpha[1]);
            gsl_matrix_set(Rcoc, 5, 2, +alpha[1]);

            gsl_matrix_set(Rcoc, 3, 3, +alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, +alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, +alpha[0]);
            break;

        case NCEM:
            //-----------------------------------------------------
            //PEM -> NCEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +1.0/gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +1.0/gamma);
            break;

        case VNCEM:
            //-----------------------------------------------------
            // PEM -> VNCEM
            //-----------------------------------------------------
            //RA = VEM_R_PEM
            rot_mat_coc(t, RA, PEM, VEM);
            //RB =  VNCEM_R_VEM
            rot_mat_coc(t, RB, VEM, VNCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //PEM -> INEM
            //-----------------------------------------------------
            //RA = VEM_R_PEM
            rot_mat_coc(t, RA, PEM, VEM);
            //RB =  INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INSEM:
            //-----------------------------------------------------
            //PEM -> INSEM
            //-----------------------------------------------------
            //RA = VEM_R_PEM
            rot_mat_coc(t, RA, PEM, VEM);
            //RB =  INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);

            //Rcoc = RB*RA= INEM_R_VEM*VEM_R_PEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INSEM_R_INEM*RB = INSEM_R_INEM*INEM_R_VEM*VEM_R_PEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case VSEM:
            //-----------------------------------------------------
            //PEM -> VSEM
            //-----------------------------------------------------
            //RA = VEM_R_PEM
            rot_mat_coc(t, RA, PEM, VEM);
            //RB =  INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);
            //RD =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RD, INSEM, VSEM);


            //Rcoc = INEM_R_VEM*VEM_R_PEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VSEM_R_INSEM*INSEM_R_INEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case PSEM:
            //-----------------------------------------------------
            //PEM -> PSEM
            //-----------------------------------------------------
            //RA = VEM_R_PEM
            rot_mat_coc(t, RA, PEM, VEM);
            //RB =  INEM_R_VEM
            rot_mat_coc(t, RB, VEM, INEM);
            //RC =  INSEM_R_INEM
            rot_mat_coc(t, RC, INEM, INSEM);
            //RD =  VSEM_R_INSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RD, INSEM, VSEM);
            //RE =  PSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RE, VSEM, PSEM);

            //Rcoc = INEM_R_VEM*VEM_R_PEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VSEM_R_INSEM*INSEM_R_INEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = VSEM_R_INSEM*INSEM_R_INEM*INEM_R_VEM*VEM_R_PEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;

        case NCSEM:
            //-----------------------------------------------------
            // PEM -> NCSEM
            //-----------------------------------------------------
            //RA = PSEM_R_PEM
            rot_mat_coc(t, RA, PEM, PSEM);
            //RB =  NCSEM_R_PSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, PSEM, NCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCSEM:
            //-----------------------------------------------------
            // PEM -> VNCSEM
            //-----------------------------------------------------
            //RA = VSEM_R_PEM
            rot_mat_coc(t, RA, PEM, VSEM);
            //RB =  VNCSEM_R_VSEM - time must be in SEM units!
            rot_mat_coc(t*SEML.us_em.ns, RB, VSEM, VNCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // For inputType = VNCSEM, NCSEM, PSEM, VSEM, the focus of SEML is
        // on the SEM SYSTEM
        //-----------------------------------------------------------------
        //-----------------------------------------------------------------
        // VSEM
        //-----------------------------------------------------------------
    case VSEM:
        switch(outputType)
        {
        case VSEM:
            //-----------------------------------------------------
            //VSEM -> VSEM
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VSEM -> PSEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, -alpha[2]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 0, +alpha[2]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 0, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 1, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 2, -alpha[1]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 3, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, 1.0/alpha[0]);
            break;

        case INSEM:
            //-----------------------------------------------------
            //VSEM -> INSEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0,  R1);
            gsl_matrix_set(Rcoc, 0, 1, -R2);
            gsl_matrix_set(Rcoc, 1, 0,  R2);
            gsl_matrix_set(Rcoc, 1, 1,  R1);

            gsl_matrix_set(Rcoc, 2, 2,  R);

            gsl_matrix_set(Rcoc, 3, 0,  R1dot);
            gsl_matrix_set(Rcoc, 3, 1, -R2dot);
            gsl_matrix_set(Rcoc, 4, 0,  R2dot);
            gsl_matrix_set(Rcoc, 4, 1,  R1dot);

            gsl_matrix_set(Rcoc, 3, 3,  R1);
            gsl_matrix_set(Rcoc, 3, 4, -R2);
            gsl_matrix_set(Rcoc, 4, 3,  R2);
            gsl_matrix_set(Rcoc, 4, 4,  R1);

            gsl_matrix_set(Rcoc, 5, 2,  Rdot);
            gsl_matrix_set(Rcoc, 5, 5,  R);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //VSEM -> VNCSEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +1.0/gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +1.0/gamma);
            break;


        case NCSEM:
            //-----------------------------------------------------
            //VSEM -> NCSEM
            //-----------------------------------------------------
            //RA = SEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, PSEM);
            //RB =  NCSEM_R_SEM
            rot_mat_coc(t, RB, PSEM, NCSEM);
            //Rcoc = NCSEM_R_SEM*SEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB, RA, 0.0, Rcoc);
            break;


        case INEM:
            //-----------------------------------------------------
            //VSEM -> INEM
            //-----------------------------------------------------
            //RA = INSEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, INSEM);
            //RB =  INEM_R_INSEM
            rot_mat_coc(t, RB, INSEM, INEM);
            //Rcoc = INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);

        case VEM:
            //-----------------------------------------------------
            //VSEM -> VEM
            //-----------------------------------------------------
            //RA = INSEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, INSEM);
            //RB =  INEM_R_INSEM
            rot_mat_coc(t, RB, INSEM, INEM);
            //RC =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RC, INEM, VEM);

            //Rcoc = INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = VEM_R_INEM*RB = VEM_R_INEM*INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case PEM:
            //-----------------------------------------------------
            //VSEM -> PEM
            //-----------------------------------------------------
            //RA = INSEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, INSEM);
            //RB =  INEM_R_INSEM
            rot_mat_coc(t, RB, INSEM, INEM);
            //RC =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RC, INEM, VEM);
            //RD =  EM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RD, VEM, PEM);

            //Rcoc = INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= EM_R_VEM*VEM_R_INEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case NCEM:
            //-----------------------------------------------------
            // VSEM -> NCEM
            //-----------------------------------------------------
            //RA = PEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, PEM);
            //RB =  NCEM_R_PEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, PEM, NCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCEM:
            //-----------------------------------------------------
            // VSEM -> VNCEM
            //-----------------------------------------------------
            //RA = VEM_R_VSEM
            rot_mat_coc(t, RA, VSEM, VEM);
            //RB =  VNCEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, VEM, VNCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;


        //-----------------------------------------------------------------
        // VNCSEM
        //-----------------------------------------------------------------
    case VNCSEM:
        switch(outputType)
        {

        case VNCSEM:
            //-----------------------------------------------------
            //VNCSEM -> VNCSEM nothing to do
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VSEM:
            //-----------------------------------------------------
            //VNCEM -> VSEM
            //-----------------------------------------------------
            //Rcoc = VSEM_R_VNCSEM
            gsl_matrix_set(Rcoc, 0, 0,  -gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +gamma);
            break;


        case NCSEM:
            //-----------------------------------------------------
            //VNCSEM -> NCSEM: equivalent to VSEM -> PSEM
            //-----------------------------------------------------
            //Rcoc = NCSEM_R_VNCSEM = PSEM_R_VSEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, -alpha[2]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 0, +alpha[2]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 0, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 1, -alpha[1]/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 2, -alpha[1]/alpha[0]);

            gsl_matrix_set(Rcoc, 3, 3, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, 1.0/alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, 1.0/alpha[0]);

            break;

        case INSEM:
            //-----------------------------------------------------
            //VNCSEM -> INSEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VSEM);
            //RB =  INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //Rcoc = INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case PSEM:
            //-----------------------------------------------------
            //VNCSEM -> PSEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VSEM);
            //RB =  PSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, PSEM);
            //Rcoc = PSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //VNCSEM -> INEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VSEM);
            //RB = INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);

            //Rcoc = RB*RA= INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INEM_R_INSEM*RB = INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case VEM:
            //-----------------------------------------------------
            //VNCSEM -> VEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VSEM);
            //RB = INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);
            //RD =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RD, INEM, VEM);


            //Rcoc = INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VEM_R_INEM*INEM_R_INSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case PEM:
            //-----------------------------------------------------
            //VNCSEM -> PEM
            //-----------------------------------------------------
            //RA = VSEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VSEM);
            //RB = INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);
            //RD =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RD, INEM, VEM);
            //RE =  PEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RE, VEM, PEM);

            //Rcoc = INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VEM_R_INEM*INEM_R_INSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = VEM_R_INEM*INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_VNCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;


        case NCEM:
            //-----------------------------------------------------
            // VNCSEM -> NCEM
            //-----------------------------------------------------
            //RA = PEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, PEM);
            //RB =  NCEM_R_PEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, PEM, NCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCEM:
            //-----------------------------------------------------
            // VNCSEM -> VNCEM
            //-----------------------------------------------------
            //RA = VEM_R_VNCSEM
            rot_mat_coc(t, RA, VNCSEM, VEM);
            //RB =  VNCEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, VEM, VNCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        }
        break;

        //-----------------------------------------------------------------
        // NCSEM
        //-----------------------------------------------------------------
    case NCSEM:
        switch(outputType)
        {

        case NCSEM:
            //-----------------------------------------------------
            //NCSEM -> NCSEM nothing to do here
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //NCSEM -> VNCSEM
            //-----------------------------------------------------
            //Rcoc = VNCSEM_R_NCSEM = VSEM_R_PSEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, +alpha[2]);
            gsl_matrix_set(Rcoc, 4, 0, -alpha[2]);

            gsl_matrix_set(Rcoc, 3, 0, +alpha[1]);
            gsl_matrix_set(Rcoc, 4, 1, +alpha[1]);
            gsl_matrix_set(Rcoc, 5, 2, +alpha[1]);

            gsl_matrix_set(Rcoc, 3, 3, +alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, +alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, +alpha[0]);
            break;


        case PSEM:
            //-----------------------------------------------------
            //NCSEM -> PSEM
            //-----------------------------------------------------
            //Rcoc = SEM_R_NCSEM
            gsl_matrix_set(Rcoc, 0, 0,  -gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +gamma);
            break;

        case VSEM:
            //-----------------------------------------------------
            //NCSEM -> VSEM
            //-----------------------------------------------------
            //RA = VNCSEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VNCSEM);
            //RB =  VSEM_R_VNCSEM
            rot_mat_coc(t, RB, VNCSEM, VSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INSEM:
            //-----------------------------------------------------
            //NCSEM -> INSEM
            //-----------------------------------------------------
            //RA = VNCSEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VNCSEM);
            //RB =  VSEM_R_VNCSEM
            rot_mat_coc(t, RB, VNCSEM, VSEM);
            //RC =  INSEM_R_VSEM
            rot_mat_coc(t, RC, VSEM, INSEM);

            //Rcoc = VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INSEM_R_VSEM*RB = INSEM_R_VSEM*VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //NCSEM -> INEM
            //-----------------------------------------------------
            //RA = VNCSEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VNCSEM);
            //RB =  VSEM_R_VNCSEM
            rot_mat_coc(t, RB, VNCSEM, VSEM);
            //RC =  INSEM_R_VSEM
            rot_mat_coc(t, RC, VSEM, INSEM);
            //RD =  INEM_R_INSEM
            rot_mat_coc(t, RD, INSEM, INEM);

            //Rcoc = VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case VEM:
            //-----------------------------------------------------
            //NCSEM -> VEM
            //-----------------------------------------------------
            //RA = VNCSEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VNCSEM);
            //RB =  VSEM_R_VNCSEM
            rot_mat_coc(t, RB, VNCSEM, VSEM);
            //RC =  INSEM_R_VSEM
            rot_mat_coc(t, RC, VSEM, INSEM);
            //RD =  INEM_R_INSEM
            rot_mat_coc(t, RD, INSEM, INEM);
            //RE =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RE, INEM, VEM);

            //Rcoc = VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;

        case PEM:
            //-----------------------------------------------------
            //NCSEM -> PEM
            //-----------------------------------------------------
            //RA = VNCSEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VNCSEM);
            //RB =  VSEM_R_VNCSEM
            rot_mat_coc(t, RB, VNCSEM, VSEM);
            //RC =  INSEM_R_VSEM
            rot_mat_coc(t, RC, VSEM, INSEM);
            //RD =  INEM_R_INSEM
            rot_mat_coc(t, RD, INSEM, INEM);
            //RE =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RE, INEM, VEM);
            //RF =  PEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RF, VEM, PEM);

            //Rcoc = VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= INEM_R_INSEM*INSEM_R_VSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_VNCSEM*VNCSEM_R_NCSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //RC = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, RC);

            //Rcoc = RF*RC = RE*RF*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RF , RC, 0.0, Rcoc);
            break;


        case NCEM:
            //-----------------------------------------------------
            // NCSEM -> NCEM
            //-----------------------------------------------------
            //RA = PEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, PEM);
            //RB =  NCEM_R_PEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, PEM, NCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCEM:
            //-----------------------------------------------------
            // NCSEM -> VNCEM
            //-----------------------------------------------------
            //RA = VEM_R_NCSEM
            rot_mat_coc(t, RA, NCSEM, VEM);
            //RB =  VNCEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, VEM, VNCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        }
        break;

        //-----------------------------------------------------------------
        // PSEM
        //-----------------------------------------------------------------
    case PSEM:
        switch(outputType)
        {

        case PSEM:
            //-----------------------------------------------------
            //PSEM -> PSEM nothing to do here
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VSEM:
            //-----------------------------------------------------
            //PSEM -> VSEM
            //-----------------------------------------------------
            //Rcoc = VSEM_R_PSEM = VNCSEM_R_NCSEM
            gsl_matrix_set(Rcoc, 0, 0,  1.0);
            gsl_matrix_set(Rcoc, 1, 1,  1.0);
            gsl_matrix_set(Rcoc, 2, 2,  1.0);

            gsl_matrix_set(Rcoc, 3, 1, +alpha[2]);
            gsl_matrix_set(Rcoc, 4, 0, -alpha[2]);

            gsl_matrix_set(Rcoc, 3, 0, +alpha[1]);
            gsl_matrix_set(Rcoc, 4, 1, +alpha[1]);
            gsl_matrix_set(Rcoc, 5, 2, +alpha[1]);

            gsl_matrix_set(Rcoc, 3, 3, +alpha[0]);
            gsl_matrix_set(Rcoc, 4, 4, +alpha[0]);
            gsl_matrix_set(Rcoc, 5, 5, +alpha[0]);
            break;


        case NCSEM:
            //-----------------------------------------------------
            //PSEM -> NCSEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 1, 1,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 2, 2,  +1.0/gamma);
            gsl_matrix_set(Rcoc, 3, 3,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 4, 4,  -1.0/gamma);
            gsl_matrix_set(Rcoc, 5, 5,  +1.0/gamma);
            break;

        case VNCSEM:
            //-----------------------------------------------------
            //PSEM -> VNCSEM
            //-----------------------------------------------------
            //RA = VSEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VSEM);
            //RB =  VNCSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, VNCSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case INSEM:
            //-----------------------------------------------------
            //PSEM -> INSEM
            //-----------------------------------------------------
            //RA = VSEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VSEM);
            //RB =  INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case INEM:
            //-----------------------------------------------------
            //SEM -> INEM
            //-----------------------------------------------------
            //RA = VSEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VSEM);
            //RB =  INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);

            //Rcoc = RB*RA= INSEM_R_VSEM*VSEM_R_PSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RB = Rcoc
            gsl_matrix_memcpy(RB, Rcoc);
            //Rcoc = INEM_R_INSEM*RB = INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_PSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RC , RB, 0.0, Rcoc);
            break;


        case VEM:
            //-----------------------------------------------------
            //PSEM -> VEM
            //-----------------------------------------------------
            //RA = VSEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VSEM);
            //RB =  INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);
            //RD =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RD, INEM, VEM);


            //Rcoc = RB*RA= INSEM_R_VSEM*VSEM_R_PSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VEM_R_INEM*INEM_R_INSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;

        case PEM:
            //-----------------------------------------------------
            //PSEM -> PEM
            //-----------------------------------------------------
            //RA = VSEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VSEM);
            //RB =  INSEM_R_VSEM
            rot_mat_coc(t, RB, VSEM, INSEM);
            //RC =  INEM_R_INSEM
            rot_mat_coc(t, RC, INSEM, INEM);
            //RD =  VEM_R_INEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RD, INEM, VEM);
            //RE =  PEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RE, VEM, PEM);

            //Rcoc = RB*RA= INSEM_R_VSEM*VSEM_R_PSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            //RA = Rcoc
            gsl_matrix_memcpy(RA, Rcoc);

            //RB= VEM_R_INEM*INEM_R_INSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RD , RC, 0.0, RB);

            //RD = RB*RA = VEM_R_INEM*INEM_R_INSEM*INSEM_R_VSEM*VSEM_R_PSEM
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, RD);

            //Rcoc = RE*RD
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RE , RD, 0.0, Rcoc);
            break;


        case NCEM:
            //-----------------------------------------------------
            // PSEM -> NCEM
            //-----------------------------------------------------
            //RA = PEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, PEM);
            //RB =  NCEM_R_PEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, PEM, NCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;


        case VNCEM:
            //-----------------------------------------------------
            // PSEM -> VNCEM
            //-----------------------------------------------------
            //RA = VEM_R_PSEM
            rot_mat_coc(t, RA, PSEM, VEM);
            //RB =  VNCEM_R_VEM - time must be in EM units!
            rot_mat_coc(t/SEML.us_em.ns, RB, VEM, VNCEM);
            //Rcoc = RB*RA
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, RB , RA, 0.0, Rcoc);
            break;
        }
        break;

        //-----------------------------------------------------------------
        // INEM
        //-----------------------------------------------------------------
    case INEM:
        switch(outputType)
        {
        case INEM:
            //-----------------------------------------------------
            //INEM -> INEM nothing to do here
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case VEM:
            //-----------------------------------------------------
            //INEM -> VEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0, +a1);
            gsl_matrix_set(Rcoc, 0, 1, +a2);
            gsl_matrix_set(Rcoc, 1, 0, -a2);
            gsl_matrix_set(Rcoc, 1, 1, +a1);

            gsl_matrix_set(Rcoc, 2, 2,  a6);

            gsl_matrix_set(Rcoc, 3, 0,  a1dot);
            gsl_matrix_set(Rcoc, 3, 1,  a2dot);
            gsl_matrix_set(Rcoc, 4, 0, -a2dot);
            gsl_matrix_set(Rcoc, 4, 1,  a1dot);

            gsl_matrix_set(Rcoc, 3, 3,  a1);
            gsl_matrix_set(Rcoc, 3, 4,  a2);
            gsl_matrix_set(Rcoc, 4, 3, -a2);
            gsl_matrix_set(Rcoc, 4, 4,  a1);

            gsl_matrix_set(Rcoc, 5, 2,  a6dot);
            gsl_matrix_set(Rcoc, 5, 5,  a6);
            break;

        case INSEM:
            //-----------------------------------------------------
            //INEM -> INSEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0,  1.0/SEML.us_em.as);
            gsl_matrix_set(Rcoc, 1, 1,  1.0/SEML.us_em.as);
            gsl_matrix_set(Rcoc, 2, 2,  1.0/SEML.us_em.as);

            gsl_matrix_set(Rcoc, 3, 3,  1.0/(SEML.us_em.as*SEML.us_em.ns));
            gsl_matrix_set(Rcoc, 4, 4,  1.0/(SEML.us_em.as*SEML.us_em.ns));
            gsl_matrix_set(Rcoc, 5, 5,  1.0/(SEML.us_em.as*SEML.us_em.ns));
            break;

        default:
            perror("rot_mat_coc. InputType = INEM. Unknown outputType.");
            break;
        }
        break;
    case INSEM:
        switch(outputType)
        {
        case INSEM:
            //-----------------------------------------------------
            //INSEM -> INSEM nothing to do here
            //-----------------------------------------------------
            gsl_matrix_set_identity(Rcoc);
            break;

        case INEM:
            //-----------------------------------------------------
            //INSEM -> INEM
            //-----------------------------------------------------
            gsl_matrix_set(Rcoc, 0, 0, SEML.us_em.as);
            gsl_matrix_set(Rcoc, 1, 1, SEML.us_em.as);
            gsl_matrix_set(Rcoc, 2, 2, SEML.us_em.as);

            gsl_matrix_set(Rcoc, 3, 3, SEML.us_em.as*SEML.us_em.ns);
            gsl_matrix_set(Rcoc, 4, 4, SEML.us_em.as*SEML.us_em.ns);
            gsl_matrix_set(Rcoc, 5, 5, SEML.us_em.as*SEML.us_em.ns);
            break;

        case VSEM:
            //-----------------------------------------------------
            //INSEM -> VSEM
            //-----------------------------------------------------
            //Set coefficient in Rcoc
            gsl_matrix_set(Rcoc, 0, 0, +A1);
            gsl_matrix_set(Rcoc, 0, 1, +A2);
            gsl_matrix_set(Rcoc, 1, 0, -A2);
            gsl_matrix_set(Rcoc, 1, 1, +A1);

            gsl_matrix_set(Rcoc, 2, 2,  A6);

            gsl_matrix_set(Rcoc, 3, 0,  A1dot);
            gsl_matrix_set(Rcoc, 3, 1,  A2dot);
            gsl_matrix_set(Rcoc, 4, 0, -A2dot);
            gsl_matrix_set(Rcoc, 4, 1,  A1dot);

            gsl_matrix_set(Rcoc, 3, 3,  A1);
            gsl_matrix_set(Rcoc, 3, 4,  A2);
            gsl_matrix_set(Rcoc, 4, 3, -A2);
            gsl_matrix_set(Rcoc, 4, 4,  A1);

            gsl_matrix_set(Rcoc, 5, 2,  A6dot);
            gsl_matrix_set(Rcoc, 5, 5,  A6);
            break;
        default:
            perror("rot_mat_coc. InputType = INSEM. Unknown outputType.");
            break;
        }
        break;

    }

    //=====================================================================
    // 4. Reset the focus in SEML, if necessary
    //=====================================================================
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);

    //=====================================================================
    // 5. Free the temporary objects
    //=====================================================================
    gsl_matrix_free(RA);
    gsl_matrix_free(RB);
    gsl_matrix_free(RC);
    gsl_matrix_free(RD);
    gsl_matrix_free(RE);
    gsl_matrix_free(RF);
}


