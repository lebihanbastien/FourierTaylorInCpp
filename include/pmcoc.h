#ifndef PMCOC_H_INCLUDED
#define PMCOC_H_INCLUDED

#include "Oftsc.h"
#include "env.h"
#include "matrix.h"
#include "gslc.h"
#include "eminsem.h"

/**
 * \file pmcoc.h
 * \brief Change of coordinates for the parameterization methods.
 * \author BLB.
 * \date May 2015
 * \version 1.0
 *
 *  Two types of change of coordinates are implemented:
 *              1. Changes between different manifold coordinates (Real Center to Complex Center...).
 *              2. Evaluations of the parameterization (Real Center to Normalized-Centered and projection the other way around).
 *
 */

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates in the manifolds
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief from CCM to RCM coordinates
 **/
void CCMtoRCM(const cdouble s1[], double si[], int nv);

/**
 *  \brief from RCM to CCM coordinates
 **/
void RCMtoCCM(const double si[], cdouble s1[], int nv);

/**
 *  \brief from RCM to CCM coordinates, with real and imag part stored separately
 **/
void RCMtoCCM8(const double si[], double s0d[], int nv);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to RCM coordinates
 **/
void CCM8toRCM(const double s0d[], double si[], int nv);

/**
 *  \brief from CCM coordinates, with real and imag part stored separately, to CCM coordinates
 **/
void CCM8toCCM(const double s0d[], cdouble s1[], int nv);

/**
 *  \brief from CCM coordinates to CCM coordinates, with real and imag part stored separately.
 **/
void CCMtoCCM8(const cdouble s1[], double s0d[], int nv);

/**
 *  \brief from TFC to TF coordinates
 **/
void TFCtoTF(const cdouble s1[6], double si[6]);


//---------------------------------------------------------------------------------------------------------------------------------------
// Change of coordinates between the systems
//---------------------------------------------------------------------------------------------------------------------------------------
///**
// *  \brief COC: Normalized-Centered coordinates to system coordinates. Use in priority instead of NCEMmtoEMm or NCSEMmtoSEMm.
// **/
//void NCtoSYS(double t, const double yNC[], double yEM[], void *params_void);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluation of the pm
//
//---------------------------------------------------------------------------------------------------------------------------------------
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
             bool isGS);

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
                  bool isGS);

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
                  bool isGS);

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
              bool isGS);



/**
 *  \brief Apply the change of variables in zIN/zOut. (OFS version).
 *       The change of variables is of the form: zOut = PC * zIN + V
 **/
void applyCOC(const matrix<Ofsc> &PC,
              const vector<Ofsc> &V,
              const vector<Ofsc> &zIn,
              vector<Ofsc> &zOut);


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
              bool isGS);

/**
 *  \brief Projection of the current NC state on the central manifold, via CCM coordinates
 **/
void NCprojCCM(const double z[], const double t, const double n, const int ofs_order,
               const matrix<Ofsc> &CQ, const vector<Ofsc> &V,
               double omega1, double omega3, cdouble sc[], int nv);

/**
 *  \brief Projection of the current TFC state on the central manifold, via CCM coordinates
 **/
void TFCprojCCM(const cdouble zh[], double omega1, double omega3, cdouble sc[], int nv);

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
                  bool isGS);

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
                  bool isGS);
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
                  bool isGS);

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
                double f8[]);


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Evaluation at time t
//
//---------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//Evaluation of parts or all the alpha/beta routines at a given time
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluation of Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoef(double *alpha, double t, double omega, int order, double *params, int number);

/**
 *  \brief Evaluation of the time derivatives Fourier series given as an array of coefficients, at a given time t.
 */
void evaluateCoefDerivatives(double *alpha, double t, double omega, int order, double *params, int number);

//-----------------------------------------------------------------------------
//Evalution of individual alpha/beta
//or derivative of alpha/beta of a given type (Odd or Even)
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) cos(k \omega t)  \f$.
 */
double evaluateEven(int order, double *coef, double *cR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N - k \omega coef(k) sin(k \omega t)  \f$.
 */
double evaluateEvenDerivative(double omega, int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N coef(k) sin(k \omega t)  \f$.
 */
double evaluateOdd(int order, double *coef, double *sR);

/**
 *  \brief Evaluate the sum \f$ \sum_{k = 0}^N  k \omega coef(k) cos(k \omega t)  \f$.
 */
double evaluateOddDerivative(double omega, int order, double *coef, double *cR);

//-----------------------------------------------------------------------------
// Evaluate the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Evaluate z(t), with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evz(const Ofsc& zr, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate dz(t)/dt, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzdot(const Ofsc& zr, const Ofsc& ztdot, double t, double n, double ni, double ai);

/**
 *  \brief Evaluate d2z(t)/dt2, with \f$ z(t) = e^{it} z_r(t) \f$ in Earth-Moon units.
 */
cdouble evzddot(const Ofsc& zr, const Ofsc& ztdot, const Ofsc& ztddot, double t, double n, double ni, double ai);


/**
 *  \brief RCOC: from inputType to outputType. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void rot_mat_coc(double t, gsl_matrix* Rcoc, int inputType, int outputType);

/**
 *  \brief RCOC: from inputType to outputType.
 *         Used ONLY in rot_mat_coc, since specific checks
 *         are made in this routine prior to any computations.
 **/
void qbcp_coc_fwrk(double t, gsl_matrix* Rcoc, int inputType, int outputType);

#endif // PMCOC_H_INCLUDED
