#ifndef VF_H_INCLUDED
#define VF_H_INCLUDED

#include "env.h"
#include "gslc.h"
#include "pmcoc.h"
#include "eminsem.h"
#include "ode.h"

#include <gsl/gsl_blas.h>


// Structure for the reduced vector field
typedef struct RVF RVF;
struct RVF
{
    vector<Oftsc>* fh;   //Reduced vector field
    Ofsc* ofs;           //Auxiliary Ofs object
    int order;           //Order of the evaluation     (<= OFTS_ORDER)
    int ofs_order;       //Order of the ofs evaluation (<= OFS_ORDER)
    double n;            //Pulsation of the QBCP
    int reduced_nv;      //Number of reduced variables
};


//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void);

int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);

//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations). The state is here (x, xdot), and not (x, px)
int qbfbp_vfn_novar_xv(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates. The state is here (x, xdot), and not (x, px)
 **/
int qbfbp_vfn_varnonlin_xv(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame. Variational equations included.
 **/
int qbfbp_vf_varnonlin(double t, const double y[], double f[], void *params_void);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Integration without STM, without Normalization
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbfbp_vf(double t, const double y[], double f[], void *params_void);

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized)
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm );


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Integration of the pm
//
//---------------------------------------------------------------------------------------------------------------------------------------
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
int qbfbp_fh(double t, const double y[], double f[], void *params_void);

//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Hamiltonians
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian of the QBCP with SYS units and SYS coordinates
 **/
double qbfbp_H(double t, const double y[], void *params_void);

/**
 *  \brief Hamiltonian of the QBCP with SEM units and SEM coordinates
 **/
double qbfbp_H_SEM(double t, const double y[], void *params_void);



//--------------------------------------------------------------------------------------------------------------------------------------------
//
// Testing the QBTBP
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, QBCP_L &qbcp_l);

#endif // VF_H_INCLUDED
