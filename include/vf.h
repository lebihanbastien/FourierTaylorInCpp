#ifndef VF_H_INCLUDED
#define VF_H_INCLUDED

#include "env.h"
#include "ephemerides.h"
#include "gslc.h"
#include "pmcoc.h"
#include "eminsem.h"
#include "ode.h"
#include <gsl/gsl_blas.h>

extern "C" {
#include "gnuplot_i.h"
}

//Pointer to a given vector field
typedef int (*vfptr)(double, const double*, double*, void*);

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

//========================================================================================
//
// Selection of the vector fields
//
//========================================================================================
/**
 *  \brief Select the vector field associated to dcs (default coordinate system)
 **/
vfptr ftc_select_vf(int dcs, int nvar);

//========================================================================================
//
// Vector fields linked to EPHEMERIDES dynamics
//
//========================================================================================
//----------------------------------------------------------------------------------------
// Subroutines
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes the acceleration vector of the Earth.
 **/
void acc_earth_from_vf(double et, double Ae[3],  double Rj[11][6], QBCP_L* qbp);

/**
 *   \brief Update the positions of the primaries of the Solar system at epoch et.
 **/
void usspos(double et, double Rj[11][6], QBCP_L* qbp);

/**
 *  \brief Computes the acceleration vectors A1 and A2 of the two primaries that defines the synodic coordinates.
 *         They are computed from the newtonian vector field, i.e. the gravitational influence of the massive bodies of the solar system.
 *         If the smaller primary is the Earth-Moon Barycenter, the acceleration is computed as the barycenter of the acceleration of the Earth and the Moon.
 **/
void acc_from_vf(double et, int pos1, int pos2, double A1[3], double A2[3],
                 double R1[3], double R2[3], double Rj[11][6], QBCP_L* qbp);

/**
 *  \brief Computes the jerk vectors A1 and A2 of the two primaries that defines the synodic coordinates.
 *         They are computed from the derivative of the newtonian vector field, i.e. the gravitational influence of the massive bodies of the solar system.
 *         If the smaller primary is the Earth-Moon Barycenter, the jerk is computed as the barycenter of the jerk of the Earth and the Moon.
 **/
void jerk_from_vf(double et, int pos1, int pos2, double J1[3], double J2[3],
                  double R1[3], double R2[3], double Rj[11][6], QBCP_L* qbp);

//----------------------------------------------------------------------------------------
// Ecliptic coordinates for JPL ephemerides (X, X')
//----------------------------------------------------------------------------------------
/**
 * \brief Vector field of the solar system (JPL ephemerides) in ecliptic coordinates
 *       (SI units) Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO
 **/
int jpl_vf(double et, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in ecliptic coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO
 **/
int jpl_vf_var(double et, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in
 *        earth-centered inertial coordinates (ECI). Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 **/
int jpl_vf_eci(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in earth-centered inertial
 *        coordinates with variationnal equations. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 **/
int jpl_vf_eci_var(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in
 *        normalized earth-centered inertial coordinates (NECI). Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 **/
int jpl_vf_neci(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in NECI
 *        coordinates with variationnal equations. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 *        Note: for now the normalization is NOT taken into account because
 *        the stepper is going wrong when the state is normalized.
 **/
int jpl_vf_neci_var(double t, const double y[], double f[], void* params_void);

//----------------------------------------------------------------------------------------
// Synodic coordinates for JPL ephemerides (x, x')
//----------------------------------------------------------------------------------------
/**
 * \brief Vector field of the solar system (JPL ephemerides) in synodical coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO (11 bodies)
 **/
int jpl_vf_syn(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the solar system (JPL ephemerides) in synodical coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO (11 bodies)
 **/
int jpl_vf_syn_var(double t, const double y[], double f[], void* params_void);

//========================================================================================
//
// Vector fields linked to QBCP dynamics
//
//========================================================================================
//----------------------------------------------------------------------------------------
// Inertial SEM coordinates, (X, V)
//----------------------------------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in SEM inertial coordinates
 **/
int qbcp_vf_insem(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the QBCP in SEM inertial coordinates
 **/
int qbcp_vf_insem_var(double t, const double y[], double f[], void* params_void);


//----------------------------------------------------------------------------------------
// Earth-centered Inertial SEM coordinates, (X, V)
//----------------------------------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in EC SEM inertial coordinates
 **/
int qbcp_vf_ecisem(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the QBCP in EC SEM inertial coordinates
 **/
int qbcp_vf_ecisem_var(double t, const double y[], double f[], void* params_void);

//----------------------------------------------------------------------------------------
// Continuation between ECISEM and NECI (JPL)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. Variational equations included.
 **/
int qbcp_ecisem_cont_necijpl(double t, const double y[], double f[], void *params_void);

//----------------------------------------------------------------------------------------
// Continuation between INSEM and ECISEM
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field for continuation procedure. The result is f = (1.0 - epsilon)*f1 + espilon*f2. Variational equations included.
 **/
int qbcp_insem_cont_ecisem(double t, const double y[], double f[], void *params_void);

//----------------------------------------------------------------------------------------
// Inertial EM coordinates, (X, V)
//----------------------------------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in EM inertial coordinates
 **/
int qbcp_vf_inem(double t, const double y[], double f[], void* params_void);

/**
 * \brief Vector field of the QBCP in EM inertial coordinates
 **/
int qbcp_vf_inem_var(double t, const double y[], double f[], void* params_void);
;

//----------------------------------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------------------------------
//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbcp_vfn(double t, const double y[], double f[], void* params_void);

//Update the normalized vector field of the state. Note that alpha[14] (alpha15) is zero for the QBCP
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);

//----------------------------------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------------------------------
// Vector field of the QBCP with EM units and Normalized-Li centered coordinates
// (no variationnal equations). The state is here (x, xdot), and not (x, px)
int qbcp_vfn_xv(double t, const double y[], double f[], void* params_void);

// Update the normalized vector field of the state. The state is here (x, xdot),
// and not (x, px)
int vfn_state_xv(const double y[], double f[], double alpha[], double alphad[],
                 double ps[], double pe[], double pm[],
                 double qps2, double qpe2, double qpm2,
                 double ms, double me, double mm,
                 double gamma);

//----------------------------------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbcp_vf(double t, const double y[], double f[], void* params_void);

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized). Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm );

//----------------------------------------------------------------------------------------
// SYS coordinates, (X, V)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM coordinates. The state is here (x, xdot), and not (x, px)
 **/
int qbcp_vf_varnonlin_xv(double t, const double y[], double f[], void* params_void);

// Vector field of the QBCP with EM units and system units centered coordinates
// (no variationnal equations). The state is here (x, xdot), and not (x, px)
int qbcp_vf_xv(double t, const double y[], double f[], void* params_void);


// Update the vector field of the state. The state is here (x, xdot),
// and not (x, px)
int vf_state_xv(const double y[], double f[], double alpha[], double alphad[],
                 double ps[], double pe[], double pm[],
                 double qps2, double qpe2, double qpm2,
                 double ms, double me, double mm);

/**
 *  \brief Update the system variational equation matrix
 **/
int vf_stm_xv(const double y[], gsl_matrix* Q, double alpha[], double alphad[],
               double ps[], double pe[], double pm[],
               double qps2, double qpe2, double qpm2,
               double ms, double me, double mm);

//----------------------------------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbcp_vfn_varnonlin(double t, const double y[], double f[], void* params_void);

/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
void set_vareq_matrix(gsl_matrix* Q, double b[], double alpha[]);
/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm(const double y[], gsl_matrix* Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma);

//----------------------------------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates. The state is here (x, xdot), and not (x, px)
 **/
int qbcp_vfn_varnonlin_xv(double t, const double y[], double f[], void* params_void);

/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha
 **/
void set_vareq_matrix_xv(gsl_matrix* Q, double b[], double alpha[], double alphad[]);

/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm_xv(const double y[], gsl_matrix* Q, double alpha[], double alphad[],
               double ps[], double pe[], double pm[],
               double qps2, double qpe2, double qpm2,
               double ms, double me, double mm,
               double gamma);

//----------------------------------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame. Variational equations included.
 **/
int qbcp_vf_varnonlin(double t, const double y[], double f[], void* params_void);

/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double y[], gsl_matrix* Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm);

//========================================================================================
//
// SUBROUTINES
//
//========================================================================================
/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j);

//========================================================================================
//
//          Reduced vector field
//
//========================================================================================
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

//========================================================================================
//
// Hamiltonians
//
//========================================================================================
/**
 *  \brief Hamiltonian of the QBCP with SYS units and SYS coordinates. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
double qbfbp_H(double t, const double y[], void *params_void);

/**
 *  \brief Hamiltonian of the QBCP with SEM units and SEM coordinates
 **/
double qbfbp_H_SEM(double t, const double y[], void *params_void);

/**
 *  \brief Hamiltonian of the QBCP with EM units and Normalized-Centered coordinates. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
double qbfbp_Hn(double t, const double y[], void *params_void);

//========================================================================================
//
// Testing the QBTBP
//
//========================================================================================
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
int qbtbp_derivatives(double t, const double y[], double f[], void *params);

/**
 * \brief Computes the semi-analytical position of the Earth, the Moon, and the Sun from the coefficients alpha and delta of the EM and SEM vector field.
 **/
void primary_analytical_position(int prim, double tem, double yEM_to_IN[], double ySEM_to_IN[], QBCP_L &qbcp_l);

/**
 * \brief Computes the "true" integrated position of the Earth, the Moon, and the Sun from the integrated bicircular movement contained in Z(t) and z(t).
 **/
void primary_integrated_position(int prim, double yIN[], cdouble z0, cdouble Z0, cdouble zdot0, cdouble Zdot0, QBCP_L &qbcp_l);

/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, QBCP_L &qbcp_l);


#endif // VF_H_INCLUDED
