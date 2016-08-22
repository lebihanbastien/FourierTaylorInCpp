#ifndef VF_H_INCLUDED
#define VF_H_INCLUDED

#include "env.h"
#include "pmcoc.h"
#include <gsl/gsl_blas.h>

//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void *params_void);

int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma);

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void *params_void);

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
void NCtoEM(double t, const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void EMtoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------
/**
 *  \brief COC: Normalized-Centered coordinates to system coordinates. Use in priority instead of NCtoEM or NCtoSEM.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from system coordinates to Normalized-Centered coordinates
 **/
void SYStoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from Sun-Earth-Moon coordinates to Normalized-Centered coordinates
 **/
void SEMtoNC(double t, const double ySEM[], double yNC[], QBCP_L *qbp);

/**
 *  \brief COC: from Normalized-Centered coordinates to Sun-Earth-Moon coordinates
 **/
void NCtoSEM(double t, const double yNC[], double ySEM[], QBCP_L *qbp);

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
// Vector <--> matrix
//
//--------------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift);

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift);

#endif // VF_H_INCLUDED
