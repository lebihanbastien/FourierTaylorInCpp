#ifndef GSLC_H_INCLUDED
#define GSLC_H_INCLUDED

/**
 * \file gslc.h
 * \brief Additional operations on GSL objects.
 * \author BLB
 * \date July 2016
 * \version 1.0
 */

//std
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <sstream>
#include <math.h>

//GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>

//Custom
#include "parameters.h"


//----------------------------------------------------------------------------------------
//
// Vector <--> matrix
//
//----------------------------------------------------------------------------------------
/**
 * \brief Transform a vector into a matrix with a given shift in the initial vector
 **/
void gslc_vectorToMatrix(gsl_matrix *m, const double y[], int rows, int columns, int shift);

/**
 * \brief Transform a matrix into a vector with a given shift in the final vector
 **/
void gslc_matrixToVector(double y[], const gsl_matrix *m, int rows, int columns, int shift);

//----------------------------------------------------------------------------------------
// Printing a real matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a real matrix
 **/
void gslc_matrix_printf(const gsl_matrix *M);

/**
 * \brief Print an approximation of a real matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_printf_approx(const gsl_matrix *M);

/**
 * \brief Print a real matrix into a txt file.
 **/
void gslc_matrix_fprintf(const gsl_matrix *M, char* fileName);

//----------------------------------------------------------------------------------------
// Complex numbers
//----------------------------------------------------------------------------------------
/**
 * \brief Create a gsl_complex from 2 doubles
 **/
gsl_complex gslc_complex(double real, double imag);

/**
 * \brief Create a gsl_complex from a complex number
 **/
gsl_complex gslc_complex(cdouble x);

/**
 * \brief Create a complex number from a gsl_complex
 **/
cdouble gslc_complex(gsl_complex c);

//----------------------------------------------------------------------------------------
// Printing a complex matrix
//----------------------------------------------------------------------------------------
/**
 * \brief Print a complex matrix.
 **/
void gslc_matrix_complex_printf(const gsl_matrix_complex *M);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_real(const gsl_matrix_complex *M);

/**
 * \brief Print the imag part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_printf_imag(const gsl_matrix_complex *M);

/**
 * \brief Print an approximation of a complex matrix: the coefficients under 1e-5 are displayed as null.
 **/
void gslc_matrix_complex_printf_approx(const gsl_matrix_complex *M);

/**
 * \brief Print a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf(const gsl_matrix_complex *M, char* fileName);

/**
 * \brief Print a complex matrix into a txt file, separating the real and imaginary part, in order to be easier to read.
 **/
void gslc_matrix_complex_fprintf_pretty(const gsl_matrix_complex *M, char* fileName);

/**
 * \brief Print the real part of a complex matrix into a txt file.
 **/
void gslc_matrix_complex_fprintf_real(const gsl_matrix_complex *M, char* fileName);

//----------------------------------------------------------------------------------------
// Complex to real matrices
//----------------------------------------------------------------------------------------
void gslc_matrix_complex_to_matrix(const gsl_matrix_complex *Mc, gsl_matrix *M);

//----------------------------------------------------------------------------------------
// Allocation
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize and return an array of GSL matrices
 **/
gsl_matrix** gslc_matrix_array_calloc(int size1, int size2, int M);

/**
 * \brief Free an array of GSL matrices
 **/
void gslc_matrix_array_free(gsl_matrix ** P, int M);

/**
 * \brief Initialize and return an array of GSL vectors
 **/
gsl_vector** gslc_vector_array_calloc(int size1, int M);

/**
 * \brief Free an array of GSL vectors
 **/
void gslc_vector_array_free(gsl_vector ** P, int M);

//----------------------------------------------------------------------------------------
// Norm
//----------------------------------------------------------------------------------------
/**
 *  \brief L2 norm of a complex matrix: res = ||Mc||_2
 **/
double gslc_matrix_complex_L2(const gsl_matrix_complex *Mc);

/**
 *  \brief L2 norm of a difference of complex matrices: res = ||M1c - M2c ||_2
 **/
double gslc_matrix_complex_diff_L2(const gsl_matrix_complex *M1c, const gsl_matrix_complex *M2c);

//----------------------------------------------------------------------------------------
//   Operations on simple pointers (vectors of double)
//----------------------------------------------------------------------------------------
/**
 *  \brief Copy src into dest, both of size n.
 **/
void vector_memcpy(double *dest, const double *src, int n);

/**
 *  \brief Copy src into dest, both of size 6.
 **/
void state_memcpy(double *dest, const double *src);

/**
 *  \brief Copy src into dest, both of size 6 x (last_indix+1)
 **/
void vecstate_memcpy(double **dest, const double **src, int last_indix);

#endif // GSLC_H_INCLUDED
