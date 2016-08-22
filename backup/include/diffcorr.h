#ifndef DIFFCORR_H_INCLUDED
#define DIFFCORR_H_INCLUDED

#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_blas.h>
#include "env.h"
#include "vf.h"

extern "C"{
 #include "nrutil.h"
 #include "gnuplot_i.h"
}


/**
 * \brief Differential correction scheme
 **/
int differential_correction_mns(double *ystart, double *ydest, double t0, double *t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);

/**
 * \brief Differential correction scheme, with fixed time
 **/
int differential_correction_ft(double *ystart, double *ydest, double t0, double *t1, double eps_diff, gsl_odeiv2_driver *d, int N, int isPlotted);

/**
 *  \brief Same as odePlot, but with more flexibility on the parameters: can choose the title, the line style, the line type and the line color.
 **/
int odePlotGen(const double y[],
               int N, double t0, double t1,
               gsl_odeiv2_driver *d,
               gnuplot_ctrl  *h1,
               int Npoints,
               char const *title,
               char const *ls,
               char const *lt,
               char const *lw,
               int lc);
#endif // DIFFCORR_H_INCLUDED
