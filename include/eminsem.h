#ifndef EMINSEM_H_INCLUDED
#define EMINSEM_H_INCLUDED

/**
 * \file  eminsem.h
 * \brief Contains all the routines to perform changes of coordinates between the EM and SEM frameworks. Including
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "vf.h"

//=============================================================================================
//
// Super COCs
//
//=============================================================================================
/**
 *  \brief COC: from inputType to outputType. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc(double t, const double y0[], double yout[], int inputType, int outputType);

/**
 *  \brief COC: from inputType to outputType. With a time update in tout. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc(double t, const double y0[], double yout[], double *tout, int inputType, int outputType);

/**
 *  \brief COC: from inputType to outputType, in vector form. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc_vec(double **y0, double *t0, double **yout, double *tout, int N, int inputType, int outputType);

/**
 *  \brief COC: from inputType to outputType, in vector form. Version with no time vector output. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc_vec(double **y0, double *t0, double **yout, int N, int inputType, int outputType);

/**
 *  \brief COC: from inputType to outputType.
 *         Used ONLY in qbcp_coc and qbcp_coc_vec, since specific checks
 *         are made in these routines prior to any computations.
 **/
void qbcp_coc_fwrk(double t, const double y0[], double yout[], int inputType, int outputType);

//=============================================================================================
//
// Change of coordinates: NC <-> SYS
//
//=============================================================================================
//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------
/**
 *  \brief COC: NC coordinates to SYSTEM coordinates. Use in priority instead of NCEMmtoEMm or NCSEMmtoSEMm.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from SYSTEM coordinates to NC coordinates. Use in priority instead of EMtoN or SEMmtoNCSEMm.
 **/
void SYStoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> SEM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from SEM coordinates to NC coordinates
 **/
void SEMmtoNCSEMm(double t, const double ySEM[], double yNC[], QBCP_L *qbp);

/**
 *  \brief COC: from NC coordinates to SEM coordinates
 **/
void NCSEMmtoSEMm(double t, const double yNC[], double ySEM[], QBCP_L *qbp);

//-----------------------------------------------------------------------------
// COC: NC <--> EM
//-----------------------------------------------------------------------------
/**
 *  \brief COC: from NC coordinates to EM coordinates
 **/
void NCEMmtoEMm(double t, const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from EM coordinates to NC coordinates
 **/
void EMmtoNCEMm(double t, const double yEM[], double yNC[], QBCP_L *qbp);


//-----------------------------------------------------------------------------
// COC: VNC <--> NC
//-----------------------------------------------------------------------------
/**
 *  \brief COC: NC coordinates to SYSTEM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCvtoSYSv(const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from SYSTEM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void SYSvtoNCv(const double yEM[], double yNC[], QBCP_L *qbp);

/**
 *  \brief COC: NC coordinates to EM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCEMvtoEMv(const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from EM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void EMvtoNCEMv(const double yEM[], double yNC[], QBCP_L *qbp);

/**
 *  \brief COC: NC coordinates to SEM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCSEMvtoSEMv(const double yNC[], double yEM[], QBCP_L *qbp);

/**
 *  \brief COC: from SEM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void SEMvtoNCSEMv(const double yEM[], double yNC[], QBCP_L *qbp);

//=============================================================================================
//
// COC: Velocities <--> Momenta
//
//=============================================================================================

//-----------------------------------------------------------------------------
// SEM
//-----------------------------------------------------------------------------
/**
 *  \brief Change the SEM velocities into SEM momenta
 **/
void SEMvtoSEMm(double t, const double ySEv[], double ySEm[], void *params_void);

/**
 *  \brief Change the SEM momenta into SEM velocities
 **/
void SEMmtoSEMv(double t, const double ySEm[], double ySEv[], void *params_void);

//-----------------------------------------------------------------------------
// EM
//-----------------------------------------------------------------------------
/**
 *  \brief Change the EM velocities into EM momenta
 **/
void EMvtoEMm(double t, const double yEMv[], double yEMm[], void *params_void);

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void EMmtoEMv(double t, const double yEMm[], double yEMv[], void *params_void);


//-----------------------------------------------------------------------------
// SYS & NC
//-----------------------------------------------------------------------------
/**
 *  \brief Change the velocities into momenta. Works for both SYS (SEM, EM) and NC (NCEM, NCSEM) coordinates.
 *         For this reason, the routine is not called SYSvSYSm, but rather VELtoMOM.
 **/
void VELtoMOM(double t, const double ySYSv[], double ySYSm[], void *params_void);

/**
 *  \brief Change the momenta into velocities. Works for both SYS (SEM, EM) and NC (NCEM, NCSEM) coordinates.
 *         For this reason, the routine is not called SYSmSYSv, but rather MOMtoVEL.
 **/
void MOMtoVEL(double t, const double ySYSm[], double ySYSv[], void *params_void);

//========================================================================================
//
// COC: Velocities <--> Momenta + NC <-> SYS
//
//========================================================================================
/**
 *  \brief From VSEM to NCSEMm
 **/
void SEMvtoNCSEMm(double t, const double ySEMv[], double yNCSEm[], QBCP_L *qbcp_l);
/**
 *  \brief From VEM to NCEMm
 **/
void EMvtoNCEMm(double t, const double yEMv[], double yNCEMm[], QBCP_L *qbcp_l);

//=============================================================================================
//
// Change of unit SYSTEM
//
//=============================================================================================
/**
 *   \brief From SEM units to EM units for a Position/Velocity and time vector in IN coordinates and time
 **/
void ussem2usem(double *tc, double yINv[], QBCP_L *qbcp_l);

/**
 *   \brief From EM units to SEM units for a Position/Velocity and time vector in IN coordinates and time
 **/
void usem2ussem(double *tc, double yINv[], QBCP_L *qbcp_l);

//=============================================================================================
//
// COC: SEM <--> IN <--> EM
//
//=============================================================================================

//-----------------------------------------------------------------------------
// COC: IN <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void EMvtoIN(double t, const double yEM[], double yIN[], QBCP_L *qbcp_l);

/**
 * \brief From IN to EM (in EM units)
 **/
void INtoEM(double t, const double yIN[], double yEM[],
            QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: IN <--> SEM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to IN (in SEM units)
 **/
void SEMtoIN(double t, const double ySE[], double yIN[],
             QBCP_L *qbcp_l);

/**
 * \brief From IN to SEM (in SEM units)
 **/
void INtoSEM(double t, const double yIN[], double ySE[],
             QBCP_L *qbcp_l);


//========================================================================================
//
// COC: INSEM <--> EM, INEM <--> SEM
//
//========================================================================================
/**
 * \brief From INEM to SEM
 **/
void INEMvtoSEMm(double t, const double yINEMv[], double ySEMm[], QBCP_L *qbcp_l);
/**
 * \brief From INEM to VSEM
 **/
void INEMvtoSEMv(double t, const double yINEMv[], double ySEMv[], QBCP_L *qbcp_l);
/**
 * \brief From INSEM to VEM
 **/
void INSEMvtoEMv(double t, const double yINSEMv[], double yEMv[], QBCP_L *qbcp_l);
/**
 * \brief From INSEM to EM
 **/
void INSEMvtoEMm(double t, const double yINSEMv[], double yEMm[], QBCP_L *qbcp_l);

//=============================================================================================
//
// COC: SEM <--> EM
//
//=============================================================================================
//-----------------------------------------------------------------------------
// COC: SEMv <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMvtoEMm(double t, const double ySEMv[], double yEMm[], QBCP_L *qbcp_l);

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMv(double t, const double yEMm[], double ySEMv[], QBCP_L *qbcp_l);


//-----------------------------------------------------------------------------
// COC: SEMm <--> EMv
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMv(double t, const double ySEMm[], double yEMv[], QBCP_L *qbcp_l);

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMvtoSEMm(double t, const double yEMv[], double ySEMm[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: SEMv <--> EMv
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMvtoEMv(double t, const double ySEv[], double yEMv[], QBCP_L *qbcp_l);

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMvtoSEMv(double t, const double yEMv[], double ySEMv[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: SEM <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMm(double t, const double ySEm[], double yEMm[], QBCP_L *qbcp_l);

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMm(double t, const double yEMm[], double ySEMm[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: SEM <--> NCEM
//-----------------------------------------------------------------------------
/**
 * \brief From VNCEM to SEMm
 **/
void NCEMvtoSEMm(double tEM, const double yNCEMv[], double ySEMm[], QBCP_L *qbcp_l);

/**
 * \brief From VNCEM to SEMv
 **/
void NCEMvtoSEMv(double tEM, const double yNCEMv[], double ySEMv[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//-----------------------------------------------------------------------------
/**
 * \brief From VNCEM to  NCSEMm
 **/
void NCEMvtoNCSEMm(double tEM, const double yNCEMv[], double yNCSEMm[], QBCP_L *qbcp_l);

/**
 * \brief From NC EM to  NC SEM (both in position/momenta form)
 **/
void NCEMmtoNCSEMm(double tEM, const double yNCEMm[], double yNCSEM[], QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to  NC EM (both in position/momenta form)
 **/
void NCSEMmtoNCEMm(double t, const double yNCSEMm[], double yNCEM[], QBCP_L *qbcp_l);

/**
 * \brief From VNCEM to VNCSEM
 **/
void NCEMvtoNCSEMv(double tEM, const double yNCEMv[], double yNCSEMv[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: SEM <--> NCEM
//-----------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM (both in position/momenta form)
 **/
void NCEMmtoSEMm(double t, const double yNCEMm[], double ySEMm[], QBCP_L *qbcp_l);

/**
 * \brief From SEM to NCEM (both in position/momenta form)
 **/
void SEMmtoNCEMm(double t, const double ySEMm[], double yNCEMm[], QBCP_L *qbcp_l);


//-----------------------------------------------------------------------------
// COC: NCSEM <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From NCSEM to EM (both in position/momenta form)
 **/
void NCSEMmtoEMm(double t, const double yNCSEMm[], double yEMm[], QBCP_L *qbcp_l);

/**
 * \brief From EM to  NCSEM (both in position/momenta form)
 **/
void EMmtoNCSEMm(double tEM, const double yEMm[], double yNCSEM[], QBCP_L *qbcp_l);

//=============================================================================================
//
//          Change of coordinates on vectors
//
//=============================================================================================
//-----------------------------------------------------------------------------
// COC: SEM <--> NCEM
//-----------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to SEMv for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMv_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//-----------------------------------------------------------------------------
/**
 * \brief From NC EM to NC SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, double *tSEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to NC SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to NC EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to NC EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, double *tEM,  int N, QBCP_L *qbcp_l);


//-----------------------------------------------------------------------------
// COC: NCSEM <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From NC SEM to EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to EMv for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMv_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: NC <--> SYS
//-----------------------------------------------------------------------------

/**
 *  \brief From NC to SYS, in SYS units (either EM or SEM).
 **/
void NCtoSYS_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l);

/**
 *  \brief From NC to SYSv, in SYS units (either EM or SEM).
 **/
void NCtoSYSv_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l);

/**
 *  \brief From SYS to NC, in SYS units (either EM or SEM).
 **/
void SYStoNC_vec(double **ySYS, double *tSYS, double **yNC, int N, QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: NC <--> VNC
//-----------------------------------------------------------------------------
/**
 *  \brief From NC to VNC, in NC units (either EM or SEM).
 **/
void NCtoVNC_vec(double **yNC, double *tNC, double **yVNC, int N, QBCP_L *qbcp_l);

/**
 *  \brief From VNC to NC, in NC units (either EM or SEM).
 **/
void VNCtoNC_vec(double **yVNC, double *tNC, double **yNC, int N, QBCP_L *qbcp_l);

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp();

/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp();

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp();

//------------------------------------------------------------------------------------
//   Print
//------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout. More precise than vector_printf.
 **/
void vector_printf_prec(const double *y, int n);
/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(const double *y, int n);
/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(const cdouble *y, int n);
/**
 *  \brief Prints an array of complex double using cout. More precise than vector_complex_printf.
 **/
void vector_complex_printf_prec(const cdouble *y, int n);
/**
 *  \brief Prints a matrix of complex double using cout.
 **/
void matrix_complex_printf(const cdouble **y, int n, int m);

//------------------------------------------------------------------------------------
//   Norm
//------------------------------------------------------------------------------------
/**
 *  Euclidian norm computed on the first k components of a complex vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(const cdouble z0[], int k);

/**
 *  Euclidian norm computed on the first k components of a double vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(const double z0[], int k);

/**
 *  Euclidian norm of the difference of two double vectors, computed on the first k: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  (z_1[p] - z_2[p]) ^2 \right)^{-1/2} \f$
 **/
double DENorm(const double z1[], const double z2[], int k);

/**
 *  Euclidian norm of the difference of two cdouble vectors, computed on the first k: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  (z_1[p] - z_2[p]) ^2 \right)^{-1/2} \f$
 **/
double DENorm(const cdouble z1[], const cdouble z2[], int k);

/**
 *  Euclidian norm of the difference of two cdouble matrices.
 **/
double DENorm(const cdouble **z1, const cdouble **z2, int k, int l);

/**
 *  Normalize the vector z0, containing k components
 **/
void ENorm(const double z0[], double z0n[], int k);

//------------------------------------------------------------------------------------
//   Default parameterization
//------------------------------------------------------------------------------------
/**
 * \brief Get the default coordinates system for variational equations, from the coord_type
 **/
int default_coordinate_system(int coord_type);

/**
 * \brief Get the default framework, from the coord_type
 **/
int default_framework(int coord_type);

#endif // EMINSEM_H_INCLUDED
