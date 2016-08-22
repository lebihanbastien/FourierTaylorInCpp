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

//-----------------------------------------------------------------------------
// COC: Velocities <--> Momenta
//-----------------------------------------------------------------------------
/**
 *  \brief Change the SEM velocities into SEM momenta
 **/
void SEMvtoSEMm(double t, const double ySEv[], double ySEm[], void *params_void);

/**
 *  \brief Change the SEM momenta into SEM velocities
 **/
void SEMmtoSEMv(double t, const double ySEm[], double ySEv[], void *params_void);

/**
 *  \brief Change the EM velocities into EM momenta
 **/
void EMvtoEMm(double t, const double yEMv[], double yEMm[], void *params_void);

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void EMmtoEMv(double t, const double yEMm[], double yEMv[], void *params_void);

//-----------------------------------------------------------------------------
// Change of unit system
//-----------------------------------------------------------------------------
/**
 *   \brief From EM unit system to SEM unit system for a Position/Velocity and time vector in IN coordinates and time
 **/
void usem2ussem(double *tc, double yINv[], QBCP_L *qbcp_l);

/**
 *   \brief From SEM unit system to EM unit system for a Position/Velocity and time vector in IN coordinates and time
 **/
void ussem2usem(double *tc, double yINv[], QBCP_L *qbcp_l);

//-----------------------------------------------------------------------------
// COC: IN <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void EMtoIN(double t, const double yEM[], double yIN[], QBCP_L *qbcp_l);

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

//-----------------------------------------------------------------------------
// COC: SEM <--> IN <--> EM
//-----------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMm(double t, const double ySEm[], double yEMm[],
              QBCP_L *qbcp_l);

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMm(double t, const double yEMm[], double ySEMm[],
              QBCP_L *qbcp_l);

/**
 * \brief From NC EM to SEM (both in position/momenta form)
 **/
void NCEMmtoSEMm(double t, const double yNCEMm[], double ySEMm[], QBCP_L *qbcp_l);

/**
 * \brief From NC EM to  NC SEM (both in position/momenta form)
 **/
void NCEMmtoNCSEMm(double t, const double yNCEMm[], double yNCSEM[], QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to  NC EM (both in position/momenta form)
 **/
void NCSEMmtoNCEMm(double t, const double yNCSEMm[], double yNCEM[], QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to EM (both in position/momenta form)
 **/
void NCSEMmtoEMm(double t, const double yNCSEMm[], double yEMm[], QBCP_L *qbcp_l);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Change of coordinates on vectors
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to NC SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC EM to NC SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, double *tSEM, int N, QBCP_L *qbcp_l);


/**
 * \brief From NC SEM to EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, int N, QBCP_L *qbcp_l);

/**
 * \brief From NC SEM to NC EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, int N, QBCP_L *qbcp_l);

/**
 *  \brief From NC to SYS, in SYS units (either EM or SEM).
 **/
void NCtoSYS_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l);

/**
 *  \brief From NC to SEM, in SEM units.
 **/
void NCtoSEM_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l);

/**
 *  \brief From SYS to NC, in SYS units (either EM or SEM).
 **/
void SYStoNC_vec(double **ySYS, double *tSYS, double **yNC, int N, QBCP_L *qbcp_l);


#endif // EMINSEM_H_INCLUDED
