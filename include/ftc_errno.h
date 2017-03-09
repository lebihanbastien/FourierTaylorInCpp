#ifndef FTC_ERRNO_H_INCLUDED
#define FTC_ERRNO_H_INCLUDED

/**
 *  \brief Error handling. It defines macros for reporting and retrieving error conditions.
 *         inspired by errno.h of the C library, and gsl_errno.h of the GSL library.
 *         As much as possible, the error codes defined in this file are consistent with
 *         its GSL counterpart gsl_errno.h.
 **/


enum {

  //Generic CODES
  FTC_SUCCESS  =  0,            //equal to GSL_SUCCESS =  0 for consistency
  FTC_FAILURE  = -1,            //equal to GSL_FAILURE = -1 for consistency
  FTC_EDOM     =  1,            //input domain error, equal to GSL_EDOM =  1 for consistency
  FTC_ENOENT   = -33,           //No such file or directory

  //For refinement procedures
  REF_EOUTOFDPC = -2,     //out of the domain of practical convergence (DPC) of the semi-analytical tools
  REF_EMAXITREACH = -3,   //the maximum iteration number is reached.


  //For orbit procedures
  ORBIT_EINT     = -4,  //integration failed
  ORBIT_EPROJ    = -5,  //projection procedure failed
  ORBIT_ECOLL    = -6   //collision with a primary
};


const char* ref_strerror (const int ref_errno);

#endif // FTC_ERRNO_H_INCLUDED
