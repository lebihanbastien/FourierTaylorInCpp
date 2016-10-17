#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <complex.h>

//------------------------------------------------------------------------------------
// Environment
//------------------------------------------------------------------------------------
// Numerotation of the Solar System planets and objects, consistent with JPL's HORIZON numerotation (NAIF ID)
#define SUN 10
#define MERCURY 199
#define VENUS 299
#define EARTH 399
#define MARS 499
#define JUPITER 599
#define SATURN 699
#define URANUS 799
#define NEPTUNE 899
#define PLUTO 999
#define MOON 301
#define EARTH_AND_MOON 3
#define EARTH_MOON_BARYCENTER 3
#define SSB 0

//Barycenters
#define MERCURY_BARYCENTER 1
#define VENUS_BARYCENTER   2
#define EARTH_BARYCENTER   3
#define MARS_BARYCENTER    4
#define JUPITER_BARYCENTER 5
#define SATURN_BARYCENTER  6
#define URANUS_BARYCENTER  7
#define NEPTUNE_BARYCENTER 8
#define PLUTO_BARYCENTER   9

// Precision on the position of the librations points L1/L2/L3 in define_env.h
#define LIBRATION_POINT_PRECISION 1e-16


//------------------------------------------------------------------------------------
//   Global constants
//------------------------------------------------------------------------------------
extern int OFTS_ORDER;
extern int OFS_ORDER;
extern int OTS_ORDER;
extern int MODEL_TYPE;
//extern int REDUCED_NV;

//------------------------------------------------------------------------------------
//   typedef
//------------------------------------------------------------------------------------
typedef complex double cdouble;


//------------------------------------------------------------------------------------
//   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
//------------------------------------------------------------------------------------

//-------------------------------------------------------------------------
//Model
//-------------------------------------------------------------------------
#define M_RTBP  0 // RTBP model indix
#define M_QBCP  1 // QBCP model indix
#define M_BCP   2 // BCP model indix
#define M_ERTBP 3 // ERTBP model indix

//-------------------------------------------------------------------------
//Frameworks
//-------------------------------------------------------------------------
#define F_EM     0
#define F_SEM    1
#define F_VNCEM  2
#define F_VNCSEM 3
#define F_NCEM   4
#define F_NCSEM  5
#define F_ECLI   6
#define F_VSEM   7
#define F_VEM    8
#define F_INSEM  9
#define F_INEM  10
#define F_J2000 11

//-------------------------------------------------------------------------
// Available coordinates system
//-------------------------------------------------------------------------
#define NCSEM  0
#define NCEM   1
#define VNCSEM 2
#define VNCEM  3
#define PSEM   4
#define PEM    5
#define VSEM   6
#define VEM    7
#define INEM   8
#define INSEM  9
#define VECLI  10
#define J2000  11

//Number of variables
#define NV 6

// Number of variables in the OFS object (a priori always 1)
#define OFS_NV 1

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp();

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp();

/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp();

#endif // PARAMETERS_H_INCLUDED
