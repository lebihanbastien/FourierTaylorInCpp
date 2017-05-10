#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED

#include <complex.h>

//----------------------------------------------------------------------------------------
// Environment
//----------------------------------------------------------------------------------------
// Numerotation of the Solar System planets and objects, consistent with
// JPL's HORIZON numerotation (NAIF ID)
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


//----------------------------------------------------------------------------------------
//   Global constants
//----------------------------------------------------------------------------------------
extern int OFTS_ORDER;
extern int OFS_ORDER;
extern int OTS_ORDER;

//----------------------------------------------------------------------------------------
//   typedef
//----------------------------------------------------------------------------------------
typedef complex double cdouble;


//----------------------------------------------------------------------------------------
//   ORDER AND NUMBER OF VARIABLES FOR OFS AND OFTS OBJECTS
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
//Model
//----------------------------------------------------------------------------------------
#define M_RTBP  0  // RTBP model index
#define M_QBCP  1  // QBCP model index
#define M_BCP   2  // BCP model index
#define M_ERTBP 3  // ERTBP model index

//----------------------------------------------------------------------------------------
//Frameworks
//----------------------------------------------------------------------------------------
#define F_EM      0
#define F_SEM     1

//----------------------------------------------------------------------------------------
//For mex communication (MATLAB)
//----------------------------------------------------------------------------------------
#define NC   0
#define SYS  1
#define VSYS 2
#define VNC  3

#define ORDERMAX_CRTBP 40
#define ORDERMAX_QBCP  20

//----------------------------------------------------------------------------------------
//Integration systems
//----------------------------------------------------------------------------------------
#define I_PEM     0
#define I_PSEM    1

#define I_VNCEM   2
#define I_VNCSEM  3
#define I_NCEM    4
#define I_NCSEM   5

#define I_VSEM    6
#define I_VEM     7
#define I_INSEM   8
#define I_INEM    9

#define I_ECISEM  10

#define I_ECLI    11
#define I_J2000   12
#define I_NJ2000  13
#define I_VSYNEM  14
#define I_VSYNSEM 15

//----------------------------------------------------------------------------------------
// Available coordinates systems
//----------------------------------------------------------------------------------------
// Synodical QBCP
#define NCSEM   0
#define NCEM    1
#define VNCSEM  2
#define VNCEM   3
#define PSEM    4
#define PEM     5
#define VSEM    6
#define VEM     7
// Inertial QBCP
#define INEM    8
#define INSEM   9
#define ECISEM  10
// Ephemerides
#define VECLI   11
#define J2000   12
#define NJ2000  13
#define VSYNEM  14
#define VSYNSEM 15

//----------------------------------------------------------------------------------------
//Number of variables
//----------------------------------------------------------------------------------------
#define NV 6

//----------------------------------------------------------------------------------------
// Number of variables in the OFS object (a priori always 1)
//----------------------------------------------------------------------------------------
#define OFS_NV 1

//----------------------------------------------------------------------------------------
//   Precisions
//----------------------------------------------------------------------------------------
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
