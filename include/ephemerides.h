#ifndef EPHEMERIDES_H_INCLUDED
#define EPHEMERIDES_H_INCLUDED

#include <string>
#include <iostream>

#include "ode.h"
#include "env.h"
#include "eminsem.h"

extern "C" {
    #include "gnuplot_i.h"
    #include "nrutil.h"
    #include "SpiceUsr.h"
}

using namespace std;

//Default
#define DEFOBS      "SSB"
#define DEFOBSINT   0
#define DEFFRAME "ECLIPJ2000"
//#define DEFFRAME "J2000"

//========================================================================================
//                         Test of the vector field
//========================================================================================
/**
 *  \brief Test of the JPL vector field using the LUTETIA asteroid
 **/
void test_asteroid();

/**
 *  \brief Generic test to display some feature of the COC
 **/
void test_coc(int coord_eph);

//========================================================================================
//                         Benchmark of the numerical constants
//========================================================================================
void comp_num_const();


//========================================================================================
//                         Subroutines to choose between SEM and EM representation
//========================================================================================
/**
 * \brief Get the ephemerides coordinate system associated to the current coordinate system.
 *        For now, only the only possible outputs are VSEM/VEM.
 *        In the future, VNCSEM and VNCEM should be made available.
 **/
int eph_coord(int coord_type);

/**
 * \brief Get the ephemerides fwrk associated to the current coordinate system.
 *        For now, only the only possible outputs are I_VSEM/I_VEM.
 *        In the future, I_VNCSEM and I_VNCEM should be made available.
 **/
int eph_fwrk(int coord_type);

/**
 *    \brief Mean motion, in rad/s, assciated to the ephemerides coordinates coord_eph
 **/
double mean_motion(int coord_eph);


//========================================================================================
//                         Looking for the best fit wrt to a given
//                              QBCP configuration
//========================================================================================
void qbcp2jpl(double t, double *et, int coord_type);
void qbcp2jpl_disp(double tSYS, double *et, int coord_type);
void qbcp2jpl_inertial(double tSYS, double *et, int coord_type);

//========================================================================================
//                         Name of the primaries
//========================================================================================
/**
 *      \brief Return the name of the first primary associated to the coord_type.
 *             Not included (since it makes no sense):
 *             case INEM   case INSEM  case VECLI
 **/
int m1name(int coord_type);

/**
 *      \brief Return the name of the second primary associated to the coord_type.
 *             Not included (since it makes no sense):
 *             case INEM   case INSEM  case VECLI
 **/
int m2name(int coord_type);


//========================================================================================
//                          Change of coordinates:
//              Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic
//                            Only in position
//========================================================================================
/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + kC r        and    r = 1/k inv(C) (R - B)
 *
 *        The time is given as a string, recognized by SPICE
 **/
void init_ecl2synpos(string epoch, double B[3], double C[3][3], double *k, int coord_type);

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + kC r        and    r = 1/k inv(C) (R - B)
 *
 *        The time is given as a double, in seconds (SPICE ephemeris time)s
 **/
void init_ecl2synpos(double et, double B[3], double C[3][3], double *k, int coord_type);

/**
 * \brief Change of coordinates: Ecliptic coordinates  <- Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *      R = B + kC r
 **/
void syn2eclpos(double vin[3], double vout[3], double B[3], double C[3][3], double k);

/**
 * \brief Change of coordinates: Ecliptic coordinates  -> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *                  r = 1/k inv(C) (R - B)
 **/
void ecl2synpos(double vin[3], double vout[3], double B[3], double C[3][3], double k);

/**
 * \brief Change of coordinates: Ecliptic coordinates  -> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *                  r = inv(C) (R - B)
 **/
void ecl2syndpos(double vin[3], double vout[3], double B[3], double C[3][3]);

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) dimensionalized synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + C r        and    r = inv(C) (R - B)
 *
 *        The time is given as a double, in seconds (SPICE ephemeris time)s
 **/
void init_ecl2syndpos(double et, double B[3], double C[3][3], int coord_eph);

//========================================================================================
//                          Change of coordinates:
//              Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic
//                            Position + Velocity
//========================================================================================
/**
 *  \brief Init the matrices C, kCprim and the scalar k necessary to compute
 *         the cocs: ecliptic -> synodic and synodic -> ecliptic.
 **/
void init_coc(const double r21[3],  const double v21[3],
              const double a21[3],  const double j21[3],
              double C[3][3], double *k, double kCprim[3][3]);

//-----------------------------
// Subroutines
//-----------------------------
/**
 * \brief Init k, k', and k" coefficients, with k = ||r21||, r21
 *        being the relative position of the primaries.
 *        See equations (3.4) of Dei Tos (2014)
 **/
void init_k(const double r21[3], const double v21[3], const double a21[3],
            double *k, double *kp, double *kpp);

/**
 *  \brief Init the cross products r21 x v21, r21 x a21, etc.
 **/
void init_cp(const double r21[3], const double v21[3],
             const double a21[3], const double j21[3],
             double r21_v21[3], double r21_a21[3],
             double r21_j21[3], double v21_a21[3]);


/**
 * \brief Init h, h', and h" coefficients, with h = ||r21 x v21||, r21/v21
 *        being the relative position/velocity of the primaries.
 *        See equations (3.5) of Dei Tos (2014)
 **/
void init_h(const double r21_v21[3], const double r21_a21[3],
            const double r21_j21[3], const double v21_a21[3],
            double *h, double *hp, double *hpp, double etemp[3]);

/**
 *  \brief Init the orthonormal basis (e1, e2, e3)
 **/
void init_ei(const double r21[3], const double v21[3],
             double e1[3], double e2[3], double e3[3]);

/**
 *  \brief Init the derivative of the orthonormal basis (e1', e2', e3')
 **/
void init_dei(const double r21[3], const double v21[3], const double a21[3],
              const double k,      const double kp,
              const double h,      const double hp,
              const double  e1[3], const double e3[3],
              double de1[3], double de2[3], double de3[3], double etemp[3]);

/**
 *  \brief Init the double derivative of the orthonormal basis (e1", e2", e3")
 **/
void init_ddei(const double r21[3], const double v21[3], const double a21[3],
               const double r21_v21[3], const double r21_a21[3],
               const double r21_j21[3], const double v21_a21[3],
               const double k, const double kp, const double kpp,
               const double h, const double hp, const double hpp,
               const double  e1[3],   const double  e3[3],
               const double  de1[3],  const double  de3[3],
               double dde1[3], double dde2[3], double dde3[3],
               double  etemp[3]);
/**
 *  \brief Init the derivative of the orthonormal basis times k:
 *           (kC)' = ((ke1)', (ke2)', (ke3)')
 **/
void init_kdei(const double  e1[3], const double  e2[3], const double  e3[3],
               const double de1[3], const double de2[3], const double de3[3],
               const double k,      const double kp,
               double dke1[3], double dke2[3], double dke3[3]);


//-----------------------------
// Init
//-----------------------------
/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate(string epoch, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph);

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate(double et, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph);

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate_2(double et, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph);

//-----------------------------
// Ecliptic coordinates -> Sun-(Earth+Moon) synodic
//-----------------------------
/**
 * \brief Change of coordinates Ecliptic coordinates  -> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  r  = 1/k inv(C) (R - B)
 *        - For the velocities: r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 *
 *         Note: n_sem is the Sun-Bem mean motion, in rad/s.
 **/
void ecl2synstate(double vin[6], double vout[6], double et, int coord_eph);

/**
 * \brief Change of coordinates Ecliptic coordinates  -> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  r  = 1/k inv(C) (R - B)
 *        - For the velocities: r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 *
 *         Note: n_sem is the Sun-Bem mean motion, in rad/s.
 **/
void ecl2synstate(double vin[6], double vout[6], double B[3], double C[3][3], double k, double Bprim[3], double kCprim[3][3], double n_sys);

//-----------------------------
// Ecliptic coordinates  <- Sun-(Earth+Moon) synodic
//-----------------------------
/**
 * \brief Change of coordinates Ecliptic coordinates  <- Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot
 **/
void syn2eclstate(double vin[6], double vout[6], double et, int coord_eph);

/**
 * \brief Change of coordinates Ecliptic coordinates  <- Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot
 **/
void syn2eclstate(double vin[6], double vout[6], double B[3], double C[3][3], double k, double Bprim[3], double kCprim[3][3], double n_sys);

//-----------------------------
// On vectors
//-----------------------------
/**
 *  From Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void ecl2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph);

/**
 *  From Normalized-Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void eci2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph);//, SS &ss);

/**
 *  To Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2eclstate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph);

/**
 *  To Normalized-Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2ecistate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph);//, SS &ss);

/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void ecl2synstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph);

/**
 * \brief From synodic to ecliptic coordinates, vector format
 **/
void syn2eclstate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, double et0, double tsys0, int coord_eph);

/**
 * \brief From synodic to normalized-ecliptic coordinates, vector format.
 *        Note that coord_eph = VEM/VSEM can be chosen independantly from the ss structure
 *        that normalizes the state
 **/
void syn2ecistate_vec(double **yin, double *tin, double **yecl, double *tecl, int N, double et0, double tsys0, int coord_eph);//, SS &ss);

/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void eci2synstate_vec(double **yecl, double *tecl, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph);//, SS & ss);
void eci2syndpos_vec(double **yecl, double *tecl, double **yout, double *tout, int N, int coord_eph);//, SS & ss);

/**
 *  \brief From full ecliptic coordinates (directly from SPICE) to normalized coordinates.
 **/
void ecl2eci(double YECL[6],  double YEARTH[6], double yeci[6]);//, SS &ss);

/**
 *  \brief To full ecliptic coordinates (directly from SPICE) from normalized coordinates.
 **/
void eci2ecl(double yeci[6],  double YEARTH[6], double YECL[6]);//, SS &ss);

//---------------------------------
// Normalized Earth-centered inertial (NECI)  <-> any other synodical type (VEM, VSEM...)
//---------------------------------
/**
 *  \brief From full ecliptic coordinates (directly from SPICE) to Earth-centered inertial coordinates.
 **/
void ecl2neci(double YECL[6], double YEARTH[6], double yneci[6], SS &ss);
/**
 *  \brief To full ecliptic coordinates (directly from SPICE) from to earth-centered inertial coordinates.
 **/
void neci2ecl(double yneci[6], double YEARTH[6], double YECL[6], SS &ss);
/**
 * \brief From synodic to earth-centered normalized coordinates, vector format.
 *        Note that coord_eph = VEM/VSEM can be chosen independantly from the ss structure
 *        that normalizes the state.
 *        Note: for now the normalization is NOT taken into account because the stepper is going wrong when the state is normalized.
 **/
void syn2necistate_vec(double **yin, double *tin, double **yneci, double *tneci, int N, double et0, double tsys0, int coord_eph, SS &ss);
/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void neci2synstate_vec(double **yneci, double *tneci, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph, SS & ss);

//---------------------------------
// Normalized Earth-centered inertial (NECI) coordinates <-> any other native type (NCEM, NCSEM...)
//---------------------------------
/**
 *  From Normalized-Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void neci2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph, SS &ss);

/**
 *  To Normalized-Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2necistate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph, SS &ss);

//========================================================================================
//                          Compute the acceleration of the
//                                  m1-m2  line,
//                           (in ecliptic coordinates)
//========================================================================================
/**
 *  \brief Return a component (along the dimension dim) of the m1-m2 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xBSecl(double et, int dim, int coord_type);

/**
 *   \brief Returns the derivative of a function func at a point x by Ridders’ method of polynomial
 *          extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but
 *          rather should be an increment in x over which func changes substantially. An estimate of the
 *          error in the derivative is returned as err. (from Numerical Recipe in C).
 **/
double dfridr(double (*func)(double, int, int), double x, int dim, double h, double *err, int coord_type);

/**
 *  \brief Return a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xm1ecl(double et, int dim, int coord_type);

/**
 *  \brief Return a component (along the dimension dim) of the m2 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xm2ecl(double et, int dim, int coord_type);


/**
 *  \brief Return the derivative of a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double dxm1ecl(double et, int dim, int coord_type);

/**
 *  \brief Return the derivative of a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double dxm2ecl(double et, int dim, int coord_type);

//========================================================================================
//                          Display information on SPICE Kernels.
//                          Routine are based on SPICE tutorials.
//========================================================================================
/**
 *  \brief Display information on a given kernel (currently de432s). The SRC code may be changed to include a user input for the name of the kernel
 **/
void displayKernelFeatures();

/**
 *  \brief Display information on a given body accross multiple kernels loaded via a meta kernel (eg: spice/kernels/metakernel.furnsh).
 **/
void displayKernelFeaturesOneBody();

#endif // EPHEMERIDES_H_INCLUDED
