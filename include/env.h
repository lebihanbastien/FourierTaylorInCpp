#ifndef ENV_H_INCLUDED
#define ENV_H_INCLUDED


/**
 * \file  env.h
 * \brief Define the working environment (the primaries). All values taken from JPL and Goddard Space Flight Center websites.
 * \author BLB.
 * \date 2016
 * \version 1.0
 */

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <complex.h>
#include <fstream>
#include <sstream>

#include "Ofsc.h"
#include "matrix.h"

//Type of unit systems
#define USYS_EM  0 //Earth-Moon      system
#define USYS_SEM 1 //Sun-Earth+Moon  system

//Type of Param Style
#define PMS_GRAPH 0
#define PMS_NORMFORM 1
#define PMS_MIXED 2

//Type of Manifolds
#define MAN_CENTER    0
#define MAN_CENTER_S  1
#define MAN_CENTER_U  2
#define MAN_CENTER_US 3

//Projection distance
#define PROJ_EPSILON     1e-6
#define PROJ_EPSILON_SEM 1e-6


using namespace std;


//----------------------------------------------------------------------------------------
//            Structures
//----------------------------------------------------------------------------------------
/**
 * \struct LibrationPoint
 * \brief Structure to describe a Libration point (position, energy...).
 **/
typedef struct LibrationPoint LibrationPoint;
struct LibrationPoint
{
    ///number
    int number;
    ///position [adim] in the corresponding CR3BP-frame
    double position[3];
    ///distance to closest primary (only for l1,l2,l3);
    double gamma_i;
    ///energy
    double Ei;
    ///jacobi constant
    double Ci;
};

/**
 * \struct Body
 * \brief Characteristic constants of a celestial body (Mass, orbital period...)
 **/
typedef struct Body Body;
struct Body
{
    ///Mass [kg]
    double M;
    ///Gravitationnal parameter [km^3/s^2]
    double GM;
    ///Equatorial Radius [km]
    double Req;
    ///Mean radius [km]
    double Rm;
    ///Sidereal orbit period [s]
    double T;
    ///Semi-major axis [km]
    double a;
    ///Name
    char name[50];
};

/**
 * \struct CR3BP
 * \brief  Environment of a given Circular Restricted Three-Body Problem.
 **/
typedef struct CR3BP CR3BP;
struct CR3BP
{
    Body m1;   //First primary
    Body m2;   //Second primary

    double mu; // Gravitational constant [-]
    double L;  // Distance parameter [km]
    double T;  // Time parameter [s]
    double R1; // Radius of the first primary  [km]
    double R2; // Radius of the second primary [km]
    double rh; // Hill's radius [-]

    //Libration points (adim)
    LibrationPoint l1;
    LibrationPoint l2;
    LibrationPoint l3;
    LibrationPoint l4;
    LibrationPoint l5;

    //precision on L1/L2/L3 position
    double gprecision;

    //Name
    char name[50];

    //Distance to manifold approximation
    double d_man;
};

/**
 * \struct USYS
 * \brief  Unit system structure containing a given set of parameters in a given unit system.
 **/
typedef struct USYS USYS;
struct USYS
{
    //Note: a name should be added!
    int label;      //type of unit system
    double mu_EM;   //mass ratio (EM)
    double mu_SEM;  //mass ratio (SEM)
    double mu_SE;   //mass ratio (SE)
    double as;      //radius of the circular motion approximating the quasi-circular exterior (Earth+Moon-Sun) motion
    double ns;      //Mean angular motion of this circular motion
    double ai;      //radius of the circular motion approximating the quasi-circular interior (Earth-Moon) motion
    double ni;      //Mean angular motion of this circular motion
    double n;       //difference of mean angular motion
    double T;       //Period
    double ms;      //Sun mass
    double me;      //Earth mass
    double mm;      //Moon mass
    double lecc;    //lunar eccentricity
};

/**
 * \struct CSYS
 * \brief  Coordinates system structure.
 **/
typedef struct CSYS CSYS;
struct CSYS
{
    //Name and label
    string name;    //name of the unit system
    int label;      //type of unit cs

    //CR3BP
    CR3BP cr3bp; //Associated CR3BP

    //Simple coefficients
    double c1;      //c1 coefficient, defined wrt to a given libration point li
    double c2;      //c2 coefficient, defined wrt to a given libration point li
    double gamma;   //gamma coefficient, defined wrt to a given libration point li
    double mu;      //mass ratio
    USYS us;        //default unit system
    int li;         //associated libration point
    int manType;    //associated manifold type
    int pmType;     //associated parameterization of manifold
    int model;      //associated model
    int fwrk;       //associated fwrk

    //Arrays
    double *coeffs; //Default set of vector field coefficients
    double *Ps;    //Sun   position in EM coordinates
    double *Pe;    //Earth position in EM coordinates
    double *Pm;    //Moon  position in EM coordinates
    double *ps;    //Sun   position in NC coordinates
    double *pe;    //Earth position in NC coordinates
    double *pm;    //Moon  position in NC coordinates

    //QBTBP
    Ofsc zt;
    Ofsc Zt;
    Ofsc ztdot;
    Ofsc Ztdot;

    //Folders
    string F_COC;   //Change of coordinates Translation+Floquet+Complexification
    string F_PRINT; //For printing results
    string F_PLOT;  //For plotting
    string F_COEF;  //Coefficients alphas for integration in the Earth-Moon system.
    string F_QBTBP; //Solution of the QBTBP

    //Parameterization folders
    string F_PMS;   //Global Parameterization folder, chosen within the folders below:
    string F_GS;    //Graph style           |
    string F_NF;    //Normal form style     |  all for Center Manifold
    string F_MS;    //Mixed style           |

    string F_CS;    //Center-Stable (always mixed style)
    string F_CU;    //Center-Unstable (always mixed style)
    string F_CUS;   //Center-Hyperbolic (always mixed style)

};

/**
 * \struct QBCP
 * \brief  Environment of a given Quasi-Bicircular Four-Body Problem.
 **/
typedef struct QBCP QBCP;
struct QBCP
{
    Body m1;   //First  primary
    Body m2;   //Second primary
    Body m3;   //Third  primary

    CR3BP cr3bp1; //First associated CR3BP
    CR3BP cr3bp2; //Second associated CR3BP
};

/**
 *  \brief Solar system structure
 **/
typedef struct SS SS;
struct SS
{
    int *id;        //id of the bodies
    double *Gmi;    //G*mi for all the bodies
    double *mui;
    double a;
    double n;
    double mu1;
    double mu2;
    double et0;
    double t0;
    int id1;
    int id2;
    int pos1;
    int pos2;
    int coord_eph;
    int maxBodies;
};


/**
 * \struct QBCP_L
 * \brief  Environment of a given Quasi-Bicircular Four-Body Problem focused on one libration point. The libration point must be L1, L2 or L3.
 *
 * Note: the alpha routines are defined here even though they are common to all libration points
 * so that the structure QBCP can be used to compute the alphas without depending on them.
 **/
typedef struct QBCP_L QBCP_L;
struct QBCP_L
{
    int nf;             //Order of the Fourier expansions
    int eff_nf_EM;      //Effective order of the Fourier expansions, for some specific computations
    int eff_nf_SEM;     //Effective order of the Fourier expansions, for some specific computations
    int isNormalized;   //Are the equations of motion normalized?

    //Specific to one libration point
    int li_EM;            //number of the libration point considered
    int li_SEM;           //number of the libration point considered


    //6*6 matrix B used for Floquet transformation
    double *B;

    //Model
    int model;
    int fwrk;
    int li;
    int pms_EM;
    int pms_SEM;

    int numberOfCoefs; //number of Fourier coefficients

    //Unit systems
    USYS us_em;   //EM unit system
    USYS us_sem;  //SEM unit system
    USYS us;      //default unit system

    //Coordinate systems
    CSYS cs_em_l1;      //EM around L1
    CSYS cs_em_l2;      //EM around L2
    CSYS cs_em_l3;      //EM around L3

    CSYS cs_sem_l1;     //SEM around L1
    CSYS cs_sem_l2;     //SEM around L2
    CSYS cs_sem_l3;     //SEM around L3

    //Default coordinate systems
    CSYS cs_em;     //Default EM csys  (for cocs)
    CSYS cs_sem;    //Default SEM csys (for cocs)
    CSYS cs;        //Default csys     (for integration)

    //Is the Moon "on"? Allows to "wipe out" the Moon from the equations of motion
    //when using Sun-Earth or Sun-(Earth+Moon) coordinates systems.
    double epsilon;

    //For ephemerides conversion we use the mean motion of the Sun-Bem system n_sem
    double n_sem;
    //For ephemerides conversion we use the mean motion of the Earth-Moon system n_em
    double n_em;

    //Solar systems
    SS ss;      //default
    SS ss_sem;  //focus on SEM
    SS ss_em;   //focus on EM
};

/**
 * \struct QBCP_I
 * \brief  Hybrid environment for contination procedures between two models (model1 and model2).
 **/
typedef struct QBCP_I QBCP_I;
struct QBCP_I
{
    QBCP_L model1;
    QBCP_L model2;
    double epsilon;
};

//----------------------------------------------------------------------------------------
// Initialize the invariant manifolds.
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize the invariant manifolds, as a form of FT fourier series.
 *        The necessary data are retrieved from data folders  whose names are stored in
 *        the CS structure.
 *
 *        IMPORTANT: to be consistent with old code, both the TFC and NC form are initialized.
 *        However, the NC form is a priori NEVER used, since the numerical computation are
 *        usually more precise evaluating the TFC form, then going back to the NC form via the equation:
 *
 *              IM_NC = P(theta)*IM_TFC + V(theta).
 *
 *        If ones would want to use the NC form, just uncomment the proper lines inside the code.
 **/
void initOFTS(vector<Oftsc> &IM_NC, vector<Oftsc> &IM_TFC,
              int reduced_nv, int ofts_order, int ofs_order,
              CSYS &csys);

//----------------------------------------------------------------------------------------
//            Init routines
//----------------------------------------------------------------------------------------
/**
* \brief Initialize one celestial body
* \param body pointer on the current body
* \param name the name of the body in integer format (consistent with HORIZON numerotation)
**/
void init_body(Body *body, int name);

/**
* \fn void init_CR3BP(CR3BP *cr3bp, int n1, int n2)
* \brief Initialize the Circular Restricted 3-Body Problem
* \param cr3bp pointer on the CR3BP
* \param n1 name of the first primary
* \param n2 name of the second primary
**/
void init_CR3BP(CR3BP *cr3bp, int n1, int n2);

/**
* \brief Initialize a unit system in the form of a usys structure, such as Earth-Moon, Sun-Earth or Sun-(Earth+Moon) unit systems.
* \param usys pointer on the usys structure to init.
* \param label the type of unit system
* \param model type of physical model associated to this unit system
*
* NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES
**/
void init_USYS(USYS *usys, int label, int model);

/**
 *  \brief From EM to NC coordinates for the primaries. Used in qbtbp_ofs_fft_*
 */
void EMtoNC_prim(double Zc[3], double zc[3], double c1, double gamma);

/**
 * \brief Initializes a coordinate systems (CSYS structure), with associated vector field coefficients, data folder names, and unit system.
 * \param csys pointer on the CSYS structure to initialize.
 * \param qbcp_l pointer on the QBCP_L structure that contains csys.
 * \param qbcp pointer on the QBCP structure that contains parameters specific to each libration points (namely, gamma)
 * \param coordsys indix of the coordinate system to use (F_EM, F_SEM).
 * \param li number of the libration point to focus on (L1, L2).
 * \param coefNumber the number of vector field coefficients to initialize. It has been set in the QBCP_init function.
 * \param isNew boolean. if true, the qbtbp has not been computed via the qbtbp() routine, so the vector field coefficients cannot be initialized.
 *
 *   Note that the QBCP structure is used only for the initialization of the coordinate systems. More precisely, it contains some parameters
 *   specific to each libration point (gamma), via its CR3BP structures.
 **/
void init_CSYS(CSYS *csys, QBCP_L *qbcp_l, QBCP *qbcp, int coordsys, int li, int coefNumber, int isNew, int pmType, int manType);

/**
* \brief Initialize the Quasi-Bicircular Four-Body Problem in the form of a QBCP structure.
* \param qbcp pointer on the QBCP structure to init.
* \param n1 name of the first primary (e.g Sun)
* \param n2 name of the second primary (e.g. Earth)
* \param n3 name of the third primary  (e.g. Moon)
*
* NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES
**/
void init_QBCP(QBCP *qbcp, int n1, int n2, int n3);

/**
* \brief Initialize a QBCP_L structure, i.e. a QBCP focused on one libration point. The libration point must be L1 or L2 of qbcp.cr3bp1.
* \param qbcp_l a pointer on the QBCP_L structure to init.
* \param qbcp a pointer on the QBCP parent structure.
* \param isNormalized: are the equations of motion normalized?
* \param li number of the libration point considered.
* \param isNew an integer: equal to 1 if no solution has been previously computed with qbtbp(int), 0 otherwise
* \param model: QBCP, BCP, CRTBP...
*
* NEEDS TO BE MODIFIED TO GET RID OF THE HARD CODED VALUES
**/
void init_QBCP_L(QBCP_L *qbcp_l, QBCP *qbcp, int isNormalized, int li_EM, int li_SEM,
                 int isNew, int model, int coordsys, int pmType_EM, int pmType_SEM,
                 int manType_EM, int manType_SEM);


void init_QBCP_I(QBCP_I *model, QBCP_L *model1, QBCP_L *model2, int n1, int n2, int n3,
                 int isNormalized, int li_EM, int li_SEM, int isNew, int mod1, int mod2,
                 int coordsys, int pmType1, int pmType2);

/**
 * \brief Initializes a libration point.
 * \param libp a pointer towards the LibrationPoint structure to init.
 * \param cr3bp a CR3BP structure that contains useful coefficients.
 * \param number the indix of the libration point to init.
 **/
void init_libp(LibrationPoint *libp, CR3BP cr3bp, int number);

/**
* \brief Initialize a solar system, composed of 11 bodies (Sun, planets and the Moon)
**/
void init_SS(SS *solarsys, QBCP_L *qbcp_l, int coordsys);

//----------------------------------------------------------------------------------------
//QBCP
//----------------------------------------------------------------------------------------
extern QBCP SEM;               //Sun-Earth-Moon system
extern QBCP_L SEML;            //Sun-Earth-Moon system around Li
extern QBCP_L SEML_EM;         //Sun-Earth-Moon system around EMLi
extern QBCP_L SEML_SEM;        //Sun-Earth-Moon system around SEMLi

/**
 *   \brief Initialization of the environnement (Sun, Earth, Moon, Li...).
 *
 *    The global variables SEM and SEML are initialized in order to describe the Sun-(Earth+Moon) Quasi-Bicircular Four-Body Problem around the Lagrangian point LI,
 *    LI being given as an integer in the file parameters.h. The default coordinates system is the normalized-centered (NC) coordinates of the Earth-Moon system.
 *    Note that, in order for the initialization to be complete - Sun-(Earth+Moon) motion given as Fourier series within SEML -
 *    the routine qbtbp(true) must have been run *once.
 **/
void initialize_environment(int li_EM, int li_SEM, int isNormalized, int model, int coordsys,  int pmStyle_EM, int pmStyle_SEM, int manType_EM, int manType_SEM);

/**
 *  Main routine for the initialization of the COC
 **/
void initCOC(QBCP_L &qbcp);

/**
 *  Main routine for the update of the COC
 **/
void updateCOC(QBCP_L &qbcp);

//----------------------------------------------------------------------------------------
//COC
//----------------------------------------------------------------------------------------
/**
 *  \brief Lighter version of the init routine, with only P, PC, Q, QC, and V retrieved.
 **/
void initCOC(matrix<Ofsc> &P,
             matrix<Ofsc> &PC,
             matrix<Ofsc> &Q,
             matrix<Ofsc> &CQ,
             vector<Ofsc> &V,
             QBCP_L& qbcp_l);

/**
 *  \brief Lighter version of the init routine, with only PC, QC, and V retrieved.
 **/
void initCOC(matrix<Ofsc>& PC,
             matrix<Ofsc>& CQ,
             vector<Ofsc>& V,
             CSYS& csys);

//----------------------------------------------------------------------------------------
//            Change Coord. System
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the default coordinate system
 **/
void changeDCS(QBCP_L &qbcp_l, int coordsys);
/**
 *  \brief Change the default coordinate system and the libration point for this coordinate system
 **/
void changeLIDCS(QBCP_L &qbcp_l, int coordsys, int li);

//----------------------------------------------------------------------------------------
//            Subroutines
//----------------------------------------------------------------------------------------
/**
 *   \brief Initialize an evenly spaced grid of size gsize, in the interval [gmin, gmax]
 **/
void init_grid(double *grid, double gmin, double gmax, int gsize);;

/**
 *  \brief Return the string corresponding to the libration point number provided (e.g. "L1" if li == 1).
 **/
string init_F_LI(int li);

/**
 *  \brief Return the string corresponding to the model indix provided (e.g. "QBCP" if model == M_QBCP).
 **/
string init_F_MODEL(int model);

/**
 *  \brief Return the string corresponding to the framework (coord. syst.) indix provided (e.g. "EM" if coordsys == F_EM).
 **/
string init_F_FRAMEWORK(int coordsys);

/**
 *  \brief Return the folder name corresponding to the prefix/model/framework/libration point number combination provided (e.g. "prefix/QBCP/EM/L1").
 **/
string init_F_FOLDER(string prefix, int model, int coordsys, int li);

/**
 * \brief Retrieve a set of coefficients, given as Fourier series from a txt file.
 * \param filename the name of the txt file.
 * \param params a pointer toward the array to update.
 * \param nf the order of the Fourier series.
 * \param shift the indix from which to start the storage of the coefficients in params.
 * \param flag: if flag == 1, the coefficients computed via FFT are used. Otherwise, the expansions obtained through Fourier series algebraic manipulations are used.
 *
 *  Warning: As of now, FFT coefficients must be used for betas and deltas (see QBCP_L structure).
 **/
void coefRetrieving(string filename, double *params,  int nf, int shift, int flag, int number);

/**
 * \brief Compute the potential energy for the given state and CR3BP (mu).
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double energy(double y[], double mu);

/**
 * \brief Compute the Jacobi constant for the given state and CR3BP (mu)
 * \param y  the state array.
 * \param mu the mass ratio of the current CR3BP.
 **/
double jacobi(double y[], double mu);

/**
 * \brief Compute the coefficient cn for a given libration point (L1 or L2 for now)
 * \param qbcp_l a reference to the QBCP_L initialized around the selected libration point.
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 **/
double cn(QBCP_L& qbcp_l, int n);

/**
 * \brief Compute the coefficient cn for a given libration point (L1, L2, and L3 for now)
 * \param li the number of the current libration point (1,2,3)
 * \param gamma the gamma parameter associated to the current libration point
 * \param mu the mass ratio of the current TBP system
 * \param n the indix of the coefficient to compute.
 *
 * We recall that the cn coefficient are given by the following equation:
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( (\pm1)^n \mu + (-1)^n \frac{(1-\mu) \gamma_j^{(n+1)}}{(1 -\mp \gamma_j)^{(n+1)}} \right) \f$ for \f$ L_j, j = 1,2 \f$, and where:
 * - The upper (resp. lower) signs stand for \f$ L_1 \f$ (resp. \f$ L_2 \f$).
 * - \f$ \gamma_j \f$ is the distance from \f$ L_j \f$ to the smallest primary.
 *
 * The L3 value is taken from Richardson, 1980:
 *
 * \f$  c_n = \frac{1}{\gamma_j^3} \left( 1 - \mu  + \frac{\mu \gamma_j^{n+1}}{(1 + \gamma_j)^{n+1}} \right) \f$.
 **/
double cn(int li, double gamma, double mu, int n);

/**
* \brief Using the Newton-Raphson method, find the root of a function known to lie close to x1 The root
*  /rtnewt will be refined until its accuracy is known within ± xacc.
*  funcd is a user-supplied routine that returns both the function value and the first derivative of the
*  function at the point x.
**/
double rtnewt(void (*funcd)(double, int, double, double *, double *), double x1, double xacc, double mu, int number);

/**
* \brief Provides the function value and its first derivative for the newton-raphson method.
* f corresponds to the equation satisfied by the Li-m2 distance for the L1/L2 cases
* and by 1-(Li-m1 distance) for the L3 case
**/
void polynomialLi(double mu, int number, double y, double *f, double *df);


/**
 *   \brief Number to string inner routine
 **/
string numTostring(double num);

#endif // ENV_H_INCLUDED
