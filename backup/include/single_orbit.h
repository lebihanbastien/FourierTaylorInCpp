#ifndef SINGLE_ORBIT_H_INCLUDED
#define SINGLE_ORBIT_H_INCLUDED

#include "env.h"
#include "ode.h"
#include "pmcoc.h"
#include "vf.h"
#include "eminsem.h"

#define TYPE_ORBIT         0
#define TYPE_CU            1
#define TYPE_CS            2
#define TYPE_MAN           3
#define TYPE_MAN_SORT_DR   4
#define TYPE_MAN_SORT_DH   5
#define TYPE_MAN_PROJ      6
#define TYPE_MAN_SORT      7
#define TYPE_MAN_SORT_IN   8


#define GSIZE 50;

extern "C"{
 #include "nrutil.h"
 #include "gnuplot_i.h"
}

template <int d1,int d2=1,int d3=1,int d4=1>
class TensorIndex {
  public:
	enum {SIZE = d1*d2*d3*d4 };
	enum {LEN1 = d1 };
	enum {LEN2 = d2 };
	enum {LEN3 = d3 };
	enum {LEN4 = d4 };

    static int indexOf(const int i) {
      return i;
    }
    static int indexOf(const int i,const int j) {
      return j*d1 + i;
    }
    static int indexOf(const int i,const int j, const int k) {
      return (k*d2 + j)*d1 + i;
    }
    static int indexOf(const int i,const int j, const int k,const int l) {
      return ((l*d3 + k)*d2 + j)*d1 + i;
    }
};


/**
 *  \struct IdxCompare
 *  \brief  Structure for comparison of indexes (see sort_indexes);
 **/
struct IdxCompare
{
    const std::vector<double>& target;
    IdxCompare(const std::vector<double>& target): target(target) {}
    bool operator()(int a, int b) const
    {
        return target[a] < target[b];
    }
};

/**
 *  Routine for comparison of indexes
 **/
vector<size_t> sort_indexes(const vector<double> &v);


/**
 *  \struct SingleOrbit
 *  \brief  Defines a given orbit with proper arrays to store results.
 **/
typedef struct SingleOrbit SingleOrbit;
struct SingleOrbit
{
    //-----------
    //Parent
    //-----------
    //-----------
    QBCP_L   *qbcp_l;              //QBCP around a given Li point (parent)

    //-----------
    //Parameterization (common to all orbits)
    //-----------
    vector<Oftsc>*  W;             //z(t) = W(s(t), t)
    vector<Oftsc>*  Wh;            //zh(t) = Wh(s(t), t)
    Ofsc* ofs;                     //Auxiliary Ofs object
    double  n;                     //Pulsation of the QBCP
    int order;                     //order of the pm
    int ofs_order;                 //order of the Fourier coefficients
    int reduced_nv;                //reduced number of variables
    bool isGS;                     //was the pm obtained through graph style?

    //-----------
    //COC (common to all orbits)
    //-----------
    matrix<Ofsc>* PC;               //COC matrix
    matrix<Ofsc>* CQ;               //COC matrix
    vector<Ofsc>* V;                //COC vector

    //Characteristics
    //-----------
    double   *z0;                    //Initial position in NC coordinates dim = 6
    double   *si;                    //Initial RCM configuration dim = REDUCED_NV
    double   *s0d;                   //Initial position in CCM8 coordinates (real+imag part) dim = 2*REDUCED_NV
    cdouble  *s0;                    //Initial position in CCM8 coordinates (real+imag part) dim = 4
    double   *xf;                    //Final position dim = 6
    double    tf;                    //final time after computation
    double    t0;                    //initial time
    double    tproj;                 //default time between each projection
    double    tprojmin;              //minimum time between each projection
    double    ePmax;                 //maximum projection distance allowed

    //-----------
    //ODE integration
    //-----------
    OdeStruct *driver;              //NC ode struct
};


/**
    \brief Initialize one SingleOrbit structure
 **/
void init_orbit(SingleOrbit &orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Ofsc>*  PC,
                matrix<Ofsc>*  CQ,
                vector<Ofsc>*  V,
                Ofsc* orbit_ofs,
                int order,
                int ofs_order,
                int reduced_nv,
                int isGS,
                double t0,
                double tf,
                double tproj,
                OdeStruct *driver,
                QBCP_L *qbcp_l);

/**
    \brief Free one orbit
 **/
void free_orbit(SingleOrbit *orbit);

/**
    \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0);

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn);

/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initian conditions st0, on the center manifold CM.
 **/
 int gridOrbit(double st0[],
                            double t0,
                            double tf,
                            double dt,
                            vector<Oftsc> &CM,
                            vector<Oftsc> &CMh,
                            matrix<Ofsc>  &Mcoc,
                            matrix<Ofsc>  &Pcoc,
                            matrix<Ofsc>  &MIcoc,
                            matrix<Ofsc>  &PIcoc,
                            vector<Ofsc>  &Vcoc);
/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file.
 **/
int cusOrbit(int Order,
             int Size,
             double epsilon,
             vector<Oftsc> &CMh,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &MIcoc,
             vector<Ofsc>  &Vcoc,
             bool isPar);

int manOrbit(double tman,
             int Order,
             int Size,
             int Nt,
             int Nman,
             int isPar,
             vector<Oftsc> &CM,
             vector<Oftsc> &CMh,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc);

int manOrbitSEM(double tman,
             int Order,
             int Size,
             int Nt,
             int Nman,
             int isPar,
             vector<Oftsc> &CM,
             vector<Oftsc> &CMh,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc);

void manOrbitPlot(int Order,
                 int Size,
                 int N,
                 int Nman,
                 int Type);

void manOrbitPostProcess(int Order,
                         int Size,
                         int Nt,
                         int isPar);
void manOrbitProj(int order_em,
                  int size_em,
                  int type_em,
                  int order_sem,
                  int Nt,
                  int Nmant,
                  matrix<Ofsc>  &Mcoc,
                  matrix<Ofsc>  &Pcoc,
                  matrix<Ofsc>  &MIcoc,
                  matrix<Ofsc>  &PIcoc,
                  vector<Ofsc>  &Vcoc,
                  int isPar);

void manOrbitLambert(int order_em,
                     int size_em,
                     int type_em,
                     int order_sem,
                     vector<Oftsc> &CM,
                     vector<Oftsc> &CMh,
                     matrix<Ofsc>  &Mcoc,
                     matrix<Ofsc>  &Pcoc,
                     matrix<Ofsc>  &MIcoc,
                     matrix<Ofsc>  &PIcoc,
                     vector<Ofsc>  &Vcoc);
/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli2(double st0[],
                 double t0,
                 double tf,
                 double tman,
                 double tplot,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CMh,
                 matrix<Ofsc>  &Mcoc,
                 matrix<Ofsc>  &Pcoc,
                 matrix<Ofsc>  &MIcoc,
                 matrix<Ofsc>  &PIcoc,
                 vector<Ofsc>  &Vcoc);

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initian conditions st0, on the center-unstable manifold CM.
 **/
int orbit_cus(double st0[],
             double t0,
             double tf,
             double tman,
             vector<Oftsc> &CM,
             vector<Oftsc> &CMh,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc);

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initian conditions st0, on the center-unstable manifold CM (SEM framework).
 **/
int orbit_cu_sem(double st0[],
                 double t0,
                 double tf,
                 double tman,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CMh,
                 matrix<Ofsc>  &Mcoc,
                 matrix<Ofsc>  &Pcoc,
                 matrix<Ofsc>  &MIcoc,
                 matrix<Ofsc>  &PIcoc,
                 vector<Ofsc>  &Vcoc);

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli_torb_vs_tman(double st0[],
                double t0,
                double tf,
                double tman,
                vector<Oftsc> &CM,
                vector<Oftsc> &CMh,
                matrix<Ofsc>  &Mcoc,
                matrix<Ofsc>  &Pcoc,
                matrix<Ofsc>  &MIcoc,
                matrix<Ofsc>  &PIcoc,
                vector<Ofsc>  &Vcoc);

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int semliToseml2(double st0[],
                double t0,
                double tf,
                double tman,
                vector<Oftsc> &CM,
                vector<Oftsc> &CMh,
                matrix<Ofsc>  &Mcoc,
                matrix<Ofsc>  &Pcoc,
                matrix<Ofsc>  &MIcoc,
                matrix<Ofsc>  &PIcoc,
                vector<Ofsc>  &Vcoc);


int writeCU_bin(double ****yNCE, double ****sNCE, double *tGrid, int gSize, int tSize, int order, int type);

int readManifold_bin(double ****yNCE, double ****sNCE, double *tGrid, int gSize, int tSize, int order, int type);

int getLenghtCU_bin(int *gSize, int *tSize, int order, int type);

/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file (obtained from gridOrbit).
 **/
int cusMan(int order,
           double epsilon,
           double tMin,
           double tMax,
           int tSize,
           double gMin,
           double gMax,
           int gSize,
           vector<Oftsc> &CMh,
           matrix<Ofsc>  &Mcoc,
           matrix<Ofsc>  &MIcoc,
           vector<Ofsc>  &Vcoc,
           bool isPar);

/**
 *  \brief Computes the manifold branches from a discrete set of unstable directions obtained from cusMan
 **/
int intMan(double tman,
           int order,
           int tSizex,
           int gSizex,
           int mSize,
           int isPar,
           vector<Oftsc> &CM,
           vector<Oftsc> &CMh,
           matrix<Ofsc>  &Mcoc,
           matrix<Ofsc>  &Pcoc,
           matrix<Ofsc>  &MIcoc,
           matrix<Ofsc>  &PIcoc,
           vector<Ofsc>  &Vcoc);

void projMan(int order_em,
             int order_sem,
             int tSizex,
             int gSizex,
             int mSizex,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc,
             int isPar);

int intAndProjMan(double tman,
               int order,
               int tSizex,
               int gSizex,
               int mSize,
               int isPar,
               vector<Oftsc> &CMEM,
               vector<Oftsc> &CMhEM,
               matrix<Ofsc>  &McocEM,
               matrix<Ofsc>  &PcocEM,
               matrix<Ofsc>  &MIcocEM,
               matrix<Ofsc>  &PIcocEM,
               vector<Ofsc>  &VcocEM);

int intSortedMan(int order,
                 int mSize,
                 int isPar,
                 vector<Oftsc> &CMEM,
                 vector<Oftsc> &CMhEM,
                 matrix<Ofsc>  &McocEM,
                 matrix<Ofsc>  &PcocEM,
                 matrix<Ofsc>  &MIcocEM,
                 matrix<Ofsc>  &PIcocEM,
                 vector<Ofsc>  &VcocEM);
//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Integration
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(SingleOrbit &orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *ePm,
                   int *nreset,
                   int isResetOn);


/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(SingleOrbit &orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *ePm,
                     int *nreset,
                     int isResetOn);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Projection on (un)stable manifold
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CMh, MIcoc, Vcoc). Then the CCM state is extended by adding a non-null direction
 *        along the hyperbolic direction (sti[4]).
 **/
void NCprojCCMtoCUS(double *yv, double tv, double sti[5], vector<Oftsc> &CMh, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc);

/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CMh, MIcoc, Vcoc).
 **/
void NCprojCCMtoCM(double *yv, double tv, double sti[5], vector<Oftsc> &CMh, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc);

//---------------------------------------------------------------------------------------------------------------------------------------
//                  Energy on vectors
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian along one orbit in SEM units and SEM coordinates
 **/
void HSEM_vec(double *tSEM, double **ySEM, double *Hvec, int N, QBCP_L *qbcp_l);

/**
 *  \brief Hamiltonian along the dyneq of SEMLi for a given time vector in SEM units and SEM coordinates
 **/
void HSEMLi_vec(double *tSEM, double *Hvec, int N, QBCP_L *qbcp_l);

/**
 *  \brief Hamiltonian along the dyneq of EMLi for a given time vector in SEM units and SEM coordinates
 **/
void HEMLi_in_SEM_vec(double *tEM, double *Hvec, int N, QBCP_L *qbcp_l);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//         Update points
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Update some key positions of notable points in the EM system
 **/
void emPoints(double t, double **emP);

/**
 *  \brief Update some key positions of notable points in the SEM system
 **/
void semPoints(double t, double **semP);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//        Plots
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Sets notable points in EM system on the gnuplot window ctrl h1
 **/
void emPlot(gnuplot_ctrl *h1, double **emP);

/**
 *  \brief Sets notable points in SEM system on the gnuplot window ctrl h2
 **/
void semPlot(gnuplot_ctrl *h2, double **semP);

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          I/O orbit
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computes a data filename as a string, depending on the order, size, and type of the data.
 **/
string filenameOrbit(int order, int sizeOrbit, int type);

/**
 * \brief Store the orbit (tNCE, yNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, int order, int sizeOrbit, int N, int type);

/**
 * \brief Store the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, double **sNCE, int order, int sizeOrbit, int N, int type);

/**
 * \brief Get the length of the data file filenameOrbit(order, sizeOrbit, type).
 **/
int getLineNumber(int order, int sizeOrbit, int type);

/**
 * \brief Read the orbit (tNCE, yNCE) in the  data file filenameOrbit(order, size, type)
 **/
int readOrbit(double *tNCE, double **yNCE, int order, int sizeOrbit, int N, int type);

/**
 * \brief Read the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int readOrbit(double *tNCE, double **yNCE, double **sNCE, int order, int sizeOrbit, int N, int type);

#endif // SINGLE_ORBIT_H_INCLUDED
