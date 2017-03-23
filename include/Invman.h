#ifndef INVMAN_H
#define INVMAN_H

//C++
#include <vector>
//Custom
#include "Oftsc.h"
#include "matrix.h"
#include "env.h"
#include "timec.h"
#include "pmcoc.h"


class Invman
{
    private:
        //--------------------------------------------------------------------------------
        // Additionnal parameters
        //--------------------------------------------------------------------------------
        Ofsc ofs;                   //temporary Ofsc object for computation
        CSYS const *cs;             //coordinate system
        double n;                   //normalized mean motion associated to the param.
        int ofts_order;             //order of the ofts objects inside the manifold.
        int ofs_order;              //order of the ofs objects inside the manifold.
        int pmType;                 //type of parameterization.
        int manType;                //type of manifold.
        int reduced_nv;             //number of reduced variables
        int fwrk;                   //framework (EM or SEM?)
        int ncs;                    //Normalized Coordinate System (NCEM or NCSEM?)
        int scs;                    //Non-normalized Coordinate System (PEM or PSEM?)

        //--------------------------------------------------------------------------------
        // Change of coordinates TFC <-> NC
        //--------------------------------------------------------------------------------
        matrix<Ofsc>  Mcoc;         //COC matrix
        matrix<Ofsc>  MIcoc;        //COC matrix = inv(Mcoc)
        vector<Ofsc>  Vcoc;         //COC vector

        //--------------------------------------------------------------------------------
        // Jacobian of the parameterization
        //--------------------------------------------------------------------------------
        matrix<Oftsc> DWh;             //Jacobian of the parameterization in TFC coordinates
        gsl_matrix_complex *CCM_R_RCM; //= dsCCM/dsRCM

        //--------------------------------------------------------------------------------
        // Parameterization
        //--------------------------------------------------------------------------------
        vector<Oftsc> Wh;           //Parameterization in TFC coordinates
        vector<Oftsc> W;            //Parameterization in NC coordinates

        //--------------------------------------------------------------------------------
        //Additionnal objects (see constructor for details)
        //--------------------------------------------------------------------------------
        vector<Oftsc> Hy;           //Additionnal parameterization object for CS/CU.

        //--------------------------------------------------------------------------------
        //For projection purposes
        //--------------------------------------------------------------------------------
        double omega1;
        double omega3;

        //--------------------------------------------------------------------------------
        //Temporary objects (see constructor for details)
        //--------------------------------------------------------------------------------
        matrix<Ofsc> mIn;

        //--------------------------------------------------------------------------------
        //For COC: NCEM <-> NCSEM
        //--------------------------------------------------------------------------------
        vector<Ofsc> nc_B_SYS;
        matrix<Ofsc> nc_R_SYS;
        matrix<Ofsc> SYS_R_nc;

        matrix<Ofsc> VSYS_R_SYS;
        matrix<Ofsc> SYS_R_VSYS;

    public:
        Invman(int ofts_order_, int ofs_order_, CSYS& csys);
        ~Invman();

        //--------------------------------------------------------------------------------
        //Getters
        //--------------------------------------------------------------------------------
        const vector<Ofsc>& getVcoc() const;
        const matrix<Ofsc>& getMcoc() const;
        const matrix<Ofsc>& getMIcoc() const;
        const matrix<Oftsc>& getDWh() const;
        const vector<Oftsc>& getWh()  const;
        const vector<Oftsc>& getW()   const;
        const vector<Oftsc>& getHy()  const;
        const CSYS* getCS() const;
        double getOmega1() const;
        double getOmega3() const;
        double getN() const;
        int getPmType() const;
        int getManType() const;
        int getRnv() const;
        int getFwrk() const;
        int getNCS() const;
        int getSCS() const;
        int getOFTSORDER() const;

        //--------------------------------------------------------------------------------
        //Evaluate
        //--------------------------------------------------------------------------------
        void evalCCMtoTFC(cdouble const  s0[], vector<Ofsc> &zIn, const int ofts_order, const int ofs_order) const;
        void evalRCMtoTFC(double const st0[], vector<Ofsc>& zIn, const int ofts_order, const int ofs_order) const;
        void evalDCCMtoTFC(cdouble const  s0[], matrix<Ofsc> &mIn, const int ofts_order, const int ofs_order) const;
        void evalRCMtoNC(double const  st0[], double const t, double z1[], const int ofts_order, const int ofs_order) const;
        void evalCCMtoNC(cdouble const s0[], double const t, double z1[], const int ofts_order, const int ofs_order) const;
        void evalDRCMtoTFC_partial(double const  st0[], double const t, gsl_matrix_complex *m1, const int ofts_order, const int ofs_order) const;
        void evalDRCMtoNC(double const  st0[], double const t, gsl_matrix *m1, const int ofts_order, const int ofs_order) const;
        void evalDRCMtoCOORD(double const  st0[], double const t, gsl_matrix *m1, const int ofts_order, const int ofs_order, const int coord_type) const;
        void evaldotRCMtoNC(double const st0[], double const t, double z1[], const int ofts_order, const int ofs_order) const;
        void evaldotRCMEMtoNCSEM(double const st0[], double const t, gsl_vector *z1, const int ofts_order, const int ofs_order, const Invman &invman_SEM) const;

        void eval_IN_B_SYS(double const t,
                           gsl_vector_complex *IN_B_SYS,
                           gsl_vector_complex *dIN_B_SYS,
                           gsl_matrix_complex *IN_R_SYS,
                           gsl_matrix_complex *dIN_R_SYS,
                           gsl_matrix_complex *SYS_R_IN,
                           gsl_matrix_complex *dSYS_R_IN,
                           const int ofs_order) const;

        //--------------------------------------------------------------------------------
        //          Projection on (un)stable manifold
        //--------------------------------------------------------------------------------
        void NCprojCCMtoCM(double *yv, double tv, double sti[5]);

        //--------------------------------------------------------------------------------
        //          Energy
        //--------------------------------------------------------------------------------
        /**
         *   \brief Computes the hamiltonian at the position st0, in outputType coordinates and units.
         **/
        double H_SYS(double st0[], double t0);
        /**
         *   \brief Computes the hamiltonian at the position st0, in outputType coordinates and units, at order ofts_order_0
         **/
        double H_SYS(double st0[], double t0, int ofts_order_0);

        //--------------------------------------------------------------------------------
        //Static
        //--------------------------------------------------------------------------------
        static int compRNV(CSYS& csys);
    protected:

};

//========================================================================================
// Tests
//========================================================================================
void test_evalCCMtoTFC();
void test_evalRCMtoNC();
void test_evalDCCMtoTFC();
void test_evalDRCMtoTFC();

//========================================================================================
//
//          SUBROUTINES
//
//========================================================================================
/**
 *  \brief Set the matrix CCM_R_RCM_C  as the rotation matrix between CCM and RCM
 *         coordinates in a center manifold (4 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER(gsl_matrix_complex *CCM_R_RCM_C);

/**
 *  \brief Set the matrix CCM_R_RCM_CH  as the rotation matrix between CCM and RCM
 *         coordinates in a center-hyperbolic manifold (5 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER_HYP(gsl_matrix_complex *CCM_R_RCM_CH);

/**
 *  \brief Set the matrix CCM_R_RCM_CH  as the rotation matrix between CCM and RCM
 *         coordinates in the complete phase space manifold (6 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER_6(gsl_matrix_complex *CCM_R_RCM_CH);


/**
 *  \struct RefineH
 *  \brief Structure for the refinement of the energy
 **/
typedef struct RefineH RefineH;
struct RefineH
{
    Invman *invman;
    int     vdim;
    double  *st0;
    double    t0;
    double    Hv;
    int       order;
    gsl_root_fsolver *s_root;
    double eps_root;

    /**
     *  \brief Constructor for RefineH
     **/
     RefineH(Invman *invman_, int vdim_, double *st0_, double t0_, double Hv_)
     {
         invman   = invman_;
         vdim     = vdim_;
         st0      = st0_;
         t0       = t0_;
         Hv       = Hv_;
         s_root   = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
         eps_root = Config::configManager().G_PREC_ROOT();
         order    = invman->getOFTSORDER();
     }
};


/**
 *   \brief Initialize an orbit wrt a Poincare map so that H(orbit.s0) = H(Pmap)
 **/
int init_s0_energy(RefineH* refineH, double st0[], double t0);


#endif // INVMAN_H
