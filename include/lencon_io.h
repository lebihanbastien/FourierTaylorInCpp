#ifndef SINGLE_ORBIT_IO_H_INCLUDED
#define SINGLE_ORBIT_IO_H_INCLUDED

#include <iostream>     // std::cout
#include <algorithm>    // std::unique, std::distance
#include <vector>       // std::vector
#include <iterator>

#include <cstdio>
#include <cstring>
#include <cerrno>
#include <sys/stat.h>

#include "single_orbit.h"
#include "Orbit.h"
#include "tinyfiledialogs.h"

//For selection of the times
#define TIME_SELECTION_ABSOLUTE 1
#define TIME_SELECTION_RATIO    0

//Number of digits in the ratios
#define RATIO_DIGITS 7

//----------------------------------------------------------------------------------------
//Parameters for the type of refinements (used in structure RefSt)
//----------------------------------------------------------------------------------------
#define REF_PLANAR     0
#define REF_3D         1
#define REF_MIXED      101

#define REF_SINGLE     2
#define REF_ORBIT      21

#define REF_CONT             30
#define REF_CONT_D           31
#define REF_CONT_D_HARD_CASE 32
#define REF_CONT_ORBIT       33

#define REF_COMP       4

#define REF_FIXED_TIME 5
#define REF_VAR_TN     6
#define REF_VAR_TIME   7

#define REF_FIXED_GRID 8
#define REF_VAR_GRID   9
#define REF_GIVEN_GRID 10

#define REF_COND_S5    11
#define REF_COND_T     12

//----------------------------------------------------------------------------------------
//I/O handling
//----------------------------------------------------------------------------------------
#define IO_DEFAULT 0
#define IO_BASH    1
#define IO_DIALOG  2

using namespace std;


//========================================================================================
//
//          I/O (filename routines)
//
//========================================================================================
/**
 *  \brief Prefix for the data filenames.
 **/
string fileprefix(int type);

/**
 *  \brief Extension for the data filenames.
 **/
string fileext(int type);

///**
// *  \brief Computes a data filename as a string, depending on the ofts_order, type,
// *         and target of the data.
// **/
//string filenameCUM(string plot_folder, int ofts_order, int type, int target);
//
///**
// *  \brief Computes a data filename as a string, depending on the ofts_order, type,
// *         target, and initial time t0 of the data.
// **/
//string filenameCUM(string plot_folder, int ofts_order, int type, int target, double t0);
//
///**
// *  \brief Computes a data filename as a string, depending on the ofts_order, type,
// *         target, and initial energy delta dH of the data.
// **/
//string filenameCUM_dH(string plot_folder, int ofts_order, int type, int target, double dH);

/**
 *  \brief Get filename via a window file dialog, with default name filename_default.
 **/
string get_filename_dialog(string filename_default);

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, type,
 *         target, initial time t0, and initial energy delta dH of the data.
 **/
string get_filenameCUM(int IO_HANDLING, string plot_folder, string filename_bash,
                       int ofts_order, int type, int target, double t0, double dH, int mode);

/**
 * \brief Structure to store average amplitude results.
 **/
typedef struct AvgSt AvgSt;
struct AvgSt
{
    //Indices
    int eml_index;
    int seml_index;
    int eml_fperiod_index;
    int seml_lperiod_index;
    int eml_nindices;
    int seml_nindices;

    //Constants
    double eml_kappa;
    double eml_omegav;
    double eml_gamma;
    double seml_kappa;
    double seml_omegav;
    double seml_gamma;

    //For the EML orbit
    double Ax_EM_mean;   //mean of Ax = mean(sqrt(x_NCEM^2 + y_NCEM^2/kappa))*gamma
    double Az_EM_mean;   //mean of Az = mean(sqrt(z_NCEM^2 + zdot_NCEM^2/omegav_CRTBP))*gamma
    double Az_EM_lsf;    //Least Square Fitting to find (Az, omega) that matches of z_NCEM^2 + zdot_NCEM^2/omega = Az^2/gamma^2
    double Omega_EM_lsf; //The corresponding frequency

    //For the SEML orbit
    double Ax_SEM_mean;   //mean of Ax = mean(sqrt(x_NCSEM^2 + y_NCSEM^2/kappa))*gamma
    double Az_SEM_mean;   //mean of Az = mean(sqrt(z_NCSEM^2 + zdot_NCSEM^2/omegav_CRTBP))*gamma
    double Az_SEM_lsf;    //Least Square Fitting to find (Az, omega) that matches of z_NCSEM^2 + zdot_NCSEM^2/omega = Az^2/gamma^2
    double Omega_SEM_lsf; //The corresponding frequency

    AvgSt():eml_index(0), seml_index(0), eml_fperiod_index(0),
    seml_lperiod_index(0), eml_nindices(0), seml_nindices(0){};


    void copy(AvgSt &avgSt_dest)
    {
        avgSt_dest.eml_index = eml_index;
        avgSt_dest.seml_index = seml_index;

        avgSt_dest.eml_fperiod_index = eml_fperiod_index;
        avgSt_dest.seml_lperiod_index = seml_lperiod_index;

        avgSt_dest.eml_nindices = eml_nindices;
        avgSt_dest.seml_nindices = seml_nindices;

        avgSt_dest.eml_kappa  = eml_kappa;
        avgSt_dest.eml_omegav = eml_omegav;
        avgSt_dest.eml_gamma  = eml_gamma;

        avgSt_dest.seml_kappa  = seml_kappa;
        avgSt_dest.seml_omegav = seml_omegav;
        avgSt_dest.seml_gamma  = seml_gamma;
    }
};

//========================================================================================
// RefSt structure
//========================================================================================
/**
 *  \struct RefSt
 *  \brief  Define a given refinement structures, with a set of parameters
 **/
typedef struct RefSt RefSt;
struct RefSt
{
    //------------------------------------------------------------------------------------
    // ... filename handling is private
    //------------------------------------------------------------------------------------
    private:

    //filename for the data file output of routines such as compute_grid_CMU_EM_3D/int_proj_CMU_EM_on_CM_SEM_3D
    string filename_output;

    //------------------------------------------------------------------------------------
    // Most of the structure is public to ease use...
    //------------------------------------------------------------------------------------
    public:

    //------------------------------------------------------------------------------------
    // General parameters
    //------------------------------------------------------------------------------------
    int    OFTS_ORDER, LI_EM, LI_SEM, LI_START, LI_TARGET;
    int    IO_HANDLING;
    string plot_folder;
    string FILE_PCU, FILE_CONT, FILE_CONT_RES, FILE_TRAJ_FROM_W, FILE_TRAJ_FROM_C, FILE_JPL_TXT, FILE_JPL_BIN, FILE_MEAN_AMP, FILE_FOR_CELESTIA;

    //------------------------------------------------------------------------------------
    // Parameters that change often
    //------------------------------------------------------------------------------------
    int type;             //single solution or continuation procedure
    int dim;              //planar or 3d
    double t0_des;        //desired  initial time
    double t0xT_des;      //desired  initial time as a percent

    // Limits for domain of research of the first guess
    double si_CMU_EM_MIN[4];
    double si_CMU_EM_MAX[4];

    // Limits for domain of research of the first guess - seeds
    double si_SEED_EM_MIN[4];
    double si_SEED_EM_MAX[4];

    // The domain of search for first guess fixed by the user if true
    int isLimUD;

    // Direction of the continuation procedure
    int isDirUD;          //the direction of refinement is fixed by the user if true
    int Dir;              //the direction of refinement if isDirUD = false

    // Limits for the time of flight during transfers - not used if negative
    double tof_MIN;
    double tof_MAX;

    // Values for crossings
    double crossings;

    // Values for collisions
    int isCollisionOn;

    // Maximum number of steps in the continuation procedure
    int cont_step_max;    //with fixed time
    int cont_step_max_vt; //with variable time

    // Initial & current step in the continuation procedure
    double ds0;           //with fixed time
    double ds0_vt;        //with variable time
    double dsc;           //current step

    // Desired number of iterations in Newton's method in the continuation procedure
    int nu0;              //with fixed time
    int nu0_vt;           //with variable time

    // User parameters
    int isFlagOn;         //do we have steps in the procedure - asking the user to press enter to go on?
    int isPlotted;        //do we plot the results during the computation?
    int isSaved;          //do we save the results in data files?
    int isFromServer;     //does the raw data comes from server files?
    int isPar;            //is parallel computation allowed?

    // Maximum angle around SEMLi if REF_COND_T is used (in degrees)
    double thetaMax;      //should be a multiple of 90°

    //------------------------------------------------------------------------------------
    // Parameters that are stable
    //------------------------------------------------------------------------------------
    int isDebug;          //if yes, additionnal tests are performed
    int gridSize;         //number of points on the refinement grid
    int mPlot;            //number of points per plot between to pach points (e.g. total plot points is gridSize*mplot)
    float fHours;         // desired frequency of plotting (in hours)

    int time;             //type of constraints on the times in REF_CONT
    int grid;             //type of grid
    int termination;      //termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
    int coord_type;       //coordinates system in the refinement procedure (usually NCSEM)

    // Maximum/Minimum step in the continuation procedure
    double dsmin;           //with fixed time
    double dsmin_vt;        //with variable time
    double dsmax;           //with fixed time
    double dsmax_vt;        //with variable time

    double xps_NCSEM;     //position of the poincaré section in NCSEM coordinates
    double xps_NCEM;      //position of the poincaré section in NCEM coordinates
    int isJPL;            //is the JPL refinement performed when possible?
    int djplcoord;        //coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with $NJ2000
    int sidim;            //0 or 2 - component of s0 that stays constant when t0 is free.

    //Sampling frequencies in REF_COMP (complete trajectory) in days
    int sf_emli;          // orbit at EML2
    int sf_man;           // transfer leg
    int sf_semli;         // orbit at SEML2

    // Integration window for each orbit
    double tspan_EM;
    double tspan_SEM;

    // Storing the orbits at each step? DEPRECATED, kept for consistency, use comp_orb_em/comp_orb_sem instead
    int isSaved_EM;    //0: don't save, 1: save using projection method
    int isSaved_SEM;   //0: don't save, 1: save using projection method, 2: save using integration in reduced coordinates

    // Type of computation for each orbit
    int comp_orb_em;
    int comp_orb_sem;

    //Type of time selection
    int typeOfTimeSelection;

    //Check if the type is of continuation type
    bool isCont(){return (type == REF_CONT || type == REF_CONT_D || type == REF_CONT_D_HARD_CASE || type == REF_CONT_ORBIT  );}
    // Check if the trajectory are 3D
    bool is3D(){return (dim == REF_3D || dim == REF_MIXED);}

    //Maximum projection distance allowed during subselection
    double pmax_dist_SEM;

    //Last error during differential correction
    double last_error;

    //Inner precision for differential correction procedures
    double inner_prec;
    double inner_prec_vt;
    double inner_prec_ft;

    //Center
    double center[3];

    //Number of solutions to be refined
    int nref;

    //------------------------------------------------------------------------------------
    //Additional variables for specific differential correctors:
    // - goodvector is used to add the pseudo-arclength constraint
    // - dH is used to add a constraint on the (initial) energy
    // - pkpos is used to add a constraint on a certain Pk section
    //------------------------------------------------------------------------------------
    double goodvector[1000]; //vector of free variables
    double dH;               //energy value at the origin
    int pkpos;               //index of the patch point that bear the poincaré section

    double dHf_SEM_MAX;

    //------------------------------------------------------------------------------------
    //Constructor
    //------------------------------------------------------------------------------------
    //RefSt():pmax_dist_SEM(1e-3) {};
    RefSt(int OFTS_ORDER_, int LI_EM_, int LI_SEM_,
          int LI_START_, int LI_TARGET_, int IO_HANDLING_, double RPS, CSYS *cs):
          OFTS_ORDER(OFTS_ORDER_), LI_EM(LI_EM_), LI_SEM(LI_SEM_),
          LI_START(LI_START_), LI_TARGET(LI_TARGET_), IO_HANDLING(IO_HANDLING_),
          plot_folder(cs->F_PLOT),
          isCollisionOn(1), pmax_dist_SEM(1e5), last_error(0.0),
          inner_prec(5e-8), inner_prec_vt(5e-8), inner_prec_ft(5e-8), nref(-1), dHf_SEM_MAX(1.4e-4)
    {
        for(int i = 0; i <4; i++) si_SEED_EM_MIN[i] = -50;
        for(int i = 0; i <4; i++) si_SEED_EM_MAX[i] = +50;

        //Center
        center[0] = (cs->li == 1)? 1:-1;
        center[1] = 0; center[2] = 0;

        //Pk section
        xps_NCEM = RPS;
        if(xps_NCEM < 0) xps_NCEM = SEML.cs->r3BSOI;

        //Frequency by default: every 5 hours
        fHours = 5.0;
    };

    /**
     *  \brief Update the filename and return it.
     **/
    string get_and_update_filename(string filename_bash, int type, int mode)
    {
        this->filename_output = get_filenameCUM(this->IO_HANDLING, this->plot_folder,
                                                filename_bash, this->OFTS_ORDER, type,
                                                this->LI_TARGET, t0xT_des, -1, mode);
        return this->filename_output;
    }
};

//========================================================================================
// ProjResSt structure
//========================================================================================
/**
 *  \struct ProjResSt
 *  \brief  Define a given structure to store the results from a projection
 **/
typedef struct ProjResSt ProjResSt;
struct ProjResSt
{
    //------------------------------------------------------------------------------------
    // Coordinate systems:
    //
    //  IC_COORD: Initial conditions are saved in IC_COORD coordinates
    //  PR_COORD: Projection is made in PR_COORD coordinates
    //  FC_COORD: Final conditions are saved in FC_COORD & FV_COORD coordinates
    //  FV_COORD: Final conditions are saved in FC_COORD & FV_COORD coordinates
    //  DP_COORD: The distance of projection is saved in PSEM coordinates,
    //            hence DP_COORD is always equal to PSEM, for now
    //------------------------------------------------------------------------------------
    int IC_COORD, PR_COORD, FC_COORD, FV_COORD, DP_COORD;

    //------------------------------------------------------------------------------------
    //Inputs, in IC_COORD
    //------------------------------------------------------------------------------------
    double init_time, init_state_CMU_NC[6], init_state_CMU_RCM[5];
    int label;

    //------------------------------------------------------------------------------------
    //Seeds
    //------------------------------------------------------------------------------------
    int seed_label;
    double seed_time;
    double seed_state_CMU_RCM[4];

    //------------------------------------------------------------------------------------
    //Outputs, in complementary system coordinates: FC_COORD & FV_COORD
    //------------------------------------------------------------------------------------
    double init_state_CMU_FC_o[6], final_state_CMU_FC_o[6];
    double projected_state_CMU_FC_o[6], projected_state_CMU_RCM_o[4];
    double dv_at_projection_FC_o;
    double final_time_FC_o;

    //------------------------------------------------------------------------------------
    // Objects that remain always in fixed coordinates
    //------------------------------------------------------------------------------------
    double min_proj_dist_SEM_o;
    double crossings_NCSEM_o;
    int collision_NCEM_o;

    //------------------------------------------------------------------------------------
    // At the EML2 Pk section
    //------------------------------------------------------------------------------------
    double ye_NCEM[6], te_NCEM, ye_NCSEM[6], te_NCSEM;
    double ve_NCEM[3], ve_NCSEM[3];

    //------------------------------------------------------------------------------------
    //Constructor, wrt NC_COORD
    //------------------------------------------------------------------------------------
    ProjResSt(int NC_COORD):DP_COORD(PSEM)
    {
        switch(NC_COORD)
        {
        case NCSEM:
        {
            // From SEML to EML
            IC_COORD = NCSEM;
            PR_COORD =  NCEM;
            FC_COORD =   PEM;
            FV_COORD =   VEM;
            break;
        }

        case NCEM:
        {
            // From EML to SEML
            IC_COORD =  NCEM;
            PR_COORD = NCSEM;
            FC_COORD =  PSEM;
            FV_COORD =  VSEM;
            break;
        }

        default:
        {
            cout << "ProjResSt. NC_COORD must be NCSEM or NCEM." << endl;
            break;
        }

        }
    };

    // Display
    void displayProjResSt()
    {
        cout << "----------------------------" << endl;
        cout << "ProjResSt = "                 << endl;
        cout << "----------------------------" << endl;
        cout << "label     = " << seed_label   << endl;
        cout << "init_time = " << init_time    << endl;

        cout << "init_state_CMU_NC = "         << endl;
        vector_printf_prec(init_state_CMU_NC, 6);

        cout << "init_state_CMU_RCM = "         << endl;
        vector_printf_prec(init_state_CMU_RCM, 5);

        cout << "seed_time = " << seed_time     << endl;
        cout << "seed_state_CMU_RCM = "         << endl;
        vector_printf_prec(seed_state_CMU_RCM, 4);

        cout << "init_state_CMU_FC_o = "         << endl;
        vector_printf_prec(init_state_CMU_FC_o, 6);

        cout << "final_state_CMU_FC_o = "         << endl;
        vector_printf_prec(final_state_CMU_FC_o, 6);

        cout << "projected_state_CMU_FC_o = "         << endl;
        vector_printf_prec(projected_state_CMU_FC_o, 6);

        cout << "projected_state_CMU_RCM_o = "         << endl;
        vector_printf_prec(projected_state_CMU_RCM_o, 4);

        cout << "dv_at_projection_FC_o = " << dv_at_projection_FC_o  << endl;
        cout << "final_time_FC_o       = " << final_time_FC_o        << endl;

        cout << "min_proj_dist_SEM_o = " << min_proj_dist_SEM_o    << endl;
        cout << "crossings_NCSEM_o  = " << crossings_NCSEM_o    << endl;
        cout << "collision_NCEM_o   = "  << collision_NCEM_o    << endl;
    }

};

//========================================================================================
// To get current path
//========================================================================================
#include <limits.h>
#include <unistd.h>

/**
 *  \brief Return the path of the program that is currently running
 **/
string getexepath();

/**
 *  \brief Return the path of the main folder from with the program is currently running.
 *         It is assumed that the binaries are stored in "main/bin/", so that we can
 *         retrieve "main/" by cutting the result of getexepath() before the string "bin/"
 **/
string getmainpath();

//========================================================================================
// Subroutines to manipulate vectors and datafiles
//========================================================================================
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
vector<size_t> sort_indexes(const vector<double>& v);

/**
 *   \brief Resize a vector to its vector of unique elemnts
 **/
void vector_getUnique(vector<double>& v0U);

/**
 *   \brief Save in indRes the indices of the occurences of t0 in t0_CMU_EM.
 **/
void vector_getIndices(vector<size_t>& indRes, vector<double>& t0_CMU_EM, double t0);

/**
 *   \brief Save in indRes the indices of the occurences of t0 in t0_CMU_EM, within a certain range
 **/
void vector_getIndices_range(vector<size_t>& indRes, vector<double>& t0_CMU_EM, double t0);

/**
 *   Check if a file exists
 *   @param[in] filename - the name of the file to check
 *   @return    true if the file exists, else false
 *
 */
bool fileExists(const std::string& filename);

/**
 *   \brief Determine the number of columns in the data file "filename",
 *          given some possible solutions in ncol[], of size nncol.
 *          The underlying hypothesis is: the first two values in the
 *          first column are identical.
 **/
int numberOfColumns(string filename, int ncol[], int nncol);



//========================================================================================
// ProjResClass class
//========================================================================================
/**
 *  \struct ProjResClass
 *  \brief  Define a given structure to store the results from a projection
 **/
class ProjResClass
{
    //====================================================================================
    // Parameters
    //====================================================================================
private:

    //Data
    vector<double> t0_CMU_EM;
    vector<double> tf_CMU_EM;

    vector<double> s1_CMU_EM;
    vector<double> s2_CMU_EM;
    vector<double> s3_CMU_EM;
    vector<double> s4_CMU_EM;
    vector<double> s5_CMU_EM;

    vector<double> pmin_dist_SEM;
    vector<double> tf_man_SEM;
    vector<double> dHf_SEM;

    vector<double> s1_CM_SEM;
    vector<double> s2_CM_SEM;
    vector<double> s3_CM_SEM;
    vector<double> s4_CM_SEM;

    vector<double> crossings_NCSEM;
    vector<double> collision_NCEM;

    vector<double> t0_CMU_EM_seed;
    vector<double> s1_CMU_EM_seed;
    vector<double> s2_CMU_EM_seed;
    vector<double> s3_CMU_EM_seed;
    vector<double> s4_CMU_EM_seed;

    vector<int> label;

    //Sort vector
    vector<size_t> sortId;


    //Size
    int csize;

    //====================================================================================
    //Public routines
    //====================================================================================
public:
    //------------------------------------------------------------------------------------
    // Constructor
    //------------------------------------------------------------------------------------
    ProjResClass():csize(0) {}

    //------------------------------------------------------------------------------------
    // Setters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Push back in all vectors of this, the element kpor of the vectors in inSt.
     *         Update csize accordingly.
     **/
    void push_back(ProjResClass& inSt, int kpor);

    /**
     *  \brief Pop back for all the vectors. Update csize accordingly.
     **/
    void pop_back();

    /**
     *  \brief Push back in all vectors of this, the elements of the vectors in inSt.
     *         that comply with the conditions defined in refSt.
     *         Update csize accordingly.
     **/
    bool push_back_conditional(ProjResClass& inSt, RefSt & refSt);

    /**
     *  \brief Performs a subselection in terms of dHf_SEM: only a few solutions between
     *         min(projSt.dHf_SEM) and max(projSt.dHf_SEM) are selected.
     **/
    bool push_back_subselection(ProjResClass& projSt, RefSt& refSt);

    //------------------------------------------------------------------------------------
    // Getters
    //------------------------------------------------------------------------------------
    /**
     *  \brief Return csize, the common size of all the vectors.
     **/
    int size();

    /**
     *  \brief Return sortId, the vector of sorted indices.
     **/
    vector<size_t>& getSortId();

    //------------------------------------------------------------------------------------
    // Display
    //------------------------------------------------------------------------------------
    /**
     *  \brief Display some selected elements of the entry k.
     **/
    void displayEntry(int k);

    /**
     *  \brief Display some selected elements of the first entry.
     **/
    void displayFirstEntry();

    /**
     *  \brief Display some selected elements of the last entry.
     **/
    void displayLastEntry();

    //------------------------------------------------------------------------------------
    // Read & sort
    //------------------------------------------------------------------------------------
    /**
     *  \brief Read in a data file the connections between EML1,2 and SEML1,2.
     **/
    int readProjRes(string filename);

    /**
    *  \brief Read in a data file the connections between EML2 and SEML1,2.
    *         Subselection in the data set to get the right desired t0 at EML2 departures.
    **/
    int readProjRes_t0(string filename, double t0_des, int typeOfTimeSelection);

    /**
     *  \brief Update sortId, the vector of sorted indices, with respect to pmin_dist_SEM.
     *         After this routine, sortId[0] is the index for which pmin_dist_SEM is minimum.
     **/
    void sort_pmin_dist_SEM();

    /**
     *  \brief Update sortId, the vector of sorted indices, with respect to dHf_SEM.
     *         After this routine, sortId[0] is the index for which dHf_SEM is minimum.
     **/
    void sort_dHf_SEM();

    /**
     *  \brief Update st_EM, st_SEM... With the elements contained in sortId[k].
     **/
    void update_ic(double st_EM[5], double st_SEM[5], double t_EM[2],
                   double st_EM_seed[4], double *t_EM_seed,
                   double* t0_SEM, double* pmin_dist_SEM_out,
                   int *label, int k);

    /**
     *  \brief Update st_EM, st_SEM... With the elements contained in sortId[k].
     **/
    void update_ic(double st_EM[5], double st_SEM[5], double t_EM[2],
                   double* t0_SEM, double* pmin_dist_SEM_out, int k);
};


//========================================================================================
//
//          I/O (complete trajectories)
//
//========================================================================================
/**
 *  \brief Stores the final_index+1 points trajectory contained in t_traj_n (time vector)
 *         and y_traj_n (state vectors) into a data file of type TYPE_COMP_FOR_JPL.
 **/
void write_cref_for_jpl_txt(double* t_traj_n, double** y_traj_n, int final_index, string filename);

/**
 *  \brief Reads the final_index+1 points trajectory contained in a data file of type
 *         TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector) and y_traj_n
 *         (state vectors). The size of the data file is given by getl_cref_for_jpl_txt(), and
 *         the data vectors should be initialized accordingly by the user prior to the use
 *         of this routine.
 **/
int read_cref_for_jpl_txt(double* t_traj_n, double** y_traj_n, int final_index, string filename);

/**
 *  \brief Get the length of the final_index+1 points trajectory contained in a
 *         data file of type TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector)
 *         and y_traj_n (state vectors).
 **/
int getl_cref_for_jpl_txt(string filename);

//========================================================================================
//
//          I/O (CU/CS/CM)
//
//========================================================================================
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1and of size (tGrid +1)*(gSize+1)
 **/
int writeCU_bin_dH(double** * yNCE, double** * sNCE, double** dH, double* tGrid,
                   int s1_grid_size, int t_grid_size, string filename);

/**
 *  \brief Get the length of the data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1 and of size (tGrid +1)*(gSize+1)
 **/
int getLenghtCU_bin_dH(int* s1_grid_size, int* t_grid_size, string filename);

/**
 *  \brief Read in a data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1 and of size (tGrid +1)*(gSize+1)
 **/
int readCU_bin_dH(double** * yNCE, double** * sNCE, double** dH, double* tGrid,
                  int s1_grid_size, int t_grid_size, string filename);

//----------------------------------------------------------------------------------------
// CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
                int s1_grid_size, int s3_grid_size, int t_grid_size, string filename);
/**
 *  \brief Read in a data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
               int s1_grid_size, int s3_grid_size, int t_grid_size, string filename);

/**
 *  \brief Get the length of the data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int getLenghtCU_bin(int* s1_grid_size, int* s3_grid_size,
                    int* t_grid_size, string filename);



//----------------------------------------------------------------------------------------
// CU 3D
//----------------------------------------------------------------------------------------

/**
 *  \brief init the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1)
 **/
int initCU_bin_3D(int* si_grid_size, int t_grid_size, string filename);

/**
 *  \brief Write the current time in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only the time is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int appTimeCU_bin_3D(double* tGrid, int nt, string filename);

/**
 *  \brief Write the current state in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only a single loop on s4 is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int writeCU_bin_3D(double** yNCE, double** sNCE, int* si_grid_size, string filename);

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 **/
int getLenghtCU_bin_3D(int* si_grid_size, int* t_grid_size, string filename);

/**
 *  \brief Read in a data file the time at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readTCU_bin_3D(int offset, double* tGrid, int nt, string filename);

/**
 *  \brief Read in a data file the state at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin_3D(int offset, double** yNCE, double** sNCE, int* si_grid_size, string filename);

//----------------------------------------------------------------------------------------
// Int CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM_3D.
 **/
void writeIntProjCU_bin(string filename,
                        ProjResSt& projResSt);

/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_ORBIT_EM_on_CM_SEM.
 **/
void writeIntProjCUSeed_bin(string filename,
                            ProjResSt& projResSt);


//========================================================================================
//
//          Display completion
//
//========================================================================================


//========================================================================================
//
//          I/O (refinement & continuation)
//
//========================================================================================
int getl_wref_traj_bin(string filename_traj, int *last_index, int *man_grid_size);

void write_wref_traj_bin(string filename_traj, Orbit& orbit_EM, Orbit& orbit_SEM,
                         double **y_traj_NCSEM, double *t_traj_NCSEM,
                         double te_NCSEM,   double* ye_NCSEM,
                         double te_NCEM,    double* ye_NCEM,
                         double ve_NCEM[3], double ve_NCSEM[3],
                         ProjResClass& projRes,
                         int isFirst, int index, int man_grid_size);

int read_wref_traj_bin(string filename_traj, int *labels,
                       double* t0_CMU_EM, double* tf_CMU_EM,
                       double** si_CMU_EM, double** si_CMS_SEM,
                       double** t_traj_NCSEM, double*** y_traj_NCSEM,
                       int man_grid_size,
                       int last_index);

int number_of_plot_points(double deltaT, double fHours, int coord_type);


/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void write_wref_conn_txt(string filename, Orbit& orbit_EM, Orbit& orbit_SEM,
                        double te_NCSEM, double* ye_NCSEM,
                        double te_NCEM, double* ye_NCEM,
                        double ve_NCEM[3], double ve_NCSEM[3],
                        ProjResClass& projRes, int isFirst,  int k);

/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void writeCONT_txt(int label, string filename, Orbit& orbit_EM, Orbit& orbit_SEM,
                   double te_NCSEM, double* ye_NCSEM,  int isFirst);

                   /**
 *  \brief Get the length the results of the continuation procedure, in txt file.
 **/
int getLengthCONT_txt(string filename);

/**
 *  \brief Reads the results of the continuation procedure, in txt file.
 **/
int readCONT_txt(double*  t0_CMU_EM, double*   tf_CMU_EM,
                 double** si_CMU_EM, double** si_CMS_SEM,
                 double** z0_CMU_NCEM, double** z0_CMS_NCSEM,
                 double* tethae, double** ye_NCSEM,
                 double* H0_NCEM, double* He_NCEM,
                 double* H0_NCSEM, double* He_NCSEM,
                 int fsize, string filename);

/**
 *  \brief Save a given solution as a complete trajectory
 **/
int writeCONT_bin(RefSt& refSt, string filename_res, int dcs, int coord_type,
                  double** y_traj_n, double* t_traj_n, int man_index, int mPlot,
                  Orbit &orbit_EM, Orbit &orbit_SEM, int label,
                  bool isFirst, int comp_orb_eml, int comp_orb_seml);

/**
 *  \brief Save a given solution as a complete trajectory
 **/
int write_wref_res_bin(RefSt& refSt, string filename_res, double** y_traj_n,
                       double* t_traj_n, int man_index,
                       Orbit& orbit_EM, Orbit& orbit_SEM,
                       double te_NCSEM,   double* ye_NCSEM,
                       double te_NCEM,    double* ye_NCEM,
                       double ve_NCEM[3], double ve_NCSEM[3],
                       bool isFirst, int comp_orb_eml,
                       int comp_orb_seml,
                       ProjResClass& projRes, int k);


/**
 *   \brief Storing the results of the W + QBCP + JPL with average amplitudes, in txt file
 **/
void write_jplref_conn_txt(string filename,
                           Orbit& orbit_EM, Orbit& orbit_SEM,
                           AvgSt &avgSt_QBCP, AvgSt &avgSt_JPL,
                           double te_NCEM, double *ye_NCEM,
                           ProjResClass& projRes,
                           int isFirst,  int index);


#endif // SINGLE_ORBIT_IO_H_INCLUDED
