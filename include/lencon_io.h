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

using namespace std;


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
    int mplot;            //number of points per plot between to pach points (e.g. total plot points is gridSize*mplot)

    int time;             //type of constraints on the times in REF_CONT
    int grid;             //type of grid
    int termination;      //termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
    int coord_type;       //coordinates system in the refinement procedure (usually NCSEM)

    // Maximum/Minimum step in the continuation procedure
    double dsmin;           //with fixed time
    double dsmin_vt;        //with variable time
    double dsmax;           //with fixed time
    double dsmax_vt;        //with variable time

    double xps;           //position of the poincaré section in NCSEM coordinates
    int isJPL;            //is the JPL refinement performed when possible?
    int djplcoord;        //coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with $NJ2000
    int sidim;            //0 or 2 - component of s0 that stays constant when t0 is free.

    //Sampling frequencies in REF_COMP (complete trajectory) in days
    int sf_eml2;          // orbit at EML2
    int sf_man;           // transfer leg
    int sf_seml2;         // orbit at SEML2

    // Integration window for each orbit
    double tspan_EM;
    double tspan_SEM;

    // Storing the orbits at each step?
    int isSaved_EM;    //0: don't save, 1: save using projection method
    int isSaved_SEM;   //0: don't save, 1: save using projection method, 2: save using integration in reduced coordinates

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

    //------------------------------------------------------------------------------------
    //Additional variables for specific differential correctors:
    // - goodvector is used to add the pseudo-arclength constraint
    // - dH is used to add a constraint on the (initial) energy
    // - pkpos is used to add a constraint on a certain Pk section
    //------------------------------------------------------------------------------------
    double goodvector[1000]; //vector of free variables
    double dH;               //energy value at the origin
    int pkpos;               //indix of the patch point that bear the poincaré section

    //------------------------------------------------------------------------------------
    //Constructor
    //------------------------------------------------------------------------------------
    //RefSt():pmax_dist_SEM(1e-3) {};
    RefSt():isCollisionOn(1), pmax_dist_SEM(1e5), last_error(0.0), inner_prec(5e-8)
    {
        for(int i = 0; i <4; i++) si_SEED_EM_MIN[i] = -50;
        for(int i = 0; i <4; i++) si_SEED_EM_MAX[i] = +50;
    };
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
    //Inputs
    double init_time_EM, init_state_CMU_NCEM[6], init_state_CMU_RCM[5];
    int label;

    //Seeds
    double seed_time_EM;
    double seed_state_CMU_RCM[4];

    //Outputs
    double init_state_CMU_SEM_o[6], final_state_CMU_SEM_o[6];
    double projected_state_CMU_SEM_o[6], projected_state_CMU_RCM_o[4];
    double min_proj_dist_SEM_o, dv_at_projection_SEM_o;
    double final_time_SEM_o;
    double crossings_NCSEM_o;
    int collision_NCEM_o;

};

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
     *         After this routine, sortId[0] is the indix for which pmin_dist_SEM is minimum.
     **/
    void sort_pmin_dist_SEM();

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
//          I/O (filename routines)
//
//========================================================================================
/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type, int destination);

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type, int destination, double t0);

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, energy, and type of the data.
 **/
string filenameCUM_dH(int ofts_order, int type, int destination, double dH);

//========================================================================================
//
//          I/O (complete trajectories)
//
//========================================================================================
/**
 *  \brief Stores the final_index+1 points trajectory contained in t_traj_n (time vector)
 *         and y_traj_n (state vectors) into a data file of type TYPE_COMP_FOR_JPL.
 **/
void writeCOMP_txt(double* t_traj_n, double** y_traj_n, int final_index);

/**
 *  \brief Reads the final_index+1 points trajectory contained in a data file of type
 *         TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector) and y_traj_n
 *         (state vectors). The size of the data file is given by getLengthCOMP_txt(), and
 *         the data vectors should be initialized accordingly by the user prior to the use
 *         of this routine.
 **/
int readCOMP_txt(double* t_traj_n, double** y_traj_n, int final_index);

/**
 *  \brief Get the length of the final_index+1 points trajectory contained in a
 *         data file of type TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector)
 *         and y_traj_n (state vectors).
 **/
int getLengthCOMP_txt();

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
                   int s1_grid_size, int t_grid_size, int ofts_order,
                   int type, int destination, double dHd);

/**
 *  \brief Get the length of the data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1 and of size (tGrid +1)*(gSize+1)
 **/
int getLenghtCU_bin_dH(int* s1_grid_size, int* t_grid_size, int ofts_order,
                       int type, int destination, double dHd);

/**
 *  \brief Read in a data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1 and of size (tGrid +1)*(gSize+1)
 **/
int readCU_bin_dH(double** * yNCE, double** * sNCE, double** dH, double* tGrid,
                  int s1_grid_size, int t_grid_size,
                  int ofts_order, int type, int destination, double dHd);

//----------------------------------------------------------------------------------------
// CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
                int s1_grid_size, int s3_grid_size, int t_grid_size,
                int ofts_order, int type, int destination);
/**
 *  \brief Read in a data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
               int s1_grid_size, int s3_grid_size, int t_grid_size,
               int ofts_order, int type, int destination);

/**
 *  \brief Get the length of the data file the containing the Initial
 *         Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int getLenghtCU_bin(int* s1_grid_size, int* s3_grid_size,
                    int* t_grid_size, int ofts_order,
                    int type, int destination);



//----------------------------------------------------------------------------------------
// CU 3D
//----------------------------------------------------------------------------------------

/**
 *  \brief init the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1)
 **/
int initCU_bin_3D(int* si_grid_size, int t_grid_size,
                  int ofts_order, int type, int destination);
/**
 *  \brief Write the current time in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only the time is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int appTimeCU_bin_3D(double* tGrid, int nt, int ofts_order, int type, int destination);

/**
 *  \brief Write the current state in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only a single loop on s4 is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int writeCU_bin_3D(double** yNCE, double** sNCE, int* si_grid_size,
                   int ofts_order, int type, int destination);

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 **/
int getLenghtCU_bin_3D(int* si_grid_size, int* t_grid_size,
                       int ofts_order, int type, int destination);

/**
 *  \brief Read in a data file the time at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readTCU_bin_3D(int offset, double* tGrid, int nt,
                   int ofts_order, int type, int destination);

/**
 *  \brief Read in a data file the state at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin_3D(int offset, double** yNCE, double** sNCE, int* si_grid_size,
                  int ofts_order, int type, int destination);

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
/**
 *   \brief Display the current completion (percent) of a routine.
 **/
void displayCompletion(string funcname, double percent);

#endif // SINGLE_ORBIT_IO_H_INCLUDED
