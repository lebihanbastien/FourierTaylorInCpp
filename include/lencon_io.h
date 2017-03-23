#ifndef SINGLE_ORBIT_IO_H_INCLUDED
#define SINGLE_ORBIT_IO_H_INCLUDED

#include <iostream>     // std::cout
#include <algorithm>    // std::unique, std::distance
#include <vector>       // std::vector
#include <iterator>

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>

#include "single_orbit.h"

//For selection of the times
#define TIME_SELECTION_ABSOLUTE 1
#define TIME_SELECTION_RATIO    0

//Number of digits in the ratios
#define RATIO_DIGITS 7

using namespace std;

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
// Function: fileExists
//========================================================================================
/**
 *   Check if a file exists
 *   @param[in] filename - the name of the file to check
 *   @return    true if the file exists, else false
 *
 */
bool fileExists(const std::string& filename);

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
void writeCOMP_txt(double *t_traj_n, double **y_traj_n, int final_index);

/**
 *  \brief Reads the final_index+1 points trajectory contained in a data file of type
 *         TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector) and y_traj_n
 *         (state vectors). The size of the data file is given by getLengthCOMP_txt(), and
 *         the data vectors should be initialized accordingly by the user prior to the use
 *         of this routine.
 **/
int readCOMP_txt(double *t_traj_n, double **y_traj_n, int final_index);

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
int writeCU_bin_dH(double*** yNCE, double*** sNCE, double** dH, double* tGrid,
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
int readCU_bin_dH(double*** yNCE, double*** sNCE, double** dH, double* tGrid,
                  int s1_grid_size, int t_grid_size,
                  int ofts_order, int type, int destination, double dHd);

//----------------------------------------------------------------------------------------
// CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double**** yNCE, double**** sNCE, double*** dH, double* tGrid,
                int s1_grid_size, int s3_grid_size, int t_grid_size,
                int ofts_order, int type, int destination);
/**
 *  \brief Read in a data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin(double**** yNCE, double**** sNCE, double*** dH, double* tGrid,
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
//          I/O (Refinement)
//
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
vector<size_t> sort_indexes(const vector<double> &v);

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
 *  \brief Read in a data file the connections between EML2 and SEML1,2.
 *         Find the data that are the closest to the desired t0 at EML2 departures.
 **/
int readClosestIntProjCU_bin(string filename, double t0_des,
                              vector<double>& t0_CMU_EM_0, vector<double>& tf_man_EM_0,
                              vector<double>& s1_CMU_EM_0, vector<double>& s2_CMU_EM_0,
                              vector<double>& s3_CMU_EM_0, vector<double>& s4_CMU_EM_0,
                              vector<double>& s5_CMU_EM_0, vector<double>& pmin_dist_SEM_0,
                              vector<double>& s1_CM_SEM_0, vector<double>& s2_CM_SEM_0,
                              vector<double>& s3_CM_SEM_0, vector<double>& s4_CM_SEM_0,
                              vector<double>& crossings_NCSEM_0,
                              vector<size_t>& sortId, int typeOfSelection);



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
