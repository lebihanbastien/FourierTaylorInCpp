#ifndef SINGLE_ORBIT_IO_H_INCLUDED
#define SINGLE_ORBIT_IO_H_INCLUDED

//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <vector>

#include "single_orbit.h"


using namespace std;

//========================================================================================
//
//          I/O (1)
//
//========================================================================================
/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type);

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type, double t0);

//-----------------------------------------------
// CU
//-----------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double ****yNCE, double ****sNCE, double *tGrid, int s1_grid_size, int s3_grid_size, int t_grid_size, int ofts_order, int type);

/**
 *  \brief Read in a data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin(double ****yNCE, double ****sNCE, double *tGrid, int s1_grid_size, int s3_grid_size, int t_grid_size, int ofts_order, int type);

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int getLenghtCU_bin(int *s1_grid_size, int *s3_grid_size, int *t_grid_size, int ofts_order, int type);

//-----------------------------------------------
// CU 3D
//-----------------------------------------------

//----------
// IN
//----------
/**
 *  \brief init the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1)
 **/
int initCU_bin_3D(int *si_grid_size, int t_grid_size, int ofts_order, int type);

/**
 *  \brief Write the current time in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only the time is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int appTimeCU_bin_3D(double *tGrid, int nt, int ofts_order, int type);

/**
 *  \brief Write the current state in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only a single loop on s4 is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int writeCU_bin_3D(double **yNCE, double **sNCE, int *si_grid_size, int ofts_order, int type);

//----------
// OUT
//----------
/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 **/
int getLenghtCU_bin_3D(int *si_grid_size, int *t_grid_size, int ofts_order, int type);

/**
 *  \brief Read in a data file the time at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readTCU_bin_3D(int offset, double *tGrid, int nt, int ofts_order, int type);

/**
 *  \brief Read in a data file the state at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin_3D(int offset, double **yNCE, double **sNCE, int *si_grid_size, int ofts_order, int type);

//-----------------------------------------------
// Int CU
//-----------------------------------------------
/**
 *  \brief Get the length of the data file containing the best connections between EML2 and SEML1,2.
 *         Used in int_sorted_sol_CMU_EM_to_CM_SEM/ref_CMU_EM_to_CM_SEM_MSD
 **/
int getLengthIntSortedCU_bin(int *number_of_sol, int ofts_order, int type);


//-----------------------------------------------
// Int CU
//-----------------------------------------------
/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM.
 **/
void writeIntProjCU_bin(string filename,
                        double *init_time_grid_EM,            //time grid in NCEM units
                        double ****init_state_CMU_NCEM,       //initial state in NCEM coordinates
                        double ****init_state_CMU_SEM,        //initial state in SEM coordinates
                        double ****init_state_CMU_RCM,        //initial state in RCM coordinates
                        double ****final_state_CMU_SEM,       //final state in SEM coordinates
                        double ****projected_state_CMU_SEM,   //projected state in SEM coordinates
                        double ****projected_state_CMU_RCM,   //projected state in RCM coordinates
                        double min_proj_dist_SEM,             //minimum distance of projection in SEM units
                        double dv_at_projection_SEM,          //associated dv
                        double *t_man_SEM,                    //time grid on manifold leg in SEM units
                        int kmin,
                        int ks1,
                        int ks3,
                        int kt);

/**
 *  \brief Read in a data file the connections between EML2 and SEML1,2.
 **/
void readIntProjCU_bin(string filename,
                       vector<double> &t0_CMU_EM,
                       vector<double> &tf_man_EM,
                       vector<double> &s1_CMU_EM,
                       vector<double> &s2_CMU_EM,
                       vector<double> &s3_CMU_EM,
                       vector<double> &s4_CMU_EM,
                       vector<double> &s5_CMU_EM,
                       vector<double> &pmin_dist_SEM,
                       vector<double> &s1_CM_SEM,
                       vector<double> &s2_CM_SEM,
                       vector<double> &s3_CM_SEM,
                       vector<double> &s4_CM_SEM,
                       vector<size_t> &sortId);

/**
 *  \brief Store in a data file the best connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM.
 **/
void writeIntProjSortCU_bin(string filename,
                            double ****init_state_CMU_NCEM,       //initial state in NCEM coordinates
                            double ****init_state_CMU_RCM,        //initial state in RCM coordinates
                            double ****final_state_CMU_SEM,       //final state in SEM coordinates
                            double ****projected_state_CMU_SEM,   //projected state in SEM coordinates
                            double ****projected_state_CMU_RCM,   //projected state in RCM coordinates
                            double ***min_proj_dist_tens_SEM,     //minimum distance of projection in SEM coordinates
                            vector<size_t> &sortId, vector<int> &ktMin,
                            vector<int> &ks1Min, vector<int> &ks3Min,
                            vector<double> &t0_min_EM, vector<double> &tf_min_EM,
                            vector<double> &distMin, int number_of_sol);

/**
 *  \brief Read from a data file the best connections between EML2 and SEML1,2.
 *         Used in int_sorted_sol_CMU_EM_to_CM_SEM. The data file must have been build via writeIntProjSortCU_bin.
 **/
void readIntProjSortCU_bin(string filename,
                           double *label,
                           double *t0_EM,
                           double *tf_EM,
                           double *s1_CMU_EM,
                           double *s3_CMU_EM,
                           double *s1_CM_SEM,
                           double *s3_CM_SEM,
                           double **init_state_CMU_NCEM,
                           double **final_state_CMU_SEM,
                           double **projected_state_CMU_SEM,
                           double *min_proj_dist_SEM_1,
                           double *min_proj_dist_SEM_2,
                           int number_of_sol);

//-----------------------------------------------
// Int CU 3D
//-----------------------------------------------
/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM_3D.
 **/
void writeIntProjCU_bin_3D(string filename,
                           double *init_time_grid_EM,          //time grid in NCEM units
                           double **init_state_CMU_NCEM,       //initial state in NCEM coordinates
                           double **init_state_CMU_SEM,        //initial state in SEM coordinates
                           double **init_state_CMU_RCM,        //initial state in RCM coordinates
                           double **final_state_CMU_SEM,       //final state in SEM coordinates
                           double **projected_state_CMU_SEM,   //projected state in SEM coordinates
                           double **projected_state_CMU_RCM,   //projected state in RCM coordinates
                           double min_proj_dist_SEM,           //minimum distance of projection in SEM units
                           double dv_at_projection_SEM,        //associated dv
                           double *t_man_SEM,                  //time grid on manifold leg in SEM units
                           int kmin,
                           int ks4,
                           int kt);
#endif // SINGLE_ORBIT_IO_H_INCLUDED
