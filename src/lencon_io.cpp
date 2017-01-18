#include "lencon_io.h"
#include <iostream>     // std::cout
#include <algorithm>    // std::unique, std::distance
#include <vector>       // std::vector
#include <iterator>

//========================================================================================
//
//          I/O (filename routines)
//
//========================================================================================
/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type, int destination)
{
    switch(type)
    {
    case TYPE_CU:
        return SEML.cs_em.F_PLOT+"cu_order_"+numTostring(ofts_order)+".bin";
    case TYPE_MAN:
        return SEML.cs_em.F_PLOT+"intcu_order_"+numTostring(ofts_order)+".bin";
    case TYPE_MAN_PROJ:
        return SEML.cs_em.F_PLOT+"projcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_MAN_SORT:
        return SEML.cs_em.F_PLOT+"sortprojcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_MAN_SORT_IN:
        return SEML.cs_em.F_PLOT+"sortprojintcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_CONT_ATF:
        return SEML.cs_em.F_PLOT+"cont_atf_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_CU_3D:
        return SEML.cs_em.F_PLOT+"cu_3d_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_MAN_PROJ_3D:
        return SEML.cs_em.F_PLOT+"projcu_3d_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".bin";
    case TYPE_MAN_PROJ_FROM_SERVER:
        return SEML.cs_em.F_PLOT+"Serv/projcu_order_20.bin";
    case TYPE_MAN_PROJ_3D_FROM_SERVER:
        return SEML.cs_em.F_PLOT+"Serv/projcu_3d_order_20.bin";
    case TYPE_COMP_FOR_JPL:
        return SEML.cs_em.F_PLOT+"comp_for_jpl_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+".txt";
    default:
        cout << "filenameOrbit: unknown type." << endl;
        return "";
    }
}

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int ofts_order, int type, int destination, double t0)
{
    switch(type)
    {
    case TYPE_CU:
        return SEML.cs_em.F_PLOT+"cu_order_"+numTostring(ofts_order)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_MAN:
        return SEML.cs_em.F_PLOT+"intcu_order_"+numTostring(ofts_order)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_MAN_PROJ:
        return SEML.cs_em.F_PLOT+"projcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_MAN_SORT:
        return SEML.cs_em.F_PLOT+"sortprojcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_MAN_SORT_IN:
        return SEML.cs_em.F_PLOT+"sortprojintcu_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_CONT_ATF:
        return SEML.cs_em.F_PLOT+"cont_atf_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".txt";
    case TYPE_CONT_ATF_TRAJ:
        return SEML.cs_em.F_PLOT+"cont_atf_traj_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    case TYPE_CU_3D:
        return SEML.cs_em.F_PLOT+"cu_3d_order_"+numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    default:
        cout << "filenameOrbit: unknown type." << endl;
        return "";
    }
}


//========================================================================================
//
//          I/O (complete trajectories)
//
//========================================================================================
/**
 *  \brief Stores the final_index+1 points trajectory contained in t_traj_n (time vector)
 *         and y_traj_n (state vectors) into a data file of type TYPE_COMP_FOR_JPL.
 **/
void writeCOMP_txt(double *t_traj_n, double **y_traj_n, int final_index)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_COMP_FOR_JPL, SEML.li_SEM);
    fstream filestream;
    filestream.open (filename.c_str(), ios::out);
    filestream << setprecision(15) <<  setiosflags(ios::scientific) << std::showpos;

    //====================================================================================
    // First value is final_index
    //====================================================================================
    filestream << final_index << endl;

    //====================================================================================
    // Store the data in the format "t x y z px py pz" on each line
    //====================================================================================
    for(int k = 0; k <= final_index; k++)
    {
        filestream << t_traj_n[k] << "  ";
        for(int i = 0; i < 6; i++) filestream << y_traj_n[i][k] << "  ";
        filestream << endl;
    }
    filestream.close();
}

/**
 *  \brief Reads the final_index+1 points trajectory contained in a data file of type
 *         TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector) and y_traj_n
 *         (state vectors). The size of the data file is given by getLengthCOMP_txt(), and
 *         the data vectors should be initialized accordingly by the user prior to the use
 *         of this routine.
 **/
void readCOMP_txt(double *t_traj_n, double **y_traj_n, int final_index)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_COMP_FOR_JPL, SEML.li_SEM);
    fstream filestream;
    filestream.open (filename.c_str(), ios::in);

    //====================================================================================
    // First value is discarded
    //====================================================================================
    string ct;
    getline(filestream,  ct);


    //====================================================================================
    // Read the data in the format "t x y z px py pz" on each line
    //====================================================================================
    double tempd = 0.0;
    for(int k = 0; k <= final_index; k++)
    {
        filestream >> tempd;
        t_traj_n[k] = tempd;
        for(int i = 0; i < 6; i++)
        {
            filestream >> tempd;
            y_traj_n[i][k]=  tempd;
        }
    }
    filestream.close();
}

/**
 *  \brief Get the length of the final_index+1 points trajectory contained in a
 *         data file of type TYPE_COMP_FOR_JPL and stores it in t_traj_n (time vector)
 *         and y_traj_n (state vectors).
 **/
int getLengthCOMP_txt()
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_COMP_FOR_JPL, SEML.li_SEM);
    fstream filestream;
    filestream.open (filename.c_str(), ios::in);

    //====================================================================================
    // First value is final_index
    //====================================================================================
    int final_index;
    filestream >> final_index;
    filestream.close();
    return final_index;
}






//========================================================================================
//
//          I/O (CU/CS/CM)
//
//========================================================================================
//-----------------------------------------------
// CU
//-----------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double**** yNCE, double**** sNCE, double* tGrid,
                int s1_grid_size, int s3_grid_size, int t_grid_size,
                int ofts_order, int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of data on the time grid
        //---------------------
        resi = t_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        resi = s1_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        resi = s3_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int nt = 0; nt <= t_grid_size; nt++)
        {
            //Store the current time
            res = tGrid[nt];
            filestream.write((char*) &res, sizeof(double));

            //Store the data at current time
            for(int n1 = 0; n1 <= s1_grid_size; n1++)
            {
                for (int n2 = 0; n2 <= s3_grid_size; n2++)
                {
                    //NC state
                    for (int k = 0; k < 6; k++)
                    {
                        res = yNCE[k][nt][n1][n2];
                        filestream.write((char*) &res, sizeof(double));
                    }

                    //RCM state
                    for (int k = 0; k < 5; k++)
                    {
                        res = sNCE[k][nt][n1][n2];
                        filestream.write((char*) &res, sizeof(double));
                    }
                }
            }
        }
        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

/**
 *  \brief Read in a data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin(double**** yNCE, double**** sNCE, double* tGrid,
               int s1_grid_size, int s3_grid_size, int t_grid_size,
               int ofts_order, int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        int tSize0, gSize0_s1, gSize0_s3;

        //---------------------
        //Number of data on the time grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        tSize0 = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        gSize0_s1 = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        gSize0_s3 = resi;


        if(tSize0 < t_grid_size || gSize0_s1 < s1_grid_size || gSize0_s3 < s3_grid_size)
        {
            cout << "readManifold_bin: wrong inputs" << endl;
            cout << "t_grid_size    = " << t_grid_size    << ", but tSize0 = " << tSize0 << endl;
            cout << "s1_grid_size = " << s1_grid_size << ", but gSize0_s1 = " << gSize0_s1 << endl;
            cout << "s3_grid_size = " << s3_grid_size << ", but gSize0_s3 = " << gSize0_s3 << endl;
            return FTC_FAILURE;
        }

        //---------------------
        //Loop
        //---------------------
        double res;
        for(int nt = 0; nt <= t_grid_size; nt++)
        {
            //Read the current time
            filestream.read((char*) &res, sizeof(double));
            tGrid[nt]= res;

            //Read the data at current time
            for(int n1 = 0; n1 <= s1_grid_size; n1++)
            {
                for (int n2 = 0; n2 <= s3_grid_size; n2++)
                {
                    //NC state
                    for (int k = 0; k < 6; k++)
                    {
                        filestream.read((char*) &res, sizeof(double));
                        yNCE[k][nt][n1][n2] = res;
                    }

                    //RCM state
                    for (int k = 0; k < 5; k++)
                    {
                        filestream.read((char*) &res, sizeof(double));
                        sNCE[k][nt][n1][n2] = res;
                    }
                }
            }
        }
        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int getLenghtCU_bin(int* s1_grid_size, int* s3_grid_size,
                    int* t_grid_size, int ofts_order,
                    int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        //---------------------
        //Number of data on the time grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *t_grid_size = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *s1_grid_size = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *s3_grid_size = resi;

        filestream.close();

    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}


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
int initCU_bin_3D(int* si_grid_size, int t_grid_size,
                  int ofts_order, int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        int resi;
        //---------------------
        //Number of data on the time grid
        //---------------------
        resi = t_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Number of data on the manifold grid on all four dimensions
        //---------------------
        for(int i = 0; i < 4; i++)
        {
            resi = si_grid_size[i];
            filestream.write((char*) &resi, sizeof(int));
        }

        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

/**
 *  \brief Write the current time in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only the time is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int appTimeCU_bin_3D(double* tGrid, int nt, int ofts_order, int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out | ios::app);
    if (filestream.is_open())
    {
        //Store the current time
        double res = tGrid[nt];
        filestream.write((char*) &res, sizeof(double));
        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

/**
 *  \brief Write the current state in the data file of Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 *         Here, only a single loop on s4 is appended to the data file, so this routine must be used within the intricated loops. See compute_grid_CMU_EM_3D src code for details.
 **/
int writeCU_bin_3D(double** yNCE, double** sNCE, int* si_grid_size,
                   int ofts_order, int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out | ios::app);
    if (filestream.is_open())
    {
        double res;
        //---------------------
        //Loop
        //---------------------
        for(int n3 = 0; n3 <= si_grid_size[2]; n3++)
        {
            //NC state
            for (int k = 0; k < 6; k++)
            {
                res = yNCE[k][n3];
                filestream.write((char*) &res, sizeof(double));
            }

            //RCM state
            for (int k = 0; k < 5; k++)
            {
                res = sNCE[k][n3];
                filestream.write((char*) &res, sizeof(double));
            }
        }

        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

//----------
// OUT
//----------

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 **/
int getLenghtCU_bin_3D(int* si_grid_size, int* t_grid_size,
                       int ofts_order, int type, int destination)
{
    //Offset
    int offset = 0;

    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;

        //---------------------
        //Number of data on the time grid
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *t_grid_size = resi;

        //---------------------
        //Number of data on the manifold grid on all four dimensions
        //---------------------
        for(int i = 0; i < 4; i++)
        {
            filestream.read((char*) &resi, sizeof(int));
            si_grid_size[i] = resi;
        }

        //---------------------
        //Get the offset
        //---------------------
        offset = filestream.tellg();

        //---------------------
        //Close
        //---------------------
        filestream.close();

    }
    else return FTC_FAILURE;

    return offset;
}

/**
 *  \brief Read in a data file the time at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readTCU_bin_3D(int offset, double* tGrid, int nt,
                   int ofts_order, int type, int destination)
{
    //Offset
    int offset2 = 0;

    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        //---------------------
        //Set the offset
        //---------------------
        filestream.seekg (offset, ios::beg);

        //---------------------
        //Read the current time
        //---------------------
        double res;
        filestream.read((char*) &res, sizeof(double));
        tGrid[nt] = res;

        //---------------------
        //Get the offset
        //---------------------
        offset2 = filestream.tellg();

        //---------------------
        //Close
        //---------------------
        filestream.close();

    }
    else return FTC_FAILURE;

    return offset2;
}

/**
 *  \brief Read in a data file the state at Initial Conditions of a planar Center-Unstable manifold. Used in int_proj_CMU_EM_on_CM_SEM and intMan
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int readCU_bin_3D(int offset, double** yNCE, double** sNCE, int* si_grid_size,
                  int ofts_order, int type, int destination)
{
    //Offset
    int offset2 = 0;

    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        //---------------------
        //Set the offset
        //---------------------
        filestream.seekg (offset, ios::beg);

        //---------------------
        //Loop
        //---------------------
        double res;
        for(int n3 = 0; n3 <= si_grid_size[2]; n3++)
        {
            //NC state
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
                yNCE[k][n3] = res;
            }

            //RCM state
            for (int k = 0; k < 5; k++)
            {
                filestream.read((char*) &res, sizeof(double));
                sNCE[k][n3] = res;
            }
        }

        //---------------------
        //Get the offset
        //---------------------
        offset2 = filestream.tellg();

        //---------------------
        //Close
        //---------------------
        filestream.close();

    }
    else return FTC_FAILURE;

    return offset2;
}




//-----------------------------------------------
// Int CU
//-----------------------------------------------
/**
 *  \brief Get the length of the data file containing the best connections between EML2 and SEML1,2.
 *         Used in int_sorted_sol_CMU_EM_to_CM_SEM/ref_CMU_EM_to_CM_SEM_MSD
 **/
int getLengthIntSortedCU_bin(int* number_of_sol, int ofts_order,
                             int type, int destination)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        //---------------------
        //Number of data
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *number_of_sol = resi;
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}


/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM.
 **/
void writeIntProjCU_bin(string filename,
                        double* init_time_grid_EM,            //time grid in NCEM units
                        double**** init_state_CMU_NCEM,       //initial state in NCEM coordinates
                        double**** init_state_CMU_SEM,        //initial state in SEM coordinates
                        double**** init_state_CMU_RCM,        //initial state in RCM coordinates
                        double**** final_state_CMU_SEM,       //final state in SEM coordinates
                        double**** projected_state_CMU_SEM,   //projected state in SEM coordinates
                        double**** projected_state_CMU_RCM,   //projected state in RCM coordinates
                        double min_proj_dist_SEM,             //minimum distance of projection in SEM units
                        double dv_at_projection_SEM,          //associated dv
                        double* t_man_SEM,                    //time grid on manifold leg in SEM units
                        int kmin,
                        int ks1,
                        int ks3,
                        int kt)
{
    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out | ios::app);
    if (filestream.is_open())
    {
        double res;
        // 1. time grid in NCEM units
        res = init_time_grid_EM[kt];
        filestream.write((char*) &res, sizeof(double));

        // 2-7. initial state in NCEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = init_state_CMU_NCEM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 8-13. initial state in SEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = init_state_CMU_SEM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 14-18. initial state in RCM coordinates
        for (int k = 0; k < 5; k++)
        {
            res = init_state_CMU_RCM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 19. minimum distance of projection
        res = min_proj_dist_SEM;
        filestream.write((char*) &res, sizeof(double));

        // 20. associated dv
        res = dv_at_projection_SEM;
        filestream.write((char*) &res, sizeof(double));

        // 21. tvMinTensor
        res = t_man_SEM[kmin];
        filestream.write((char*) &res, sizeof(double));

        // 22-27. final_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = final_state_CMU_SEM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 28-33. projected_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projected_state_CMU_SEM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 34-37. projected_state_CMU_RCM state in SE coordinates
        for (int k = 0; k < 4; k++)
        {
            res = projected_state_CMU_RCM[k][kt][ks1][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        filestream.close();
    }
}

//-----------------------------------------------
// Int CU 3D
//-----------------------------------------------
/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM_3D.
 **/
void writeIntProjCU_bin_3D(string filename,
                           double* init_time_grid_EM,          //time grid in NCEM units
                           double** init_state_CMU_NCEM,       //initial state in NCEM coordinates
                           double** init_state_CMU_SEM,        //initial state in SEM coordinates
                           double** init_state_CMU_RCM,        //initial state in RCM coordinates
                           double** final_state_CMU_SEM,       //final state in SEM coordinates
                           double** projected_state_CMU_SEM,   //projected state in SEM coordinates
                           double** projected_state_CMU_RCM,   //projected state in RCM coordinates
                           double min_proj_dist_SEM,           //minimum distance of projection in SEM units
                           double dv_at_projection_SEM,        //associated dv
                           double* t_man_SEM,                  //time grid on manifold leg in SEM units
                           int kmin,
                           int ks3,
                           int kt)
{
    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out | ios::app);
    if (filestream.is_open())
    {
        double res;
        // 1. time grid in NCEM units
        res = init_time_grid_EM[kt];
        filestream.write((char*) &res, sizeof(double));

        // 2-7. initial state in NCEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = init_state_CMU_NCEM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 8-13. initial state in SEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = init_state_CMU_SEM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 14-18. initial state in RCM coordinates
        for (int k = 0; k < 5; k++)
        {
            res = init_state_CMU_RCM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 19. minimum distance of projection
        res = min_proj_dist_SEM;
        filestream.write((char*) &res, sizeof(double));

        // 20. associated dv
        res = dv_at_projection_SEM;
        filestream.write((char*) &res, sizeof(double));

        // 21. tvMinTensor
        res = t_man_SEM[kmin];
        filestream.write((char*) &res, sizeof(double));

        // 22-27. final_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = final_state_CMU_SEM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 28-33. projected_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projected_state_CMU_SEM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        // 34-37. projected_state_CMU_RCM state in SE coordinates
        for (int k = 0; k < 4; k++)
        {
            res = projected_state_CMU_RCM[k][ks3];
            filestream.write((char*) &res, sizeof(double));
        }

        filestream.close();
    }
}


void vector_getUnique(vector<double>& v0U)
{
    //Creating iterator on v0U
    std::vector<double>::iterator it;
    it = std::unique (v0U.begin(), v0U.end());

    //Resize
    v0U.resize(std::distance(v0U.begin(),it) );

    //Print out content:
    //std::cout << "v0U contains:";
    //for (it=v0U.begin(); it!=v0U.end(); ++it)
    //    std::cout << ' ' << *it;
    //std::cout << '\n';
}

void vector_getIndices(vector<size_t>& indRes, vector<double>& t0_CMU_EM, double t0)
{
    //Create an iterator that will contain the position of each component of t0_CMU_EM that matches t0
    //Here is the first one
    std::vector<double>::iterator itf = std::find(t0_CMU_EM.begin(), t0_CMU_EM.end(), t0);

    // Loop on all t0_CMU_EM, update of indRes
    while (itf != t0_CMU_EM.end())
    {
        indRes.push_back(std::distance(t0_CMU_EM.begin(), itf));
        std::advance(itf, 1);
        itf = std::find(itf, t0_CMU_EM.end(),  t0);
    }
}

/**
 *  \brief Read in a data file the connections between EML2 and SEML1,2.
 **/
void readIntProjCU_bin(string filename,
                       vector<double>& t0_CMU_EM_0,
                       vector<double>& tf_man_EM_0,
                       vector<double>& s1_CMU_EM_0,
                       vector<double>& s2_CMU_EM_0,
                       vector<double>& s3_CMU_EM_0,
                       vector<double>& s4_CMU_EM_0,
                       vector<double>& s5_CMU_EM_0,
                       vector<double>& pmin_dist_SEM_0,
                       vector<double>& s1_CM_SEM_0,
                       vector<double>& s2_CM_SEM_0,
                       vector<double>& s3_CM_SEM_0,
                       vector<double>& s4_CM_SEM_0,
                       vector<size_t>& sortId)
{
    //==========================================================
    //Temporary variables
    //==========================================================
    vector<double> t0_CMU_EM;
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


    //==========================================================
    //Open and read datafile
    //==========================================================
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        double res;

        do
        {
            // 1. time grid in NCEM units
            filestream.read((char*) &res, sizeof(double));
            t0_CMU_EM.push_back(res);

            // 2-7. initial state in NCEM coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));

            }

            // 8-13. initial state in SEM coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 14-18. initial state in RCM coordinates: SAVED
            filestream.read((char*) &res, sizeof(double));
            s1_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s2_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s3_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s4_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s5_CMU_EM.push_back(res);

            // 19. minimum distance of projection: SAVED
            filestream.read((char*) &res, sizeof(double));
            pmin_dist_SEM.push_back(res);

            // 20. associated dv: NOT SAVED
            filestream.read((char*) &res, sizeof(double));

            // 21. tf at SEM: SAVED
            filestream.read((char*) &res, sizeof(double));
            tf_man_SEM.push_back(res);


            // 22-27. final_state_CMU_SEM state in SE coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 28-33. projected_state_CMU_SEM state in SE coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 34-37. projected_state_CMU_RCM state in SE coordinates: SAVED
            filestream.read((char*) &res, sizeof(double));
            s1_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s2_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s3_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s4_CM_SEM.push_back(res);
        }
        while(!filestream.eof());

        filestream.close();
    }


    //==========================================================
    //Delete last element that is not a real value
    //==========================================================
    t0_CMU_EM.pop_back();
    s1_CMU_EM.pop_back();
    s2_CMU_EM.pop_back();
    s3_CMU_EM.pop_back();
    s4_CMU_EM.pop_back();
    s5_CMU_EM.pop_back();
    pmin_dist_SEM.pop_back();
    tf_man_SEM.pop_back();
    s1_CM_SEM.pop_back();
    s2_CM_SEM.pop_back();
    s3_CM_SEM.pop_back();
    s4_CM_SEM.pop_back();


    //==========================================================
    //Get the unique elements in t0_CMU_EM
    //==========================================================
    //Copy t0_CMU_EM into t0_CMU_EM_UNIQUE
    vector<double> t0_CMU_EM_UNIQUE(t0_CMU_EM);
    //Get unique elements
    vector_getUnique(t0_CMU_EM_UNIQUE);

    //==========================================================
    // Get all the indices that match t0_CMU_EM_UNIQUE[xxx]
    //==========================================================
    int ti = 0;
    cout << "--------------------------------------" << endl;
    cout << "There is " << t0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
    cout << "[" << t0_CMU_EM_UNIQUE[0] << ", " << t0_CMU_EM_UNIQUE[t0_CMU_EM_UNIQUE.size()-1] << "]" << endl;
    do
    {
        cout << "Please enter a number between " << 0 << " and " << t0_CMU_EM_UNIQUE.size()-1;
        cout << " to select a specific starting time" << endl;
        scanf("%d", &ti);
    }
    while(ti < 0 || ti > (int) t0_CMU_EM_UNIQUE.size()-1);
    cout << "You have selected " << ti << ", which corresponds to t0_EM = " << t0_CMU_EM_UNIQUE[ti] << endl;

    std::vector<size_t> indRes;
    vector_getIndices(indRes, t0_CMU_EM, t0_CMU_EM_UNIQUE[ti]);

    //==========================================================
    // Display info on the indices that match t0_CMU_EM_UNIQUE[xxx]
    //==========================================================
    coutmp();
    cout << "--------------------------------------" << endl;
    cout << "The first entry for this time is:" << endl;
    int ind = indRes[0];
    cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
    cout << "tf_EM    = " << tf_man_SEM[ind]/SEML.us_em.ns << endl;
    cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
    cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
    cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;

    cout << "--------------------------------------" << endl;
    cout << "The last entry for this time is:" << endl;
    ind = indRes[indRes.size()-1];
    cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
    cout << "tf_EM    = " << tf_man_SEM[ind]/SEML.us_em.ns << endl;
    cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
    cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
    cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;
    coutlp();

    //==========================================================
    // Copy the selected results in the inputs
    //==========================================================
    for(int ind = 0; ind < (int) indRes.size(); ind++)
    {
        //Times
        t0_CMU_EM_0.push_back(t0_CMU_EM[indRes[ind]]);
        tf_man_EM_0.push_back(tf_man_SEM[indRes[ind]]/SEML.us_em.ns);
        //CMU of EM
        s1_CMU_EM_0.push_back(s1_CMU_EM[indRes[ind]]);
        s2_CMU_EM_0.push_back(s2_CMU_EM[indRes[ind]]);
        s3_CMU_EM_0.push_back(s3_CMU_EM[indRes[ind]]);
        s4_CMU_EM_0.push_back(s4_CMU_EM[indRes[ind]]);
        s5_CMU_EM_0.push_back(s5_CMU_EM[indRes[ind]]);
        //Projection distance
        pmin_dist_SEM_0.push_back(pmin_dist_SEM[indRes[ind]]);
        //CM of SEM
        s1_CM_SEM_0.push_back(s1_CM_SEM[indRes[ind]]);
        s2_CM_SEM_0.push_back(s2_CM_SEM[indRes[ind]]);
        s3_CM_SEM_0.push_back(s3_CM_SEM[indRes[ind]]);
        s4_CM_SEM_0.push_back(s4_CM_SEM[indRes[ind]]);
    }

    //==========================================================
    // Sort data wrt to the projection distance
    //==========================================================
    sortId = sort_indexes(pmin_dist_SEM_0);

}



/**
 *  \brief Read in a data file the connections between EML2 and SEML1,2.
 *         Interpolate in the data set to get the right desired t0 at EML2 departures.
 **/
void readAndInterpolateIntProjCU_bin(string filename,
                                     double t0_des,
                                     vector<double>& t0_CMU_EM_0,
                                     vector<double>& tf_man_EM_0,
                                     vector<double>& s1_CMU_EM_0,
                                     vector<double>& s2_CMU_EM_0,
                                     vector<double>& s3_CMU_EM_0,
                                     vector<double>& s4_CMU_EM_0,
                                     vector<double>& s5_CMU_EM_0,
                                     vector<double>& pmin_dist_SEM_0,
                                     vector<double>& s1_CM_SEM_0,
                                     vector<double>& s2_CM_SEM_0,
                                     vector<double>& s3_CM_SEM_0,
                                     vector<double>& s4_CM_SEM_0,
                                     vector<size_t>& sortId)
{
    //==========================================================
    //Temporary variables
    //==========================================================
    vector<double> t0_CMU_EM;
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


    //==========================================================
    //Open and read datafile
    //==========================================================
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        double res;

        do
        {
            // 1. time grid in NCEM units
            filestream.read((char*) &res, sizeof(double));
            t0_CMU_EM.push_back(res);

            // 2-7. initial state in NCEM coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));

            }

            // 8-13. initial state in SEM coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 14-18. initial state in RCM coordinates: SAVED
            filestream.read((char*) &res, sizeof(double));
            s1_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s2_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s3_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s4_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s5_CMU_EM.push_back(res);

            // 19. minimum distance of projection: SAVED
            filestream.read((char*) &res, sizeof(double));
            pmin_dist_SEM.push_back(res);

            // 20. associated dv: NOT SAVED
            filestream.read((char*) &res, sizeof(double));

            // 21. tf at SEM: SAVED
            filestream.read((char*) &res, sizeof(double));
            tf_man_SEM.push_back(res);


            // 22-27. final_state_CMU_SEM state in SE coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 28-33. projected_state_CMU_SEM state in SE coordinates: NOT SAVED
            for (int k = 0; k < 6; k++)
            {
                filestream.read((char*) &res, sizeof(double));
            }

            // 34-37. projected_state_CMU_RCM state in SE coordinates: SAVED
            filestream.read((char*) &res, sizeof(double));
            s1_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s2_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s3_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            s4_CM_SEM.push_back(res);
        }
        while(!filestream.eof());

        filestream.close();
    }else{
        cout << "readAndInterpolateIntProjCU_bin. Unable to open the file. Check file name. " << endl;
        return;
    }


    //==========================================================
    //Delete last element that is not a real value
    //==========================================================
    t0_CMU_EM.pop_back();
    s1_CMU_EM.pop_back();
    s2_CMU_EM.pop_back();
    s3_CMU_EM.pop_back();
    s4_CMU_EM.pop_back();
    s5_CMU_EM.pop_back();
    pmin_dist_SEM.pop_back();
    tf_man_SEM.pop_back();
    s1_CM_SEM.pop_back();
    s2_CM_SEM.pop_back();
    s3_CM_SEM.pop_back();
    s4_CM_SEM.pop_back();


    //==========================================================
    //Get the unique elements in t0_CMU_EM
    //==========================================================
    //Copy t0_CMU_EM into t0_CMU_EM_UNIQUE
    vector<double> t0_CMU_EM_UNIQUE(t0_CMU_EM);
    //Get unique elements
    vector_getUnique(t0_CMU_EM_UNIQUE);

    //==========================================================
    // Print all the indices that match t0_CMU_EM_UNIQUE[xxx]
    //==========================================================
    cout << "--------------------------------------" << endl;
    cout << "There is " << t0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
    cout << "[" << t0_CMU_EM_UNIQUE[0]/SEML.us_em.T << ", " << t0_CMU_EM_UNIQUE[t0_CMU_EM_UNIQUE.size()-1]/SEML.us_em.T << "]x SEML.us_em.T" << endl;

    //==========================================================
    // Find the nearest t0 value
    //==========================================================
    double dmin = fabs(t0_CMU_EM_UNIQUE[0] - t0_des);
    int ti = 0;
    for(int i = 1; i < t0_CMU_EM_UNIQUE.size(); i++)
    {
        if(fabs(t0_CMU_EM_UNIQUE[i] - t0_des) < dmin)
        {
            dmin = fabs(t0_CMU_EM_UNIQUE[i] - t0_des);
            ti = i;
        }
    }
    cout << "The index of the time closer to the desired one is " << ti;
    cout << ", which corresponds to;" << endl;
    cout << "t0_EM = " << t0_CMU_EM_UNIQUE[ti]/SEML.us_em.T;
    cout << " x SEML.us_em.T" <<  endl;

    std::vector<size_t> indRes;
    vector_getIndices(indRes, t0_CMU_EM, t0_CMU_EM_UNIQUE[ti]);

    //==========================================================
    // Display info on the indices that match t0_CMU_EM_UNIQUE[xxx]
    //==========================================================
    coutmp();
    cout << "--------------------------------------" << endl;
    cout << "The first entry for this time is:" << endl;
    int ind = indRes[0];
    cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
    cout << "tf_EM    = " << tf_man_SEM[ind]/SEML.us_em.ns << endl;
    cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
    cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
    cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;

    cout << "--------------------------------------" << endl;
    cout << "The last entry for this time is:" << endl;
    ind = indRes[indRes.size()-1];
    cout << "t0_EM    = " << t0_CMU_EM[ind]  << endl;
    cout << "tf_EM    = " << tf_man_SEM[ind]/SEML.us_em.ns << endl;
    cout << "pmin_SEM = " << pmin_dist_SEM[ind] << endl;
    cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
    cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;
    coutlp();

    //==========================================================
    // Copy the selected results in the inputs
    //==========================================================
    for(int ind = 0; ind < (int) indRes.size(); ind++)
    {
        //Times
        t0_CMU_EM_0.push_back(t0_CMU_EM[indRes[ind]]);
        tf_man_EM_0.push_back(tf_man_SEM[indRes[ind]]/SEML.us_em.ns);
        //CMU of EM
        s1_CMU_EM_0.push_back(s1_CMU_EM[indRes[ind]]);
        s2_CMU_EM_0.push_back(s2_CMU_EM[indRes[ind]]);
        s3_CMU_EM_0.push_back(s3_CMU_EM[indRes[ind]]);
        s4_CMU_EM_0.push_back(s4_CMU_EM[indRes[ind]]);
        s5_CMU_EM_0.push_back(s5_CMU_EM[indRes[ind]]);
        //Projection distance
        pmin_dist_SEM_0.push_back(pmin_dist_SEM[indRes[ind]]);
        //CM of SEM
        s1_CM_SEM_0.push_back(s1_CM_SEM[indRes[ind]]);
        s2_CM_SEM_0.push_back(s2_CM_SEM[indRes[ind]]);
        s3_CM_SEM_0.push_back(s3_CM_SEM[indRes[ind]]);
        s4_CM_SEM_0.push_back(s4_CM_SEM[indRes[ind]]);
    }

    //==========================================================
    // Sort data wrt to the projection distance
    //==========================================================
    sortId = sort_indexes(pmin_dist_SEM_0);

}



/**
 *  \brief Store in a data file the best connections between EML2 and SEML1,2.
 *         Used in int_proj_CMU_EM_on_CM_SEM.
 **/
void writeIntProjSortCU_bin(string filename,
                            double**** init_state_CMU_NCEM,      //initial state in NCEM coordinates
                            double**** init_state_CMU_RCM,       //initial state in RCM coordinates
                            double**** final_state_CMU_SEM,      //final state in SEM coordinates
                            double**** projected_state_CMU_SEM,  //projected state in SEM coordinates
                            double**** projected_state_CMU_RCM,  //projected state in RCM coordinates
                            double** *min_proj_dist_tens_SEM,     //minimum distance of projection in SEM coordinates
                            vector<size_t>& sortId, vector<int>& ktMin,
                            vector<int>& ks1Min, vector<int>& ks3Min,
                            vector<double>& t0_min_EM, vector<double>& tf_min_EM,
                            vector<double>& distMin, int number_of_sol)
{
    int ksortpos, ks1pos, ks3pos, ktpos;
    //----------------------------------------------------------
    //Open datafile
    //----------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of stored solutions
        //---------------------
        resi = number_of_sol;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int kpos = 0; kpos <= number_of_sol; kpos++)
        {
            //Get sorted indices
            ksortpos = sortId[kpos];
            ks1pos   = ks1Min[ksortpos];
            ks3pos   = ks3Min[ksortpos];
            ktpos    = ktMin[ksortpos];

            //1. label
            res = kpos;
            filestream.write((char*) &res, sizeof(double));

            //2. t0 in EM coordinates
            res = t0_min_EM[ksortpos];
            filestream.write((char*) &res, sizeof(double));

            // 3-8. init_state_CMU_NCEM state in NCEM coordinates again
            for (int i = 0; i < 6; i++)
            {
                res = init_state_CMU_NCEM[i][ktpos][ks1pos][ks3pos];
                filestream.write((char*) &res, sizeof(double));
            }

            //9. s1 (EM)
            res = init_state_CMU_RCM[0][ktpos][ks1pos][ks3pos];
            filestream.write((char*) &res, sizeof(double));

            //10. s3 (EM)
            res = init_state_CMU_RCM[2][ktpos][ks1pos][ks3pos];
            filestream.write((char*) &res, sizeof(double));

            //11. tf in EM coordinates
            res = tf_min_EM[ksortpos];
            filestream.write((char*) &res, sizeof(double));

            //12-17. yf in SEM coordinates
            for(int i = 0; i <6; i++)
            {
                res = final_state_CMU_SEM[i][ktpos][ks1pos][ks3pos];
                filestream.write((char*) &res, sizeof(double));
            }

            //18-23. yp in SEM coordinates
            for(int i = 0; i <6; i++)
            {
                res = projected_state_CMU_SEM[i][ktpos][ks1pos][ks3pos];
                filestream.write((char*) &res, sizeof(double));
            }

            //24. s1 (SEM)
            res = projected_state_CMU_RCM[0][ktpos][ks1pos][ks3pos];
            filestream.write((char*) &res, sizeof(double));

            //25. s3 (SEM)
            res = projected_state_CMU_RCM[2][ktpos][ks1pos][ks3pos];
            filestream.write((char*) &res, sizeof(double));

            //26. min_proj_dist_SEM (1)
            res = min_proj_dist_tens_SEM[ktpos][ks1pos][ks3pos];
            filestream.write((char*) &res, sizeof(double));

            //27. min_proj_dist_SEM (2)
            res = distMin[ksortpos];
            filestream.write((char*) &res, sizeof(double));
        }
        filestream.close();
    }
}

/**
 *  \brief Read from a data file the best connections between EML2 and SEML1,2.
 *         Used in int_sorted_sol_CMU_EM_to_CM_SEM. The data file must have been build via writeIntProjSortCU_bin.
 **/
void readIntProjSortCU_bin(string filename,
                           double* label,
                           double* t0_EM,
                           double* tf_EM,
                           double* s1_CMU_EM,
                           double* s3_CMU_EM,
                           double* s1_CM_SEM,
                           double* s3_CM_SEM,
                           double** init_state_CMU_NCEM,
                           double** final_state_CMU_SEM,
                           double** projected_state_CMU_SEM,
                           double* min_proj_dist_SEM_1,
                           double* min_proj_dist_SEM_2,
                           int number_of_sol)
{
    //---------------------
    //Open datafile
    //---------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of stored solutions
        //---------------------
        filestream.read((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int kpos = 0; kpos <= number_of_sol; kpos++)
        {
            //1. label
            filestream.read((char*) &res, sizeof(double));
            label[kpos] = res;

            //2. t0 in EM coordinates
            filestream.read((char*) &res, sizeof(double));
            t0_EM[kpos] = res;

            // 3-8. init_state_CMU_NCEM state in NCEM coordinates again
            for (int i = 0; i < 6; i++)
            {
                filestream.read((char*) &res, sizeof(double));
                init_state_CMU_NCEM[i][kpos] = res;

            }

            //9. s1 (EM)
            filestream.read((char*) &res, sizeof(double));
            s1_CMU_EM[kpos] = res;

            //10. s3 (EM)
            filestream.read((char*) &res, sizeof(double));
            s3_CMU_EM[kpos] = res;

            //11. tf in EM coordinates
            filestream.read((char*) &res, sizeof(double));
            tf_EM[kpos] = res;

            //12-17. yf in SEM coordinates
            for(int i = 0; i <6; i++)
            {
                filestream.read((char*) &res, sizeof(double));
                final_state_CMU_SEM[i][kpos] = res;
            }

            //18-23. yp in SEM coordinates
            for(int i = 0; i <6; i++)
            {
                filestream.read((char*) &res, sizeof(double));
                projected_state_CMU_SEM[i][kpos] = res;
            }

            //24. s1 (SEM)
            filestream.read((char*) &res, sizeof(double));
            s1_CM_SEM[kpos] = res;

            //25. s3 (SEM)
            filestream.read((char*) &res, sizeof(double));
            s3_CM_SEM[kpos] = res;

            //26. min_proj_dist_SEM (1)
            filestream.read((char*) &res, sizeof(double));
            min_proj_dist_SEM_1[kpos] = res;

            //27. min_proj_dist_SEM (2)
            filestream.read((char*) &res, sizeof(double));
            min_proj_dist_SEM_2[kpos] = res;

        }
        filestream.close();
    }
}


