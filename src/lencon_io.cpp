#include "lencon_io.h"



//========================================================================================
// Function: fileExists
//========================================================================================
/**
 *   Check if a file exists
 *   @param[in] filename - the name of the file to check
 *   @return    true if the file exists, else false
 *
 */
bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}


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
    string order_t0_bin      = numTostring(ofts_order)+"_t0_"+numTostring(t0)+".bin";
    string order_dest_t0_bin = numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".bin";
    string order_dest_t0_txt = numTostring(ofts_order)+"_dest_L"+numTostring(destination)+"_t0_"+numTostring(t0)+".txt";
    switch(type)
    {
    case TYPE_CU:
        return SEML.cs_em.F_PLOT+"cu_order_"+order_t0_bin;
    case TYPE_CU_3D:
        return SEML.cs_em.F_PLOT+"cu_3d_order_"+order_t0_bin;
    case TYPE_MAN:
        return SEML.cs_em.F_PLOT+"intcu_order_"+order_t0_bin;

    case TYPE_MAN_PROJ:
        return SEML.cs_em.F_PLOT+"projcu_order_"+order_dest_t0_bin;
    case TYPE_MAN_SORT:
        return SEML.cs_em.F_PLOT+"sortprojcu_order_"+order_dest_t0_bin;
    case TYPE_MAN_SORT_IN:
        return SEML.cs_em.F_PLOT+"sortprojintcu_order_"+order_dest_t0_bin;
    case TYPE_CONT_ATF_TRAJ:
        return SEML.cs_em.F_PLOT+"cont_atf_traj_order_"+order_dest_t0_bin;
    case TYPE_CONT_JPL_TRAJ:
        return SEML.cs_em.F_PLOT+"cont_jpl_order_"+order_dest_t0_bin;

    case TYPE_CONT_ATF:
        return SEML.cs_em.F_PLOT+"cont_atf_order_"+order_dest_t0_txt;

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
int readCOMP_txt(double *t_traj_n, double **y_traj_n, int final_index)
{

    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_COMP_FOR_JPL, SEML.li_SEM);

    //Check the existence of the file
    if(!fileExists(filename))
    {
        cout << "readCOMP_txt. " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }


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

    return FTC_SUCCESS;
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

    //Check the existence of the file
    if(!fileExists(filename))
    {
        cout << "getLengthCOMP_txt. " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }

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
//----------------------------------------------------------------------------------------
// CU
//----------------------------------------------------------------------------------------
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
    //------------------------------------------------------------------------------------
    //Filename
    //------------------------------------------------------------------------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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


//----------------------------------------------------------------------------------------
// CU 3D
//----------------------------------------------------------------------------------------

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

/**
 *  \brief Get the length of the data file the containing the Initial Conditions of a 3D Center-Unstable manifold. Used in compute_grid_CMU_EM_3D.
 *         The data are of type t0*s1*s2*s3*s4 and of size (t_grid_size +1)*(si_grid_size[0]+1)*(si_grid_size[1]+1)*(si_grid_size[2]+1)*(si_grid_size[3]+1).
 **/
int getLenghtCU_bin_3D(int* si_grid_size, int* t_grid_size,
                       int ofts_order, int type, int destination)
{
    //Offset
    int offset = 0;

    //------------------------------------------------------------------------------------
    //Filename
    //------------------------------------------------------------------------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;

        //--------------------------------------------------------------------------------
        //Number of data on the time grid
        //--------------------------------------------------------------------------------
        filestream.read((char*) &resi, sizeof(int));
        *t_grid_size = resi;

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid on all four dimensions
        //--------------------------------------------------------------------------------
        for(int i = 0; i < 4; i++)
        {
            filestream.read((char*) &resi, sizeof(int));
            si_grid_size[i] = resi;
        }

        //--------------------------------------------------------------------------------
        //Get the offset
        //--------------------------------------------------------------------------------
        offset = filestream.tellg();

        //--------------------------------------------------------------------------------
        //Close
        //--------------------------------------------------------------------------------
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

    //------------------------------------------------------------------------------------
    //Filename
    //------------------------------------------------------------------------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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

    //------------------------------------------------------------------------------------
    //Filename
    //------------------------------------------------------------------------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        //--------------------------------------------------------------------------------
        //Set the offset
        //--------------------------------------------------------------------------------
        filestream.seekg (offset, ios::beg);

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
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

        //--------------------------------------------------------------------------------
        //Get the offset
        //--------------------------------------------------------------------------------
        offset2 = filestream.tellg();

        //--------------------------------------------------------------------------------
        //Close
        //--------------------------------------------------------------------------------
        filestream.close();

    }
    else return FTC_FAILURE;

    return offset2;
}


//----------------------------------------------------------------------------------------
// Int CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Get the length of the data file containing the best connections between EML2 and SEML1,2.
 *         Used in int_sorted_sol_CMU_EM_to_CM_SEM/ref_CMU_EM_to_CM_SEM_MSD
 **/
int getLengthIntSortedCU_bin(int* number_of_sol, int ofts_order,
                             int type, int destination)
{
    //------------------------------------------------------------------------------------
    //Filename
    //------------------------------------------------------------------------------------
    string filename = filenameCUM(ofts_order, type, destination);

    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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
                        double crossings_NCSEM,               //number of crossings of the x = -1 line (clock/counterclockwise)
                        int collision_NCEM,                   //collision flag, from NCEM flow
                        int kmin,
                        int ks1,
                        int ks3,
                        int kt)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out | ios::app);
    if (filestream.is_open())
    {
        double res;

        double H0_NCEM, H0_NCSEM, H0_EM, H0_SEM;
        double H0_emli_NCEM, H0_emli_NCSEM, H0_emli_EM, H0_emli_SEM;
        double Hf_NCEM, Hf_NCSEM, Hf_EM, Hf_SEM;
        double Hf_semli_NCEM, Hf_semli_NCSEM, Hf_semli_EM, Hf_semli_SEM;

        double yv_NCEM[6], yv_SEM[6];
        double yv_emli_NCEM[6], yv_semli_NCSEM[6];
        double tv_EM, tv_SEM;

        //Origins at both ends
        for(int i = 0; i <6; i++)
        {
            yv_emli_NCEM[i]   = 0.0;
            yv_semli_NCSEM[i] = 0.0;
        }


        //--------------------------------------------------------------------------------
        // Store data column by column
        //--------------------------------------------------------------------------------
        // 1. time grid in NCEM units
        res   = init_time_grid_EM[kt];
        tv_EM = init_time_grid_EM[kt];
        filestream.write((char*) &res, sizeof(double));

        // 2-7. initial state in NCEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res        = init_state_CMU_NCEM[k][kt][ks1][ks3];
            yv_NCEM[k] = init_state_CMU_NCEM[k][kt][ks1][ks3];
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

        // 21. t_man_SEM
        res    = t_man_SEM[kmin];
        tv_SEM = t_man_SEM[kmin];
        filestream.write((char*) &res, sizeof(double));

        // 22-27. final_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res       = final_state_CMU_SEM[k][kt][ks1][ks3];
            yv_SEM[k] = final_state_CMU_SEM[k][kt][ks1][ks3];
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

        //38. Number of crossings of the x = -1 line (clock/counterclockwise)
        res = crossings_NCSEM;
        filestream.write((char*) &res, sizeof(double));

        //39. Collision flag, from NCEM flow
        res = collision_NCEM;
        filestream.write((char*) &res, sizeof(double));

        //--------------------------------------------------------------------------------
        // Computing the energies before storing
        //--------------------------------------------------------------------------------
        // H0 at IC
        H0_NCEM  = qbcp_H_complete(tv_EM, yv_NCEM, NCEM, NCEM);
        H0_NCSEM = qbcp_H_complete(tv_EM, yv_NCEM, NCEM, NCSEM);
        H0_EM    = qbcp_H_complete(tv_EM, yv_NCEM, NCEM, PEM);
        H0_SEM   = qbcp_H_complete(tv_EM, yv_NCEM, NCEM, PSEM);

        // H0 at emli
        H0_emli_NCEM  = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, NCEM);
        H0_emli_NCSEM = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, NCSEM);
        H0_emli_EM    = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, PEM);
        H0_emli_SEM   = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, PSEM);

        // Hf - careful, the state is given in PSEM
        Hf_NCEM  = qbcp_H_complete(tv_SEM, yv_SEM, PSEM, NCEM);
        Hf_NCSEM = qbcp_H_complete(tv_SEM, yv_SEM, PSEM, NCSEM);
        Hf_EM    = qbcp_H_complete(tv_SEM, yv_SEM, PSEM, PEM);
        Hf_SEM   = qbcp_H_complete(tv_SEM, yv_SEM, PSEM, PSEM);

        // Hf at semli - careful, the state is given in NCSEM
        Hf_semli_NCEM  = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCEM);
        Hf_semli_NCSEM = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCSEM);
        Hf_semli_EM    = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PEM);
        Hf_semli_SEM   = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PSEM);


        //--------------------------------------------------------------------------------
        // Then store
        //--------------------------------------------------------------------------------
        //40-43: H0 at IC
        res = H0_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_SEM;
        filestream.write((char*) &res, sizeof(double));

        //44-47: H0 at emli
        res = H0_emli_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_SEM;
        filestream.write((char*) &res, sizeof(double));

        //48-51: Hf
        res = Hf_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_EM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_SEM;
        filestream.write((char*) &res, sizeof(double));

        //52-55: Hf at semli
        res = Hf_semli_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_EM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_SEM;
        filestream.write((char*) &res, sizeof(double));

        filestream.close();
    }
}

//----------------------------------------------------------------------------------------
// Int CU 3D
//----------------------------------------------------------------------------------------
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
                           double crossings_NCSEM,             //number of crossings of the x = -1 line (clock/counterclockwise)
                           int collision_NCEM,                 //collision flag, from NCEM flow
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

        //38. Number of crossings of the x = -1 line (clock/counterclockwise)
        res = crossings_NCSEM;
        filestream.write((char*) &res, sizeof(double));

        //39. Collision flag, from NCEM flow
        res = collision_NCEM;
        filestream.write((char*) &res, sizeof(double));

        filestream.close();
    }
}

//========================================================================================
//
//          I/O (Refinement)
//
//========================================================================================
/**
 *  Routine for comparison of indexes
 **/
vector<size_t> sort_indexes(const vector<double> &v)
{

    // initialize original index locations
    vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(), IdxCompare(v));

    return idx;
}

/**
 *   \brief Resize a vector to its vector of unique elemnts
 **/
void vector_getUnique(vector<double>& v0U)
{
    //Creating iterator on v0U
    std::vector<double>::iterator it;
    it = std::unique (v0U.begin(), v0U.end());

    //Resize
    v0U.resize(std::distance(v0U.begin(),it) );
}

/**
 *   \brief Save in indRes the indices of the occurences of t0 in t0_CMU_EM.
 **/
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
 *   \brief Determine the number of columns in the data file "filename",
 *          given some possible solutions in ncol[], of size nncol.
 *          The underlying hypothesis is: the first two values in the
 *          first column are identical.
 **/
int numberOfColumns(string filename, int ncol[], int nncol)
{
    if(!fileExists(filename))
    {
        cout << "numberOfColumns. " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }


    fstream filestream;
    double c1, c2, res;
    for(int n = 0; n < nncol; n++)
    {
        filestream.open (filename.c_str(), ios::binary | ios::in);

        //First value of first column
        filestream.read((char*) &c1, sizeof(double));
        //Garbage
        for(int k = 2; k <= ncol[n]; k++) filestream.read((char*) &res, sizeof(double));
        //Possible second value of first column
        filestream.read((char*) &c2, sizeof(double));

        cout << "c1 = " << c1 << endl;
        cout << "c2 = " << c2 << endl;

        filestream.close();

        if(c1 == c2) return ncol[n];
    }


    return FTC_FAILURE;
}


/**
 *  \brief Read in a data file the connections between EML2 and SEML1,2.
 *         Interpolate in the data set to get the right desired t0 at EML2 departures.
 **/
int readClosestIntProjCU_bin(string filename, double t0_des,
                              vector<double>& t0_CMU_EM_0, vector<double>& tf_man_EM_0,
                              vector<double>& s1_CMU_EM_0, vector<double>& s2_CMU_EM_0,
                              vector<double>& s3_CMU_EM_0, vector<double>& s4_CMU_EM_0,
                              vector<double>& s5_CMU_EM_0, vector<double>& pmin_dist_SEM_0,
                              vector<double>& s1_CM_SEM_0, vector<double>& s2_CM_SEM_0,
                              vector<double>& s3_CM_SEM_0, vector<double>& s4_CM_SEM_0,
                              vector<size_t>& sortId)
{
    //====================================================================================
    //Check the existence of the filename
    //====================================================================================
    if(!fileExists(filename))
    {
        cout << "readClosestIntProjCU_bin. " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }


    //====================================================================================
    //Temporary variables
    //====================================================================================
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
    vector<double> crossings_NCSEM;
    vector<double> collision_NCEM;

    //====================================================================================
    //Get the number of columns
    //====================================================================================
    int ncol[3] = {37, 39, 55};
    int ncol0   = numberOfColumns(filename, ncol, 3);

    cout << "Number of columns in " << filename << " is: " << ncol0 << endl;
    pressEnter(true);

    //====================================================================================
    //Open and read datafile (the existence) has already been checked)
    //====================================================================================
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        double res;

        //--------------------------------------------------------------------------------
        // Read data
        //--------------------------------------------------------------------------------
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

            //----------------------------------------------------------------------------
            // Last columns depend on the number ncol0
            //----------------------------------------------------------------------------
            switch(ncol0)
            {
            case 37:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                res = -1.0;
                crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                res = -1.0;
                collision_NCEM.push_back(res);

                break;
            }
            case 39:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                collision_NCEM.push_back(res);

                break;
            }

            case 55:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                collision_NCEM.push_back(res);

                //40-43: H0 at IC - NOT SAVED FOR NOW
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));

                //44-47: H0 at emli - NOT SAVED FOR NOW
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));

                //48-51: Hf - NOT SAVED FOR NOW
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));

                //52-55: Hf at semli - NOT SAVED FOR NOW
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));
                filestream.read((char*) &res, sizeof(double));

                break;
            }
            }

        }
        while(!filestream.eof());

        filestream.close();
    }
    else
    {

        cout << "readClosestIntProjCU_bin. Unable to open the file. Check file name: " << endl;
        cout << filename << endl;
        return FTC_FAILURE;
    }


    //====================================================================================
    //Delete last element that is not a real value
    //====================================================================================
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
    crossings_NCSEM.pop_back();
    collision_NCEM.pop_back();


    //====================================================================================
    //Get the unique elements in t0_CMU_EM
    //====================================================================================
    //Copy t0_CMU_EM into t0_CMU_EM_UNIQUE
    vector<double> t0_CMU_EM_UNIQUE(t0_CMU_EM);
    //Get unique elements
    vector_getUnique(t0_CMU_EM_UNIQUE);

    //====================================================================================
    // Print all the indices that match t0_CMU_EM_UNIQUE[xxx]
    //====================================================================================
    cout << "--------------------------------------" << endl;
    cout << "There is " << t0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
    cout << "[" << t0_CMU_EM_UNIQUE[0]/SEML.us_em.T << ", " << t0_CMU_EM_UNIQUE[t0_CMU_EM_UNIQUE.size()-1]/SEML.us_em.T << "]x SEML.us_em.T" << endl;

    //====================================================================================
    // Find the nearest t0 value, or let the user choose if t0_des < 0
    //====================================================================================
    int ti = 0;

    if(t0_des >= 0)
    {
        double dmin = fabs(t0_CMU_EM_UNIQUE[0] - t0_des);

        for(int i = 1; i < (int) t0_CMU_EM_UNIQUE.size(); i++)
        {
            if(fabs(t0_CMU_EM_UNIQUE[i] - t0_des) < dmin)
            {
                dmin = fabs(t0_CMU_EM_UNIQUE[i] - t0_des);
                ti = i;
            }
        }
        cout << "The index of the time closer to the desired one is " << ti;
        cout << ", which corresponds to:" << endl;
    }
    else
    {
        do
        {
            cout << "Please enter a number between " << 0 << " and " << t0_CMU_EM_UNIQUE.size()-1;
            cout << " to select a specific starting time" << endl;
            scanf("%d", &ti);
        }
        while(ti < 0 || ti > (int) t0_CMU_EM_UNIQUE.size()-1);
        cout << "You have selected " << ti << ", which corresponds to:" << endl;
    }

    cout << "t0_EM = " << t0_CMU_EM_UNIQUE[ti]/SEML.us_em.T;
    cout << " x SEML.us_em.T" <<  endl;

    //------------------------------------------------------------------------------------
    // Select the values that matches the desired time
    //------------------------------------------------------------------------------------
    std::vector<size_t> indRes;
    vector_getIndices(indRes, t0_CMU_EM, t0_CMU_EM_UNIQUE[ti]);

    //====================================================================================
    // Copy the selected results in the inputs
    //====================================================================================
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

    //====================================================================================
    // Sort data wrt to the tof (this is done here using the tf, since t0 is constant)
    //====================================================================================
    sortId = sort_indexes(tf_man_EM_0);
    int ind;

    cout << "--------------------------------------" << endl;
    cout << "The min and max time of flights are:" << endl;
    ind = sortId[0];
    cout << "min(tof_EM)  = " << tf_man_EM_0[ind]/SEML.us->T - t0_CMU_EM_0[ind]/SEML.us->T  << " x T" << endl;
    ind = sortId[tf_man_EM_0.size() -1];
    cout << "max(tof_EM)  = " << tf_man_EM_0[ind]/SEML.us->T - t0_CMU_EM_0[ind]/SEML.us->T  << " x T" << endl;

    //====================================================================================
    // Sort data wrt to the projection distance
    //====================================================================================
    sortId = sort_indexes(pmin_dist_SEM_0);
    coutlp();

     return FTC_SUCCESS;
}


//========================================================================================
//
//          I/O (sorted solutions, deprecated)
//
//========================================================================================
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


//========================================================================================
//
//          Display completion
//
//========================================================================================
/**
 *   \brief Display the current completion (percent) of a routine.
 **/
void displayCompletion(string funcname, double percent)
{
    if(floor(percent*0.1) > COMPLETION)
    {
        cout << resetiosflags(ios::floatfield) << resetiosflags(ios::showpos);
        cout << cout <<  setw(5) << setprecision(5);
        cout << "\r" << funcname << ": " << percent << "% completed: ";
        cout << string((int)floor(0.1*percent), '|') << endl;
        cout.flush();
        cout << std::showpos << setiosflags(ios::scientific);
        COMPLETION++;
    }
}
