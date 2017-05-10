#include "lencon_io.h"


//========================================================================================
// To get current path
//========================================================================================
/**
 *  \brief Return the path of the program that is currently running
 **/
string getexepath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

/**
 *  \brief Return the path of the main folder from with the program is currently running.
 *         It is assumed that the binaries are stored in "main/bin/", so that we can
 *         retrieve "main/" by cutting the result of getexepath() before the string "bin/"
 **/
string getmainpath()
{
    string exepath   = getexepath();
    string delimiter = "bin/";
    string mainpath  = exepath.substr(0, exepath.find(delimiter));

    return mainpath;
}


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
 *  \brief Prefix for the data filenames.
 **/
string fileprefix(int type)
{
     switch(type)
    {
    case TYPE_CU:
        return "cu_order_";
    case TYPE_CU_3D:
        return "cu_3d_order_";
    case TYPE_MAN_PROJ:
        return "projcu_order_";
    case TYPE_MAN_PROJ_ORBIT:
        return "projcu_orbit_order_";
    case TYPE_MAN_PROJ_3D:
        return "projcu_3d_order_";
    case TYPE_CONT_ATF_TRAJ:
        return "cont_atf_traj_order_";
    case TYPE_COMP_FOR_JPL:
        return "comp_for_jpl_order_";
    case TYPE_CONT_JPL_TRAJ:
        return "cont_jpl_order_";
    case TYPE_CONT_ATF:
        return "cont_atf_order_";
    case TYPE_TRAJ_FROM_W:
        return "traj_from_w";
    case TYPE_TRAJ_FROM_C:
        return "traj_from_c";
    case TYPE_TRAJ_CELESTIA:
        return "traj_for_celestia";

    default:
        cout << "fileprefix: unknown type. An empty string is returned." << endl;
        return "";
    }
}

/**
 *  \brief Extension for the data filenames.
 **/
string fileext(int type)
{
     switch(type)
    {
    case TYPE_CU:
    case TYPE_MAN_PROJ:
    case TYPE_MAN_PROJ_ORBIT:
    case TYPE_CU_3D:
    case TYPE_MAN_PROJ_3D:
    case TYPE_CONT_ATF_TRAJ:
    case TYPE_CONT_JPL_TRAJ:
    case TYPE_TRAJ_FROM_W:
    case TYPE_TRAJ_FROM_C:
        return ".bin";
    case TYPE_COMP_FOR_JPL:
    case TYPE_CONT_ATF:
        return ".txt";
    case TYPE_TRAJ_CELESTIA: //no ext because it is just a prefix
        return "";
    default:
        cout << "fileext: unknown type. An empty .txt is returned." << endl;
        return ".txt";
    }
}

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, type,
 *         target, and initial energy delta dH of the data.
 **/
string get_filenameCUM_default(string plot_folder, int ofts_order, int type, int target, double t0, double dH)
{
    string order_dt0dH      = numTostring(ofts_order);
    string order_dest_dt0dH = numTostring(ofts_order)+"_dest_L"+numTostring(target);


    if(t0 != -1)
    {
        order_dt0dH      += "_t0_"+numTostring(t0);
        order_dest_dt0dH += "_t0_"+numTostring(t0);
    }

    if(dH != -1)
    {
        order_dt0dH      += "_dH_"+numTostring(dH);
        order_dest_dt0dH += "_dH_"+numTostring(dH);
    }

    switch(type)
    {
    case TYPE_CU:
        return plot_folder+fileprefix(type)+order_dt0dH+fileext(type);
    case TYPE_MAN_PROJ:
    case TYPE_MAN_PROJ_ORBIT:
    case TYPE_CU_3D:
    case TYPE_MAN_PROJ_3D:
    case TYPE_CONT_ATF_TRAJ:
    case TYPE_COMP_FOR_JPL:
    case TYPE_CONT_JPL_TRAJ:
    case TYPE_CONT_ATF:
    case TYPE_TRAJ_FROM_W:
    case TYPE_TRAJ_FROM_C:
    case TYPE_TRAJ_CELESTIA:
        return plot_folder+fileprefix(type)+order_dest_dt0dH+fileext(type);
    default:
        cout << "get_filenameCUM_default: unknown type." << endl;
        return "";
    }
}

/**
 *  \brief Get filename via a window file dialog, with default name filename_default.
 **/
string get_filename_dialog(string filename_default, int mode)
{
#ifdef _MSC_VER
#pragma warning(disable:4996) /* silences warning about strcpy strcat fopen*/
#endif
        string filedef = getmainpath()+filename_default;
        char const* lFilterPatterns[3] = { "*.txt", "*.text", "*.bin" };

        char const *test;
        switch(mode)
        {
            case ios::out:
            {
                test  = tinyfd_saveFileDialog("Save data", filedef.c_str(), 3, lFilterPatterns, NULL);
                break;
            }

            case ios::in:
            default:
            {
                test  = tinyfd_openFileDialog("Read data", filedef.c_str(), 3, lFilterPatterns, NULL, 0);
                break;
            }
        }

        string filename_output = "";
        if(test == NULL) filename_output = filename_default;
        else filename_output = test;

#ifdef _MSC_VER
#pragma warning(default:4996)
#endif

return filename_output;

}

/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, type,
 *         target, initial time t0, and initial energy delta dH of the data.
 **/
string get_filenameCUM(int IO_HANDLING, string plot_folder, string filename_bash,
                       int ofts_order, int type, int target, double t0, double dH, int mode)
{
    string filename_output = "";

    switch(IO_HANDLING)
    {
    case IO_DEFAULT:
    {
        filename_output = get_filenameCUM_default(plot_folder, OFTS_ORDER, type, target, t0, dH);
        break;
    }

    case IO_BASH:
    {
        filename_output = plot_folder+filename_bash;
        break;
    }

    case IO_DIALOG:
    {
        string filename_default = get_filenameCUM_default(plot_folder, OFTS_ORDER, type, target, t0, dH);
        filename_output = get_filename_dialog(filename_default, mode);
        break;
    }
    }

    return filename_output;
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
void write_cref_for_jpl_txt(double* t_traj_n, double** y_traj_n, int final_index, string filename)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    //    string filename = get_filenameCUM(IO_DEFAULT, SEML.cs->F_PLOT, "", OFTS_ORDER,
    //                                      TYPE_COMP_FOR_JPL, SEML.li_SEM, -1, -1, ios::in);
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
 *         (state vectors). The size of the data file is given by getl_cref_for_jpl_txt(), and
 *         the data vectors should be initialized accordingly by the user prior to the use
 *         of this routine.
 **/
int read_cref_for_jpl_txt(double* t_traj_n, double** y_traj_n, int final_index, string filename)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    //    string filename = get_filenameCUM(IO_DEFAULT, SEML.cs->F_PLOT, "", OFTS_ORDER,
    //                                      TYPE_COMP_FOR_JPL, SEML.li_SEM, -1, -1, ios::in);

    //Check the existence of the file
    if(!fileExists(filename))
    {
        cout << "read_cref_for_jpl_txt. " << filename << ": " << strerror(errno) << endl;
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
int getl_cref_for_jpl_txt(string filename)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    //    string filename = get_filenameCUM(IO_DEFAULT, SEML.cs->F_PLOT, "", OFTS_ORDER,
    //                                      TYPE_COMP_FOR_JPL, SEML.li_SEM, -1, -1, ios::in);

    //Check the existence of the file
    if(!fileExists(filename))
    {
        cout << "getl_cref_for_jpl_txt. " << filename << ": " << strerror(errno) << endl;
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
int writeCU_bin_dH(double** * yNCE, double** * sNCE, double** dH, double* tGrid,
                   int s1_grid_size, int t_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //--------------------------------------------------------------------------------
        //Number of data on the time grid
        //--------------------------------------------------------------------------------
        resi = t_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid
        //--------------------------------------------------------------------------------
        resi = s1_grid_size;
        filestream.write((char*) &resi, sizeof(int));


        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        for(int nt = 0; nt <= t_grid_size; nt++)
        {
            //Store the current time
            res = tGrid[nt];
            filestream.write((char*) &res, sizeof(double));

            //Store the data at current time
            for(int n1 = 0; n1 <= s1_grid_size; n1++)
            {
                //NC state
                for (int k = 0; k < 6; k++)
                {
                    res = yNCE[k][nt][n1];
                    filestream.write((char*) &res, sizeof(double));
                }

                //RCM state
                for (int k = 0; k < 5; k++)
                {
                    res = sNCE[k][nt][n1];
                    filestream.write((char*) &res, sizeof(double));
                }

                //Is dH valid
                res = dH[nt][n1];
                filestream.write((char*) &res, sizeof(double));
            }
        }
        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

int getLenghtCU_bin_dH(int* s1_grid_size, int* t_grid_size, string filename)
{

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

        filestream.close();

    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

int readCU_bin_dH(double** * yNCE, double** * sNCE, double** dH, double* tGrid,
                  int s1_grid_size, int t_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        int tSize0, gSize0_s1;

        //--------------------------------------------------------------------------------
        //Number of data on the time grid
        //--------------------------------------------------------------------------------
        filestream.read((char*) &resi, sizeof(int));
        tSize0 = resi;

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid
        //--------------------------------------------------------------------------------
        filestream.read((char*) &resi, sizeof(int));
        gSize0_s1 = resi;


        if(tSize0 < t_grid_size || gSize0_s1 < s1_grid_size)
        {
            cout << "readManifold_bin: wrong inputs" << endl;
            cout << "t_grid_size    = " << t_grid_size    << ", but tSize0 = " << tSize0 << endl;
            cout << "s1_grid_size = " << s1_grid_size << ", but gSize0_s1 = " << gSize0_s1 << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        double res;
        for(int nt = 0; nt <= t_grid_size; nt++)
        {
            //Read the current time
            filestream.read((char*) &res, sizeof(double));
            tGrid[nt]= res;

            //Read the data at current time
            for(int n1 = 0; n1 <= s1_grid_size; n1++)
            {
                //NC state
                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                    yNCE[k][nt][n1] = res;
                }

                //RCM state
                for (int k = 0; k < 5; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                    sNCE[k][nt][n1] = res;
                }

                //Is dH valid
                filestream.read((char*) &res, sizeof(double));
                dH[nt][n1] = res;
            }
        }
        filestream.close();
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

//----------------------------------------------------------------------------------------
// CU
//----------------------------------------------------------------------------------------
/**
 *  \brief Store in a data file the Initial Conditions of a planar Center-Unstable manifold. Used in compute_grid_CMU_EM.
 *         The data are of type t0*s1*s3 and of size (tGrid +1)*(gSize+1)*(gSize+1)
 **/
int writeCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
                int s1_grid_size, int s3_grid_size, int t_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //--------------------------------------------------------------------------------
        //Number of data on the time grid
        //--------------------------------------------------------------------------------
        resi = t_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid
        //--------------------------------------------------------------------------------
        resi = s1_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid
        //--------------------------------------------------------------------------------
        resi = s3_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
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

                    //Is dH valid
                    res = dH[nt][n1][n2];
                    filestream.write((char*) &res, sizeof(double));

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
int readCU_bin(double**** yNCE, double**** sNCE, double** * dH, double* tGrid,
               int s1_grid_size, int s3_grid_size, int t_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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

                    //Is dH valid
                    filestream.read((char*) &res, sizeof(double));
                    dH[nt][n1][n2] = res;
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
                    int* t_grid_size, string filename)
{
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
int initCU_bin_3D(int* si_grid_size, int t_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        int resi;
        //--------------------------------------------------------------------------------
        //Number of data on the time grid
        //--------------------------------------------------------------------------------
        resi = t_grid_size;
        filestream.write((char*) &resi, sizeof(int));

        //--------------------------------------------------------------------------------
        //Number of data on the manifold grid on all four dimensions
        //--------------------------------------------------------------------------------
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
int appTimeCU_bin_3D(double* tGrid, int nt, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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
int writeCU_bin_3D(double** yNCE, double** sNCE, int* si_grid_size, string filename)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
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
int getLenghtCU_bin_3D(int* si_grid_size, int* t_grid_size, string filename)
{
    //Offset
    int offset = 0;

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
int readTCU_bin_3D(int offset, double* tGrid, int nt, string filename)
{
    //Offset
    int offset2 = 0;

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
int readCU_bin_3D(int offset, double** yNCE, double** sNCE, int* si_grid_size, string filename)
{
    //Offset
    int offset2 = 0;

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
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in e.g. int_proj_CMU_EM_on_CM_SEM.
 **/
void writeIntProjCU_bin(string filename,
                        ProjResSt& projResSt)
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

        double yv_emli_NCEM[6], yv_semli_NCSEM[6];
        double yv_IC[6], yv_FC[6];
        double tv_IC, tv_FC;
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
        // 1. label
        res  = projResSt.label;
        filestream.write((char*) &res, sizeof(double));

        // 2. time grid in IC_COORD units
        res  = projResSt.init_time;
        filestream.write((char*) &res, sizeof(double));

        // 3-8. initial state in IC_COORD coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.init_state_CMU_NC[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 9-14. initial state in FC_COORD coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.init_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 15-19. initial state in RCM coordinates
        for (int k = 0; k < 5; k++)
        {
            res = projResSt.init_state_CMU_RCM[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 20. minimum distance of projection, in DP_COORD
        res = projResSt.min_proj_dist_SEM_o;
        filestream.write((char*) &res, sizeof(double));

        // 21. associated dv, in FV_COORD
        res = projResSt.dv_at_projection_FC_o;
        filestream.write((char*) &res, sizeof(double));

        // 22. Final time, in FC_COORD
        res    = projResSt.final_time_FC_o;
        filestream.write((char*) &res, sizeof(double));

        // 23-28. final_state_CMU_SEM state in FC_COORD coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.final_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 29-34. projected_state_CMU_SEM state in FC_COORD coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.projected_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 35-38. projected_state_CMU_RCM state in RCM coordinates
        for (int k = 0; k < 4; k++)
        {
            res = projResSt.projected_state_CMU_RCM_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        //39. Number of crossings of the x = -1 line (clock/counterclockwise), in NCSEM
        res = projResSt.crossings_NCSEM_o;
        filestream.write((char*) &res, sizeof(double));

        //40. Collision flag, from NCEM flow
        res = projResSt.collision_NCEM_o;
        filestream.write((char*) &res, sizeof(double));

        //--------------------------------------------------------------------------------
        // Computing the energies before storing
        //--------------------------------------------------------------------------------
        // States and time, in IC_COORD
        tv_IC = projResSt.init_time;
        state_memcpy(yv_IC, projResSt.init_state_CMU_NC);
        // States and time, in FC_COORD
        tv_FC = projResSt.final_time_FC_o;
        state_memcpy(yv_FC, projResSt.final_state_CMU_FC_o);
        //Time, in EM and SEM units
        tv_EM  = qbcp_coc_time(tv_IC, projResSt.IC_COORD, NCEM);
        tv_SEM = qbcp_coc_time(tv_IC, projResSt.IC_COORD, NCSEM);

        // H0 at IC - careful, the state is given in IC_COORD
        H0_NCEM  = qbcp_H_complete(tv_IC, yv_IC, projResSt.IC_COORD, NCEM);
        H0_NCSEM = qbcp_H_complete(tv_IC, yv_IC, projResSt.IC_COORD, NCSEM);
        H0_EM    = qbcp_H_complete(tv_IC, yv_IC, projResSt.IC_COORD, PEM);
        H0_SEM   = qbcp_H_complete(tv_IC, yv_IC, projResSt.IC_COORD, PSEM);

        // H0 at emli
        H0_emli_NCEM  = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, NCEM);
        H0_emli_NCSEM = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, NCSEM);
        H0_emli_EM    = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, PEM);
        H0_emli_SEM   = qbcp_H_complete(tv_EM, yv_emli_NCEM, NCEM, PSEM);

        // Hf - careful, the state is given in FC_COORD
        Hf_NCEM  = qbcp_H_complete(tv_FC, yv_FC, projResSt.FC_COORD, NCEM);
        Hf_NCSEM = qbcp_H_complete(tv_FC, yv_FC, projResSt.FC_COORD, NCSEM);
        Hf_EM    = qbcp_H_complete(tv_FC, yv_FC, projResSt.FC_COORD, PEM);
        Hf_SEM   = qbcp_H_complete(tv_FC, yv_FC, projResSt.FC_COORD, PSEM);

        // Hf at semli - careful, the state is given in NCSEM
        Hf_semli_NCEM  = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCEM);
        Hf_semli_NCSEM = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCSEM);
        Hf_semli_EM    = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PEM);
        Hf_semli_SEM   = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PSEM);


        //--------------------------------------------------------------------------------
        // Then store
        //--------------------------------------------------------------------------------
        //41-44: H0 at IC
        res = H0_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_SEM;
        filestream.write((char*) &res, sizeof(double));

        //45-48: H0 at emli
        res = H0_emli_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_SEM;
        filestream.write((char*) &res, sizeof(double));

        //49-52: Hf
        res = Hf_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_EM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_SEM;
        filestream.write((char*) &res, sizeof(double));

        //53-56: Hf at semli
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

/**
 *  \brief Store in a data file the connections between EML2 and SEML1,2.
 *         Used in int_proj_ORBIT_EM_on_CM_SEM.
 **/
void writeIntProjCUSeed_bin(string filename,
                            ProjResSt& projResSt)
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
        // 1. label
        res  = projResSt.label;
        filestream.write((char*) &res, sizeof(double));

        // 2. time grid in NCEM units
        res  = projResSt.init_time;
        filestream.write((char*) &res, sizeof(double));

        // 3-8. initial state in NCEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.init_state_CMU_NC[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 9-14. initial state in SEM coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.init_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 15-19. initial state in RCM coordinates
        for (int k = 0; k < 5; k++)
        {
            res = projResSt.init_state_CMU_RCM[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 20. minimum distance of projection
        res = projResSt.min_proj_dist_SEM_o;
        filestream.write((char*) &res, sizeof(double));

        // 21. associated dv
        res = projResSt.dv_at_projection_FC_o;
        filestream.write((char*) &res, sizeof(double));

        // 22. t_man_SEM
        res    = projResSt.final_time_FC_o;
        filestream.write((char*) &res, sizeof(double));

        // 23-28. final_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.final_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 29-34. projected_state_CMU_SEM state in SE coordinates
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.projected_state_CMU_FC_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        // 35-38. projected_state_CMU_RCM state in SE coordinates
        for (int k = 0; k < 4; k++)
        {
            res = projResSt.projected_state_CMU_RCM_o[k];
            filestream.write((char*) &res, sizeof(double));
        }

        //39. Number of crossings of the x = -1 line (clock/counterclockwise)
        res = projResSt.crossings_NCSEM_o;
        filestream.write((char*) &res, sizeof(double));

        //40. Collision flag, from NCEM flow
        res = projResSt.collision_NCEM_o;
        filestream.write((char*) &res, sizeof(double));

        //--------------------------------------------------------------------------------
        // Computing the energies before storing
        //--------------------------------------------------------------------------------
        // States and time
        tv_EM = projResSt.init_time;
        state_memcpy(yv_NCEM, projResSt.init_state_CMU_NC);

        tv_SEM = projResSt.final_time_FC_o;
        state_memcpy(yv_SEM, projResSt.final_state_CMU_FC_o);

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
        //41-44: H0 at IC
        res = H0_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_SEM;
        filestream.write((char*) &res, sizeof(double));

        //45-48: H0 at emli
        res = H0_emli_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_EM;
        filestream.write((char*) &res, sizeof(double));
        res = H0_emli_SEM;
        filestream.write((char*) &res, sizeof(double));

        //49-52: Hf
        res = Hf_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_EM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_SEM;
        filestream.write((char*) &res, sizeof(double));

        //53-56: Hf at semli
        res = Hf_semli_NCEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_NCSEM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_EM;
        filestream.write((char*) &res, sizeof(double));
        res = Hf_semli_SEM;
        filestream.write((char*) &res, sizeof(double));


        //--------------------------------------------------------------------------------
        // Seeds!
        //--------------------------------------------------------------------------------
        //57. seed_time_EM
        res  = projResSt.seed_time;
        filestream.write((char*) &res, sizeof(double));

        //58-61. initial seed in RCM coordinates
        for (int k = 0; k < 4; k++)
        {
            res = projResSt.seed_state_CMU_RCM[k];
            filestream.write((char*) &res, sizeof(double));
        }


        //--------------------------------------------------------------------------------
        // States @ the EML2 pk section
        //--------------------------------------------------------------------------------
        //62. NCEM time at the Pk section
        res  = projResSt.te_NCEM;
        filestream.write((char*) &res, sizeof(double));


        //63-68. NCEM coordinates at the Pk section
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.ye_NCEM[k];
            filestream.write((char*) &res, sizeof(double));
        }

        //69-71. NCEM velocity at the Pk section
        for (int k = 0; k < 3; k++)
        {
            res = projResSt.ve_NCEM[k];
            filestream.write((char*) &res, sizeof(double));
        }

        //72. NCSEM time at the Pk section
        res  = projResSt.te_NCSEM;
        filestream.write((char*) &res, sizeof(double));

        //73-78. NCSEM coordinates at the Pk section
        for (int k = 0; k < 6; k++)
        {
            res = projResSt.ye_NCSEM[k];
            filestream.write((char*) &res, sizeof(double));
        }

        //79-81. NCSEM velocity at the Pk section
        for (int k = 0; k < 3; k++)
        {
            res = projResSt.ve_NCSEM[k];
            filestream.write((char*) &res, sizeof(double));
        }

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
vector<size_t> sort_indexes(const vector<double>& v)
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

//functor for vector_getIndices_range
struct is_within_range
{
public:
    is_within_range(const double& t0) : t0_desired(t0) {}   //constructor
    bool operator()(const double& t0)                 //operator
    {
        return fabs(t0 - t0_desired) < 1e-10;
    }

private:
    double t0_desired;
};

/**
 *   \brief Save in indRes the indices of the occurences of t0 in t0_CMU_EM, within a certain range
 **/
void vector_getIndices_range(vector<size_t>& indRes, vector<double>& t0_CMU_EM, double t0)
{
    //Create an iterator that will contain the position of each component of t0_CMU_EM that matches t0 within 1e-10
    //Here is the first one
    std::vector<double>::iterator itf = std::find_if(t0_CMU_EM.begin(), t0_CMU_EM.end(), is_within_range(t0));

    // Loop on all t0_CMU_EM, update of indRes
    while (itf != t0_CMU_EM.end())
    {
        indRes.push_back(std::distance(t0_CMU_EM.begin(), itf));
        std::advance(itf, 1);
        itf = std::find_if(itf, t0_CMU_EM.end(),  is_within_range(t0));
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

        filestream.close();

        if(c1 == c2) return ncol[n];
    }


    return FTC_FAILURE;
}

//========================================================================================
//
//          I/O (continuation procedures)
//
//========================================================================================
/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void writeCONT_txt(int label, string filename, Orbit& orbit_EM, Orbit& orbit_SEM,
                   double te_NCSEM, double* ye_NCSEM,  int isFirst)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    fstream filestream;

    if(isFirst)
    {
        //================================================================================
        // If it is the first entry, the title of the columns are written
        //================================================================================
        filestream.open (filename.c_str(), ios::out);
        //Title
        filestream << "label t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
        filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
        filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
        filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
        filestream << "te_NCSEM xe_CMS_NCSEM ye_CMS_NCSEM ze_CMS_NCSEM pxe_CMS_NCSEM pye_CMS_NCSEM pze_CMS_NCSEM ";
        filestream << "H0_NCEM H0_NCSEM H0_EM H0_SEM ";
        filestream << "H0_emli_NCEM H0_emli_NCSEM H0_emli_EM H0_emli_SEM ";
        filestream << "Hf_NCEM Hf_NCSEM Hf_EM Hf_SEM ";
        filestream << "Hf_semli_NCEM Hf_semli_NCSEM Hf_semli_EM Hf_semli_SEM ";
        filestream << endl;
    }
    else
    {
        //================================================================================
        // Else, we append
        //================================================================================
        filestream.open (filename.c_str(), ios::out | ios::app);
    }

    //====================================================================================
    //Data storage
    //====================================================================================
    filestream << setprecision(15) <<  setiosflags(ios::scientific) << std::showpos;

    //------------------------------------------------------------------------------------
    // Store from 1 to 8
    //------------------------------------------------------------------------------------
    //1. label
    filestream << label << "  ";
    //2. t0 in EM units
    filestream << orbit_EM.getT0() << "  ";
    //3-7. s0 in RCM coordinates
    for(int i = 0; i < 5; i++) filestream << orbit_EM.getSi()[i]  << "  ";
    //8. tf in EM units
    filestream << orbit_EM.getTf() << "  ";
    //9-13. sf in RCM coordinates
    for(int i = 0; i <5; i++) filestream << orbit_SEM.getSi()[i] << "  ";
    //14-19. z0 in NCEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_EM.getZ0()[i]  << "  ";
    //20-25. zf in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_SEM.getZ0()[i] << "  ";
    //26. thetae_NCSEM
    filestream << te_NCSEM* SEML.us_sem.n << "  ";
    //27-32. ze in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << ye_NCSEM[i] << "  ";


    //------------------------------------------------------------------------------------
    // Computing Energies at t = t0 & t = tf
    //------------------------------------------------------------------------------------
    double H0_NCEM, H0_NCSEM, H0_EM, H0_SEM;
    double H0_emli_NCEM, H0_emli_NCSEM, H0_emli_EM, H0_emli_SEM;
    double Hf_NCEM, Hf_NCSEM, Hf_EM, Hf_SEM;
    double Hf_semli_NCEM, Hf_semli_NCSEM, Hf_semli_EM, Hf_semli_SEM;

    double yv_NCEM[6], yv_NCSEM[6];
    double yv_emli_NCEM[6], yv_semli_NCSEM[6];
    double tv_EM, tv_SEM;

    //Origins at both ends
    for(int i = 0; i <6; i++)
    {
        yv_emli_NCEM[i]   = 0.0;
        yv_semli_NCSEM[i] = 0.0;
    }

    //t0, z0 in NCEM coordinates
    tv_EM = orbit_EM.getT0();
    for(int i = 0; i <6; i++) yv_NCEM[i] = orbit_EM.getZ0()[i];

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

    //tf, yf in NCSEM coordinates
    tv_SEM = orbit_EM.getTf()*SEML.us_em.ns;
    for(int i = 0; i <6; i++) yv_NCSEM[i] = orbit_SEM.getZ0()[i];

    // Hf
    Hf_NCEM  = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, NCEM);
    Hf_NCSEM = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, NCSEM);
    Hf_EM    = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, PEM);
    Hf_SEM   = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, PSEM);

    // Hf at semli
    Hf_semli_NCEM  = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCEM);
    Hf_semli_NCSEM = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCSEM);
    Hf_semli_EM    = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PEM);
    Hf_semli_SEM   = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PSEM);


    //------------------------------------------------------------------------------------
    // Storing Energies
    //------------------------------------------------------------------------------------
    filestream << H0_NCEM  << "  "; //33
    filestream << H0_NCSEM << "  "; //34
    filestream << H0_EM    << "  "; //35
    filestream << H0_SEM   << "  "; //36

    filestream << H0_emli_NCEM  << "  "; //37
    filestream << H0_emli_NCSEM << "  "; //38
    filestream << H0_emli_EM    << "  "; //39
    filestream << H0_emli_SEM   << "  "; //40

    filestream << Hf_NCEM  << "  "; //41
    filestream << Hf_NCSEM << "  "; //42
    filestream << Hf_EM    << "  "; //43
    filestream << Hf_SEM   << "  "; //44

    filestream << Hf_semli_NCEM  << "  "; //45
    filestream << Hf_semli_NCSEM << "  "; //46
    filestream << Hf_semli_EM    << "  "; //47
    filestream << Hf_semli_SEM   << "  "; //48

    filestream << endl;

    filestream.close();
}

/**
 *  \brief Get the length the results of the continuation procedure, in txt file.
 **/
int getLengthCONT_txt(string filename)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    //string filename  = filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, t0xT);

    //Check the existence of the file
    if(!fileExists(filename))
    {
        cout << "getLengthCONT_txt. " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }


    fstream filestream;
    filestream.open (filename.c_str(), ios::in);

    //====================================================================================
    // Check the opening
    //====================================================================================
    if (!filestream.is_open())
    {
        cerr << "getLengthCONT_txt. Cannot open file." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // First value is discarded: these are the column titles
    //====================================================================================
    string ct;
    getline(filestream,  ct);

    //====================================================================================
    // Get the size of the file
    //====================================================================================
    int fsize = 0;
    while (!filestream.eof())
    {
        getline(filestream,  ct);
        fsize++;
    }
    filestream.close();

    return fsize-1;
}


/**
 *  \brief Reads the results of the continuation procedure, in txt file.
 **/
int readCONT_txt(double*  t0_CMU_EM, double*   tf_CMU_EM,
                 double** si_CMU_EM, double** si_CMS_SEM,
                 double** z0_CMU_NCEM, double** z0_CMS_NCSEM,
                 double* tethae, double** ye_NCSEM,
                 double* H0_NCEM, double* He_NCEM,
                 double* H0_NCSEM, double* He_NCSEM,
                 int fsize, string filename)
{
    //====================================================================================
    // Get the size of the file
    //====================================================================================
    int fsize0 = getLengthCONT_txt(filename);
    if(fsize0 != fsize)
    {
        cerr << "readCONT_txt. The user-defined file size mismatch the true size." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    //string filename  = filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, tr0);
    fstream filestream;
    filestream.open (filename.c_str(), ios::in);

    //====================================================================================
    // Check the opening
    //====================================================================================
    if (!filestream.is_open())
    {
        cerr << "readCONT_txt. Cannot open file." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // First value: these are the column titles.
    //  We can then compute the number of columns.
    //====================================================================================
    string ct, s;
    getline(filestream,  ct);
    istringstream iss(ct);

    int ncolumns = 0;
    while ( getline( iss, s, ' ' ) ) {
        if(s.size() > 0)
        {
            ncolumns++;
            //printf( "`%s', %d, \n", s.c_str(), ncolumns );
        }
    }

    //====================================================================================
    // Read the data on each line and close
    //====================================================================================
    switch(ncolumns)
    {
        case 35:
        {

        for(int k = 0; k < fsize; k++)
        {
            //1. t0 in EM units
            filestream >> t0_CMU_EM[k];
            //2-6. s0 in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMU_EM[i][k];
            //7. tf in EM units
            filestream >> tf_CMU_EM[k];
            //8-12. sf in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMS_SEM[i][k];
            //13-18. z0 in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMU_NCEM[i][k];
            //19-24. zf in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMS_NCSEM[i][k];
            //25. tethae
            filestream >> tethae[k];
            //26-31. ze in NCSEM coordinates
            for(int i = 0; i <6; i++) filestream >> ye_NCSEM[i][k];
            // 32.  Initial energy in NCEM coordinates
            filestream >> H0_NCEM[k];
            // 33. Final energy in NCEM coordinates
            filestream >> He_NCEM[k];
            // 34. Initial energy in NCSEM coordinates
            filestream >> H0_NCSEM[k];
            // 35. Final energy in NCSEM coordinates
            filestream >> He_NCSEM[k];
        }

        break;

        }
        case 48:
        {
        for(int k = 0; k < fsize; k++)
        {
            //1. Label
            filestream >> s;

            //2. t0 in EM units
            filestream >> t0_CMU_EM[k];
            //3-7. s0 in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMU_EM[i][k];

            //8. tf in EM units
            filestream >> tf_CMU_EM[k];
            //9-13. sf in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMS_SEM[i][k];

            //14-19. z0 in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMU_NCEM[i][k];

            //20-25. zf in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMS_NCSEM[i][k];

            //26. tethae
            filestream >> tethae[k];
            //27-32. ze in NCSEM coordinates
            for(int i = 0; i <6; i++) filestream >> ye_NCSEM[i][k];

            // 33.  Initial energy in NCEM coordinates
            filestream >> H0_NCEM[k];
            // 34. Initial energy in NCSEM coordinates
            filestream >> H0_NCSEM[k];
            // 35. Initial energy in EM coordinates
            filestream >> s;
            // 36. Initial energy in SEM coordinates
            filestream >> s;

            // 37-40: energy of emli in NCEM, NCSEM, EM, and SEM coordinates
            filestream >> s;
            filestream >> s;
            filestream >> s;
            filestream >> s;

            // 41. Final energy in NCEM coordinates
            filestream >> He_NCEM[k];
            // 42. Final energy in NCSEM coordinates
            filestream >> He_NCSEM[k];
            // 43. Final energy in EM coordinates
            filestream >> s;
            // 44. Final energy in SEM coordinates
            filestream >> s;

            // 45-48: energy of emli in NCEM, NCSEM, EM, and SEM coordinates
            filestream >> s;
            filestream >> s;
            filestream >> s;
            filestream >> s;
        }

        break;


        }

        case 53:
        {
        for(int k = 0; k < fsize; k++)
        {
            //1. Label
            filestream >> s;

            //2. t0 in EM units
            filestream >> t0_CMU_EM[k];
            //3-7. s0 in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMU_EM[i][k];

            //8. t0_seed in EM units
            filestream >> s;
            //9-12. s0_seed in RCM coordinates
            for(int i = 0; i <4; i++) filestream >> s;

            //13. tf in EM units
            filestream >> tf_CMU_EM[k];
            //14-18. sf in RCM coordinates
            for(int i = 0; i <5; i++) filestream >> si_CMS_SEM[i][k];

            //19-24. z0 in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMU_NCEM[i][k];

            //25-30. zf in NCEM coordinates
            for(int i = 0; i <6; i++) filestream >> z0_CMS_NCSEM[i][k];

            //31. tethae
            filestream >> tethae[k];
            //32-37. ze in NCSEM coordinates
            for(int i = 0; i <6; i++) filestream >> ye_NCSEM[i][k];

            // 38.  Initial energy in NCEM coordinates
            filestream >> H0_NCEM[k];
            // 39. Initial energy in NCSEM coordinates
            filestream >> H0_NCSEM[k];
            // 40. Initial energy in EM coordinates
            filestream >> s;
            // 41. Initial energy in SEM coordinates
            filestream >> s;

            // 42-45: energy of emli in NCEM, NCSEM, EM, and SEM coordinates
            filestream >> s;
            filestream >> s;
            filestream >> s;
            filestream >> s;

            // 46. Final energy in NCEM coordinates
            filestream >> He_NCEM[k];
            // 47. Final energy in NCSEM coordinates
            filestream >> He_NCSEM[k];
            // 48. Final energy in EM coordinates
            filestream >> s;
            // 49. Final energy in SEM coordinates
            filestream >> s;

            // 50-53: energy of emli in NCEM, NCSEM, EM, and SEM coordinates
            filestream >> s;
            filestream >> s;
            filestream >> s;
            filestream >> s;
        }

        break;
        }


        default:
        {
            cout << "Error in readCONT_txt: unknown number of columns." << endl;
            return FTC_FAILURE;
        }

    }


    filestream.close();

    return FTC_SUCCESS;
}

/**
 *  \brief Save a given solution as a complete trajectory
 **/
int writeCONT_bin(RefSt& refSt, string filename_res, int dcs, int coord_type,
                  double** y_traj_n, double* t_traj_n, int man_index, int mPlot,
                  Orbit& orbit_EM, Orbit& orbit_SEM, int label,
                  bool isFirst, int comp_orb_eml, int comp_orb_seml)
{
    string fname = "writeCONT_bin";

    //====================================================================================
    // Initialization of temporary variables
    //====================================================================================
    fstream filestream;
    double yv[42], res = 0.0;
    int ode78coll  = 0, status;

    double** ymc_NCSEM  = dmatrix(0, 5, 0, mPlot);
    double* tmc_SEM     = dvector(0, mPlot);
    double** ymc_NCEM   = dmatrix(0, 5, 0, mPlot);
    double* tmc_EM      = dvector(0, mPlot);

    //====================================================================================
    // Transfer leg
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Open the stream. If it is not the first element of the continuation procedure,
    // we append results to the preexisting file
    //------------------------------------------------------------------------------------
    if(isFirst) filestream.open (filename_res.c_str(), ios::out | ios::binary);
    else filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

    //------------------------------------------------------------------------------------
    //Final trajectory on lines, segment by segment + saving
    //------------------------------------------------------------------------------------
    for(int k = 0; k < man_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        status = ode78(ymc_NCSEM, tmc_SEM, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Checks and warnings, if necessary
        if(status != FTC_SUCCESS)
        {
            //At this step (final plot), a simple warning is issued
            cout << fname << ". Warning: ode78 returned a flag." << endl;
            cout << "A collision may have occured during the transfer." << endl;
        }

        if(ode78coll)
        {
            //At this step (plot), a simple warning is issued
            cout << "plottrajsegbyseg. Warning: a collision with a primary has occurred." << endl;
        }

        //--------------------------------------------------------------------------------
        //To NCEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(ymc_NCSEM, tmc_SEM, ymc_NCEM, tmc_EM, mPlot, NCSEM, NCEM);

        //--------------------------------------------------------------------------------
        // Save to data file
        // We save the following outputs:
        //  1. Label of the solution (number)
        //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //  3. Current time in NCSEM coordinates
        //  4. Current state in NCSEM coordinates
        //  5. Current state in NCEM coordinates
        //  6. Hamiltonian in NCSEM coordinates
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= mPlot; p++)
        {
            //1. Label of the solution
            res = label;
            filestream.write((char*) &res, sizeof(double));

            //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
            res = 2;
            filestream.write((char*) &res, sizeof(double));

            //3. Current time in NCSEM coordinates
            res = tmc_SEM[p];
            filestream.write((char*) &res, sizeof(double));

            //4. Current state in NCSEM coordinates
            for(int i = 0; i < 6; i++)
            {
                res  = ymc_NCSEM[i][p];
                yv[i] = ymc_NCSEM[i][p];
                filestream.write((char*) &res, sizeof(double));
            }

            //5. Current state in NCEM coordinates
            for(int i = 0; i < 6; i++)
            {
                res = ymc_NCEM[i][p];
                filestream.write((char*) &res, sizeof(double));
            }

            //6. Hamiltonian in NCSEM coordinates
            res = qbcp_Hn_SEM(tmc_SEM[p], yv, &SEML);
            filestream.write((char*) &res, sizeof(double));
        }
    }
    filestream.close();


    //====================================================================================
    // Compute the initial orbit if necessary
    //====================================================================================
    if(comp_orb_eml)
    {
        //--------------------------------------------------------------------------------
        // Initialize
        //--------------------------------------------------------------------------------
        //We plot every 0.5 day
        int oPlot = 1.0/0.5*refSt.tspan_EM*SEML.cs_em.cr3bp.T/(86400*2*M_PI);
        int nPlot = oPlot;
        double** yorb_NCEM  = dmatrix(0, 5, 0, oPlot);
        double*  torb_EM    = dvector(0, oPlot);
        double** yorb_NCSEM = dmatrix(0, 5, 0, oPlot);
        double*  torb_SEM   = dvector(0, oPlot);


        //Save & Reset the unstable direction
        double hyp_back = orbit_EM.getSi()[4];
        orbit_EM.setSi(0, 4);

        //Set the final time
        orbit_EM.setTf(orbit_EM.getT0()-refSt.tspan_EM);

        //Update the initial state in the orbit, with the RCM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //--------------------------------------------------------------------------------
        //Integration on oPlot+1 fixed grid
        //--------------------------------------------------------------------------------
        nPlot = orbit_EM.traj_int_main(orbit_EM.getTf(), yorb_NCEM, torb_EM, oPlot, INT_PROJ_CHECK);

        //--------------------------------------------------------------------------------
        // old implemenation, for refererence (equivalent to last line)
        //--------------------------------------------------------------------------------
        //        int output = orbit_EM.traj_int_grid(orbit_EM.getTf(), yorb_NCEM, torb_EM, oPlot, 1);
        //
        //        //If output is strictly greater than 0, then the projection procedure inside
        //        //the integration went wrong, and the new index is output.
        //        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //        //is no point stored in the data, except the first one, and the new index is zero
        //        if(output == ORBIT_EPROJ) nPlot = 0;
        //        else if(output > 0) nPlot = output;


        //--------------------------------------------------------------------------------
        //To NCSEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(yorb_NCEM, torb_EM, yorb_NCSEM, torb_SEM, nPlot, NCEM, NCSEM);

        //--------------------------------------------------------------------------------
        // Save to data file
        // We save the following outputs:
        //  1. Label of the solution (number)
        //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //  3. Current time in NCSEM coordinates
        //  4. Current state in NCSEM coordinates
        //  5. Current state in NCEM coordinates
        //  6. Hamiltonian in NCSEM coordinates
        //--------------------------------------------------------------------------------
        if(nPlot > 0)
        {
            //Reopen stream
            filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

            //Storing all points between 0 and nPlot
            for(int p = 0; p <= nPlot; p++)
            {
                //1. Label of the solution
                res = label;
                filestream.write((char*) &res, sizeof(double));

                //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                res = 1;
                filestream.write((char*) &res, sizeof(double));

                //3. Current time in NCSEM coordinates
                res = torb_SEM[p];
                filestream.write((char*) &res, sizeof(double));

                //4. Current state in NCSEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = yorb_NCSEM[i][p];
                    yv[i] = yorb_NCSEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //5. Current state in NCEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = yorb_NCEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //6. Hamiltonian in NCSEM coordinates
                res = qbcp_Hn_SEM(torb_SEM[p], yv, &SEML);
                filestream.write((char*) &res, sizeof(double));
            }
            filestream.close();
        }
        else
        {
            cout << fname << ". Warning: computation of the orbit at EMLi went wrong." << endl;
            cout << "The orbit is not stored." << endl;
        }

        //================================================================================
        // At EML2. The first 4 RCM components are good, as well as the initial time.
        // Hence, we need to update:
        // 1. The last RCM component (unstable part),
        // 2. The final time.
        //================================================================================
        orbit_EM.setSi(hyp_back, 4);
        orbit_EM.setTf(t_traj_n[man_index]/SEML.us_em.ns);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());


        //--------------------------------------------------------------------------------
        //Free variables
        //--------------------------------------------------------------------------------
        free_dmatrix(yorb_NCEM, 0, 5, 0, oPlot);
        free_dvector(torb_EM, 0, oPlot);
        free_dmatrix(yorb_NCSEM, 0, 5, 0, oPlot);
        free_dvector(torb_SEM, 0, oPlot);
    }


    //====================================================================================
    // Compute the final orbit if necessary
    //====================================================================================
    if(comp_orb_seml)
    {
        switch(comp_orb_seml)
        {
        case 1: //computation using projection method
        {
            //----------------------------------------------------------------------------
            // Initialize
            //----------------------------------------------------------------------------
            //We plot every 0.5 days
            int oPlot = 1.0/0.5*refSt.tspan_SEM*SEML.cs_sem.cr3bp.T/(86400*2*M_PI);
            int nPlot = oPlot;
            double** yorb_NCEM  = dmatrix(0, 5, 0, oPlot);
            double*  torb_EM    = dvector(0, oPlot);
            double** yorb_NCSEM = dmatrix(0, 5, 0, oPlot);
            double*  torb_SEM   = dvector(0, oPlot);

            //Save unstable value
            double s5 = orbit_SEM.getSi(4);

            //Reset the unstable direction
            orbit_SEM.setSi(0, 4);

            //Set the final time
            orbit_SEM.setTf(orbit_SEM.getT0()+refSt.tspan_SEM);

            //Update the initial state in the orbit, with the RCM coordinates
            orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());

            //----------------------------------------------------------------------------
            //Integration on oPlot+1 fixed grid
            //----------------------------------------------------------------------------
            int output = orbit_SEM.traj_int_grid(orbit_SEM.getTf(), yorb_NCSEM, torb_SEM, oPlot, 1);

            //If output is strictly greater than 0, then the projection procedure inside
            //the integration went wrong, and the new index is output.
            //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
            //is no point stored in the data, except the first one, and the new index is zero
            if(output == ORBIT_EPROJ) nPlot = 0;
            else if(output > 0) nPlot = output;

            //----------------------------------------------------------------------------
            //To NCEM coordinates
            //----------------------------------------------------------------------------
            qbcp_coc_vec(yorb_NCSEM, torb_SEM, yorb_NCEM, torb_EM, nPlot, NCSEM, NCEM);

            //----------------------------------------------------------------------------
            // Save to data file
            // We save the following outputs:
            //  1. Label of the solution (number)
            //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
            //  3. Current time in NCSEM coordinates
            //  4. Current state in NCSEM coordinates
            //  5. Current state in NCEM coordinates
            //  6. Hamiltonian in NCSEM coordinates
            //----------------------------------------------------------------------------
            if(nPlot > 0)
            {
                //Reopen stream
                filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

                //Storing all points between 0 and nPlot
                for(int p = 0; p <= nPlot; p++)
                {
                    //1. Label of the solution
                    res = label;
                    filestream.write((char*) &res, sizeof(double));

                    //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                    res = 3;
                    filestream.write((char*) &res, sizeof(double));

                    //3. Current time in NCSEM coordinates
                    res = torb_SEM[p];
                    filestream.write((char*) &res, sizeof(double));

                    //4. Current state in NCSEM coordinates
                    for(int i = 0; i < 6; i++)
                    {
                        res   = yorb_NCSEM[i][p];
                        yv[i] = yorb_NCSEM[i][p];
                        filestream.write((char*) &res, sizeof(double));
                    }

                    //4. Current state in NCEM coordinates
                    for(int i = 0; i < 6; i++)
                    {
                        res = yorb_NCEM[i][p];
                        filestream.write((char*) &res, sizeof(double));
                    }

                    //6. Hamiltonian in NCSEM coordinates
                    res = qbcp_Hn_SEM(torb_SEM[p], yv, &SEML);
                    filestream.write((char*) &res, sizeof(double));
                }
                filestream.close();
            }
            else
            {
                cout << fname << ". Warning: computation of the orbit at EMLi went wrong." << endl;
                cout << "The orbit is not stored." << endl;
            }

            //============================================================================
            // At SEMLi, we need to update:
            // 1. The unstable direction, that has been previously erased.
            // 2. The initial time.
            //============================================================================
            orbit_SEM.setSi(s5, 4);
            orbit_SEM.setT0(t_traj_n[man_index]);

            //----------------------------------------------------------------------------
            // Update the initial state in the orbit, with the RCM coordinates
            //----------------------------------------------------------------------------
            orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());

            //----------------------------------------------------------------------------
            //Free variables
            //----------------------------------------------------------------------------
            free_dmatrix(yorb_NCEM, 0, 5, 0, oPlot);
            free_dvector(torb_EM, 0, oPlot);
            free_dmatrix(yorb_NCSEM, 0, 5, 0, oPlot);
            free_dvector(torb_SEM, 0, oPlot);

            break;
        }
        case 2: //computation using reduced coordinates
        {
            //----------------------------------------------------------------------------
            //Read the reduced vector field
            //----------------------------------------------------------------------------
            vector<Oftsc> Fh;
            Fh.reserve(5);
            for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
            readVOFTS_bin(Fh, SEML_SEM.cs->F_PMS+"rvf/fh");

            //----------------------------------------------------------------------------
            //For dot(s) = fh(s)
            //----------------------------------------------------------------------------
            RVF rvf;
            rvf.ofs_order  = SEML.eff_nf_SEM;
            Ofsc AUX(rvf.ofs_order);
            rvf.fh         = &Fh;
            rvf.ofs        = &AUX;
            rvf.order      = OFTS_ORDER;
            rvf.n          = orbit_SEM.getN();
            rvf.reduced_nv = 5;

            gsl_odeiv2_system sys_fh;
            sys_fh.function  = qbfbp_fh;
            sys_fh.jacobian  = NULL;
            sys_fh.dimension = 2*rvf.reduced_nv;
            sys_fh.params    = &rvf;
            const gsl_odeiv2_step_type* T_fh = gsl_odeiv2_step_rk8pd;

            gsl_odeiv2_driver* d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh,
                                      Config::configManager().G_PREC_HSTART(),
                                      Config::configManager().G_PREC_ABS(),
                                      Config::configManager().G_PREC_REL());


            //----------------------------------------------------------------------------
            // Temp variables
            //----------------------------------------------------------------------------
            double t0_SEM = orbit_SEM.getT0();
            double t1_SEM = t0_SEM+refSt.tspan_SEM;
            double z[6], z_EM[6];
            double t2 = t0_SEM;
            int k  = 0;
            double  s1ccm8[2*rvf.reduced_nv]; //CCM8

            //----------------------------------------------------------------------------
            // Initial state in CCM8 form
            //----------------------------------------------------------------------------
            RCMtoCCM8(orbit_SEM.getSi(), s1ccm8, 5);

            //----------------------------------------------------------------------------
            // Reopen  stream
            //----------------------------------------------------------------------------
            filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

            //----------------------------------------------------------------------------
            // Loop
            //----------------------------------------------------------------------------
            while(t2 < t1_SEM && k <= mPlot)
            {
                cout << "t1_SEM -  t2 = " << t1_SEM -  t2 << endl;
                cout << "z[0]         = " << z[0]         << endl;
                cout << "mPlot - k    = " << mPlot - k    << endl;

                //------------------------------------------------------------------------
                //------------------------------------------------------------------------
                //To NCSEM coordinates
                //------------------------------------------------------------------------
                //------------------------------------------------------------------------
                orbit_SEM.evalRCMtoNC(t2, z);

                //------------------------------------------------------------------------
                //To NCEM coordinates
                //------------------------------------------------------------------------
                //------------------------------------------------------------------------
                qbcp_coc(t2, z, z_EM, NCSEM, NCEM);


                //------------------------------------------------------------------------
                // Save to data file
                // We save the following outputs:
                //  1. Label of the solution (number)
                //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                //  3. Current time in NCSEM coordinates
                //  4. Current state in NCSEM coordinates
                //  5. Current state in NCEM coordinates
                //  6. Hamiltonian in NCSEM coordinates
                //------------------------------------------------------------------------
                //1. Label of the solution
                res = label;
                filestream.write((char*) &res, sizeof(double));

                //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                res = 3;
                filestream.write((char*) &res, sizeof(double));

                //3. Current time in NCSEM coordinates
                res = t2;
                filestream.write((char*) &res, sizeof(double));

                //4. Current state in NCSEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res   = z[i];
                    filestream.write((char*) &res, sizeof(double));
                }

                //4. Current state in NCEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = z_EM[i];
                    filestream.write((char*) &res, sizeof(double));
                }

                //6. Hamiltonian in NCSEM coordinates
                res = qbcp_Hn_SEM(t2, z, &SEML);
                filestream.write((char*) &res, sizeof(double));


                //------------------------------------------------------------------------
                //Advance one step
                //------------------------------------------------------------------------
                gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
                orbit_SEM.ccm8torcm(s1ccm8);
                k++;
            }
            filestream.close();

            break;
        }
        }

    }


    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc_NCSEM, 0, 5, 0, mPlot);
    free_dvector(tmc_SEM, 0, mPlot);
    free_dmatrix(ymc_NCEM, 0, 5, 0, mPlot);
    free_dvector(tmc_EM, 0, mPlot);



    return FTC_SUCCESS;
}


//========================================================================================
//
//          I/O (continuation procedures on one orbit)
//
//========================================================================================
/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void write_wref_conn_txt(string filename, Orbit& orbit_EM, Orbit& orbit_SEM,
                        double te_NCSEM,   double* ye_NCSEM,
                        double te_NCEM,    double* ye_NCEM,
                        double ve_NCEM[3], double ve_NCSEM[3],
                        ProjResClass& projRes,
                        int isFirst,  int index)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    fstream filestream;

    if(isFirst)
    {
        //================================================================================
        // If it is the first entry, the title of the columns are written
        //================================================================================
        filestream.open (filename.c_str(), ios::out);
        //Title
        filestream << "label t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
        filestream << "t0_CMU_EM_seed  s1_CMU_EM_seed  s2_CMU_EM_seed  s3_CMU_EM_seed  s4_CMU_EM_seed ";
        filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
        filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
        filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
        filestream << "te_NCSEM xe_CMS_NCSEM ye_CMS_NCSEM ze_CMS_NCSEM pxe_CMS_NCSEM pye_CMS_NCSEM pze_CMS_NCSEM ";
        filestream << "vxe_CMS_NCSEM vye_CMS_NCSEM vze_CMS_NCSEM ";
        filestream << "te_NCEM xe_CMS_NCEM ye_CMS_NCEM ze_CMS_NCEM pxe_CMS_NCEM pye_CMS_NCEM pze_CMS_NCEM ";
        filestream << "vxe_CMS_NCEM vye_CMS_NCEM vze_CMS_NCEM ";
        filestream << "H0_NCEM H0_NCSEM H0_EM H0_SEM ";
        filestream << "H0_emli_NCEM H0_emli_NCSEM H0_emli_EM H0_emli_SEM ";
        filestream << "Hf_NCEM Hf_NCSEM Hf_EM Hf_SEM ";
        filestream << "Hf_semli_NCEM Hf_semli_NCSEM Hf_semli_EM Hf_semli_SEM ";
        filestream << endl;
    }
    else
    {
        //================================================================================
        // Else, we append
        //================================================================================
        filestream.open (filename.c_str(), ios::out | ios::app);
    }

    //====================================================================================
    //Update the seed via projRes
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    double st_EM_seed[4], t_EM_seed = 0.0; int label = 0;
    projRes.update_ic(st_EM, st_SEM, t_EM, st_EM_seed, &t_EM_seed, &t0_SEM, &pmin_dist_SEM_out, &label, index);

    //====================================================================================
    //Data storage
    //====================================================================================
    filestream << setprecision(15) <<  setiosflags(ios::scientific) << std::showpos;

    //------------------------------------------------------------------------------------
    // Store from 1 to 8
    //------------------------------------------------------------------------------------
    //1. label
    filestream << label << "  ";
    //2. t0 in EM units
    filestream << orbit_EM.getT0() << "  ";
    //3-7. s0 in RCM coordinates
    for(int i = 0; i <5; i++) filestream << orbit_EM.getSi()[i]  << "  ";
    //8. t0_seed in EM units
    filestream << t_EM_seed << "  ";
    //9-12. s0_seed in RCM coordinates
    for(int i = 0; i < 4; i++) filestream << st_EM_seed[i]  << "  ";
    //13. tf in EM units
    filestream << orbit_EM.getTf() << "  ";
    //14-18. sf in RCM coordinates
    for(int i = 0; i <5; i++) filestream << orbit_SEM.getSi()[i] << "  ";
    //19-24. z0 in NCEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_EM.getZ0()[i]  << "  ";
    //25-30. zf in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_SEM.getZ0()[i] << "  ";
    //31. te_NCSEM
    filestream << te_NCSEM << "  ";
    //32-37. ze in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << ye_NCSEM[i] << "  ";
    //38-40. ve in NCSEM coordinates
    for(int i = 0; i <3; i++) filestream << ve_NCSEM[i] << "  ";
    //41. te_NCEM
    filestream << te_NCEM << "  ";
    //42-47. ze in NCEM coordinates
    for(int i = 0; i <6; i++) filestream << ye_NCEM[i] << "  ";
    //48-50. ve in NCEM coordinates
    for(int i = 0; i <3; i++) filestream << ve_NCEM[i] << "  ";

    //------------------------------------------------------------------------------------
    // Computing Energies at t = t0 & t = tf
    //------------------------------------------------------------------------------------
    double H0_NCEM, H0_NCSEM, H0_EM, H0_SEM;
    double H0_emli_NCEM, H0_emli_NCSEM, H0_emli_EM, H0_emli_SEM;
    double Hf_NCEM, Hf_NCSEM, Hf_EM, Hf_SEM;
    double Hf_semli_NCEM, Hf_semli_NCSEM, Hf_semli_EM, Hf_semli_SEM;

    double yv_NCEM[6], yv_NCSEM[6];
    double yv_emli_NCEM[6], yv_semli_NCSEM[6];
    double tv_EM, tv_SEM;

    //Origins at both ends
    for(int i = 0; i <6; i++)
    {
        yv_emli_NCEM[i]   = 0.0;
        yv_semli_NCSEM[i] = 0.0;
    }

    //t0, z0 in NCEM coordinates
    tv_EM = orbit_EM.getT0();
    for(int i = 0; i <6; i++) yv_NCEM[i] = orbit_EM.getZ0()[i];

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

    //tf, yf in NCSEM coordinates
    tv_SEM = orbit_EM.getTf()*SEML.us_em.ns;
    for(int i = 0; i <6; i++) yv_NCSEM[i] = orbit_SEM.getZ0()[i];

    // Hf
    Hf_NCEM  = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, NCEM);
    Hf_NCSEM = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, NCSEM);
    Hf_EM    = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, PEM);
    Hf_SEM   = qbcp_H_complete(tv_SEM, yv_NCSEM, NCSEM, PSEM);

    // Hf at semli
    Hf_semli_NCEM  = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCEM);
    Hf_semli_NCSEM = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, NCSEM);
    Hf_semli_EM    = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PEM);
    Hf_semli_SEM   = qbcp_H_complete(tv_SEM, yv_semli_NCSEM, NCSEM, PSEM);


    //------------------------------------------------------------------------------------
    // Storing Energies
    //------------------------------------------------------------------------------------
    filestream << H0_NCEM  << "  "; //51
    filestream << H0_NCSEM << "  "; //52
    filestream << H0_EM    << "  "; //53
    filestream << H0_SEM   << "  "; //54

    filestream << H0_emli_NCEM  << "  "; //55
    filestream << H0_emli_NCSEM << "  "; //56
    filestream << H0_emli_EM    << "  "; //57
    filestream << H0_emli_SEM   << "  "; //58

    filestream << Hf_NCEM  << "  "; //59
    filestream << Hf_NCSEM << "  "; //60
    filestream << Hf_EM    << "  "; //61
    filestream << Hf_SEM   << "  "; //62

    filestream << Hf_semli_NCEM  << "  "; //63
    filestream << Hf_semli_NCSEM << "  "; //64
    filestream << Hf_semli_EM    << "  "; //65
    filestream << Hf_semli_SEM   << "  "; //66

    filestream << endl;

    filestream.close();
}


void write_wref_traj_bin(string filename_traj, Orbit& orbit_EM, Orbit& orbit_SEM,
                         double **y_traj_NCSEM, double *t_traj_NCSEM,
                         double te_NCSEM,   double* ye_NCSEM,
                         double te_NCEM,    double* ye_NCEM,
                         double ve_NCEM[3], double ve_NCSEM[3],
                         ProjResClass& projRes,
                         int isFirst, int index, int man_grid_size)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    fstream filestream;
    double res;
    int resi;

    //====================================================================================
    //Update the seed via projRes
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    double st_EM_seed[4], t_EM_seed = 0.0; int label = 0;
    projRes.update_ic(st_EM, st_SEM, t_EM, st_EM_seed, &t_EM_seed, &t0_SEM, &pmin_dist_SEM_out, &label, index);

    //------------------------------------------------------------------------------------
    // Open the stream. If it is not the first element of the continuation procedure,
    // we append results to the preexisting file
    //------------------------------------------------------------------------------------
    if(isFirst) filestream.open (filename_traj.c_str(), ios::out | ios::binary);
    else filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);


    //====================================================================================
    //Data storage
    //====================================================================================
    // 1. label
    resi = label;
    filestream.write((char*) &resi, sizeof(int));

    //2. t0 in EM units
    res = orbit_EM.getT0();
    filestream.write((char*) &res, sizeof(double));


    //3-8. s0 in RCM coordinates
    for(int i = 0; i <5; i++)
    {
        res = orbit_EM.getSi()[i] ;
        filestream.write((char*) &res, sizeof(double));
    }

    //9. tf in EM units
    res = orbit_EM.getTf();
    filestream.write((char*) &res, sizeof(double));


    //10-14. sf in RCM coordinates
    for(int i = 0; i <5; i++)
    {
        res = orbit_SEM.getSi()[i] ;
        filestream.write((char*) &res, sizeof(double));
    }

    //15. Number of points (man_grid_size+1) on the grid
    resi = man_grid_size;
    filestream.write((char*) &resi, sizeof(int));

    //The rest: the whole grid y_traj_NCSEM/t_traj_NCSEM in NCSEM coordinates
    for(int k = 0; k <= man_grid_size; k++)
    {
        //Time
        res = t_traj_NCSEM[k];
        filestream.write((char*) &res, sizeof(double));


        //State
        for(int i = 0; i <6; i++)
        {
            res = y_traj_NCSEM[i][k];
            filestream.write((char*) &res, sizeof(double));
        }
    }

    //Again, number of points (man_grid_size+1) on the grid
    resi = man_grid_size;
    filestream.write((char*) &resi, sizeof(int));


    //Finally: the index
    resi = index;
    filestream.write((char*) &resi, sizeof(int));


    filestream.close();
}

int getl_wref_traj_bin(string filename_traj, int *last_index, int *man_grid_size)
{
    //------------------------------------------------------------------------------------
    //Open datafile
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename_traj.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        filestream.seekg(-2*sizeof(int), ios::end);
        filestream.read((char*) man_grid_size, sizeof(int));
        filestream.read((char*) last_index, sizeof(int));
        filestream.close();
    }
    else
    {
        cerr << "getl_wref_traj_bin. Cannot open file." << endl;
        return FTC_FAILURE;
    }

   return FTC_SUCCESS;
}


int read_wref_traj_bin(string filename_traj, int *labels,
                       double* t0_CMU_EM, double* tf_CMU_EM,
                       double** si_CMU_EM, double** si_CMS_SEM,
                       double** t_traj_NCSEM, double*** y_traj_NCSEM,
                       int man_grid_size,
                       int last_index)
{
    //------------------------------------------------------------------------------------
    //Check the last label consistency
    //------------------------------------------------------------------------------------
    int last_index_i = -1, man_grid_size_i = -1;
    getl_wref_traj_bin(filename_traj, &last_index_i, &man_grid_size_i);

    if(last_index_i != last_index || man_grid_size_i != man_grid_size )
    {
            cout << "read_wref_traj_bin: wrong inputs" << endl;
            cout << "last_index     = " << last_index    << ", but last_index_i    = " << last_index_i << endl;
            cout << "man_grid_size  = " << man_grid_size << ", but man_grid_size_i = " << man_grid_size_i << endl;
            return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //Open datafile.
    //------------------------------------------------------------------------------------
    fstream filestream;
    filestream.open (filename_traj.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi, man_grid_size = 0;
        double res;

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        for(int label = 0; label <= last_index; label++)
        {
            // 1. label
            filestream.read((char*) &resi, sizeof(int));
            labels[label] = resi;

            //2. t0 in EM units
            filestream.read((char*) &res, sizeof(double));
            t0_CMU_EM[label] = res;


            //3-8. s0 in RCM coordinates
            for(int i = 0; i <5; i++)
            {
                filestream.read((char*) &res, sizeof(double));
                si_CMU_EM[i][label] = res;
            }

            //9. tf in EM units
            filestream.read((char*) &res, sizeof(double));
            tf_CMU_EM[label] = res;


            //10-14. sf in RCM coordinates
            for(int i = 0; i <5; i++)
            {
                filestream.read((char*) &res, sizeof(double));
                si_CMS_SEM[i][label] = res;
            }

            //15. Number of points (man_grid_size+1) on the grid
            filestream.read((char*) &resi, sizeof(int));
            man_grid_size = resi;

            //The rest: the whole grid y_traj_NCSEM/t_traj_NCSEM in NCSEM coordinates
            for(int k = 0; k <= man_grid_size; k++)
            {
                //Time
                filestream.read((char*) &res, sizeof(double));
                t_traj_NCSEM[label][k] = res;


                //State
                for(int i = 0; i < 6; i++)
                {
                    filestream.read((char*) &res, sizeof(double));
                    y_traj_NCSEM[i][label][k] = res;
                }
            }

            //Finally: the man_grid_size is discarded
            filestream.read((char*) &resi, sizeof(int));


            //Finally: the index is discarded
            filestream.read((char*) &resi, sizeof(int));

        }
        filestream.close();
    }
    else return FTC_FAILURE;

   return FTC_SUCCESS;
}

/**
 *  \brief Save a given solution as a complete trajectory
 **/
int write_wref_res_bin(RefSt& refSt, string filename_res, double** y_traj_n, double* t_traj_n,
                        int man_index, Orbit& orbit_EM, Orbit& orbit_SEM, bool isFirst,
                        int comp_orb_eml, int comp_orb_seml, ProjResClass& projRes, int k)
{
    string fname = "writeCONT_bin";

    //====================================================================================
    //Framework
    //====================================================================================
    int coord_type    = refSt.coord_type;
    int dcs           = default_coordinate_system(coord_type);

    //====================================================================================
    //Update the seed via projRes
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    double st_EM_seed[4], t_EM_seed = 0.0; int label = 0;
    projRes.update_ic(st_EM, st_SEM, t_EM, st_EM_seed, &t_EM_seed, &t0_SEM, &pmin_dist_SEM_out, &label, k);
    double r0_CMU_EMT = t_EM[0]/SEML.us_em.T;

    //====================================================================================
    // Initialization of temporary variables
    //====================================================================================
    fstream filestream;
    double yv[42], res = 0.0;
    int ode78coll  = 0, status;

    double** ymc_NCSEM  = dmatrix(0, 5, 0, refSt.mPlot);
    double* tmc_SEM     = dvector(0, refSt.mPlot);
    double** ymc_NCEM   = dmatrix(0, 5, 0, refSt.mPlot);
    double* tmc_EM      = dvector(0, refSt.mPlot);

    //====================================================================================
    // Transfer leg
    //====================================================================================
    cout << fname << ". computing the Transfer leg... " << endl;
    //------------------------------------------------------------------------------------
    // Open the stream. If it is not the first element of the continuation procedure,
    // we append results to the preexisting file
    //------------------------------------------------------------------------------------
    if(isFirst) filestream.open (filename_res.c_str(), ios::out | ios::binary);
    else filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

    //------------------------------------------------------------------------------------
    //Final trajectory on lines, segment by segment + saving
    //------------------------------------------------------------------------------------
    for(int k = 0; k < man_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        status = ode78(ymc_NCSEM, tmc_SEM, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, refSt.mPlot, dcs, coord_type, coord_type);

        //Checks and warnings, if necessary
        if(status != FTC_SUCCESS)
        {
            //At this step (final plot), a simple warning is issued
            cout << fname << ". Warning: ode78 returned a flag." << endl;
            cout << "A collision may have occured during the transfer." << endl;
        }

        if(ode78coll)
        {
            //At this step (plot), a simple warning is issued
            cout << "plottrajsegbyseg. Warning: a collision with a primary has occurred." << endl;
        }

        //--------------------------------------------------------------------------------
        //To NCEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(ymc_NCSEM, tmc_SEM, ymc_NCEM, tmc_EM, refSt.mPlot, NCSEM, NCEM);

        //--------------------------------------------------------------------------------
        // Save to data file
        // We save the following outputs:
        //  1. Label of the solution (number)
        //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //  3. Current time in NCSEM coordinates
        //  4. Current state in NCSEM coordinates
        //  5. Current state in NCEM coordinates
        //  6. Hamiltonian in NCSEM coordinates
        //  7. t_EM_0 as a ratio
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= refSt.mPlot; p++)
        {
            //1. Label of the solution
            res = label;
            filestream.write((char*) &res, sizeof(double));

            //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
            res = 2;
            filestream.write((char*) &res, sizeof(double));

            //3. Current time in NCSEM coordinates
            res = tmc_SEM[p];
            filestream.write((char*) &res, sizeof(double));

            //4. Current state in NCSEM coordinates
            for(int i = 0; i < 6; i++)
            {
                res  = ymc_NCSEM[i][p];
                yv[i] = ymc_NCSEM[i][p];
                filestream.write((char*) &res, sizeof(double));
            }

            //5. Current state in NCEM coordinates
            for(int i = 0; i < 6; i++)
            {
                res = ymc_NCEM[i][p];
                filestream.write((char*) &res, sizeof(double));
            }

            //6. Hamiltonian in NCSEM coordinates
            res = qbcp_Hn_SEM(tmc_SEM[p], yv, &SEML);
            filestream.write((char*) &res, sizeof(double));

            //7. t_EM_0 as a ratio
            res = r0_CMU_EMT;
            filestream.write((char*) &res, sizeof(double));
        }
    }
    filestream.close();


    //====================================================================================
    // Compute the initial orbit if necessary
    //====================================================================================
    if(comp_orb_eml)
    {
        cout << fname << ". computing the EMLi orbit... " << endl;
        //--------------------------------------------------------------------------------
        // Initialize
        //--------------------------------------------------------------------------------
        //We plot every 0.5 day
        int oPlot = 1.0/0.5*refSt.tspan_EM*SEML.cs_em.cr3bp.T/(86400*2*M_PI);
        int nPlot = oPlot;
        double** yorb_NCEM  = dmatrix(0, 5, 0, oPlot);
        double*  torb_EM    = dvector(0, oPlot);
        double** yorb_NCSEM = dmatrix(0, 5, 0, oPlot);
        double*  torb_SEM   = dvector(0, oPlot);


        //Save & Reset the unstable direction
        double hyp_back = orbit_EM.getSi()[4];
        orbit_EM.setSi(0, 4);

        //Set the final time
        orbit_EM.setTf(orbit_EM.getT0()-refSt.tspan_EM);

        //Update the initial state in the orbit, with the RCM coordinates
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //--------------------------------------------------------------------------------
        //Integration on oPlot+1 fixed grid
        //--------------------------------------------------------------------------------
        nPlot = orbit_EM.traj_int_main(orbit_EM.getTf(), yorb_NCEM, torb_EM, oPlot, INT_PROJ_CHECK);

        //--------------------------------------------------------------------------------
        // old implemenation, for refererence (equivalent to last line)
        //--------------------------------------------------------------------------------
        //        int output = orbit_EM.traj_int_grid(orbit_EM.getTf(), yorb_NCEM, torb_EM, oPlot, 1);
        //
        //        //If output is strictly greater than 0, then the projection procedure inside
        //        //the integration went wrong, and the new index is output.
        //        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //        //is no point stored in the data, except the first one, and the new index is zero
        //        if(output == ORBIT_EPROJ) nPlot = 0;
        //        else if(output > 0) nPlot = output;


        //--------------------------------------------------------------------------------
        //To NCSEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(yorb_NCEM, torb_EM, yorb_NCSEM, torb_SEM, nPlot, NCEM, NCSEM);

        //--------------------------------------------------------------------------------
        // Save to data file
        // We save the following outputs:
        //  1. Label of the solution (number)
        //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //  3. Current time in NCSEM coordinates
        //  4. Current state in NCSEM coordinates
        //  5. Current state in NCEM coordinates
        //  6. Hamiltonian in NCSEM coordinates
        //  7. t_EM_0 as a ratio
        //--------------------------------------------------------------------------------
        if(nPlot > 0)
        {
            //Reopen stream
            filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

            //Storing all points between 0 and nPlot
            for(int p = 0; p <= nPlot; p++)
            {
                //1. Label of the solution
                res = label;
                filestream.write((char*) &res, sizeof(double));

                //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                res = 1;
                filestream.write((char*) &res, sizeof(double));

                //3. Current time in NCSEM coordinates
                res = torb_SEM[p];
                filestream.write((char*) &res, sizeof(double));

                //4. Current state in NCSEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = yorb_NCSEM[i][p];
                    yv[i] = yorb_NCSEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //5. Current state in NCEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = yorb_NCEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //6. Hamiltonian in NCSEM coordinates
                res = qbcp_Hn_SEM(torb_SEM[p], yv, &SEML);
                filestream.write((char*) &res, sizeof(double));

                //7. t_EM_0 as a ratio
                res = r0_CMU_EMT;
                filestream.write((char*) &res, sizeof(double));
            }
            filestream.close();
        }
        else
        {
            cout << fname << ". Warning: computation of the orbit at EMLi went wrong." << endl;
            cout << "The orbit is not stored." << endl;
        }

        //================================================================================
        // At EML2. The first 4 RCM components are good, as well as the initial time.
        // Hence, we need to update:
        // 1. The last RCM component (unstable part),
        // 2. The final time.
        //================================================================================
        orbit_EM.setSi(hyp_back, 4);
        orbit_EM.setTf(t_traj_n[man_index]/SEML.us_em.ns);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());


        //--------------------------------------------------------------------------------
        //Free variables
        //--------------------------------------------------------------------------------
        free_dmatrix(yorb_NCEM, 0, 5, 0, oPlot);
        free_dvector(torb_EM, 0, oPlot);
        free_dmatrix(yorb_NCSEM, 0, 5, 0, oPlot);
        free_dvector(torb_SEM, 0, oPlot);
    }

    //====================================================================================
    // Compute the final orbit if necessary
    //====================================================================================
    if(comp_orb_seml)
    {
        //--------------------------------------------------------------------------------
        // Initialize
        //--------------------------------------------------------------------------------
        //We plot every 0.5 days
        int oPlot = 1.0/0.5*refSt.tspan_SEM*SEML.cs_sem.cr3bp.T/(86400*2*M_PI);
        int nPlot = oPlot;
        double** yorb_NCEM  = dmatrix(0, 5, 0, oPlot);
        double*  torb_EM    = dvector(0, oPlot);
        double** yorb_NCSEM = dmatrix(0, 5, 0, oPlot);
        double*  torb_SEM   = dvector(0, oPlot);

        //Save unstable value
        double s5 = orbit_SEM.getSi(4);

        //Reset the unstable direction
        //orbit_SEM.setSi(0, 4);

        //Set the final time
        orbit_SEM.setTf(orbit_SEM.getT0()+refSt.tspan_SEM);

        //Update the initial state in the orbit, with the RCM coordinates
        orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());

        //--------------------------------------------------------------------------------
        // Switch on the integration method, in comp_orb_seml
        //--------------------------------------------------------------------------------
        nPlot = orbit_SEM.traj_int_main(orbit_SEM.getTf(), yorb_NCSEM, torb_SEM, oPlot, comp_orb_seml);

        //--------------------------------------------------------------------------------
        // OLD CODE, FOR REFERENCE (EQUIVALENT TO LAST LINE)
        //--------------------------------------------------------------------------------
        //        switch(comp_orb_seml)
        //        {
        //        case 1: //computation using projection method
        //        {
        //            cout << fname << ". computing the SEMLi orbit using projection method... " << endl;
        //            //----------------------------------------------------------------------------
        //            //Integration on oPlot+1 fixed grid
        //            //----------------------------------------------------------------------------
        //            int output = orbit_SEM.traj_int_grid(orbit_SEM.getTf(), yorb_NCSEM, torb_SEM, oPlot, 1);
        //
        //            //If output is strictly greater than 0, then the projection procedure inside
        //            //the integration went wrong, and the new index is output.
        //            //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //            //is no point stored in the data, except the first one, and the new index is zero
        //            if(output == ORBIT_EPROJ) nPlot = 0;
        //            else if(output > 0) nPlot = output;
        //            break;
        //        }
        //        case 2: //computation using reduced coordinates
        //        {
        //            /*
        //            // ### Old version ###
        //            //----------------------------------------------------------------------------
        //            //Read the reduced vector field
        //            //----------------------------------------------------------------------------
        //            vector<Oftsc> Fh;
        //            Fh.reserve(4);
        //            for(int i = 0; i < 4; i++) Fh.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
        //            readVOFTS_bin(Fh, SEML_SEM.cs->F_GS+"rvf/fh");
        //
        //            //----------------------------------------------------------------------------
        //            //For dot(s) = fh(s)
        //            //----------------------------------------------------------------------------
        //            RVF rvf;
        //            rvf.ofs_order  = SEML.eff_nf_SEM;
        //            Ofsc AUX(rvf.ofs_order);
        //            rvf.fh         = &Fh;
        //            rvf.ofs        = &AUX;
        //            rvf.order      = OFTS_ORDER;
        //            rvf.n          = orbit_SEM.getN();
        //            rvf.reduced_nv = 4;
        //
        //            gsl_odeiv2_system sys_fh;
        //            sys_fh.function  = qbfbp_fh;
        //            sys_fh.jacobian  = NULL;
        //            sys_fh.dimension = 2*rvf.reduced_nv;
        //            sys_fh.params    = &rvf;
        //            const gsl_odeiv2_step_type* T_fh = gsl_odeiv2_step_rk8pd;
        //
        //            gsl_odeiv2_driver* d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh,
        //                                      Config::configManager().G_PREC_HSTART(),
        //                                      Config::configManager().G_PREC_ABS(),
        //                                      Config::configManager().G_PREC_REL());
        //
        //
        //            //----------------------------------------------------------------------------
        //            // Temp variables
        //            //----------------------------------------------------------------------------
        //            double t0_SEM = orbit_SEM.getT0();
        //            double t1_SEM = t0_SEM+refSt.tspan_SEM;
        //            double z[6], z_EM[6];
        //            double t2 = t0_SEM;
        //            int k  = 0;
        //            double  s1ccm8[2*rvf.reduced_nv]; //CCM8
        //
        //            //----------------------------------------------------------------------------
        //            // Initial state in CCM8 form
        //            //----------------------------------------------------------------------------
        //            RCMtoCCM8(orbit_SEM.getSi(), s1ccm8, 4);
        //
        //            cout << "s1ccm8 = " << endl;
        //            vector_printf_prec(s1ccm8, 8);
        //            pressEnter(true);
        //
        //            //----------------------------------------------------------------------------
        //            // Reopen  stream
        //            //----------------------------------------------------------------------------
        //            filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
        //
        //            //----------------------------------------------------------------------------
        //            // Loop
        //            //----------------------------------------------------------------------------
        //            while(t2 < t1_SEM && k <= refSt.mPlot)
        //            {
        //                //------------------------------------------------------------------------
        //                //------------------------------------------------------------------------
        //                //To NCSEM coordinates
        //                //------------------------------------------------------------------------
        //                //------------------------------------------------------------------------
        //                orbit_SEM.evalRCMtoNC(t2, z);
        //
        //                if(k == 0)
        //                {
        //                    cout << "t0 = " << t0_SEM << endl;
        //                    cout << "z0 = " << endl;
        //                    vector_printf_prec(z, 6);
        //                }
        //
        //                //------------------------------------------------------------------------
        //                //To NCEM coordinates
        //                //------------------------------------------------------------------------
        //                //------------------------------------------------------------------------
        //                qbcp_coc(t2, z, z_EM, NCSEM, NCEM);
        //
        //
        //                //------------------------------------------------------------------------
        //                // Save to data file
        //                // We save the following outputs:
        //                //  1. Label of the solution (number)
        //                //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //                //  3. Current time in NCSEM coordinates
        //                //  4. Current state in NCSEM coordinates
        //                //  5. Current state in NCEM coordinates
        //                //  6. Hamiltonian in NCSEM coordinates
        //                //  7. t_EM_0 as a ratio
        //                //------------------------------------------------------------------------
        //                //1. Label of the solution
        //                res = label;
        //                filestream.write((char*) &res, sizeof(double));
        //
        //                //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //                res = 3;
        //                filestream.write((char*) &res, sizeof(double));
        //
        //                //3. Current time in NCSEM coordinates
        //                res = t2;
        //                filestream.write((char*) &res, sizeof(double));
        //
        //                //4. Current state in NCSEM coordinates
        //                for(int i = 0; i < 6; i++)
        //                {
        //                    res   = z[i];
        //                    filestream.write((char*) &res, sizeof(double));
        //                }
        //
        //                //4. Current state in NCEM coordinates
        //                for(int i = 0; i < 6; i++)
        //                {
        //                    res = z_EM[i];
        //                    filestream.write((char*) &res, sizeof(double));
        //                }
        //
        //                //6. Hamiltonian in NCSEM coordinates
        //                res = qbcp_Hn_SEM(t2, z, &SEML);
        //                filestream.write((char*) &res, sizeof(double));
        //
        //                //7. t_EM_0 as a ratio
        //                res = r0_CMU_EMT;
        //                filestream.write((char*) &res, sizeof(double));
        //
        //
        //                //------------------------------------------------------------------------
        //                //Advance one step
        //                //------------------------------------------------------------------------
        //                gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
        //                orbit_SEM.ccm8torcm(s1ccm8);
        //                k++;
        //            }
        //            filestream.close();
        //            //*/
        //
        //            cout << fname << ". computing the SEMLi orbit using reduced vector field... " << endl;
        //            //----------------------------------------------------------------------------
        //            //Integration on oPlot+1 fixed grid
        //            //----------------------------------------------------------------------------
        //            orbit_SEM.traj_red_grid(orbit_SEM.getTf(), yorb_NCSEM, torb_SEM, oPlot);
        //            break;
        //        }
        //        }

        //--------------------------------------------------------------------------------
        //To NCEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(yorb_NCSEM, torb_SEM, yorb_NCEM, torb_EM, nPlot, NCSEM, NCEM);

        //--------------------------------------------------------------------------------
        // Save to data file
        // We save the following outputs:
        //  1. Label of the solution (number)
        //  2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
        //  3. Current time in NCSEM coordinates
        //  4. Current state in NCSEM coordinates
        //  5. Current state in NCEM coordinates
        //  6. Hamiltonian in NCSEM coordinates
        //  7. t_EM_0 as a ratio
        //--------------------------------------------------------------------------------
        if(nPlot > 0)
        {
            //Reopen stream
            filestream.open (filename_res.c_str(), ios::out | ios::binary | ios::app);

            //Storing all points between 0 and nPlot
            for(int p = 0; p <= nPlot; p++)
            {
                //1. Label of the solution
                res = label;
                filestream.write((char*) &res, sizeof(double));

                //2. Type of leg: either emli orbit (1), transfer leg (2) or semli orbit (3)
                res = 3;
                filestream.write((char*) &res, sizeof(double));

                //3. Current time in NCSEM coordinates
                res = torb_SEM[p];
                filestream.write((char*) &res, sizeof(double));

                //4. Current state in NCSEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res   = yorb_NCSEM[i][p];
                    yv[i] = yorb_NCSEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //4. Current state in NCEM coordinates
                for(int i = 0; i < 6; i++)
                {
                    res = yorb_NCEM[i][p];
                    filestream.write((char*) &res, sizeof(double));
                }

                //6. Hamiltonian in NCSEM coordinates
                res = qbcp_Hn_SEM(torb_SEM[p], yv, &SEML);
                filestream.write((char*) &res, sizeof(double));

                //7. t_EM_0 as a ratio
                res = r0_CMU_EMT;
                filestream.write((char*) &res, sizeof(double));
            }
            filestream.close();
        }
        else
        {
            cout << fname << ". Warning: computation of the orbit at EMLi went wrong." << endl;
            cout << "The orbit is not stored." << endl;
        }

        //================================================================================
        // At SEMLi, we need to update:
        // 1. The unstable direction, that has been previously erased.
        // 2. The initial time.
        //================================================================================
        orbit_SEM.setSi(s5, 4);
        orbit_SEM.setT0(t_traj_n[man_index]);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());

        //--------------------------------------------------------------------------------
        //Free variables
        //--------------------------------------------------------------------------------
        free_dmatrix(yorb_NCEM, 0, 5, 0, oPlot);
        free_dvector(torb_EM, 0, oPlot);
        free_dmatrix(yorb_NCSEM, 0, 5, 0, oPlot);
        free_dvector(torb_SEM, 0, oPlot);
    }


    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc_NCSEM, 0, 5, 0, refSt.mPlot);
    free_dvector(tmc_SEM, 0, refSt.mPlot);
    free_dmatrix(ymc_NCEM, 0, 5, 0, refSt.mPlot);
    free_dvector(tmc_EM, 0, refSt.mPlot);



    return FTC_SUCCESS;
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
        cout << cout <<  setw(2) << setprecision(2);
        cout << "\r" << funcname << ": " << percent << "% completed: ";
        cout << string((int)floor(0.1*percent), '|') << endl;
        cout.flush();
        cout << std::showpos << setiosflags(ios::scientific);
        COMPLETION++;
    }
}


//========================================================================================
//
//          ProjResClass class
//
//========================================================================================
//----------------------------------------------------------------------------------------
// Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Push back in all vectors of this, the element kpor of the vectors in inSt.
 *         Update csize accordingly.
 **/
void ProjResClass::push_back(ProjResClass& inSt, int kpor)
{
    //------------------------------------------------------------------------------------
    // Push back in all the vectors
    //------------------------------------------------------------------------------------
    this->t0_CMU_EM.push_back(inSt.t0_CMU_EM[kpor]);
    this->tf_CMU_EM.push_back(inSt.tf_CMU_EM[kpor]);

    this->s1_CMU_EM.push_back(inSt.s1_CMU_EM[kpor]);
    this->s2_CMU_EM.push_back(inSt.s2_CMU_EM[kpor]);
    this->s3_CMU_EM.push_back(inSt.s3_CMU_EM[kpor]);
    this->s4_CMU_EM.push_back(inSt.s4_CMU_EM[kpor]);
    this->s5_CMU_EM.push_back(inSt.s5_CMU_EM[kpor]);

    this->pmin_dist_SEM.push_back(inSt.pmin_dist_SEM[kpor]);
    this->tf_man_SEM.push_back(inSt.tf_man_SEM[kpor]);

    this->s1_CM_SEM.push_back(inSt.s1_CM_SEM[kpor]);
    this->s2_CM_SEM.push_back(inSt.s2_CM_SEM[kpor]);
    this->s3_CM_SEM.push_back(inSt.s3_CM_SEM[kpor]);
    this->s4_CM_SEM.push_back(inSt.s4_CM_SEM[kpor]);

    this->crossings_NCSEM.push_back(inSt.crossings_NCSEM[kpor]);
    this->collision_NCEM.push_back(inSt.collision_NCEM[kpor]);

    this->t0_CMU_EM_seed.push_back(inSt.t0_CMU_EM_seed[kpor]);
    this->s1_CMU_EM_seed.push_back(inSt.s1_CMU_EM_seed[kpor]);
    this->s2_CMU_EM_seed.push_back(inSt.s2_CMU_EM_seed[kpor]);
    this->s3_CMU_EM_seed.push_back(inSt.s3_CMU_EM_seed[kpor]);
    this->s4_CMU_EM_seed.push_back(inSt.s4_CMU_EM_seed[kpor]);

    this->label.push_back(inSt.label[kpor]);

    //------------------------------------------------------------------------------------
    // Update csize accordingly
    //------------------------------------------------------------------------------------
    this->csize = this->t0_CMU_EM.size();
}

/**
 *  \brief Pop back for all the vectors. Update csize accordingly.
 **/
void ProjResClass::pop_back()
{
    //------------------------------------------------------------------------------------
    // Pop back for all the vectors
    //------------------------------------------------------------------------------------
    t0_CMU_EM.pop_back();
    tf_CMU_EM.pop_back();
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
    t0_CMU_EM_seed.pop_back();
    s1_CMU_EM_seed.pop_back();
    s2_CMU_EM_seed.pop_back();
    s3_CMU_EM_seed.pop_back();
    s4_CMU_EM_seed.pop_back();
    label.pop_back();

    //------------------------------------------------------------------------------------
    // Update csize accordingly
    //------------------------------------------------------------------------------------
    csize = t0_CMU_EM.size();
}

/**
 *  \brief Push back in all vectors of this, the elements of the vectors in inSt.
 *         that comply with the conditions defined in refSt.
 *         Update csize accordingly.
 *
 *         Returns a boolean (flag). If true, the subselection that comply with the
 *         conditions defined in refSt is not empty.
 **/
bool ProjResClass::push_back_conditional(ProjResClass& inSt, RefSt& refSt)
{
    //------------------------------------------------------------------------------------
    // Select given intervals for the inputs, via info in refSt
    //------------------------------------------------------------------------------------
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;
    double s2_CMU_EM_MIN, s2_CMU_EM_MAX;
    double s4_CMU_EM_MIN, s4_CMU_EM_MAX;
    double tof_MIN, tof_MAX;

    if(refSt.isLimUD)
    {
        cout << "-------------------------------------------------------------" << endl;
        cout << "   push_back_conditional. User-defined interval of research  " << endl;
        cout << "-------------------------------------------------------------" << endl;
        cout << "Enter a value for s1_CMU_EM_MIN: ";
        cin >> s1_CMU_EM_MIN;
        cout << "Enter a value for s1_CMU_EM_MAX: ";
        cin >> s1_CMU_EM_MAX;
        cout << "Enter a value for s3_CMU_EM_MIN: ";
        cin >> s3_CMU_EM_MIN;
        cout << "Enter a value for s3_CMU_EM_MAX: ";
        cin >> s3_CMU_EM_MAX;

        cout << "Enter a value for tof_MIN (% of T, -1 if no use): ";
        cin >> tof_MIN;
        tof_MIN *= SEML.us->T;
        cout << "Enter a value for tof_MAX (% of T, -1 if no use): ";
        cin >> tof_MAX;
        tof_MAX *= SEML.us->T;
    }
    else
    {
        s1_CMU_EM_MIN = refSt.si_CMU_EM_MIN[0];
        s1_CMU_EM_MAX = refSt.si_CMU_EM_MAX[0];
        s3_CMU_EM_MIN = refSt.si_CMU_EM_MIN[2];
        s3_CMU_EM_MAX = refSt.si_CMU_EM_MAX[2];
        tof_MIN = refSt.tof_MIN;
        tof_MAX = refSt.tof_MAX;
    }

    s2_CMU_EM_MIN = refSt.si_CMU_EM_MIN[1];
    s2_CMU_EM_MAX = refSt.si_CMU_EM_MAX[1];
    s4_CMU_EM_MIN = refSt.si_CMU_EM_MIN[3];
    s4_CMU_EM_MAX = refSt.si_CMU_EM_MAX[3];


    //------------------------------------------------------------------------------------
    // Subselection via loop on all the data in inSt
    //------------------------------------------------------------------------------------
    int kpor  = 0, kpos = 0;
    bool cst  = 0;
    bool flag = 0;

    // Note: we can add the condition on the pmin directly on the while condition
    // because we suppose that inSt.sortId has been updated
    // (hence, once the pmin condition is broken, all subsequent values will be broken)
    while(kpos < inSt.size() && (inSt.pmin_dist_SEM[inSt.sortId[kpos]] <= refSt.pmax_dist_SEM))
    {
        kpor = inSt.sortId[kpos];

        //--------------------------------------------------------------------------------
        // Limit in the reduced coordinates (s1, s2, s3, s4)
        //--------------------------------------------------------------------------------
        cst  = (inSt.s1_CMU_EM[kpor] >= s1_CMU_EM_MIN) & (inSt.s1_CMU_EM[kpor] <= s1_CMU_EM_MAX);
        cst  = cst & (inSt.s3_CMU_EM[kpor] >= s3_CMU_EM_MIN) & (inSt.s3_CMU_EM[kpor] <= s3_CMU_EM_MAX);

        if(refSt.is3D())
        {
            cst  = cst & (inSt.s2_CMU_EM[kpor] >= s2_CMU_EM_MIN) & (inSt.s2_CMU_EM[kpor] <= s2_CMU_EM_MAX);
            cst  = cst & (inSt.s4_CMU_EM[kpor] >= s4_CMU_EM_MIN) & (inSt.s4_CMU_EM[kpor] <= s4_CMU_EM_MAX);
        }

        //--------------------------------------------------------------------------------
        // Limit in the reduced coordinates (s1, s2, s3, s4) - SEEDS
        //--------------------------------------------------------------------------------
        cst  = cst & (inSt.s1_CMU_EM_seed[kpor] >= refSt.si_SEED_EM_MIN[0]) & (inSt.s1_CMU_EM_seed[kpor] <= refSt.si_SEED_EM_MAX[0]);
        cst  = cst & (inSt.s3_CMU_EM_seed[kpor] >= refSt.si_SEED_EM_MIN[2]) & (inSt.s3_CMU_EM_seed[kpor] <= refSt.si_SEED_EM_MAX[2]);

        if(refSt.is3D())
        {
            cst  = cst & (inSt.s2_CMU_EM_seed[kpor] >= refSt.si_SEED_EM_MIN[1]) & (inSt.s2_CMU_EM_seed[kpor] <= refSt.si_SEED_EM_MAX[1]);
            cst  = cst & (inSt.s4_CMU_EM_seed[kpor] >= refSt.si_SEED_EM_MIN[3]) & (inSt.s4_CMU_EM_seed[kpor] <= refSt.si_SEED_EM_MAX[3]);
        }

        //--------------------------------------------------------------------------------
        // Limits in the time of flight (TOF) if necessary
        //--------------------------------------------------------------------------------
        if(tof_MIN > 0) cst = cst & ( (inSt.tf_CMU_EM[kpor] - inSt.t0_CMU_EM[kpor]) > tof_MIN );
        if(tof_MAX > 0) cst = cst & ( (inSt.tf_CMU_EM[kpor] - inSt.t0_CMU_EM[kpor]) < tof_MAX );

        //--------------------------------------------------------------------------------
        // Limits in the crossings if necessary
        //--------------------------------------------------------------------------------
        if(refSt.crossings > 0)
        {
            cst = cst & (inSt.crossings_NCSEM[kpor] == refSt.crossings);
        }

        //--------------------------------------------------------------------------------
        // Getting rid of the collisions, if necessary
        //--------------------------------------------------------------------------------
        if(refSt.isCollisionOn)
        {
            cst = cst & (inSt.collision_NCEM[kpor] == 0);
        }

        if(cst)
        {
            flag = 1;
            this->push_back(inSt, kpor);
        }

        //--------------------------------------------------------------------------------
        // One step
        //--------------------------------------------------------------------------------
        kpos++;
    }


    //------------------------------------------------------------------------------------
    // Update the csize
    //------------------------------------------------------------------------------------
    this->csize = this->t0_CMU_EM.size();

    //------------------------------------------------------------------------------------
    // Sort
    //------------------------------------------------------------------------------------
    this->sort_pmin_dist_SEM();

    return flag;
}


//----------------------------------------------------------------------------------------
// Getters
//----------------------------------------------------------------------------------------
/**
 *  \brief Return csize, the common size of all the vectors.
 **/
int ProjResClass::size()
{
    return csize;
}

/**
 *  \brief Return sortId, the vector of sorted indices.
 **/
vector<size_t>& ProjResClass::getSortId()
{
    return sortId;
}

//----------------------------------------------------------------------------------------
// Display
//----------------------------------------------------------------------------------------
/**
 *  \brief Display some selected elements of the entry k.
 **/
void ProjResClass::displayEntry(int k)
{
    int ind = sortId[k];
    cout << "t0_EM    = "  << t0_CMU_EM[ind]  << endl;
    cout << "tf_EM    = "  << tf_CMU_EM[ind] << endl;
    cout << "pmin_SEM = "  << pmin_dist_SEM[ind] << endl;
    cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind];
    cout << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ", " << s5_CMU_EM[ind] << ")" << endl;
    cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind];
    cout << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;
}

/**
 *  \brief Display some selected elements of the first entry.
 **/
void ProjResClass::displayFirstEntry()
{
    displayEntry(0);
}

/**
 *  \brief Display some selected elements of the last entry.
 **/
void ProjResClass::displayLastEntry()
{
    displayEntry(csize-1);
}

//----------------------------------------------------------------------------------------
// Read & sort
//----------------------------------------------------------------------------------------
/**
 *  \brief Read in a data file the connections between EML1,2 and SEML1,2.
 **/
int ProjResClass::readProjRes(string filename)
{
    string fname = "readProjRes";

    //====================================================================================
    //Check the existence of the filename
    //====================================================================================
    if(!fileExists(filename))
    {
        cout << fname << ". " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }
    else cout << "Getting back data from " << filename << endl;

    //====================================================================================
    //Get the number of columns
    //====================================================================================
    int ncol[6] = {61, 81, 56, 55, 39, 37};
    int ncol0   = numberOfColumns(filename, ncol, 5);

    //cout << "Number of columns in " << filename << " is: " << ncol0 << endl;
    //pressEnter(true);

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
            //0. Label if ncol >= 56
            if(ncol0 >= 56)
            {
                filestream.read((char*) &res, sizeof(double));
                this->label.push_back( (int) res);
            }else this->label.push_back(0);

            // 1. time grid in NCEM units
            filestream.read((char*) &res, sizeof(double));
            this->t0_CMU_EM.push_back(res);

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
            this->s1_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s2_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s3_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s4_CMU_EM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s5_CMU_EM.push_back(res);

            // 19. minimum distance of projection: SAVED
            filestream.read((char*) &res, sizeof(double));
            this->pmin_dist_SEM.push_back(res);

            // 20. associated dv: NOT SAVED
            filestream.read((char*) &res, sizeof(double));

            // 21. tf at SEM: SAVED
            filestream.read((char*) &res, sizeof(double));
            this->tf_man_SEM.push_back(res);

            //Moreover, we save the value in tf_CMU_EM
            this->tf_CMU_EM.push_back(res/SEML.us_em.ns);


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
            this->s1_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s2_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s3_CM_SEM.push_back(res);

            filestream.read((char*) &res, sizeof(double));
            this->s4_CM_SEM.push_back(res);

            //----------------------------------------------------------------------------
            // Last columns depend on the number ncol0
            //----------------------------------------------------------------------------
            switch(ncol0)
            {
            case 37:
            {
                res = -1.0;
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                this->crossings_NCSEM.push_back(res);

                res = 0.0;
                //39. Collision flag, from NCEM flow
                this->collision_NCEM.push_back(res);

                //57. t0_CMU_EM_seed - SET to -1
                this->t0_CMU_EM_seed.push_back(res);

                //58-61. initial seed in RCM coordinates - SET to -1
                this->s1_CMU_EM_seed.push_back(res);
                this->s2_CMU_EM_seed.push_back(res);
                this->s3_CMU_EM_seed.push_back(res);
                this->s4_CMU_EM_seed.push_back(res);

                break;
            }
            case 39:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                this->crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                this->collision_NCEM.push_back(res);

                //57. t0_CMU_EM_seed - SET to -1
                this->t0_CMU_EM_seed.push_back(res);

                //58-61. initial seed in RCM coordinates - SET to -1
                this->s1_CMU_EM_seed.push_back(res);
                this->s2_CMU_EM_seed.push_back(res);
                this->s3_CMU_EM_seed.push_back(res);
                this->s4_CMU_EM_seed.push_back(res);

                break;
            }
            case 55:
            case 56:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                this->crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                this->collision_NCEM.push_back(res);

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

                //57. t0_CMU_EM_seed - SET to -1
                this->t0_CMU_EM_seed.push_back(res);

                //58-61. initial seed in RCM coordinates - SET to -1
                this->s1_CMU_EM_seed.push_back(res);
                this->s2_CMU_EM_seed.push_back(res);
                this->s3_CMU_EM_seed.push_back(res);
                this->s4_CMU_EM_seed.push_back(res);

                break;
            }
            case 61:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                this->crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                this->collision_NCEM.push_back(res);

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

                //57. init_time_EM_seed
                filestream.read((char*) &res, sizeof(double));
                this->t0_CMU_EM_seed.push_back(res);

                //58-61. initial seed in RCM coordinates
                filestream.read((char*) &res, sizeof(double));
                this->s1_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s2_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s3_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s4_CMU_EM_seed.push_back(res);

                break;
            }

            case 81:
            {
                //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                filestream.read((char*) &res, sizeof(double));
                this->crossings_NCSEM.push_back(res);

                //39. Collision flag, from NCEM flow
                filestream.read((char*) &res, sizeof(double));
                this->collision_NCEM.push_back(res);

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

                //57. init_time_EM_seed
                filestream.read((char*) &res, sizeof(double));
                this->t0_CMU_EM_seed.push_back(res);

                //58-61. initial seed in RCM coordinates
                filestream.read((char*) &res, sizeof(double));
                this->s1_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s2_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s3_CMU_EM_seed.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s4_CMU_EM_seed.push_back(res);


                //------------------------------------------------------------------------
                // States @ the EML2 pk section
                //------------------------------------------------------------------------
                //62. NCEM time at the Pk section
                filestream.read((char*) &res, sizeof(double));

                //63-68. NCEM coordinates at the Pk section
                for (int k = 0; k < 6; k++) filestream.read((char*) &res, sizeof(double));

                //69-71. NCEM velocity at the Pk section
                for (int k = 0; k < 3; k++) filestream.read((char*) &res, sizeof(double));

                //72. NCSEM time at the Pk section
                filestream.read((char*) &res, sizeof(double));

                //73-78. NCSEM coordinates at the Pk section
                for (int k = 0; k < 6; k++) filestream.read((char*) &res, sizeof(double));

                //79-81. NCSEM velocity at the Pk section
                for (int k = 0; k < 3; k++) filestream.read((char*) &res, sizeof(double));

                break;
            }
            }
        }
        while(!filestream.eof());

        filestream.close();
    }
    else
    {

        cout << "readProjRes_t0. Unable to open the file. Check file name: " << endl;
        cout << filename << endl;
        return FTC_FAILURE;
    }


    //====================================================================================
    //Delete last element that is not a real value
    //====================================================================================
    this->pop_back();

    //====================================================================================
    //Update the size
    //====================================================================================
    this->csize = this->t0_CMU_EM.size();

    //====================================================================================
    //Sort
    //====================================================================================
    this->sort_pmin_dist_SEM();

    return FTC_SUCCESS;
}

/**
*  \brief Read in a data file the connections between EML2 and SEML1,2.
*         Subselection in the data set to get the right desired t0 at EML2 departures.
**/
int ProjResClass::readProjRes_t0(string filename, double t0_des, int typeOfTimeSelection)
{
    string fname = "readProjRes_t0";

    //====================================================================================
    //Check the existence of the filename
    //====================================================================================
    if(!fileExists(filename))
    {
        cout << fname << ". " << filename << ": " << strerror(errno) << endl;
        return FTC_ENOENT;
    }


    //====================================================================================
    // Open and read datafile
    //====================================================================================
    ProjResClass readSt;
    readSt.readProjRes(filename);

    //====================================================================================
    //Type of selection
    //====================================================================================
    coutmp();
    cout << fname << ". ";

    if(typeOfTimeSelection == TIME_SELECTION_ABSOLUTE)
        cout << " the initial time is searched in [0, Inf[ " << endl;
    else
        cout << " the initial time is searched in [0, T[ " << endl;

    int ti = 0;
    std::vector<size_t> indRes;
    switch(typeOfTimeSelection)
    {
    case TIME_SELECTION_ABSOLUTE:
    {
        //--------------------------------------------------------------------------------
        //Get the unique elements in t0_CMU_EM
        //--------------------------------------------------------------------------------
        //Copy t0_CMU_EM into t0_CMU_EM_UNIQUE
        vector<double> t0_CMU_EM_UNIQUE(readSt.t0_CMU_EM);
        //Get unique elements
        vector_getUnique(t0_CMU_EM_UNIQUE);

        //--------------------------------------------------------------------------------
        // Print the range in t0_CMU_EM_UNIQUE
        //--------------------------------------------------------------------------------
        this->sortId = sort_indexes(t0_CMU_EM_UNIQUE);
        cout << "--------------------------------------" << endl;
        cout << "There is " << t0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
        cout << "[" << t0_CMU_EM_UNIQUE[this->sortId[0]]/SEML.us_em.T << ", ";
        cout << t0_CMU_EM_UNIQUE[this->sortId[t0_CMU_EM_UNIQUE.size()-1]]/SEML.us_em.T;
        cout << "] x T." << endl;

        //--------------------------------------------------------------------------------
        // Find the nearest t0 value, or let the user choose if t0_des < 0
        //--------------------------------------------------------------------------------
        cout << "Desired time is  = " << t0_des/SEML.us_em.T;
        cout << " x T." <<  endl;
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

        cout << "t0 = " << t0_CMU_EM_UNIQUE[ti]/SEML.us_em.T;
        cout << " x T." <<  endl;

        //--------------------------------------------------------------------------------
        // Select the values that matches the desired time
        //--------------------------------------------------------------------------------
        vector_getIndices_range(indRes, readSt.t0_CMU_EM, t0_CMU_EM_UNIQUE[ti]);

        break;
    }

    case TIME_SELECTION_RATIO:
    default:
    {
        //----------------------------------------------------------------------------
        //We get all the unique times, as %T, between [0 and T]
        //----------------------------------------------------------------------------
        vector<double> r0_CMU_EM(readSt.t0_CMU_EM);

        //Get the modulo[T], up to RATIO_DIGITS digits
        double factor = 1.0*pow(10, RATIO_DIGITS);
        for(int k = 0; k < (int) r0_CMU_EM.size(); k++)
        {
            r0_CMU_EM[k] = round(factor*fmod(r0_CMU_EM[k]/SEML.us_em.T, 1.0))/factor;
        }

        //Get unique elements
        vector<double> r0_CMU_EM_UNIQUE(r0_CMU_EM);
        vector_getUnique(r0_CMU_EM_UNIQUE);

        //----------------------------------------------------------------------------
        // Print the range in r0_CMU_EM_UNIQUE
        //----------------------------------------------------------------------------
        this->sortId = sort_indexes(r0_CMU_EM_UNIQUE);
        cout << "--------------------------------------" << endl;
        cout << "There is " << r0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
        cout << "[" << r0_CMU_EM_UNIQUE[this->sortId[0]] << ", ";
        cout << r0_CMU_EM_UNIQUE[this->sortId[r0_CMU_EM_UNIQUE.size()-1]];
        cout << "] x T." << endl;

        //----------------------------------------------------------------------------
        // Find the nearest t0 value, or let the user choose if t0_des < 0
        //----------------------------------------------------------------------------
        double r0_des = round(factor*fmod(t0_des/SEML.us_em.T, 1.0))/factor;

        cout << "Desired time is  = " << r0_des;
        cout << " x T." <<  endl;

        if(t0_des >= 0)
        {
            double dmin = fabs(r0_CMU_EM_UNIQUE[0] - r0_des);

            for(int i = 1; i < (int) r0_CMU_EM_UNIQUE.size(); i++)
            {
                if(fabs(r0_CMU_EM_UNIQUE[i] - r0_des) < dmin)
                {
                    dmin = fabs(r0_CMU_EM_UNIQUE[i] - r0_des);
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
                cout << "Please enter a number between " << 0 << " and " << r0_CMU_EM_UNIQUE.size()-1;
                cout << " to select a specific starting time" << endl;
                scanf("%d", &ti);
            }
            while(ti < 0 || ti > (int) r0_CMU_EM_UNIQUE.size()-1);
            cout << "You have selected " << ti << ", which corresponds to:" << endl;
        }

        cout << "t0 = " << r0_CMU_EM_UNIQUE[ti];
        cout << " x T." <<  endl;

        //----------------------------------------------------------------------------
        // Select the values that matches the desired time
        //----------------------------------------------------------------------------
        vector_getIndices_range(indRes, r0_CMU_EM, r0_CMU_EM_UNIQUE[ti]);
        break;
    }


    }
    cout << "--------------------------------------" << endl;
    coutlp();


    //====================================================================================
    // Copy the selected results in the inputs
    //====================================================================================
    for(int ind = 0; ind < (int) indRes.size(); ind++)
    {
        //Times
        this->t0_CMU_EM.push_back(readSt.t0_CMU_EM[indRes[ind]]);
        this->tf_CMU_EM.push_back(readSt.tf_CMU_EM[indRes[ind]]);
        //CMU of EM
        this->s1_CMU_EM.push_back(readSt.s1_CMU_EM[indRes[ind]]);
        this->s2_CMU_EM.push_back(readSt.s2_CMU_EM[indRes[ind]]);
        this->s3_CMU_EM.push_back(readSt.s3_CMU_EM[indRes[ind]]);
        this->s4_CMU_EM.push_back(readSt.s4_CMU_EM[indRes[ind]]);
        this->s5_CMU_EM.push_back(readSt.s5_CMU_EM[indRes[ind]]);
        //Projection distance
        this->pmin_dist_SEM.push_back(readSt.pmin_dist_SEM[indRes[ind]]);
        //tf_man_SEM
        this->tf_man_SEM.push_back(readSt.tf_man_SEM[indRes[ind]]);
        //CM of SEM
        this->s1_CM_SEM.push_back(readSt.s1_CM_SEM[indRes[ind]]);
        this->s2_CM_SEM.push_back(readSt.s2_CM_SEM[indRes[ind]]);
        this->s3_CM_SEM.push_back(readSt.s3_CM_SEM[indRes[ind]]);
        this->s4_CM_SEM.push_back(readSt.s4_CM_SEM[indRes[ind]]);
        //Crossings
        this->crossings_NCSEM.push_back(readSt.crossings_NCSEM[indRes[ind]]);
        this->collision_NCEM.push_back(readSt.collision_NCEM[indRes[ind]]);
        //Seeds
        this->t0_CMU_EM_seed.push_back(readSt.t0_CMU_EM_seed[indRes[ind]]);
        this->s1_CMU_EM_seed.push_back(readSt.s1_CMU_EM_seed[indRes[ind]]);
        this->s2_CMU_EM_seed.push_back(readSt.s2_CMU_EM_seed[indRes[ind]]);
        this->s3_CMU_EM_seed.push_back(readSt.s3_CMU_EM_seed[indRes[ind]]);
        this->s4_CMU_EM_seed.push_back(readSt.s4_CMU_EM_seed[indRes[ind]]);
        //Label
        this->label.push_back(readSt.label[indRes[ind]]);
    }

    //====================================================================================
    //Update the size
    //====================================================================================
    this->csize = this->t0_CMU_EM.size();

    //====================================================================================
    // Sort data wrt to the tof (this is done here using the tf, since t0 is constant)
    //====================================================================================
    this->sortId = sort_indexes(this->tf_CMU_EM);
    int ind;

    cout << "--------------------------------------" << endl;
    cout << "The min and max time of flights are:" << endl;
    ind = this->sortId[0];
    cout << "min(tof_EM)  = " << this->tf_CMU_EM[ind]/SEML.us->T - this->t0_CMU_EM[ind]/SEML.us->T  << " x T" << endl;
    ind = this->sortId[this->csize -1];
    cout << "max(tof_EM)  = " << this->tf_CMU_EM[ind]/SEML.us->T - this->t0_CMU_EM[ind]/SEML.us->T  << " x T" << endl;

    //====================================================================================
    // Sort data wrt to the projection distance
    //====================================================================================
    this->sortId = sort_indexes(this->pmin_dist_SEM);
    coutlp();

    return FTC_SUCCESS;
}

/**
 *  \brief Update sortId, the vector of sorted indices, with respect to pmin_dist_SEM.
 *         After this routine, sortId[0] is the index for which pmin_dist_SEM is minimum.
 **/
void ProjResClass::sort_pmin_dist_SEM()
{
    if(this->csize > 0) this->sortId = sort_indexes(this->pmin_dist_SEM);
    else cout << "sort_pmin_dist_SEM. Warning: data is empty. No sorting is performed." << endl;
}

/**
 *  \brief Update st_EM, st_SEM... With the elements contained in sortId[k].
 **/
void ProjResClass::update_ic(double st_EM[5], double st_SEM[5], double t_EM[2],
                             double st_EM_seed[4], double *t_EM_seed,
                             double* t0_SEM, double* pmin_dist_SEM_out,
                             int *label, int k)
{
    //====================================================================================
    //Check
    //====================================================================================
    if(k > this->csize)
    {
        cout << "update_ic. Data length is smaller than the desired index." << endl;
        return;
    }

    int kpos = this->sortId[k];

    //====================================================================================
    // Initialize local variables: EM
    //====================================================================================
    //RCM coordinates
    st_EM[0] = this->s1_CMU_EM[kpos];
    st_EM[1] = this->s2_CMU_EM[kpos];
    st_EM[2] = this->s3_CMU_EM[kpos];
    st_EM[3] = this->s4_CMU_EM[kpos];
    st_EM[4] = this->s5_CMU_EM[kpos];
    //Time
    t_EM[0]  = this->t0_CMU_EM[kpos];
    t_EM[1]  = this->tf_CMU_EM[kpos];


    //====================================================================================
    // Initialize local variables: SEM
    //====================================================================================
    //RCM coordinates
    st_SEM[0] = this->s1_CM_SEM[kpos];
    st_SEM[1] = this->s2_CM_SEM[kpos];
    st_SEM[2] = this->s3_CM_SEM[kpos];
    st_SEM[3] = this->s4_CM_SEM[kpos];
    st_SEM[4] = 0.0;
    //Initial time in SEM units
    *t0_SEM = this->tf_CMU_EM[kpos]*SEML.us_em.ns;

    //====================================================================================
    // Minimum projection distance
    //====================================================================================
    *pmin_dist_SEM_out = this->pmin_dist_SEM[kpos];

    //====================================================================================
    // Seed
    //====================================================================================
    //RCM coordinates
    st_EM_seed[0] = this->s1_CMU_EM_seed[kpos];
    st_EM_seed[1] = this->s2_CMU_EM_seed[kpos];
    st_EM_seed[2] = this->s3_CMU_EM_seed[kpos];
    st_EM_seed[3] = this->s4_CMU_EM_seed[kpos];
    //Time
    *t_EM_seed    = this->t0_CMU_EM_seed[kpos];

    //====================================================================================
    // Label
    //====================================================================================
    *label = this->label[kpos];
}

/**
 *  \brief Update st_EM, st_SEM... With the elements contained in sortId[k].
 **/
void ProjResClass::update_ic(double st_EM[5], double st_SEM[5], double t_EM[2],
                             double* t0_SEM, double* pmin_dist_SEM_out,
                             int k)
{
    //====================================================================================
    //Check
    //====================================================================================
    if(k > this->csize)
    {
        cout << "update_ic. Data length is smaller than the desired index." << endl;
        return;
    }

    int kpos = this->sortId[k];

    //====================================================================================
    // Initialize local variables: EM
    //====================================================================================
    //RCM coordinates
    st_EM[0] = this->s1_CMU_EM[kpos];
    st_EM[1] = this->s2_CMU_EM[kpos];
    st_EM[2] = this->s3_CMU_EM[kpos];
    st_EM[3] = this->s4_CMU_EM[kpos];
    st_EM[4] = this->s5_CMU_EM[kpos];
    //Time
    t_EM[0]  = this->t0_CMU_EM[kpos];
    t_EM[1]  = this->tf_CMU_EM[kpos];


    //====================================================================================
    // Initialize local variables: SEM
    //====================================================================================
    //RCM coordinates
    st_SEM[0] = this->s1_CM_SEM[kpos];
    st_SEM[1] = this->s2_CM_SEM[kpos];
    st_SEM[2] = this->s3_CM_SEM[kpos];
    st_SEM[3] = this->s4_CM_SEM[kpos];
    st_SEM[4] = 0.0;
    //Initial time in SEM units
    *t0_SEM = this->tf_CMU_EM[kpos]*SEML.us_em.ns;

    //====================================================================================
    // Minimum projection distance
    //====================================================================================
    *pmin_dist_SEM_out = this->pmin_dist_SEM[kpos];
}


