#include "oolencon.h"


//========================================================================================
//
//          I/O (to set in oolencon_io.cpp)
//
//========================================================================================
/**
 *   \brief Storing the results of the continuation procedure, in txt file.
 **/
void writeCONT_txt(Orbit& orbit_EM, Orbit& orbit_SEM, double* te_NCSEM, double* ye_NCSEM,  int isFirst)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename  = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, orbit_EM.getT0()/SEML.us_em.T);
    fstream filestream;

    if(isFirst)
    {
        //====================================================================================
        // If it is the first entry, the title of the columns are written
        //====================================================================================
        filestream.open (filename.c_str(), ios::out);
        //Title
        filestream << "t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
        filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
        filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
        filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
        filestream << "te_NCSEM xe_CMS_NCSEM ye_CMS_NCSEM ze_CMS_NCSEM pxe_CMS_NCSEM pye_CMS_NCSEM pze_CMS_NCSEM ";
        filestream << endl;
    }
    else
    {
        //====================================================================================
        // Else, we append
        //====================================================================================
        filestream.open (filename.c_str(), ios::out | ios::app);
    }

    //====================================================================================
    //Data storage
    //====================================================================================
    filestream << setprecision(15) <<  setiosflags(ios::scientific) << std::showpos;

    //1. t0 in EM units
    filestream << orbit_EM.getT0() << "  ";
    //2. s0 in RCM coordinates
    for(int i = 0; i <5; i++) filestream << orbit_EM.getSi()[i]  << "  ";
    //3. tf in EM units
    filestream << orbit_EM.getTf() << "  ";
    //4. sf in RCM coordinates
    for(int i = 0; i <5; i++) filestream << orbit_SEM.getSi()[i] << "  ";
    //5. z0 in NCEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_EM.getZ0()[i]  << "  ";
    //6. zf in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << orbit_SEM.getZ0()[i] << "  ";
    //7. thetae_NCSEM
    filestream << *te_NCSEM* SEML.us_sem.n << "  ";
    //8. ze in NCSEM coordinates
    for(int i = 0; i <6; i++) filestream << ye_NCSEM[i] << "  ";
    filestream << endl;

    filestream.close();
}

/**
 *  \brief Get the length the results of the continuation procedure, in txt file.
 **/
int getLengthCONT_txt(double t0)
{
    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename  = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, t0);
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
                 double tr0, int fsize)
{
    //====================================================================================
    // Get the size of the file
    //====================================================================================
    int fsize0 = getLengthCONT_txt(tr0);
    if(fsize0 != fsize)
    {
        cerr << "readCONT_txt. The user-defined file size mismatch the true size." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Initialize the I/O objects
    //====================================================================================
    string filename  = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, tr0);
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
    // First value is discarded: these are the column titles
    //====================================================================================
    string ct;
    getline(filestream,  ct);

    //====================================================================================
    // Read the data on each line and close
    //====================================================================================
    for(int k = 0; k < fsize; k++)
    {
        //1. t0 in EM units
        filestream >> t0_CMU_EM[k];
        //2. s0 in RCM coordinates
        for(int i = 0; i <5; i++) filestream >> si_CMU_EM[i][k];
        //3. tf in EM units
        filestream >> tf_CMU_EM[k];
        //4. sf in RCM coordinates
        for(int i = 0; i <5; i++) filestream >> si_CMS_SEM[i][k];
        //5. z0 in NCEM coordinates
        for(int i = 0; i <6; i++) filestream >> z0_CMU_NCEM[i][k];
        //6. zf in NCEM coordinates
        for(int i = 0; i <6; i++) filestream >> z0_CMS_NCSEM[i][k];
        //7. tethae
        filestream >> tethae[k];
        //8. ze in NCSEM coordinates
        for(int i = 0; i <6; i++) filestream >> ye_NCSEM[i][k];
    }

    filestream.close();

    return FTC_SUCCESS;
}


//========================================================================================
//
//          Computation of the CMU about EML2
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a 5-dimensional
 *         box: one dimension for the starting time, four dimensions for the
 *         parameterization of the Center Manifold (s1 to s4 coordinates). The RCM
 *         coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:      the value in RCM coordinates on the unstable direction s5.
 *  \param projSt.TLIM:     the min/max starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:    the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI:  the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                          in the IC box.
 *  \param projSt.GSIZE_SI: the number of points on the  s1, s2, s3, s4 values  grids
 *                          in the IC box.
 *  \param projSt.ISPAR;    if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/EM/L2/cu_3d_order_16.bin"
 **/
int oo_compute_grid_CMU_EM_3D(double dist_to_cm, ProjSt &projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Get the invariant manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs);

    //====================================================================================
    // Splash screen
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);
    cout << resetiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;
    cout << "       Computation of center-unstable 3D IC at EML2            " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << projSt.GSIZE_SI[0]+1  << " value(s) of s1 in [" << projSt.GLIM_SI[0][0];
    cout << ", " << projSt.GLIM_SI[0][1] << "]" << endl;
    cout << "  - " << projSt.GSIZE_SI[1]+1  << " value(s) of s2 in [" << projSt.GLIM_SI[1][0];
    cout << ", " << projSt.GLIM_SI[1][1] << "]" << endl;
    cout << "  - " << projSt.GSIZE_SI[2]+1  << " value(s) of s3 in [" << projSt.GLIM_SI[2][0];
    cout << ", " << projSt.GLIM_SI[2][1] << "]" << endl;
    cout << "  - " << projSt.GSIZE_SI[3]+1  << " value(s) of s4 in [" << projSt.GLIM_SI[3][0];
    cout << ", " << projSt.GLIM_SI[3][1] << "]" << endl;
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us.T;
    cout << ", " << projSt.TLIM[1]/SEML.us.T << "] x T" << endl;
    cout << " The data will be stored in " << filename << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // Check that the invman is an unstable-manifold
    //====================================================================================
    if(invman.getManType() != MAN_CENTER_U)
    {
        cout << "oo_compute_grid_CMU_EM_3D. The invariant manifold must be of center-unstable type. return." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------
    //Building the working grids
    //------------------------------------------
    double** grid_si_CMU_RCM = (double**) calloc(4, sizeof(double*));
    for(int i = 0; i <4; i++)
    {
        grid_si_CMU_RCM[i] = (double*) calloc(projSt.GSIZE_SI[i]+1, sizeof(double));
        init_grid(grid_si_CMU_RCM[i], projSt.GLIM_SI[i][0], projSt.GLIM_SI[i][1], projSt.GSIZE_SI[i]);
    }


    //------------------------------------------
    //Building the time grid
    //------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.GSIZE_SI[1])*(1+projSt.GSIZE_SI[2])*(1+projSt.GSIZE_SI[3])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------
    // Data structures
    //------------------------------------------
    double** init_state_CMU_NCEM = dmatrix(0, 5, 0, projSt.GSIZE_SI[2]);
    double** init_state_CMU_RCM  = dmatrix(0, 4, 0, projSt.GSIZE_SI[2]);

    //------------------------------------------
    // Reset the data file
    //------------------------------------------
    initCU_bin_3D(projSt.GSIZE_SI, projSt.TSIZE, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);

    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= projSt.TSIZE; kt++)
    {
        //----------------------
        //Append the time in data file
        //----------------------
        appTimeCU_bin_3D(grid_t_EM, kt, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);


        for(int ks2 = 0; ks2 <= projSt.GSIZE_SI[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= projSt.GSIZE_SI[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= projSt.GSIZE_SI[0]; ks1++)
                {
                    #pragma omp parallel for if(isPar)  shared(iter)
                    for(int ks3 = 0; ks3 <= projSt.GSIZE_SI[2]; ks3++)
                    {
                        Ofsc ofs(OFS_ORDER);
                        double* yvu = dvector(0,5);
                        double* sti = dvector(0,4);

                        //----------------------
                        // Initialization on the center-unstable manifold
                        //----------------------
                        //Init sti
                        sti[0] = grid_si_CMU_RCM[0][ks1];
                        sti[1] = grid_si_CMU_RCM[1][ks2];
                        sti[2] = grid_si_CMU_RCM[2][ks3];
                        sti[3] = grid_si_CMU_RCM[3][ks4];
                        sti[4] = dist_to_cm;

                        //Equivalent state
                        invman.evalRCMtoNC(sti, grid_t_EM[kt], yvu, OFTS_ORDER, OFS_ORDER);

                        //Save
                        #pragma omp critical
                        {

                            for(int i = 0; i < 6; i++) init_state_CMU_NCEM[i][ks3] = yvu[i];
                            for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][ks3] = sti[i];
                            //Display
                            displayCompletion("oo_compute_grid_CMU_EM", (double) iter++/noe*100);
                        }

                        free_dvector(yvu, 0, 5);
                        free_dvector(sti, 0, 4);
                    }

                    //------------------------------------------
                    //Store values
                    //------------------------------------------
                    writeCU_bin_3D(init_state_CMU_NCEM, init_state_CMU_RCM, projSt.GSIZE_SI,
                                   OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);
                }
            }
        }
    }

    //------------------------------------------
    //Free
    //------------------------------------------
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, projSt.GSIZE_SI[2]);
    free_dmatrix(init_state_CMU_RCM,  0, 4, 0, projSt.GSIZE_SI[2]);

    return FTC_SUCCESS;
}

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         three-dimensional box: one dimension for the starting time, two dimensions for
 *         the parameterization of the Center Manifold (s1 and s3 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:           the value in RCM coordinates on the unst. direction s5.
 *  \param projSt.TLIM[0]:       the minimum starting time (in EM units) in the IC box.
 *  \param projSt.TLIM[1]:       the maximum starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:         the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI[0][0]: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[0][1]: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[2][0]: the minimum s3 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[2][1]: the maximum s3 value (in RCM coordinates) in the IC box.
 *  \param projSt.GSIZE_SI[0]:   the number of points on the s1 grid in the IC box.
 *  \param projSt.GSIZE_SI[2]:   the number of points on the s3 grid in the IC box.
 *  \param projSt.ISPAR:         if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int oo_compute_grid_CMU_EM(double dist_to_cm, ProjSt &projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Splash screen
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_CU, SEML.li_SEM);
    cout << resetiosflags(ios::scientific) << setprecision(15);

    cout << "===================================================================" << endl;
    cout << "       Computation of center-unstable planar IC at EML2            " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << projSt.GSIZE_SI[0]+1  << " value(s) of s1 in [" << projSt.GLIM_SI[0][0];
    cout << ", " << projSt.GLIM_SI[0][1] << "]" << endl;
    cout << "  - " << projSt.GSIZE_SI[2]+1  << " value(s) of s3 in [" << projSt.GLIM_SI[2][0];
    cout << ", " << projSt.GLIM_SI[2][1] << "]" << endl;
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us.T;
    cout << ", " << projSt.TLIM[1]/SEML.us.T << "] x T" << endl;
    cout << " The data will be stored in " << filename << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // Get the invariant manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs);

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------
    //Building the working grids
    //------------------------------------------
    double* grid_s1_CMU_RCM = dvector(0,  projSt.GSIZE_SI[0]);
    double* grid_s3_CMU_RCM = dvector(0,  projSt.GSIZE_SI[2]);
    init_grid(grid_s1_CMU_RCM, projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1], projSt.GSIZE_SI[0]);
    init_grid(grid_s3_CMU_RCM, projSt.GLIM_SI[2][0], projSt.GLIM_SI[2][1], projSt.GSIZE_SI[2]);

    //------------------------------------------
    //Building the time grid
    //------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.GSIZE_SI[2])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------
    // Data structures
    //------------------------------------------
    double**** init_state_CMU_NCEM = d4tensor(0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    double**** init_state_CMU_RCM  = d4tensor(0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);

    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= projSt.TSIZE; kt++)
    {
        for(int ks1 = 0; ks1 <= projSt.GSIZE_SI[0]; ks1++)
        {
            #pragma omp parallel for if(isPar)  shared(iter)
            for(int ks3 = 0; ks3 <= projSt.GSIZE_SI[2]; ks3++)
            {
                Ofsc ofs(OFS_ORDER);
                double* yvu = dvector(0,5);
                double* sti = dvector(0,4);

                //----------------------
                // Initialization on the center-unstable manifold
                //----------------------
                //Init sti
                sti[0] = grid_s1_CMU_RCM[ks1];
                sti[1] = 0.0;
                sti[2] = grid_s3_CMU_RCM[ks3];
                sti[3] = 0.0;
                sti[4] = dist_to_cm;

                //Equivalent state
                invman.evalRCMtoNC(sti, grid_t_EM[kt], yvu, OFTS_ORDER, OFS_ORDER);

                //Save
                #pragma omp critical
                {

                    for(int i = 0; i < 6; i++) init_state_CMU_NCEM[i][kt][ks1][ks3] = yvu[i];
                    for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][kt][ks1][ks3] = sti[i];
                    //Display
                    displayCompletion("oo_compute_grid_CMU_EM", (double) iter++/noe*100);
                }

                free_dvector(yvu, 0, 5);
                free_dvector(sti, 0, 4);
            }
        }
    }

    //------------------------------------------
    //Store values
    //------------------------------------------
    writeCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, grid_t_EM,
                projSt.GSIZE_SI[0], projSt.GSIZE_SI[2], projSt.TSIZE, OFTS_ORDER,
                TYPE_CU, SEML.li_SEM);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    free_d4tensor(init_state_CMU_RCM,  0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);

    return FTC_SUCCESS;
}


//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEMLi
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM. Then, each point on the
 *         integration grid is projected on the Center Manifold CM_SEM_NC about SEMLi.
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM: the maximum integration time on each leg, in EM units.
 *
 *  \param projSt.MSIZE:       the number of points on each manifold leg.
 *
 *  \param projSt.NOD:         the number of dimensions on which the distance of
 *                             projection is computed (usually either 3 (the physical
 *                             distance) or 6 (the whole phase space)).
 *
 *  \param projSt.ISPAR:       if TRUE, the computation is parallelized.
 *
 *  \param projSt.YNMAX:       the maximum norm in NCSEM coordinates for which a given
 *                             state on the integration grid is projected on CM_SEM_NC
 *                             More precisely: for a given state y along the manifold leg,
 *                             if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                             Otherwise, it is considered too far away from SEMLi to be
 *                             a good candidate for projection.
 *
 *  \param projSt.SNMAX:       the maximum norm in RCM SEM coordinates for which a given
 *                             projection state on the CM of SEMLi (CM_SEM_NC) is
 *                             computed back in NCSEM coordinates. More precisely, for a
 *                             given state y in NCSEM coordinates, the result of the
 *                             projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                             coordinates.
 *                             if norm(sproj, 4) < projSt.SNMAX, the computation
 *                             yproj = CM_SEM_NC(sproj, t) is performed. Otherwise, the
 *                             state sproj is considered too far away from the RCM origin
 *                             to be a good candidate - it is out of the domain of
 *                             practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form:
 *          "plot/QBCP/EM/L2/projcu_3d_order_16.bin".
 **/
int oo_int_proj_CMU_EM_on_CM_SEM_3D(ProjSt &projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, si_grid_size[4], offset;
    offset = getLenghtCU_bin_3D(si_grid_size, &t_grid_size,
                                OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);


    //====================================================================================
    // Splash screen
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_CU, SEML.li_SEM);
    cout << resetiosflags(ios::scientific) << setprecision(15);

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double** init_state_CMU_NCEM = dmatrix(0, 5, 0, si_grid_size[2]);
    double** init_state_CMU_RCM  = dmatrix(0, 4, 0, si_grid_size[2]);
    double* init_time_grid_EM    = dvector(0, t_grid_size);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = (1+si_grid_size[0])*(1+si_grid_size[1])*(1+si_grid_size[2])*(1+si_grid_size[3])*(1+t_grid_size);


    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs.manType != MAN_CENTER)
    {
        cout << "oo_int_proj_CMU_EM_on_CM_SEM_3D. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);


    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    double** final_state_CMU_SEM      = dmatrix(0, 5, 0, si_grid_size[2]);
    double** projected_state_CMU_SEM  = dmatrix(0, 5, 0, si_grid_size[2]);
    double** projected_state_CMU_RCM  = dmatrix(0, 3, 0, si_grid_size[2]);
    double** init_state_CMU_SEM       = dmatrix(0, 5, 0, si_grid_size[2]);
    double** init_state_CMU_NCEM_0    = dmatrix(0, 5, 0, si_grid_size[2]);
    double* min_proj_dist_tens_SEM    = dvector(0, si_grid_size[2]);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================
    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //----------------------------------------------------------
    // projection by default is arbitrary big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Filename
    filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ_3D, SEML.li_SEM);
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    filestream.close();


    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    int index  = 0;
    for(int kt = 0; kt <= t_grid_size; kt++)
    {
        //--------------------------------------------------------------------------------
        //Read time from file
        //--------------------------------------------------------------------------------
        offset = readTCU_bin_3D(offset, init_time_grid_EM, kt,
                                OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);

        for(int ks2 = 0; ks2 <= si_grid_size[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= si_grid_size[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= si_grid_size[0]; ks1++)
                {
                    //--------------------------------------------------------------------
                    //Read data from file
                    //--------------------------------------------------------------------
                    offset = readCU_bin_3D(offset, init_state_CMU_NCEM,
                                           init_state_CMU_RCM, si_grid_size, OFTS_ORDER,
                                           TYPE_CU_3D, SEML.li_SEM);

                    //--------------------------------------------------------------------
                    //Most inner loop is parallelized
                    //--------------------------------------------------------------------
                    #pragma omp parallel for if(isPar)  shared(index)
                    for(int ks3 = 0; ks3 <= si_grid_size[2]; ks3++)
                    {

                        //================================================================
                        // 3.1. Integration of the manifold leg
                        //================================================================
                        //----------------------------------------------------------------
                        //Initialize the initial conditions (both NC and RCM coordinates)
                        //----------------------------------------------------------------
                        double yv[6], tv;
                        tv  = init_time_grid_EM[kt];
                        for(int i = 0; i < 6; i++) yv[i] = init_state_CMU_NCEM[i][ks3];

                        //----------------------------------------------------------------
                        //Local variables to store the manifold leg
                        //----------------------------------------------------------------
                        double** y_man_NCSEM = dmatrix(0, 5, 0, projSt.MSIZE);
                        double* t_man_SEM    = dvector(0, projSt.MSIZE);


                        //----------------------------------------------------------------
                        //Integration on projSt.MSIZE+1 fixed grid
                        // PB: when Release + I_NCSEM!!
                        //----------------------------------------------------------------
                        int ode78coll;
                        int status = ode78(y_man_NCSEM, t_man_SEM, &ode78coll,
                                           tv, tv+projSt.TM, yv, 6,
                                           projSt.MSIZE, I_NCEM, NCEM, NCSEM);

                        //================================================================
                        // 3.2. Projection on the center manifold of SEMLi.
                        //      No use of SEML_EM or SEML after this point!
                        //      We need first to check that the integration went well
                        //================================================================
                        //----------------------------------------------------------
                        //Temp variables
                        //----------------------------------------------------------
                        Ofsc ofs(OFS_ORDER);
                        double yvproj_NCSEM[6], sproj[4], yv_SEM[6], yvproj_SEM[6], yvEM[6], yv_VSEM[6], yvproj_VSEM[6];
                        double proj_dist_SEM, min_proj_dist_SEM = ePdef, dv_at_projection_SEM = 0.0, y_man_norm_NCSEM = 0.0;
                        int kmin = 0;

                        //----------------------------------------------------------------
                        // We need first to check that the integration went well
                        //---------------------------------------------------------------
                        if(status)
                        {
                            // If status != 0, something went wrong during the integration
                            // in ode78 and we use the maximum value ePdef (the solution
                            // is basically discarded
                            proj_dist_SEM = ePdef;
                        }
                        else
                        {
                            //------------------------------------------------------------
                            //Loop on trajectory
                            //------------------------------------------------------------
                            for(int kman = 0; kman <= projSt.MSIZE; kman++)
                            {
                                //Current state
                                for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
                                tv = t_man_SEM[kman];

                                //Current distance from SEMLi in NCSEM units
                                y_man_norm_NCSEM = 0.0;
                                for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
                                y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

                                //--------------------------------------------------------
                                //Check 1: the current state is close enough to SEMLi
                                //--------------------------------------------------------
                                if(y_man_norm_NCSEM < projSt.YNMAX)
                                {
                                    // Projection on the center manifold
                                    invman_SEM.NCprojCCMtoCM(yv, tv, sproj);

                                    //----------------------------------------------------
                                    //Check 2: the projection is close enough to SEMLi
                                    //----------------------------------------------------
                                    if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< projSt.SNMAX)
                                    {
                                        //yvproj_NCSEM = W(sproj, tv)
                                        invman_SEM.evalRCMtoNC(sproj, tv, yvproj_NCSEM, OFTS_ORDER, OFS_ORDER);

                                        //Back to SEM coordinates
                                        NCSEMmtoSEMm(tv, yv, yv_SEM, &SEML_SEM);
                                        NCSEMmtoSEMm(tv, yvproj_NCSEM, yvproj_SEM, &SEML_SEM);

                                        //Back to SEM coordinates with velocity
                                        SEMmtoSEMv(tv, yv_SEM, yv_VSEM, &SEML_SEM);
                                        SEMmtoSEMv(tv, yvproj_SEM, yvproj_VSEM, &SEML_SEM);

                                        //Distance of projection in SEM coordinates
                                        proj_dist_SEM = 0.0;
                                        for(int i = 0; i < projSt.NOD; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
                                        proj_dist_SEM = sqrt(proj_dist_SEM);
                                    }
                                    else proj_dist_SEM = ePdef;
                                }
                                else proj_dist_SEM = ePdef;

                                //Update distance min if necessary
                                if(proj_dist_SEM < min_proj_dist_SEM)
                                {
                                    min_proj_dist_SEM = proj_dist_SEM;
                                    kmin = kman;
                                    for(int i = 0; i < 6; i++) final_state_CMU_SEM[i][ks3]     = yv_VSEM[i];
                                    for(int i = 0; i < 6; i++) projected_state_CMU_SEM[i][ks3] = yvproj_VSEM[i];
                                    for(int i = 0; i < 4; i++) projected_state_CMU_RCM[i][ks3] = sproj[i];

                                    //Associated DV
                                    dv_at_projection_SEM = 0.0;
                                    for(int i = 3; i < 6; i++) dv_at_projection_SEM += (yvproj_VSEM[i] - yv_VSEM[i])*(yvproj_VSEM[i] - yv_VSEM[i]);
                                    dv_at_projection_SEM = sqrt(dv_at_projection_SEM);
                                }
                            }

                            //------------------------------------------------------------
                            //We check for collisions
                            //------------------------------------------------------------
                            if(ode78coll)
                            {
                                //If ode78coll !=  0 a collision has occured, and we
                                //send back the code of the associated primary in the
                                //dv_at_projection_SEM.
                                //cout << "Collision! with " << ode78coll << endl;
                                dv_at_projection_SEM = ode78coll;
                            }


                        }


                        //================================================================
                        // 3.3. Save outputs
                        //================================================================
                        if(min_proj_dist_SEM < ePdef)
                        {
                            #pragma omp critical
                            {
                                //Initial position in SEM coordinates
                                for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][0];
                                NCSEMmtoSEMm(t_man_SEM[0], yv, yv_SEM, &SEML_SEM);
                                for(int i = 0; i < 6; i++) init_state_CMU_SEM[i][ks3] = yv_SEM[i];

                                //Same in NCEM coordinates
                                NCSEMmtoNCEMm(t_man_SEM[0], yv, yvEM, &SEML_SEM);
                                for(int i = 0; i < 6; i++) init_state_CMU_NCEM_0[i][ks3] = yvEM[i];

                                //Minimum projection distance
                                min_proj_dist_tens_SEM[ks3] = min_proj_dist_SEM;

                                //----------------------------------------------------------
                                //Open datafile
                                //----------------------------------------------------------
                                writeIntProjCU_bin_3D(filename, init_time_grid_EM, init_state_CMU_NCEM,
                                init_state_CMU_SEM, init_state_CMU_RCM, final_state_CMU_SEM, projected_state_CMU_SEM,
                                projected_state_CMU_RCM, min_proj_dist_SEM, dv_at_projection_SEM, t_man_SEM, kmin, ks3, kt);
                            }
                        }

                        //----------------------------------------------------------
                        //Display completion
                        //----------------------------------------------------------
                        #pragma omp critical
                        {
                            displayCompletion("oo_int_proj_CMU_EM_on_CM_SEM_3D", 100.0*index++/noe);
                        }

                        //----------------------------------------------------------------
                        //Free
                        //----------------------------------------------------------------
                        //free_dmatrix(y_man_NCEM, 0, 5, 0, projSt.MSIZE);
                        free_dmatrix(y_man_NCSEM, 0, 5, 0, projSt.MSIZE);
                        //free_dvector(t_man_EM, 0, projSt.MSIZE);
                        free_dvector(t_man_SEM, 0, projSt.MSIZE);


                    }
                }
            }
        }
    }
    return FTC_SUCCESS;
}


/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM.
 *         Then, each point on the integration grid is projected on the center Manifold
 *         about SEMLi (denoted here CM_SEM_NC).
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM:......the maximum integration time on each leg, in EM units.
 *  \param projSt.MSIZE:...the number of points on each manifold leg.
 *  \param projSt.NSMIN:...the number of best solutions that are kept
 *  \param projSt.NOD:.....the number of dimensions on which the distance of
 *                         projection is computed - usually either 3
 *                         (the physical distance) or 6 (the whole phase space).
 *  \param projSt.ISPAR:...if TRUE, the computation is parallelized.
 *  \param projSt.YNMAX:...the maximum norm in NCSEM coordinates for which a given
 *                         state on the integration grid is projected on CM_SEM_NC.
 *                         More precisely: for a given state y along the manifold leg,
 *                         if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                         Otherwise, it is considered too far away from SEMLi to
 *                         be a good candidate for projection.
 *  \param projSt.SNMAX:...the maximum norm in RCM SEM coordinates for which a given
 *                         projection state on the CM of SEMLi (CM_SEM_NC) is
 *                         computed back in NCSEM coordinates. More precisely, for a
 *                         given state y in NCSEM coordinates, the result of the
 *                         projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                         coordinates. If norm(sproj, 4) < projSt.SNMAX, the computation
 *                         yproj = CM_SEM_NC(sproj, t) is performed. Otherwise,
 *                         the state sproj is considered too far away from the RCM
 *                         origin to be a good candidate - it is out of the domain of
 *                         practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form
 * "plot/QBCP/EM/L2/projcu_order_16.bin", and
 * "plot/QBCP/EM/L2/sortprojcu_order_16.bin" for the projSt.NSMIN best solutions.
 **/
int oo_int_proj_CMU_EM_on_CM_SEM(ProjSt &projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, s1_grid_size, s3_grid_size;
    getLenghtCU_bin(&s1_grid_size, &s3_grid_size, &t_grid_size,
                    OFTS_ORDER, TYPE_CU, SEML.li_SEM);

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double**** init_state_CMU_NCEM = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** init_state_CMU_RCM  = d4tensor(0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double* init_time_grid_EM      = dvector(0, t_grid_size);

    //----------------------------------------------------------
    //Read data from file
    //----------------------------------------------------------
    readCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, init_time_grid_EM,
               s1_grid_size, s3_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU, SEML.li_SEM);

    //====================================================================================
    // Splash screen
    //====================================================================================
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ, SEML.li_SEM);
    cout << "===================================================================" << endl;
    cout << "              Computation of the connections between:              " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << s1_grid_size+1  << " value(s) of s1" << endl;
    cout << "  - " << s3_grid_size+1  << " value(s) of s3" << endl;
    cout << "  - " << t_grid_size+1   << " value(s) of t"  << endl;
    cout << " The data will be stored in " << filename << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // 2.1. Initialize tools for the sorting phase
    //====================================================================================
    int Nsort = (s1_grid_size+1)*(s3_grid_size+1)*(t_grid_size+1);
    vector<int>    indexMin(Nsort);
    vector<int>    ks1Min(Nsort);
    vector<int>    ks3Min(Nsort);
    vector<int>    ktMin(Nsort);
    vector<double> t0_min_EM(Nsort);
    vector<double> tf_min_EM(Nsort);
    vector<double> distMin(Nsort);

    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs.manType != MAN_CENTER)
    {
        cout << "oo_int_proj_CMU_EM_on_CM_SEM. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);


    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    double**** final_state_CMU_SEM      = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** projected_state_CMU_SEM  = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** projected_state_CMU_RCM  = d4tensor(0, 3, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** init_state_CMU_SEM       = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** init_state_CMU_NCEM_0    = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double** *min_proj_dist_tens_SEM    = d3tensor(0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //----------------------------------------------------------
    // projection by default is arbitrary big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Filename
    filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ, SEML.li_SEM);
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    filestream.close();

    //====================================================================================
    // 3. Loop: only the inner loop is parallelized, since
    //    open_mp doest not allow nested // loops by default
    //====================================================================================
    COMPLETION = 0;
    int index  = 0;
    for(int kt = 0; kt <= t_grid_size; kt++)
    {
        for(int ks1 = 0; ks1 <= s1_grid_size; ks1++)
        {
            #pragma omp parallel for if(isPar) shared(index)
            for(int ks3 = 0; ks3 <= s3_grid_size; ks3++)
            {
                //========================================================================
                // 3.1. Integration of the manifold leg
                //========================================================================
                //------------------------------------------------------------------------
                //Initialize the initial conditions (both NC and RCM coordinates)
                //------------------------------------------------------------------------
                double yv[6], tv;
                tv  = init_time_grid_EM[kt];
                for(int i = 0; i < 6; i++) yv[i] = init_state_CMU_NCEM[i][kt][ks1][ks3];

                //------------------------------------------------------------------------
                //Local variables to store the manifold leg
                //------------------------------------------------------------------------
                double** y_man_NCSEM = dmatrix(0, 5, 0, projSt.MSIZE);
                double* t_man_SEM    = dvector(0, projSt.MSIZE);


                //------------------------------------------------------------------------
                //Integration on projSt.MSIZE+1 fixed grid
                // PB: when Release + I_NCSEM!! not used for now
                //------------------------------------------------------------------------
                int ode78coll;
                int status = ode78(y_man_NCSEM, t_man_SEM, &ode78coll, tv, tv+projSt.TM, yv, 6, projSt.MSIZE, I_NCEM, NCEM, NCSEM);

                //========================================================================
                // 3.2. Projection on the center manifold of SEMLi. No use of SEML_EM or SEML after this point!
                // We need first to check that the integration went well
                //========================================================================
                //------------------------------------------------------------------------
                //Temp variables
                //------------------------------------------------------------------------
                Ofsc ofs(OFS_ORDER);
                double yvproj_NCSEM[6], sproj[4], yv_SEM[6], yvproj_SEM[6], yvEM[6], yv_VSEM[6], yvproj_VSEM[6];
                double proj_dist_SEM, min_proj_dist_SEM = ePdef, dv_at_projection_SEM = 0.0, y_man_norm_NCSEM = 0.0;
                int ksort = (s1_grid_size+1)*(s3_grid_size+1)*kt + (s3_grid_size+1)*ks1 + ks3;
                int kmin = 0;

                //------------------------------------------------------------------------
                // We need first to check that the integration went well
                //------------------------------------------------------------------------
                if(status)
                {
                    // If status != 0, something went wrong during the integration
                    // in ode78 and we use the maximum value ePdef (the solution
                    // is basically discarded
                    proj_dist_SEM = ePdef;
                }
                else
                {

                    //--------------------------------------------------------------------
                    //Loop on trajectory
                    //--------------------------------------------------------------------
                    for(int kman = 0; kman <= projSt.MSIZE; kman++)
                    {
                        //Current state
                        for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
                        tv = t_man_SEM[kman];

                        //Current distance from SEMLi in NCSEM units
                        y_man_norm_NCSEM = 0.0;
                        for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
                        y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

                        //----------------------------------------------------------------
                        //Check n°1: the current state is close enough to SEMLi
                        //----------------------------------------------------------------
                        if(y_man_norm_NCSEM < projSt.YNMAX)
                        {
                            // Projection on the center manifold
                            invman_SEM.NCprojCCMtoCM(yv, tv, sproj);

                            //------------------------------------------------------------
                            //Check n°2: the projection is close enough to SEMLi
                            //------------------------------------------------------------
                            if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< projSt.SNMAX)
                            {
                                //yvproj_NCSEM = W(sproj, tv)
                                invman_SEM.evalRCMtoNC(sproj, tv, yvproj_NCSEM, OFTS_ORDER, OFS_ORDER);

                                //Back to SEM coordinates
                                NCSEMmtoSEMm(tv, yv, yv_SEM, &SEML_SEM);
                                NCSEMmtoSEMm(tv, yvproj_NCSEM, yvproj_SEM, &SEML_SEM);

                                //Back to SEM coordinates with velocity
                                SEMmtoSEMv(tv, yv_SEM, yv_VSEM, &SEML_SEM);
                                SEMmtoSEMv(tv, yvproj_SEM, yvproj_VSEM, &SEML_SEM);

                                //Distance of projection in SEM coordinates
                                proj_dist_SEM = 0.0;
                                for(int i = 0; i < projSt.NOD; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
                                proj_dist_SEM = sqrt(proj_dist_SEM);
                            }
                            else proj_dist_SEM = ePdef;
                        }
                        else proj_dist_SEM = ePdef;

                        //Update distance min if necessary
                        if(proj_dist_SEM < min_proj_dist_SEM)
                        {
                            min_proj_dist_SEM = proj_dist_SEM;
                            kmin = kman;
                            for(int i = 0; i < 6; i++) final_state_CMU_SEM[i][kt][ks1][ks3]     = yv_VSEM[i];
                            for(int i = 0; i < 6; i++) projected_state_CMU_SEM[i][kt][ks1][ks3] = yvproj_VSEM[i];
                            for(int i = 0; i < 4; i++) projected_state_CMU_RCM[i][kt][ks1][ks3] = sproj[i];

                            //Associated DV
                            dv_at_projection_SEM = 0.0;
                            for(int i = 3; i < 6; i++) dv_at_projection_SEM += (yvproj_VSEM[i] - yv_VSEM[i])*(yvproj_VSEM[i] - yv_VSEM[i]);
                            dv_at_projection_SEM = sqrt(dv_at_projection_SEM);
                        }

                    }

                    //--------------------------------------------------------------------
                    //We check for collisions
                    //--------------------------------------------------------------------
                    if(ode78coll)
                    {
                        //If ode78coll !=  0 a collision has occured, and we
                        //send back the code of the associated primary in the
                        //dv_at_projection_SEM.
                        //cout << "Collision! with " << ode78coll << endl;
                        dv_at_projection_SEM = ode78coll;
                    }
                }


                //========================================================================
                // 3.3. Save outputs
                //========================================================================
                if(min_proj_dist_SEM < ePdef)
                {
                    #pragma omp critical
                    {
                        //Initial position in SEM coordinates
                        for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][0];
                        NCSEMmtoSEMm(t_man_SEM[0], yv, yv_SEM, &SEML_SEM);
                        for(int i = 0; i < 6; i++) init_state_CMU_SEM[i][kt][ks1][ks3] = yv_SEM[i];

                        //Same in NCEM coordinates
                        NCSEMmtoNCEMm(t_man_SEM[0], yv, yvEM, &SEML_SEM);
                        for(int i = 0; i < 6; i++) init_state_CMU_NCEM_0[i][kt][ks1][ks3] = yvEM[i];

                        //Initial & final time in EM coordinates
                        t0_min_EM[ksort]= t_man_SEM[0]/SEML.us_em.ns;       //in EM coordinates
                        tf_min_EM[ksort]= t_man_SEM[kmin]/SEML.us_em.ns;    //in EM coordinates


                        //Minimum projection distance
                        min_proj_dist_tens_SEM[kt][ks1][ks3] = min_proj_dist_SEM;

                        //For sorting
                        indexMin[ksort] = (int) kmin;
                        ks1Min[ksort]   = (int) ks1;
                        ks3Min[ksort]   = (int) ks3;
                        ktMin[ksort]    = (int) kt;
                        distMin[ksort]  =  min_proj_dist_SEM;

                        //----------------------------------------------------------
                        //Open datafile
                        //----------------------------------------------------------
                        writeIntProjCU_bin(filename, init_time_grid_EM,
                        init_state_CMU_NCEM, init_state_CMU_SEM,
                        init_state_CMU_RCM, final_state_CMU_SEM,
                        projected_state_CMU_SEM,
                        projected_state_CMU_RCM,
                        min_proj_dist_SEM, dv_at_projection_SEM,
                        t_man_SEM, kmin, ks1, ks3, kt);
                    }
                }

                //----------------------------------------------------------
                //Display completion
                //----------------------------------------------------------
                #pragma omp critical
                {
                    displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/((1+s1_grid_size)*(1+s3_grid_size)*(1+t_grid_size)));
                }

                //------------------------------------------------------------------------
                //Free
                //------------------------------------------------------------------------
                free_dmatrix(y_man_NCSEM, 0, 5, 0, projSt.MSIZE);
                free_dvector(t_man_SEM, 0, projSt.MSIZE);
            }
        }
    }

    //====================================================================================
    // 4. Sorting the best solutions
    //====================================================================================
    vector<size_t> sortId = sort_indexes(distMin);

    //----------------------------------------------------------
    //Saving the projSt.NSMIN best results or all results if less than 50 have been computed
    //----------------------------------------------------------
    int number_of_sol  = min(projSt.NSMIN, Nsort-2);


    //---------------------
    //Filename
    //---------------------
    filename = filenameCUM(OFTS_ORDER, TYPE_MAN_SORT, SEML.li_SEM);

    //---------------------
    //Write sorted solutions
    //---------------------
    writeIntProjSortCU_bin(filename, init_state_CMU_NCEM, init_state_CMU_RCM,
                           final_state_CMU_SEM,  projected_state_CMU_SEM, projected_state_CMU_RCM,
                           min_proj_dist_tens_SEM, sortId, ktMin, ks1Min, ks3Min, t0_min_EM, tf_min_EM, distMin, number_of_sol);


    //----------------------------------------------------------
    //Free
    //----------------------------------------------------------
    free_d4tensor(final_state_CMU_SEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(projected_state_CMU_SEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(projected_state_CMU_RCM, 0, 3, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_SEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_NCEM_0, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d3tensor(min_proj_dist_tens_SEM, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);


    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_RCM, 0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_dvector(init_time_grid_EM, 0, t_grid_size);


    return FTC_SUCCESS;
}

//========================================================================================
//
//         PLOTTING TRAJECTORIES
//
//========================================================================================
/**
 *  \brief Plotting the notable points in the SEM system, in coord_type coordinates
 **/
int notablePoints_sem(gnuplot_ctrl* h2, int coord_type)
{
    //------------------------------------------------------------------------------------
    //Container
    //------------------------------------------------------------------------------------
    double** semP_coord = dmatrix(0, 6, 0, 2);

    //------------------------------------------------------------------------------------
    //Switch on the coord type
    //------------------------------------------------------------------------------------
    switch(coord_type)
    {
    case PSEM:
    case VSEM:
        semPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        break;
    case PEM:
    case VEM:
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        break;
    default:
        cerr << "notablePoints_sem. Unknown coord_type. " << endl;
        free_dmatrix(semP_coord, 0, 6, 0, 2);
        return FTC_EDOM;
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_dmatrix(semP_coord, 0, 6, 0, 2);

    return FTC_SUCCESS;
}

/**
 *  \brief Plotting the notable points in the SEM & EM systems in coord_type coordinates,
 *         as well as complementary coordinates.
 **/
int notablePoints(gnuplot_ctrl* h2, gnuplot_ctrl* h3, int coord_type)
{
    //------------------------------------------------------------------------------------
    //Container
    //------------------------------------------------------------------------------------
    double** semP_coord = dmatrix(0, 6, 0, 2);
    double** semP_comp  = dmatrix(0, 6, 0, 2);

    //------------------------------------------------------------------------------------
    //Switch on the coord type
    //------------------------------------------------------------------------------------
    switch(coord_type)
    {
    case PSEM:
    case VSEM:
    {
        semPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        //Comp
        emPoints(0.0, semP_comp);
        emPlot(h3, semP_comp);

        //h2 displays the system in SEM coordinates
        gnuplot_set_xlabel(h2, (char*) "Xsem");
        gnuplot_set_ylabel(h2, (char*) "Ysem");
        gnuplot_set_zlabel(h2, (char*) "Zsem");

        //h3 displays the system in EM coordinates
        gnuplot_set_xlabel(h3, (char*) "Xem");
        gnuplot_set_ylabel(h3, (char*) "Yem");
        gnuplot_set_zlabel(h3, (char*) "Zem");
        break;
    }
    case NCSEM:
    case VNCSEM:
    {
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        //Comp
        emNCPoints(0.0, semP_comp);
        emPlot(h3, semP_comp);

        //h2 displays the system in NCSEM coordinates
        gnuplot_set_xlabel(h2, (char*) "xsem");
        gnuplot_set_ylabel(h2, (char*) "ysem");
        gnuplot_set_zlabel(h2, (char*) "zsem");

        //h3 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h3, (char*) "xem");
        gnuplot_set_ylabel(h3, (char*) "yem");
        gnuplot_set_zlabel(h3, (char*) "zem");
        break;
    }
    case PEM:
    case VEM:
    {
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);

        //h3 displays the system in SEM coordinates
        gnuplot_set_xlabel(h3, (char*) "Xsem");
        gnuplot_set_ylabel(h3, (char*) "Ysem");
        gnuplot_set_zlabel(h3, (char*) "Zsem");

        //h2 displays the system in EM coordinates
        gnuplot_set_xlabel(h2, (char*) "Xem");
        gnuplot_set_ylabel(h2, (char*) "Yem");
        gnuplot_set_zlabel(h2, (char*) "Zem");
        break;
    }
    case NCEM:
    case VNCEM:
    {
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semNCPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);

        //h3 displays the system in NCSEM coordinates
        gnuplot_set_xlabel(h3, (char*) "xsem");
        gnuplot_set_ylabel(h3, (char*) "ysem");
        gnuplot_set_zlabel(h3, (char*) "zsem");

        //h2 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h2, (char*) "xem");
        gnuplot_set_ylabel(h2, (char*) "yem");
        gnuplot_set_zlabel(h2, (char*) "zem");
        break;
    }
    default:
        cerr << "notablePoints. Unknown coord_type. " << endl;
        free_dmatrix(semP_coord, 0, 6, 0, 2);
        free_dmatrix(semP_comp, 0, 6, 0, 2);
        return FTC_EDOM;

    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);

    return FTC_SUCCESS;
}


/**
 *  \brief Plot a trajectory, segment by segment, in to complementary coordinate systems.
 **/
int plottrajsegbyseg(double** y_traj, double* t_traj,
                     int final_index, int mPlot, int coord_int,
                     double et0,     double tsys0,
                     int coordsys1,  gnuplot_ctrl* h2,
                     int coordsys2,  gnuplot_ctrl* h3,
                     int color, string title)
{
    //------------------------------------------------------------------------------------
    //Initialize local variables
    //------------------------------------------------------------------------------------
    double yv[6];
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);
    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //Choosing the dcs
    int dcs_int  = default_coordinate_system(coord_int);

    //tsys0 in the second coordinate system
    double tsys0_comp = tsys0;
    if(eph_coord(coordsys1) != eph_coord(coordsys2))
    {
        switch(eph_coord(coordsys2))
        {
        case VEM:
            tsys0_comp /= SEML.us_em.ns;
            break;

        case VSEM:
            tsys0_comp *= SEML.us_em.ns;
            break;
        }
    }

    //------------------------------------------------------------------------------------
    //Trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    int status;
    int ode78coll;
    for(int k = 0; k < final_index; k++)
    {
        //--------------------------------------------------------------------------------
        //Integration segment by segment
        //--------------------------------------------------------------------------------
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        status = ode78(ymc_comp, tmc_comp, &ode78coll, t_traj[k], t_traj[k+1], yv, 6, mPlot, dcs_int, coord_int, coord_int);

        //--------------------------------------------------------------------------------
        //Checks and warnings, if necessary
        //--------------------------------------------------------------------------------
        if(status != FTC_SUCCESS)
        {
            //At this step (plot), a simple warning is issued
            cout << "plottrajsegbyseg. Warning: ode78 returned a flag." << endl;
        }

        if(ode78coll)
        {
            //At this step (plot), a simple warning is issued
            cout << "plottrajsegbyseg. Warning: a collision with a primary has occurred." << endl;
        }

        //--------------------------------------------------------------------------------
        //Back to coordsys1 coordinates.
        //--------------------------------------------------------------------------------
        switch(coord_int)
        {
            //----------------------------------------------------------------------------
            //If the computation was performed using JPL ephemerides, some specific COC
            //is needed.
            //----------------------------------------------------------------------------
        case VECLI:
            ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys1, et0, tsys0, eph_coord(coordsys1));
            break;

        case J2000:
            eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys1, et0, tsys0, eph_coord(coordsys1));
            break;

        case NJ2000:
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys1, et0, tsys0, eph_coord(coordsys1), SEML.ss);
            break;

        default:
            //----------------------------------------------------------------------------
            //The default case is QBCP to QBCP, and we can use the usual COC
            //----------------------------------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_int, coordsys1);
        }


        //--------------------------------------------------------------------------------
        //Store
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }

        //--------------------------------------------------------------------------------
        //Back to coordsys2 coordinates
        //--------------------------------------------------------------------------------
        switch(coord_int)
        {
            //----------------------------------------------------------------------------
            //If the computation was performed using JPL ephemerides, some specific COC
            //is needed.
            //----------------------------------------------------------------------------
        case VECLI:
            ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys2, et0, tsys0_comp, eph_coord(coordsys2));
            break;

        case J2000:
            eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys2, et0, tsys0_comp, eph_coord(coordsys2));
            break;

        case NJ2000:
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys2, et0, tsys0_comp, eph_coord(coordsys2), SEML.ss);
            break;

        default:
            //----------------------------------------------------------------------------
            //The default case is QBCP to QBCP, and we can use the usual COC
            //----------------------------------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_int, coordsys2);
        }

        //--------------------------------------------------------------------------------
        //Store
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_comp_v[i][k*mPlot + p] = ymc[i][p];
            tmc_comp_v[k*mPlot + p] = tmc[p];
        }
    }

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)title.c_str(), "lines", "1", "2", color);

    //Plot on h3
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)title.c_str(), "lines", "1", "2", color);

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);
    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v, 0, mPlot*final_index);

    return FTC_SUCCESS;
}

//========================================================================================
//
//         REFINING TRAJECTORIES
//
//========================================================================================
int ooconteml2seml(RefSt& refst)
{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int coord_type = refst.coord_type;
    int status = 0;

    //====================================================================================
    //         Initialization of the containers
    //====================================================================================
    //Name
    string fname = "ooconteml2seml";

    // Get size of the data file
    double tr0 = refst.t0_des/SEML.us_em.T;
    int fsize = getLengthCONT_txt(tr0);

    // Init the vectors
    double* t0_CMU_EM     = dvector(0, fsize-1);
    double* tf_CMU_EM     = dvector(0, fsize-1);
    double** si_CMU_EM    = dmatrix(0, 4, 0, fsize-1);
    double** si_CMS_SEM   = dmatrix(0, 4, 0, fsize-1);
    double** z0_CMU_NCEM  = dmatrix(0, 5, 0, fsize-1);
    double** z0_CMS_NCSEM = dmatrix(0, 5, 0, fsize-1);
    double* thetae        = dvector(0, fsize-1);
    double** ze_NCSEM     = dmatrix(0, 5, 0, fsize-1);

    //====================================================================================
    //         Read
    //====================================================================================
    status = readCONT_txt(t0_CMU_EM, tf_CMU_EM, si_CMU_EM, si_CMS_SEM,
                          z0_CMU_NCEM, z0_CMS_NCSEM, thetae, ze_NCSEM, tr0, fsize);

    if(status)
    {
        cerr << fname << ". Data reading went wrong." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    //         Display
    //====================================================================================
    //for(int k = 0; k < fsize; k++) cout << ze_NCSEM[0][k] << endl;

    //====================================================================================
    // 2. Structures to compute the invariant manifolds
    //====================================================================================
    Invman invman_EM(OFTS_ORDER, OFS_ORDER, SEML.cs);
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);

    //====================================================================================
    // Select the parameters
    //====================================================================================
    for(int k = 0; k < 1 /*fsize*/; k+=50)
    {
        double t_EM[2] = {t0_CMU_EM[k], tf_CMU_EM[k]};
        double st_EM[5]  = {si_CMU_EM[0][k], si_CMU_EM[1][k], si_CMU_EM[2][k], si_CMU_EM[3][k], si_CMU_EM[4][k]};

        double t0_SEM = tf_CMU_EM[k]*SEML.us_em.ns;
        double st_SEM[5] = {si_CMS_SEM[0][k], si_CMS_SEM[1][k], si_CMS_SEM[2][k], si_CMS_SEM[3][k], si_CMS_SEM[4][k]};

        //================================================================================
        // Initialize local variables: EM
        //================================================================================
        //--------------------------------------------------------------------------------
        // Initialisation of the orbit structure
        //--------------------------------------------------------------------------------
        OdeStruct driver_EM;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        int coll_EM;
        OdeParams odeParams_EM(&coll_EM, &SEML_EM);
        //Init ode structure
        init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &odeParams_EM);
        //Init routine
        Orbit orbit_EM(&invman_EM, &SEML_EM, &driver_EM, OFTS_ORDER, OFS_ORDER, t_EM[0], t_EM[1]);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit_EM, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_EM.update_ic(st_EM,  t_EM[0]);

        //================================================================================
        // Initialize local variables: SEM
        //================================================================================

        //--------------------------------------------------------------------------------
        // Initialisation of the orbit structure
        //--------------------------------------------------------------------------------
        OdeStruct driver_SEM;
        //Parameters
        int coll_SEM;
        OdeParams odeParams_SEM(&coll_SEM, &SEML_SEM);
        //Init ode structure
        init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &odeParams_SEM);
        //Init routine
        Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &driver_SEM, OFTS_ORDER, OFS_ORDER,
                        t0_SEM, t0_SEM+5*SEML.us_sem.T);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit_SEM, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_SEM.update_ic(st_SEM, t0_SEM);

        //--------------------------------------------------------------------------------
        // We can advance to REF_COMP for the rest of the computation
        //--------------------------------------------------------------------------------
        refst.type = REF_COMP;
        int sampfreq[3] = {refst.sf_eml2, refst.sf_man, refst.sf_seml2};
        status = oocomprefft3d(sampfreq, coord_type, orbit_EM, orbit_SEM, refst);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure on the whole trajectory. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

    }

    return FTC_SUCCESS;
}


//========================================================================================
//         Refinement of solutions: CMU to CMS - general routines
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2.
 *         The final conditions vary in the paramerization of the CMS of SEMLi.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,
 *         the refst.coord_type must be NCSEM. However, the user can put
 *         other coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines
 *         are not used.
 **/
int oorefeml2seml(RefSt& refst)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    string fname = "oorefeml2seml";
    cout << "===================================================================" << endl;
    cout << "   oorefeml2seml. Refinement of EML2-SEMLi arc                     " << endl;
    cout << "===================================================================" << endl;

    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int coord_type = refst.coord_type;
    int man_grid_size = refst.gridSize;
    int dcs  = default_coordinate_system(coord_type);
    int status = 0;

    //====================================================================================
    // 2. Structures to compute the invariant manifolds
    //====================================================================================
    Invman invman_EM(OFTS_ORDER, OFS_ORDER, SEML.cs);
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);

    cout << "SEML_SEM.li_SEM = " << SEML_SEM.li_SEM << endl;
    cout << "SEML_SEM.li_SEM = " << SEML_SEM.li << endl;

    //====================================================================================
    // 3. User-defined number of continuation steps, if necessary
    //====================================================================================
    int cont_steps_MAX = 0;
    if(refst.type == REF_CONT || refst.type == REF_CONT_D ||
            refst.type == REF_CONT_D_HARD_CASE ||  refst.type == REF_COMP)
    {
        if(refst.isLimUD)
        {
            cout << "oorefeml2seml. User-defined number of continuation steps" << endl;
            cout << "Enter a value for cont_steps_MAX: ";
            cin >> cont_steps_MAX;
            refst.cont_step_max    = cont_steps_MAX;

            cout << "oorefeml2seml. User-defined number of continuation steps" << endl;
            cout << "Enter a value for cont_steps_MAX_vt: ";
            cin >> cont_steps_MAX;
            refst.cont_step_max_vt = cont_steps_MAX;
        }
    }
    else
    {
        cout << "oorefeml2seml. No continuation will be performed." << endl;
        refst.cont_step_max    = 0;
        refst.cont_step_max_vt = 0;
    }

    //====================================================================================
    // 4. Select the good IC for EML2-to-SEMLi connections in data files
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    ooseleml2seml(refst, st_EM, st_SEM, t_EM, &t0_SEM, &pmin_dist_SEM_out);

    cout << "===================================================================" << endl;
    cout << "Estimated error at patch point (km):                               " << endl;
    cout <<  pmin_dist_SEM_out* SEML.cs_sem.cr3bp.L                                << endl;
    cout << "Estimated error at patch point (SEMSU):                            " << endl;
    cout <<  pmin_dist_SEM_out                                                    << endl;
    cout << "===================================================================" << endl;
    pressEnter(refst.isFlagOn);

    //====================================================================================
    // 5.1 Initialize local variables: EM
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Initialisation of the orbit structure
    //------------------------------------------------------------------------------------
    OdeStruct driver_EM;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    int coll_EM;
    OdeParams odeParams_EM(&coll_EM, &SEML_EM);
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &odeParams_EM);
    //Init routine
    Orbit orbit_EM(&invman_EM, &SEML_EM, &driver_EM, OFTS_ORDER, OFS_ORDER, t_EM[0], t_EM[1]);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ic(st_EM,  t_EM[0]);

    //====================================================================================
    // 5.2 Initialize local variables: SEM
    //====================================================================================

    //------------------------------------------------------------------------------------
    // Initialisation of the orbit structure
    //------------------------------------------------------------------------------------
    OdeStruct driver_SEM;
    //Parameters
    int coll_SEM;
    OdeParams odeParams_SEM(&coll_SEM, &SEML_SEM);
    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &odeParams_SEM);
    //Init routine
    Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &driver_SEM, OFTS_ORDER, OFS_ORDER,
                    t0_SEM, t0_SEM+5*SEML.us_sem.T);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_SEM.update_ic(st_SEM, t0_SEM);

    //------------------------------------------------------------------------------------
    // Time and state vectors (twice the size)
    //------------------------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, 2*man_grid_size);
    double*  t_traj  = dvector(0, 2*man_grid_size);

    //====================================================================================
    // 6 Multiple shooting procedure
    //====================================================================================
    gnuplot_ctrl* h2;
    switch(refst.type)
    {
    case REF_CONT:
    case REF_SINGLE:
    {
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type);

        //--------------------------------------------------------------------------------
        //Multiple shooting procedure
        //--------------------------------------------------------------------------------
        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        gnuplot_close(h2);
        if(status)
        {
            cerr << fname << ". Error during the refinement procedure. ref_errno = " << ref_strerror (status) << endl;
            return FTC_FAILURE;
        }
        break;
    }

    case REF_CONT_D:
    {
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type);

        //--------------------------------------------------------------------------------
        //First step: variable tn
        //--------------------------------------------------------------------------------
        refst.time = REF_VAR_TN;

        //There is no need to save at this point,
        //since everything is gonna be erased by the second step
        int isSaved = refst.isSaved;
        refst.isSaved = 0;

        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }


        //--------------------------------------------------------------------------------
        //Second step: fixed time + saved, if desired by the user
        //--------------------------------------------------------------------------------
        refst.time    = REF_FIXED_TIME;
        refst.isSaved = isSaved;

        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with fixed time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        gnuplot_close(h2);
        break;
    }

    case REF_CONT_D_HARD_CASE:
    {
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type);

        //--------------------------------------------------------------------------------
        //First step: variable tn
        //--------------------------------------------------------------------------------
        refst.time = REF_VAR_TN;

        //There is no need to save at this point,
        //since everything is gonna be erased by the second step
        int isSaved = refst.isSaved;
        refst.isSaved = 0;

        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }


        //--------------------------------------------------------------------------------
        //Second step: fixed time + given grid + saved, if desired by the user
        // The main change with respect to REF_CONT_D is REF_GIVEN_GRID: this means
        // that the time grid used for the first phase is kept for the next phase (fixed time)
        //--------------------------------------------------------------------------------
        refst.time    = REF_FIXED_TIME;
        refst.grid    = REF_GIVEN_GRID;
        refst.isSaved = isSaved;

        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with fixed time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        gnuplot_close(h2);
        break;
    }


    case REF_COMP:
    {
        cout << "===============================================================" << endl;
        cout << "oorefeml2seml.                                                 " << endl;
        cout << "The computation of the entire trajectory has been selected.    " << endl;
        cout << "===============================================================" << endl;
        cout << "oorefeml2seml. First part: "                                     << endl;
        cout << "Refinement of the connection leg in the parameterization space." << endl;
        //--------------------------------------------------------------------------------
        //First part: REF_CONT with variable tn
        //--------------------------------------------------------------------------------
        refst.type = REF_CONT;
        refst.time = REF_VAR_TN;

        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type);

        //Multiple shooting procedure
        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        //Second part: REF_SINGLE
        //--------------------------------------------------------------------------------
        refst.type = REF_SINGLE;
        refst.time = REF_FIXED_TIME;


        //Multiple shooting procedure
        status = oosrefeml2seml(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refst, h2);
        gnuplot_close(h2);

        if(status)
        {
            cerr << fname << ". Error during the refinement procedure. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        // After this point, the orbits are updated so that the connection is
        // truly natural in the parameterization space.
        // We can then advance to REF_COMP for the rest of the computation
        //--------------------------------------------------------------------------------
        refst.type = REF_COMP;
        int sf_eml2, sf_man, sf_seml2;
        cout << "===============================================================" << endl;
        cout << "oorefeml2seml. Second part: "                                    << endl;
        cout << "Refinement of the entire trajectory."                            << endl;
        cout << "This is a 3-legged trajectory: "                                 << endl;
        cout << "EML2 orbit + manifold leg + SEMLi orbit."                        << endl;
        cout << "Enter a value for sf_eml2, "                                     << endl;
        cout << "the sampling frequeny (in days) during the EML2 leg: "           << endl;
        cin >> sf_eml2;
        cout << "Enter a value for sf_man, "                                      << endl;
        cout << "the sampling frequeny (in days) during the man leg: "            << endl;
        cin >> sf_man;
        cout << "Enter a value for sf_seml2, "                                    << endl;
        cout << "the sampling frequeny (in days) during the SEMLi leg: "          << endl;
        cin >> sf_seml2;
        char ch;
        scanf("%c",&ch);

        int sampfreq[3] = {sf_eml2, sf_man, sf_seml2};
        status = oocomprefft3d(sampfreq, coord_type, orbit_EM, orbit_SEM, refst);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure on the whole trajectory. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        //CAREFUL!! the rest of the routines are still using msd_grid_size
        //(a number of points) and NOT the sampling frequencies!! To be adapted...
        //--------------------------------------------------------------------------------
        //oocomprefft3d_test_seml_synjpl(msd_grid_size, coord_type, orbit_SEM, refst);
        //oocomprefft3d_test_eml_synjpl(msd_grid_size, coord_type, orbit_EM, refst);
        //comprefft3d_test_eml2seml_synjpl(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        //comprefft3d_test_eml2seml_insem(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        break;
    }
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    pressEnter(refst.isFlagOn);


    return FTC_SUCCESS;
}

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - subroutines
//
//========================================================================================
/**
 *  \brief Find the intersection of a EML2-SEMLi connection contained in
 *         y_traj_n/t_traj_n with a certain Pk section x = cst defined by refst.
 **/
int ooxpseml2seml(double ye[6], double* te, double* t_traj_n, double** y_traj_n, int man_index, RefSt& refst)
{
    double yv[42];
    //Get the default coordinates system from the coord_type
    int dcs  = default_coordinate_system(NCSEM);

    //====================================================================================
    //1. Find, from the end, the couple of patch points before and after the PS
    //====================================================================================
    int k2 = man_index;
    int k1 = man_index-1;
    bool test = false;
    do
    {
        k1--;
        k2--;
        test = (y_traj_n[0][k2] > refst.xps && y_traj_n[0][k1] < refst.xps) ||
               (y_traj_n[0][k2] < refst.xps && y_traj_n[0][k1] > refst.xps);
    }
    while(!test && k1 >= 0);


    //====================================================================================
    //2. If a solution has been found:
    //====================================================================================
    if(test && k1 >= 0)
    {
        cout << "ooxpseml2seml. Range accross the PS found:" << endl;
        cout << y_traj_n[0][k1] << " < " << refst.xps << " < " << y_traj_n[0][k2] << endl;


        //====================================================================================
        //2. Integrate until x = xps
        //====================================================================================
        double center[3];
        struct value_params val_par;
        val_par.max_events = 1;
        val_par.direction  = 0;
        val_par.dim        = 0;
        val_par.value      = refst.xps;
        val_par.center     = center;
        val_par.type       = 'X';

        double** ye_NCSEM = dmatrix(0, 5, 0, 1);
        double* te_NCSEM  = dvector(0, 1);

        //====================================================================================
        //3. After this step: ye_NCSEM[*][0] & te_NCSEM[0] contains the intersection
        //   Note: the (possible) collisions with the primaries are not taken into account here:
        //   We suppose that for this particular application (intersection with a Pk section
        //   far away from any primary), such collision are very unlikely!
        //====================================================================================
        int ode78coll;
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k1];
        ode78_qbcp_event(ye_NCSEM, te_NCSEM, &ode78coll, t_traj_n[k1], t_traj_n[k2], yv, 6, dcs,
                         NCSEM, NCSEM, &val_par);

        //====================================================================================
        //4. Store
        //====================================================================================
        *te = te_NCSEM[0];
        for(int i = 0; i < 6; i++) ye[i] = ye_NCSEM[i][0];

        //====================================================================================
        //5. Free
        //====================================================================================
        free_dmatrix(ye_NCSEM, 0, 5, 0, 1);
        free_dvector(te_NCSEM, 0, 1);

        return FTC_SUCCESS;
    }
    else
    {
        *te = -1;
        for(int i = 0; i < 6; i++) ye[i] = 0.0;
        return FTC_FAILURE;
    }
}

/**
 *  \brief Continuation of a single of EML2-to-SEMLi connection, between orbit_EM and
 *         orbit_SEM.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,the
 *         coord_type must be NCSEM. However, the user can put other
 *         coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines are
 *         not used.
 *
 *         In general, we recommend to use NCSEM at least for precision
 *         reasons. In this way, a warning is issued when coord_type is different
 *         from NCSEM.
 **/
int oosrefeml2seml(Orbit& orbit_EM, Orbit& orbit_SEM, double** y_traj, double* t_traj,
                   int dcs, int coord_type, int* man_grid_size_t,
                   RefSt& refst, gnuplot_ctrl* h2)
{
    //====================================================================================
    // 0. Check on  coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cout << "oosrefeml2seml. WARNING: the coord_type inside refst is different  " << endl;
        cout << "from NCSEM, which means that the variable-time capability cannot   " << endl;
        cout << "be used. If you still wish to continue, press enter.               " << endl;
        cout << "===================================================================" << endl;
        char ch;
        scanf("%c",&ch);
    }

    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "oosrefeml2seml";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = 2;//floor((double)refst.cont_step_max/floor(0.1*GP_MAX_TMP_FILES))+1;

    //------------------------------------------------------------------------------------
    //If the grid is fixed, then we use the user-provided size. Else, we use a
    //big value (max_grid_size) so that the arrays would not be saturated.
    //------------------------------------------------------------------------------------
    int max_grid_size = 1000;
    int man_grid_size = (refst.grid == REF_VAR_GRID) ? max_grid_size:*man_grid_size_t;

    //------------------------------------------------------------------------------------
    //Local variables for plotting
    //------------------------------------------------------------------------------------
    int mPlot = 500;
    double** y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double** y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double*  t_man_EM          = dvector(0, mPlot);
    double*  t_man_SEM         = dvector(0, mPlot);
    double** y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double*  t_man_coord_plot  = dvector(0, mPlot);

    //Number of potential new points, during the continuation procedure with variable time
    int nnew = 4;
    int ode78coll = 0;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    double** ymt   = dmatrix(0, 5, 0, nnew);
    double* tmt    = dvector(0, nnew);
    //double** y_traj       = dmatrix(0, 41, 0, man_grid_size+2*nnew);
    //double*  t_traj       = dvector(0, man_grid_size+2*nnew);

    //------------------------------------------------------------------------------------
    //Type of corrector & predictor
    //------------------------------------------------------------------------------------
    diffcorrptr  diffcorr  = ftc_select_diffcorr(refst);
    predictorptr predictor = ftc_select_predictor(refst);

    //------------------------------------------------------------------------------------
    //To store final data
    //------------------------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, man_grid_size+2*nnew);
    double* t_traj_n    = dvector(0, man_grid_size+2*nnew);

    //------------------------------------------------------------------------------------
    // Initialize the COC matrices
    //------------------------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex* CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEMLi
    gsl_matrix_complex* CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //====================================================================================
    // 2.  Compute first guess
    //====================================================================================
    int man_index = oomanfgeml2seml(y_traj, t_traj, orbit_EM, orbit_SEM, dcs, coord_type, man_grid_size, refst);

    // Plot the resulting trajectory
    gnuplot_plot_X(h2, y_traj, man_index,   (char*)"", "lines", "1", "1", 4);
    gnuplot_plot_X(h2, y_traj, man_index+1, (char*)"", "points", "1", "1", 4);


    //====================================================================================
    // 3.  Differential correction
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    pressEnter(refst.isFlagOn);

    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    int nfv = nfreevariables(refst, man_index);
    double* nullvector = dvector(0, nfv-1);

    //------------------------------------------------------------------------------------
    //GNUPLOT:
    //  1. if there is a continuation procedure based on variable time,
    //     the component s5 at SEML is plotted.
    //  2. if there is a continuation procedure based on fixed time,
    //     several computations are produced.
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h3 = 0, *h4 = 0, *h5 = 0;
    if(refst.isCont())
    {
        switch(refst.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:

            h5 = gnuplot_init();
            gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
            gnuplot_cmd(h5, "set grid");
            break;


        case REF_FIXED_TIME:

            h3 = gnuplot_init();
            gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");
            gnuplot_cmd(h3, "set grid");

            h4 = gnuplot_init();
            gnuplot_cmd(h4,  "set title \"t0 vs s1_EM\" ");
            gnuplot_cmd(h4, "set grid");

            h5 = gnuplot_init();
            gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
            gnuplot_cmd(h5, "set grid");
            break;


        }
    }


    //====================================================================================
    //6.1. First step of the continuation procedure.
    //====================================================================================
    int status = 0;
    double inner_prec = 5e-8;
    int niter = 1;

    status = diffcorr(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                      man_index, coord_type, inner_prec, true, orbit_EM, orbit_SEM,
                      h2, refst, &niter);

    //====================================================================================
    //6.2. If it is a sucess, we go on.
    //====================================================================================
    //------------------------------------------------------------------------------------
    // If it is not a success, we print the value, and we return
    //------------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". Error during the first refinement procedure. ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }
    else
    {
        double yv[42];
        //================================================================================
        // Find the intersection with a certain Poincaré Section (PS)
        //================================================================================
        double ye[6], te = 0.0;
        ooxpseml2seml(ye, &te, t_traj_n, y_traj_n, man_index, refst);

        //================================================================================
        // Save first entry if the continuation process is on
        //================================================================================
        double** ymc   = dmatrix(0, 5, 0, mPlot);
        double* tmc    = dvector(0, mPlot);
        string filename_traj = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF_TRAJ, SEML.li_SEM, orbit_EM.getT0()/SEML.us_em.T);

        cout << "filename_traj = " << filename_traj << endl;
        pressEnter(true);

        fstream filestream;
        if(refst.isSaved && refst.isCont())
        {
            cout << "-----------------------------------------------------------"  << endl;
            cout << fname << ". Save first entry in " << filename_traj             << endl;
            cout << "-----------------------------------------------------------"  << endl;

            //============================================================================
            // Main parameters
            //============================================================================
            writeCONT_txt(orbit_EM, orbit_SEM, &te, ye, true);

            //============================================================================
            // Entire trajectory
            //============================================================================
            filestream.open (filename_traj.c_str(), ios::out | ios::binary);
            double res;
            //Final trajectory on lines, segment by segment
            int status;
            for(int k = 0; k < man_index; k++)
            {
                //Integration segment by segment
                for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                status = ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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

                //Storing
                for(int p = 0; p < mPlot; p++)
                {
                    res = 0.0;
                    filestream.write((char*) &res, sizeof(double));

                    res = tmc[p];
                    filestream.write((char*) &res, sizeof(double));

                    for(int i = 0; i < 6; i++)
                    {
                        res = ymc[i][p];
                        filestream.write((char*) &res, sizeof(double));
                    }
                }
            }
            filestream.close();


            //============================================================================
            // Compute final orbit in reduced coordinates
            //============================================================================
            /*
            //            cout << "-------------------------------------------"  << endl;
            //            cout << "Compute the initial SEMLi orbit            "  << endl;
            //            cout << "-------------------------------------------"  << endl;
            //            printf("Press ENTER to go on\n");
            //            scanf("%c",&ch);
            //
            //
            //            vector<Oftsc> Fh;
            //            Fh.reserve(5);
            //            for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
            //            readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");
            //
            //            //------------------------------------------------------------------------------------
            //            //For dot(s) = fh(s)
            //            //------------------------------------------------------------------------------------
            //            RVF rvf;
            //            rvf.ofs_order = SEML.eff_nf_SEM;
            //            Ofsc AUX(rvf.ofs_order);
            //            rvf.fh         = &Fh;
            //            rvf.ofs        = &AUX;
            //            rvf.order      = OFTS_ORDER;
            //            rvf.n          = orbit_SEM.getN();
            //            rvf.reduced_nv = 5;
            //
            //            gsl_odeiv2_system sys_fh;
            //            sys_fh.function  = qbfbp_fh;
            //            sys_fh.jacobian  = NULL;
            //            sys_fh.dimension = 2*rvf.reduced_nv;
            //            sys_fh.params    = &rvf;
            //
            //            const gsl_odeiv2_step_type* T_fh = gsl_odeiv2_step_rk8pd;
            //            gsl_odeiv2_driver* d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);
            //
            //
            //            //------------------------------------------------------------------------------------
            //            // Temp variables
            //            //------------------------------------------------------------------------------------
            //            double t0_SEM = t_traj_n[man_index];
            //            double t1_SEM = t0_SEM+tof_seml_SEM;
            //            double z[6];
            //            double t2 = t0_SEM;
            //            int k  = 0;
            //            double  s1ccm8[2*rvf.reduced_nv]; //CCM8
            //
            //            //------------------------------------------------------------------------------------
            //            // Initial state in CCM8 form
            //            //------------------------------------------------------------------------------------
            //            RCMtoCCM8(orbit_SEM.getSi(), s1ccm8, 5);
            //
            //            //------------------------------------------------------------------------------------
            //            // Reopen the file
            //            //------------------------------------------------------------------------------------
            //            filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
            //
            //            //------------------------------------------------------------------------------------
            //            // Loop
            //            //------------------------------------------------------------------------------------
            //            while(k <= mPlot && orbit_SEM.getSi()[4] > 1e-5)
            //            {
            //                cout << "orbit_SEM.si[4] = " << orbit_SEM.getSi()[4] << endl;
            //                //To NC coordinates
            //                orbit_SEM.evalRCMtoNC(t2, z);
            //
            //                //Save
            //                for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
            //                t_man_SEM[k] = t2;
            //
            //                //Plot
            //                gnuplot_plot_X(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);
            //
            //                //Save in file
            //                res = -1.0;
            //                filestream.write((char*) &res, sizeof(double));
            //
            //                res = t2;
            //                filestream.write((char*) &res, sizeof(double));
            //
            //                for(int i = 0; i < 6; i++)
            //                {
            //                    res = z[i];
            //                    filestream.write((char*) &res, sizeof(double));
            //                }
            //
            //                //Advance one step
            //                gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
            //                orbit_SEM.ccm8torcm(s1ccm8);
            //                k++;
            //            }
            //            filestream.close()
            */

        }

        //================================================================================
        //Continuation procedure
        //================================================================================
        if(refst.isCont())
        {
            cout << fname << ". Continuation procedure.                         "  << endl;
            cout << "-----------------------------------------------------------"  << endl;
            double ds0 = 1e-1;
            double niter0 = 4;
            switch(refst.time)
            {
            case REF_VAR_TIME:
            case REF_VAR_TN:
                ds0 = 8e-2;
                break;
            case REF_FIXED_TIME:
                ds0 = 5e-1;
                break;
            }
            double ds  = ds0;
            double dkn = 0.0;
            int kn = 0;

            //----------------------------------------------------------------------------
            //An additional condition can be set:
            //
            // REF_COND_S5: if the time is not fixed, we want to decrease the
            //s5 component at SEMLi. Therefore, we add a stopping condition on the loop.
            //
            // REF_COND_T: we want to make enough "turns" about SEMLi
            //----------------------------------------------------------------------------
            bool addCondition = true;
            double theta_acc  = 0.0, theta_old = 0.0, theta_max  = 360;
            int mani_old = man_index;
            switch(refst.termination)
            {
            case REF_COND_S5:
            {
                switch(refst.time)
                {
                case REF_VAR_TIME:
                case REF_VAR_TN:
                    addCondition = fabs(orbit_SEM.getSi()[4]) > ORBIT_SEM_UNSTABLE_MIN;
                    break;


                case REF_FIXED_TIME:
                    addCondition =  true;
                    break;
                }
                break;
            }

            case REF_COND_T:
            {
                addCondition =  true;
                break;
            }
            }

            //============================================================================
            // Continuation loop
            //============================================================================

            //----------------------------------------------------------------------------
            // The maximum number of steps depends on the strategy
            //----------------------------------------------------------------------------
            int cont_step_max_local = refst.cont_step_max;
            switch(refst.time)
            {
            case REF_VAR_TIME:
            case REF_VAR_TN:
                cont_step_max_local = refst.cont_step_max_vt;
                break;


            case REF_FIXED_TIME:
                cont_step_max_local = refst.cont_step_max;
                break;
            }

            //----------------------------------------------------------------------------
            //Loop
            //----------------------------------------------------------------------------
            while(kn < cont_step_max_local && status == GSL_SUCCESS && addCondition)
            {
                //========================================================================
                // 5.2.1. Updating the free variables via the predictor
                //========================================================================
                predictor(y_traj_n, t_traj_n, &ds, ds0, nullvector,
                          orbit_EM, orbit_SEM, man_index, coord_type, refst);

                //========================================================================
                // 5.2.2. Diff correction
                //========================================================================
                status = diffcorr(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                  nullvector, 42, man_index, coord_type, inner_prec,
                                  false, orbit_EM, orbit_SEM, h2, refst, &niter);

                //========================================================================
                // Find the intersection with a certain Poincaré Section (PS)
                //========================================================================
                double ye[6], te = 0.0;
                if(status == GSL_SUCCESS) ooxpseml2seml(ye, &te, t_traj_n, y_traj_n, man_index, refst);

                //========================================================================
                // 5.2.3. Save
                //========================================================================
                if(refst.isSaved && status == GSL_SUCCESS)
                {

                    //====================================================================
                    // Main parameters
                    //====================================================================
                    writeCONT_txt(orbit_EM, orbit_SEM, &te, ye, false);

                    //====================================================================
                    // Entire trajectory
                    //====================================================================
                    filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
                    double res;
                    //Final trajectory on lines, segment by segment
                    for(int k = 0; k < man_index; k++)
                    {
                        //Integration segment by segment
                        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                        status = ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
                            cout << fname << ". Warning: a collision with a primary has occurred." << endl;
                        }

                        //Storing
                        for(int p = 0; p < mPlot; p++)
                        {
                            res = 1.0*(kn+1);
                            filestream.write((char*) &res, sizeof(double));

                            res = tmc[p];
                            filestream.write((char*) &res, sizeof(double));

                            for(int i = 0; i < 6; i++)
                            {
                                res = ymc[i][p];
                                filestream.write((char*) &res, sizeof(double));
                            }
                        }
                    }
                    filestream.close();
                }


                //========================================================================
                // 5.2.4. Display
                //========================================================================
                dkn = (double) kn;
                double t0 = orbit_EM.getT0()/SEML.us_sem.T;//t_traj_n[0]/SEML.us_sem.T;
                if(refst.isCont() && status == GSL_SUCCESS && kn % plotfreq == 0)
                {
                    switch(refst.time)
                    {
                    case REF_VAR_TIME:
                    case REF_VAR_TN:
                    {
                        gnuplot_plotc_xy(h5, &dkn, &orbit_SEM.getSi()[4], 1, (char*)"", "points", "1", "2", 0);
                        break;
                    }

                    case REF_FIXED_TIME:
                    {
                        gnuplot_plotc_xy(h3, &orbit_EM.getSi()[0],  &orbit_EM.getSi()[2], 1, (char*)"", "points", "1", "2", 0);
                        gnuplot_plotc_xy(h4, &t0, &orbit_EM.getSi()[0], 1, (char*)"", "points", "1", "2", 0);
                        gnuplot_plotc_xy(h5, &dkn, &orbit_SEM.getSi()[4], 1, (char*)"", "points", "1", "2", 0);
                        break;
                    }

                    }
                }

                //========================================================================
                // 5.2.5. update the additionnal condition, if necessary
                //========================================================================
                switch(refst.termination)
                {
                case REF_COND_S5:
                {
                    //--------------------------------------------------------------------
                    // First type of condition: stop when we are close enough to the
                    // center manifold (unstable component is small enough)
                    // note: this condition is not used when all the times are fixed.
                    //--------------------------------------------------------------------
                    switch(refst.time)
                    {
                    case REF_VAR_TIME:
                    case REF_VAR_TN:
                        addCondition = fabs(orbit_SEM.getSi()[4]) > ORBIT_SEM_UNSTABLE_MIN;
                        break;
                    case REF_FIXED_TIME:
                        addCondition = true;
                        break;
                    }
                    break;
                }
                case REF_COND_T:
                {
                    //--------------------------------------------------------------------
                    // Another possible condition: enough turns around SEMLi.
                    //--------------------------------------------------------------------
                    double det, dot, theta;

                    //Compute the angle between the last point (entry in the manifold)
                    // and the second last, at the beginning of the procedure
                    dot = y_traj_n[0][man_index]*y_traj_n[0][mani_old-1] +
                          y_traj_n[1][man_index]*y_traj_n[1][mani_old-1];  // dot product
                    det = y_traj_n[0][man_index]*y_traj_n[1][mani_old-1] -
                          y_traj_n[1][man_index]*y_traj_n[0][mani_old-1];  // determinant

                    theta = 180/M_PI*atan(det/dot);//180/M_PI*atan2(det, dot);  // atan2(sin, cos) or tan(sin/cos)

                    //We save the value
                    if(kn == 0)
                    {
                        theta_old = theta;
                    }
                    else
                    {
                        //----------------------------------------------------------------
                        // Each time that theta_old*theta < 0, it means that we have made
                        // about 180° around the last point. We then do 2 actions:
                        //  1. We update the accumulated angle (+90°)
                        //  2. We add a few patch points, to ensure numerical stability.
                        //----------------------------------------------------------------
                        if(theta_old*theta < 0)
                        {
                            //  1. We update the accumulated angle (+90°)
                            theta_acc += 90;
                        }

                        theta_old = theta;
                    }



                    //--------------------------------------------------------------------
                    // The condition is only applied when some times are left free to vary
                    //--------------------------------------------------------------------
                    switch(refst.time)
                    {
                    case REF_VAR_TIME:
                    case REF_VAR_TN:
                        cout << fname << ". theta     = " << theta << endl;
                        cout << fname << ". theta_acc = " << theta_acc << endl;
                        cout << fname << ". theta_max = " << theta_max << endl;
                        addCondition = theta_acc < theta_max;
                        break;
                    case REF_FIXED_TIME:
                        addCondition = true;
                        break;
                    }

                    break;
                }
                }


                //========================================================================
                // 5.2.4. Advance
                //========================================================================
                kn++;
                cout << "Step n°" << kn << "/" << cont_step_max_local << " complete"  << endl;
                cout << "----------------------------------------------------------"  << endl;


                //========================================================================
                // 5.2.6. update the stepper. For now, only in the REF_VAR_TIME case !
                //========================================================================
                switch(refst.time)
                {
                case REF_VAR_TIME:
                case REF_VAR_TN:
                    ds *= niter0/niter;
                    break;
                case REF_FIXED_TIME:
                    ds *= niter0/niter;
                    break;
                }
            }

            //============================================================================
            // If the trajectory has been extended... We had a few points.
            // Careful with this piece of code: if nnew is big, then a big number of points
            // is added to the state vectors, which can make the Newton procedure unstable.
            //
            // Consequently, this option is taken only if refst.type = REF_CONT_D_HARD_CASE
            // i.e. basically if the user has tried refst.type = REF_CONT_D before
            // but it did not work...
            //
            // IMPORTANT: after this step, the number of points man_grid_size_t is changed!
            // this means that any further use of this routine with man_grid_size_t will
            // take into account this change.
            //
            //============================================================================
            if(theta_acc >= theta_max && refst.type == REF_CONT_D_HARD_CASE)
            {
                //  Add a few patch points
                for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][man_index-1];
                status = ode78(ymt, tmt, &ode78coll, t_traj_n[man_index-1],
                               t_traj_n[man_index], yv, 6, nnew, dcs, coord_type,
                               coord_type);

                //Copy y_traj_n/t_traj_n in y_traj/t_traj
                for(int k = 0; k <= man_index; k++)
                {
                    for(int i = 0; i < 6; i++) y_traj[i][k] = y_traj_n[i][k];
                    t_traj[k] = t_traj_n[k];
                }

                //Copy y_traj/t_traj back in y_traj_n/t_traj_n, except the last point
                for(int k = 0; k < man_index; k++)
                {
                    for(int i = 0; i < 6; i++) y_traj_n[i][k] = y_traj[i][k];
                    t_traj_n[k] = t_traj[k];
                }

                //Copy additionnal points, except the first one
                for(int k = 1; k <= nnew; k++)
                {
                    for(int i = 0; i < 6; i++) y_traj_n[i][k+man_index-1] = ymt[i][k];
                    t_traj_n[k+man_index-1] = tmt[k];
                }

                //////////////////////////////////////////////////////////////////////////
                //New size is man_index+nnew-1 +
                //IMPORTANT: we also update *man_grid_size_t for
                //next continuation !!!
                //////////////////////////////////////////////////////////////////////////
                man_index += nnew-1;
                *man_grid_size_t = man_index;

                //Plot the new distribution
                gnuplot_plot_X(h2, y_traj_n, man_index+1, (char*)"", "points", "1", "1", 8);
            }
        }

    }

    //====================================================================================
    // Final correction with better precision
    // CAREFUL: it is NOT done for now because some elements depends on man_index for their size,
    // and man_index CAN move if REF_CONT_D_HARD_CASE is used !!
    // TODO: make a check on who is using man_index...
    //====================================================================================
    /*
    //    inner_prec = 1e-10;                    //harder precision in diff corr procedures
    //    Config::configManager().C_PREC_HARD(); //harder precision in numerical integration
    //    switch(refst.dim)
    //    {
    //    case REF_3D:
    //        switch(refst.time)
    //        {
    //        case REF_FIXED_TIME:
    //            status = msft3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                            nullvector, 42,
    //                            man_index, coord_type, inner_prec, false,
    //                            orbit_EM, orbit_SEM,
    //                            h2, refst, &niter);
    //            break;
    //
    //        case REF_VAR_TIME:
    //        case REF_VAR_TN:
    //            status = msvt3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                            nullvector, 42,
    //                            man_index, coord_type, inner_prec, false,
    //                            orbit_EM, orbit_SEM,
    //                            h2, refst, &niter);
    //            break;
    //        }
    //        break;
    //    case REF_PLANAR:
    //        switch(refst.time)
    //        {
    //        case REF_FIXED_TIME:
    //            status = msftplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                              nullvector, 42,
    //                              man_index, coord_type, inner_prec, false,
    //                              orbit_EM, orbit_SEM,
    //                              h2, refst, &niter);
    //            break;
    //
    //        case REF_VAR_TIME:
    //            status = msvtplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                              nullvector, 42,
    //                              man_index, coord_type, inner_prec, false,
    //                              orbit_EM, orbit_SEM,
    //                              h2, refst, &niter);
    //            break;
    //
    //        case REF_VAR_TN:
    //            status = msvltplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                               nullvector, 42,
    //                               man_index, coord_type, inner_prec, false,
    //                               orbit_EM, orbit_SEM,
    //                               h2, refst, &niter);
    //            break;
    //        }
    //        break;
    //    }
    //
    //    Config::configManager().C_PREC_BACK(); //back to initial precision
    */


    //====================================================================================
    // 7. Display final solution
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Error during the continuation procedure. ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }
    else
    {
        //================================================================================
        // 7.1. Compute the initial orbit
        //================================================================================
        //TOF on EML2 orbit
        double tof_eml_EM   = 5*SEML.us.T;
        //TOF on SEMLi orbit
        //double tof_seml_SEM = 30*SEML.us_sem.T;
        cout << "-----------------------------------------------------------"  << endl;
        cout << fname << ". Computation of the initial EML2 orbit           "  << endl;
        cout << "-----------------------------------------------------------"  << endl;

        //--------------------------------------------------------------------------------
        // Reset the unstable direction
        //--------------------------------------------------------------------------------
        orbit_EM.setSi(0, 4);
        orbit_EM.setTf(orbit_EM.getT0()-tof_eml_EM);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //--------------------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //--------------------------------------------------------------------------------
        int output = orbit_EM.traj_int_grid(orbit_EM.getTf(), y_man_NCEM, t_man_EM, mPlot, 1);

        //If output is strictly greater than 0, then the projection procedure inside
        //the integration went wrong, and the new indix is output.
        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //is no point stored in the data, except the first one, and the new indix is zero
        if(output == ORBIT_EPROJ) mPlot = 0;
        else if(output > 0) mPlot = output;

        //--------------------------------------------------------------------------------
        //To SEM coordinates
        //--------------------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //--------------------------------------------------------------------------------
        //Plot
        //--------------------------------------------------------------------------------
        gnuplot_plot_X(h2, y_man_coord_plot, mPlot+1, (char*)"", "lines", "1", "1", 2);

        //================================================================================
        // 7.2. Final trajectory, on a grid
        //================================================================================
        cout << fname << ". Final trajectory, on a grid.                    "  << endl;
        cout << "-----------------------------------------------------------"  << endl;
        double** ymc   = dmatrix(0, 5, 0, mPlot);
        double* tmc    = dvector(0, mPlot);
        double yv[42];
        int ode78coll;
        //Final trajectory on lines, segment by segment
        for(int k = 0; k < man_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            status = ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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

            //Plot on h2
            if(k == 0)
                gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Final.", "lines", "1", "2", 0);
            else
                gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 0);
        }

        //================================================================================
        // Free.
        //================================================================================
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
    }


    //====================================================================================
    // 8. Update the orbits for next step
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Error during the plotting procedure. ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }
    else
    {
        //================================================================================
        // Finally, copy y_traj/t_traj if necessary for the next step.
        //================================================================================
        for(int k = 0; k <= man_index; k++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][k] = y_traj_n[i][k];
            t_traj[k] = t_traj_n[k];
        }

        //================================================================================
        // At EML2. The first 4 RCM components are good, as well as the initial time.
        // Hence, we need to update:
        // 1. The last RCM component (unstable part),
        // 2. The final time.
        //================================================================================
        orbit_EM.setSi(PROJ_EPSILON, 4);
        orbit_EM.setTf(t_traj_n[man_index]/SEML.us_em.ns);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //================================================================================
        // At SEMLi. The first 5 RCM components are good.
        // Hence, we need to update:
        // 1. The initial time.
        //================================================================================
        orbit_SEM.setT0(t_traj_n[man_index]);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());
    }

    cout << fname << ". "                                    << endl;
    cout << "          End of computation.                "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    pressEnter(refst.isFlagOn);

    //====================================================================================
    // 9. Free
    //====================================================================================
    if(refst.isCont())
    {
        switch(refst.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:

            gnuplot_close(h5);
            break;

        case REF_FIXED_TIME:

            gnuplot_close(h3);
            gnuplot_close(h4);
            gnuplot_close(h5);
            break;
        }
    }

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    //free_dmatrix(y_traj, 0, 41, 0, man_grid_size+2*nnew);
    //free_dvector(t_traj, 0, man_grid_size+2*nnew);
    free_dmatrix(y_traj_n, 0, 41, 0, man_grid_size+2*nnew);
    free_dvector(t_traj_n, 0, man_grid_size+2*nnew);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);


    free_dmatrix(ymt, 0, 5, 0, nnew);
    free_dvector(tmt, 0, nnew);

    return FTC_SUCCESS;
}

//----------------------------------------------------------------------------------------
//         Brick A: select the good IC for EML2-to-SEMli connections in data files
//----------------------------------------------------------------------------------------
/**
 *  \brief Selects good initial conditions for EML2-to-SEMLi connections, searching
 *         through data files produced by oo_int_proj_CMU_EM_on_CM_SEM_3D or
 *         oo_int_proj_CMU_EM_on_CM_SEM.
 **/
int ooseleml2seml(RefSt& refst, double st_EM[5], double st_SEM[5], double t_EM[2],
                  double* t0_SEM, double* pmin_dist_SEM_out)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    cout << "-------------------------------------------------------------------" << endl;
    cout << "   ooseleml2seml. Selects good IC for EML2-SEMLi arc"                << endl;
    cout << "-------------------------------------------------------------------" << endl;
    //====================================================================================
    // 1. Init the data containers
    //====================================================================================
    //------------------------------------------------------------------------------------
    //To store data from the data file
    //------------------------------------------------------------------------------------
    vector<double> t0_EM;
    vector<double> tf_EM;
    vector<double> s1_CMU_EM;
    vector<double> s2_CMU_EM;
    vector<double> s3_CMU_EM;
    vector<double> s4_CMU_EM;
    vector<double> s5_CMU_EM;
    vector<double> pmin_dist_SEM;
    vector<double> s1_CM_SEM;
    vector<double> s2_CM_SEM;
    vector<double> s3_CM_SEM;
    vector<double> s4_CM_SEM;
    vector<size_t> sortId;

    //------------------------------------------------------------------------------------
    //To store data from a sub selection of the previous elements
    //------------------------------------------------------------------------------------
    vector<double> t0_EM_R;
    vector<double> tf_EM_R;
    vector<double> s1_CMU_EM_R;
    vector<double> s2_CMU_EM_R;
    vector<double> s3_CMU_EM_R;
    vector<double> s4_CMU_EM_R;
    vector<double> s5_CMU_EM_R;
    vector<double> pmin_dist_SEM_R;
    vector<double> s1_CM_SEM_R;
    vector<double> s2_CM_SEM_R;
    vector<double> s3_CM_SEM_R;
    vector<double> s4_CM_SEM_R;
    vector<size_t> sortId_R;


    //====================================================================================
    // 2. Getting back the data
    //====================================================================================
    string filename;
    if(refst.isFromServer)
    {

        double t0_des_mod = fmod(refst.t0_des/SEML.us_em.T, 1.0);

        cout << "refst.t0_des/SEML.us_em.T = " << refst.t0_des/SEML.us_em.T << endl;
        cout << "t0_des_mod                = " << t0_des_mod                << endl;

        filename = SEML.cs_em.F_PLOT+"Serv/projcu_order_20_tspan_";


        if(t0_des_mod >= 0 && t0_des_mod < 0.25)
        {
            if(SEML.li_SEM == 2) filename += "0T_025T_FINAL.bin";
            else filename += "0T_025T_FINAL_dest_L1.bin";
        }
        else if(t0_des_mod >= 0.25 && t0_des_mod < 0.5)
        {
            filename += "025T_05T_FINAL.bin";
        }
        else if(t0_des_mod >= 0.5 && t0_des_mod < 0.75)
        {
            filename += "05T_075T_FINAL.bin";
        }
        else
        {
            filename += "075T_T_FINAL.bin";
        }

        cout << "ooseleml2seml. The data will be retrieved from the server-computed file named:" << endl;
        cout << filename << endl;

        readAndInterpolateIntProjCU_bin(filename, refst.t0_des, t0_EM, tf_EM,
                                        s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM, s5_CMU_EM,
                                        pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM,
                                        sortId);
    }
    else
    {
        int type = TYPE_MAN_PROJ;
        switch(refst.dim)
        {
        case REF_PLANAR:
            //type = refst.isFromServer? TYPE_MAN_PROJ_FROM_SERVER:TYPE_MAN_PROJ;
            type = TYPE_MAN_PROJ;
            break;

        case REF_3D:
            //type = refst.isFromServer? TYPE_MAN_PROJ_3D_FROM_SERVER:TYPE_MAN_PROJ_3D;
            type = TYPE_MAN_PROJ_3D;
            break;
        }
        filename = filenameCUM(OFTS_ORDER, type, SEML.li_SEM);
        readIntProjCU_bin(filename, t0_EM, tf_EM,
                          s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM, s5_CMU_EM,
                          pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM,
                          sortId);
    }



    //====================================================================================
    // 3. Select given intervals for the inputs
    //====================================================================================
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;

    if(refst.isLimUD)
    {
        cout << "-----------------------------------------------------" << endl;
        cout << "   ooseleml2seml. User-defined interval of research  " << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "Enter a value for s1_CMU_EM_MIN: ";
        cin >> s1_CMU_EM_MIN;
        cout << "Enter a value for s1_CMU_EM_MAX: ";
        cin >> s1_CMU_EM_MAX;
        cout << "Enter a value for s3_CMU_EM_MIN: ";
        cin >> s3_CMU_EM_MIN;
        cout << "Enter a value for s3_CMU_EM_MAX: ";
        cin >> s3_CMU_EM_MAX;
    }
    else
    {
        s1_CMU_EM_MIN = refst.s1_CMU_EM_MIN;
        s1_CMU_EM_MAX = refst.s1_CMU_EM_MAX;
        s3_CMU_EM_MIN = refst.s3_CMU_EM_MIN;
        s3_CMU_EM_MAX = refst.s3_CMU_EM_MAX;
    }

    //--------------------------------------
    // Subselection
    //--------------------------------------
    int kpor = 0;
    bool cst = 0;
    for(int kpos = 0; kpos < (int) sortId.size(); kpos++)
    {
        kpor = sortId[kpos];
        cst  = (s1_CMU_EM[kpor] > s1_CMU_EM_MIN) & (s1_CMU_EM[kpor] < s1_CMU_EM_MAX);
        cst  = cst & (s3_CMU_EM[kpor] > s3_CMU_EM_MIN) & (s3_CMU_EM[kpor] < s3_CMU_EM_MAX);
        if(cst)
        {
            t0_EM_R.push_back(t0_EM[kpor]);
            tf_EM_R.push_back(tf_EM[kpor]);
            s1_CMU_EM_R.push_back(s1_CMU_EM[kpor]);
            s2_CMU_EM_R.push_back(s2_CMU_EM[kpor]);
            s3_CMU_EM_R.push_back(s3_CMU_EM[kpor]);
            s4_CMU_EM_R.push_back(s4_CMU_EM[kpor]);
            s5_CMU_EM_R.push_back(s5_CMU_EM[kpor]);
            pmin_dist_SEM_R.push_back(pmin_dist_SEM[kpor]);
            s1_CM_SEM_R.push_back(s1_CM_SEM[kpor]);
            s2_CM_SEM_R.push_back(s2_CM_SEM[kpor]);
            s3_CM_SEM_R.push_back(s3_CM_SEM[kpor]);
            s4_CM_SEM_R.push_back(s4_CM_SEM[kpor]);
        }
    }
    //Sort the subselection in sortId_R
    sortId_R = sort_indexes(pmin_dist_SEM_R);

    //====================================================================================
    // 5. The best solution in the subselection will serve as the first guess.
    //====================================================================================
    int kpos = sortId_R[0];

    //====================================================================================
    // 6. Initialize local variables: EM
    //====================================================================================
    //RCM coordinates
    st_EM[0] = s1_CMU_EM_R[kpos];
    st_EM[1] = s2_CMU_EM_R[kpos];
    st_EM[2] = s3_CMU_EM_R[kpos];
    st_EM[3] = s4_CMU_EM_R[kpos];
    st_EM[4] = PROJ_EPSILON;
    //Time
    t_EM[0] = t0_EM_R[kpos];
    t_EM[1] = tf_EM_R[kpos];


    //====================================================================================
    // 7. Initialize local variables: SEM
    //====================================================================================
    //RCM coordinates
    st_SEM[0] = s1_CM_SEM_R[kpos];
    st_SEM[1] = s2_CM_SEM_R[kpos];
    st_SEM[2] = s3_CM_SEM_R[kpos];
    st_SEM[3] = s4_CM_SEM_R[kpos];
    st_SEM[4] = 0.0;
    //Initial time in SEM units
    *t0_SEM = tf_EM_R[kpos]*SEML.us_em.ns;

    //Minimum projection distance
    *pmin_dist_SEM_out = pmin_dist_SEM_R[kpos];


    return FTC_SUCCESS;
}

//----------------------------------------------------------------------------------------
//         Brick B: generate a first guess (either a single unstable manifold leg or
//              a complete trajectory EML2 orbit + man leg + SEMLi orbit).
//----------------------------------------------------------------------------------------
/**
 *  \brief Computing the first guess for the connection leg between the orbit orbit_EM and
 *         and the orbit orbit_SEM. Only the manifold leg is computed.
 *         The grid size is returned.
 **/
int oomanfgeml2seml(double** y_traj, double* t_traj,
                    Orbit& orbit_EM, Orbit& orbit_SEM,
                    int dcs, int coord_type, int man_grid_size,
                    RefSt& refst)
{
    //====================================================================================
    // 1.  Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Local variables to store the manifold leg
    //------------------------------------------------------------------------------------
    double** y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double*  t_man_coord  = dvector(0, man_grid_size);

    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    cout << " oomanfgeml2seml. Compute the manifold leg..."  << endl;
    //------------------------------------------------------------------------------------
    // Integration: from orbit_EM.t0 to orbit_EM.tf
    // Note that the collisionner is not used at this step. It is used later, after
    // the refinement procedures
    //------------------------------------------------------------------------------------
    int man_index = man_grid_size;
    int ode78coll;

    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    {
        //----------------------------------------------------------------------------
        // Integration using ode78
        //----------------------------------------------------------------------------
        ode78(y_man_coord, t_man_coord, &ode78coll, orbit_EM.getT0(), orbit_EM.getTf(),
              orbit_EM.getZ0(), 6, man_grid_size, dcs, NCEM, coord_type);

        //----------------------------------------------------------------------------
        // Save All BUT the last point
        //----------------------------------------------------------------------------
        for(int kman = 0; kman < man_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
            t_traj[kman] = t_man_coord[kman];
        }
        break;
    }

    case REF_VAR_GRID:
    {
        //----------------------------------------------------------------------------
        // Integration using ode78_qbcp_vg (variable grid)
        //----------------------------------------------------------------------------
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, &ode78coll, orbit_EM.getT0(),
                                  orbit_EM.getTf(), orbit_EM.getZ0(), 6, man_grid_size,
                                  dcs, NCEM, coord_type, man_grid_size);
        //----------------------------------------------------------------------------
        // Save All BUT the last point
        //----------------------------------------------------------------------------
        for(int kman = 0; kman < man_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
            t_traj[kman] = t_man_coord[kman];
        }
        break;
    }

    case REF_GIVEN_GRID:
    {
        //----------------------------------------------------------------------------
        // Nothing to do here, since a previous computation has been performed
        //----------------------------------------------------------------------------
        //            double *t_traj_loc = dvector(0, man_grid_size);
        //            qbcp_coc_time(t_traj, t_traj_loc, man_grid_size, coord_type, NCEM);
        //
        //            cout << "orbit_EM.getT0() = " << orbit_EM.getT0() << endl;
        //            cout << "t_traj[0] = " <<  t_traj_loc[0]  << endl;
        //
        //            ode78_qbcp_gg(y_man_coord, t_man_coord, &ode78coll, t_traj_loc, orbit_EM.getZ0(),
        //                          6, man_grid_size, dcs, NCEM, coord_type);
        //
        //            free_dvector(t_traj_loc, 0, man_grid_size);
        break;
    }
    }


    //====================================================================================
    // 5. Compute the first point of the final SEMLi orbit and add it to the manifold leg
    //====================================================================================
    cout << " oomanfgeml2seml. Compute the first point of the final SEMLi orbit..."<< endl;

    //------------------------------------------------------------------------------------
    // Save the last point at position man_index
    //------------------------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.getT0(), orbit_SEM.getZ0(), yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_index] = yout[i];
    t_traj[man_index] = tout;


    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    //The (true) grid size is returned
    return man_index;
}

/**
 *  \brief Computing the first guess for the connection leg between the orbit orbit_EM and
 *         and the orbit orbit_SEM. The complete trajectory (EML2 + man leg + SEMLi)
 *         is computed.
 *         The grid size is returned.
 **/
int oofgeml2seml(double** y_traj, double* t_traj,
                 double** y_traj_comp, double* t_traj_comp,
                 Orbit& orbit_EM, Orbit& orbit_SEM,
                 int dcs, int coord_type, int grid_points_des[3], int grid_points_eff[3], int max_grid,
                 RefSt& refst, gnuplot_ctrl* h2, gnuplot_ctrl* h3)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "oofgeml2seml";

    //Maximum number of points on ONE leg
    int grid_points_eff_max = max(grid_points_eff[0], max(grid_points_eff[1], grid_points_eff[2]));

    //------------------------------------------------------------------------------------
    // Local variables to store the manifold leg
    //------------------------------------------------------------------------------------
    double** y_man_NCEM  = dmatrix(0, 5, 0, grid_points_eff[0]);
    double* t_man_EM     = dvector(0, grid_points_eff[0]);

    double** y_man_NCSEM = dmatrix(0, 5, 0, grid_points_eff[2]);
    double* t_man_SEM    = dvector(0, grid_points_eff[2]);

    double** y_man_coord = dmatrix(0, 5, 0, grid_points_eff_max);
    double** y_man_comp  = dmatrix(0, 5, 0, grid_points_eff_max);
    double* t_man_coord  = dvector(0, grid_points_eff_max);
    double* t_man_comp   = dvector(0, grid_points_eff_max);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_eml_EM   = refst.tspan_EM;    //TOF on EML2 orbit
    double tof_seml_SEM = refst.tspan_SEM;   //TOF on SEMLi orbit

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 2. Compute the initial orbit
    //====================================================================================
    cout << fname << ". Compute the initial orbit..."  << endl;
    //------------------------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //------------------------------------------------------------------------------------
    orbit_EM.setSi(0.0, 4);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //------------------------------------------------------------------------------------
    //Integration on des_grid_size+1 fixed grid
    //------------------------------------------------------------------------------------
    int em_index = grid_points_eff[0];
    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    case REF_GIVEN_GRID:
    {
        int output = orbit_EM.traj_int_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, grid_points_eff[0], true);
        //--------------------------------------------------------------------------------
        //If output is strictly greater than 0, then the projection procedure inside
        //the integration went wrong, and the new indix is output.
        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //is no point stored in the data, except the first one, and the new indix is zero
        //--------------------------------------------------------------------------------
        if(output == ORBIT_EPROJ) em_index = 0;
        else if(output > 0) em_index = output;
        break;
    }

    case REF_VAR_GRID:
    {
        em_index = orbit_EM.traj_int_var_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, grid_points_eff[0], true);
        break;
    }

    }

    //------------------------------------------------------------------------------------
    //To coord_type coordinates
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, em_index, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp,  t_man_comp,  em_index, NCEM, comp_type);

    //------------------------------------------------------------------------------------
    // Save All BUT the last point
    //------------------------------------------------------------------------------------
    for(int kman = 0; kman < em_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][em_index-kman];
        t_traj[kman] = t_man_coord[em_index-kman];
    }

    //------------------------------------------------------------------------------------
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(h2, y_man_coord, em_index+1, (char*)"", "points", "5", "3", 5);
    gnuplot_plot_X(h3, y_man_comp, em_index+1, (char*)"", "points", "5", "3", 5);

    //====================================================================================
    // 4.3 Compute the manifold leg
    //====================================================================================
    cout << fname << ". Compute the manifold leg..."  << endl;

    double tf_EM = orbit_SEM.getT0()/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //------------------------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //------------------------------------------------------------------------------------
    orbit_EM.setSi(PROJ_EPSILON, 4);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //------------------------------------------------------------------------------------
    // Integration
    // Note that the collisionner is not used at this step. It is used later, after
    // the refinement procedures
    //------------------------------------------------------------------------------------
    int ode78coll;
    int man_index = grid_points_eff[1];

    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    case REF_GIVEN_GRID:
    {
        ode78(y_man_coord, t_man_coord, &ode78coll, orbit_EM.getT0(), tf_EM, orbit_EM.getZ0(), 6, grid_points_eff[1], dcs, NCEM, coord_type);
        break;
    }

    case REF_VAR_GRID:
    {
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, &ode78coll, orbit_EM.getT0(), tf_EM, orbit_EM.getZ0(), 6, grid_points_eff[1], dcs, NCEM, coord_type, -1);
        break;
    }

    }

    //------------------------------------------------------------------------------------
    //To comp_type coordinates
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_man_coord, t_man_coord, y_man_comp, t_man_comp, man_index, coord_type, comp_type);


    //------------------------------------------------------------------------------------
    // Save All BUT the last point
    //------------------------------------------------------------------------------------
    for(int kman = 0; kman < man_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][em_index + kman] = y_man_coord[i][kman];
        t_traj[em_index + kman] = t_man_coord[kman];
    }

    //------------------------------------------------------------------------------------
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(h2, y_man_coord, man_index+1, (char*)"Man", "points", "1", "3", 4);
    gnuplot_plot_X(h3, y_man_comp, man_index+1, (char*)"", "points", "1", "3", 4);

    //====================================================================================
    // 4.4 Compute the final SEMLi orbit
    //====================================================================================
    cout << fname << ". Compute the final SEMLi orbit..."  << endl;
    //------------------------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //------------------------------------------------------------------------------------
    orbit_SEM.setSi(0.0, 4);
    orbit_SEM.update_ic(orbit_SEM.getSi());


    //------------------------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //------------------------------------------------------------------------------------
    int sem_index = grid_points_eff[2];

    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    case REF_GIVEN_GRID:
    {
        int output = orbit_SEM.traj_int_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, grid_points_eff[2], true);
        //--------------------------------------------------------------------------------
        //If output is strictly greater than 0, then the projection procedure inside
        //the integration went wrong, and the new indix is output.
        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //is no point stored in the data, except the first one, and the new indix is zero
        //--------------------------------------------------------------------------------
        if(output == ORBIT_EPROJ) sem_index = 0;
        else if(output > 0) sem_index = output;
        break;
    }

    case REF_VAR_GRID:
    {
        sem_index = orbit_SEM.traj_int_var_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, grid_points_eff[2], 1);
        break;
    }
    }


    //------------------------------------------------------------------------------------
    //To coord_type coordinates
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, sem_index, NCSEM, coord_type);
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_comp, t_man_comp, sem_index, NCSEM, comp_type);


    //------------------------------------------------------------------------------------
    // Save ALL points
    //------------------------------------------------------------------------------------
    for(int kman = 0; kman <= sem_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][em_index + man_index + kman] = y_man_coord[i][kman];
        t_traj[em_index + man_index + kman] = t_man_coord[kman];
    }

    //------------------------------------------------------------------------------------
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(h2, y_man_coord, sem_index+1, (char*)"SEMLi", "points", "1", "3", 6);
    gnuplot_plot_X(h3, y_man_comp, sem_index+1, (char*)"", "points", "1", "3", 6);

    //====================================================================================
    // Entire size is: no more than sum(grid_points_des) points or so.
    //====================================================================================
    int final_index = 0;
    if(refst.grid == REF_VAR_GRID)
    {
        final_index = em_index + man_index + sem_index;
        int freq = final_index/(grid_points_des[0]+grid_points_des[1]+grid_points_des[2]);
        cout << "----------------------------------"<< endl;
        cout << "Number of MSD points:        "     << endl;
        cout << "At EML:           " << em_index    << endl;
        cout << "During coast arc: " << man_index   << endl;
        cout << "At SEML:          " << sem_index   << endl;
        cout << "Total:            " << final_index << endl;
        cout << "Selecting 1/" << freq << " points" << endl;
        cout << "Subtotal:         " << (int) floor(final_index/freq) << endl;
        cout << "----------------------------------"<< endl;

        //------------------------------------------------------------------------------------
        // Subset of points in y_traj_comp/t_traj_comp
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= floor(final_index/freq); kman++)
        {
            for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
            t_traj_comp[kman] = t_traj[freq*kman];
        }

        //New final index
        final_index = floor(final_index/freq);

        //------------------------------------------------------------------------------------
        // Subset of points back in y_traj/t_traj
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= final_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
            t_traj[kman] = t_traj_comp[kman];
        }
    }
    else
    {
        //If the grid size is fixed, we still return the sum of the variable indices, because
        //some integration procedure could have gone wrong so that the returned indices
        //are smaller than their optimal desired value
        final_index = em_index + man_index + sem_index;
    }

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(y_man_NCEM, 0, 5, 0, grid_points_eff[0]);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, grid_points_eff[2]);
    free_dmatrix(y_man_coord, 0, 5, 0, grid_points_eff_max);
    free_dmatrix(y_man_comp, 0, 5, 0, grid_points_eff_max);
    free_dvector(t_man_EM, 0, grid_points_eff[0]);
    free_dvector(t_man_SEM, 0, grid_points_eff[2]);
    free_dvector(t_man_coord, 0, grid_points_eff_max);
    free_dvector(t_man_comp, 0, grid_points_eff_max);

    return final_index;
}

/**
 *  \brief Get the complementary coordinates associated to the coordinates coord_type.
 **/
int comp_coord_typ(int coord_type)
{
    int comp_type = 0;
    switch(coord_type)
    {
    case PSEM:
        comp_type = PEM;
        break;
    case NCSEM:
    case VNCSEM:
        comp_type = NCEM;
        break;
    case PEM:
        comp_type = PSEM;
        break;
    case NCEM:
    case VNCEM:
        comp_type = NCSEM;
        break;
    }
    return comp_type;

}


//========================================================================================
//
//         Refinement of solutions: Complete trajectory
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the WHOLE trajectory
 *         (EML2 orbit + manifold leg + SEMLi orbit).
 **/
int oocomprefft3d(int grid_freq_days[3], int coord_type,
                  Orbit& orbit_EM, Orbit& orbit_SEM,
                  RefSt& refst)

{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "oocomprefft3d";
    cout << fname << ". Initialization of the local variables..."  << endl;

    // Status along the computation
    int status = 0;

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //------------------------------------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //------------------------------------------------------------------------------------
    //Complementary coordinate type
    //------------------------------------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //------------------------------------------------------------------------------------
    //Notable points in SEM & EM systems
    //------------------------------------------------------------------------------------
    notablePoints(h2, h3, coord_type);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Initialize the number of points on the grids. If the number of points is fixed,
    // the solution is straightforward. If the number of points is desired variable,
    // the default number are set to an arbitrarily high number of points.
    //------------------------------------------------------------------------------------
    cout << " oocomprefft3d. Initialize the number of points on the grids..."  << endl;
    int max_grid    = 50000;
    int grid_points_des[3];

    //Desired number of points
    grid_points_des[0] = refst.tspan_EM/(2*M_PI*86400*grid_freq_days[0]/SEML.cs_em.cr3bp.T);
    grid_points_des[1] = (orbit_SEM.getT0()/SEML.us_em.ns - orbit_EM.getT0())/(2*M_PI*86400*grid_freq_days[1]/SEML.cs_em.cr3bp.T);
    grid_points_des[2] = refst.tspan_SEM/(2*M_PI*86400*grid_freq_days[2]/SEML.cs_sem.cr3bp.T);

    cout << "Desired frequency on leg 1: "  << grid_freq_days[0]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[0] << endl    << endl;
    cout << "Desired frequency on leg 2: "  << grid_freq_days[1]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[1] << endl    << endl;
    cout << "Desired frequency on leg 3: "  << grid_freq_days[2]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[2] << endl    << endl;


    //Effective number of points for the following computation: equal to grid_points_des
    //if the number of points if fixed, equal to max_grid otherwise
    int grid_points_eff[3];
    grid_points_eff[0]  = (refst.grid == REF_VAR_GRID) ? max_grid:grid_points_des[0];
    grid_points_eff[1]  = (refst.grid == REF_VAR_GRID) ? max_grid:grid_points_des[1];
    grid_points_eff[2]  = (refst.grid == REF_VAR_GRID) ? max_grid:grid_points_des[2];

    //Number of points on the entire trajectory
    int traj_grid_eff   = grid_points_eff[0]+grid_points_eff[1]+grid_points_eff[2];

    //------------------------------------------------------------------------------------
    //To store final data
    //------------------------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, traj_grid_eff);
    double*  t_traj  = dvector(0, traj_grid_eff);

    double** y_traj_comp  = dmatrix(0, 41, 0, traj_grid_eff);
    double*  t_traj_comp  = dvector(0, traj_grid_eff);

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, traj_grid_eff);
    double* t_traj_n    = dvector(0, traj_grid_eff);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================
    int final_index = oofgeml2seml(y_traj, t_traj, y_traj_comp, t_traj_comp,
                                   orbit_EM, orbit_SEM, dcs, coord_type,
                                   grid_points_des,
                                   grid_points_eff,
                                   max_grid,
                                   refst, h2, h3);


    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    cout << " oocomprefft3d. Initial trajectory, on a grid..."  << endl;
    int mPlot = 100;

    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);

    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj, t_traj, final_index, mPlot, coord_type, 0.0, 0.0,
                     coord_type, h2, comp_type,  h3, 7, "Initial guess");


    //====================================================================================
    // 5. Differential correction & final trajectory
    //====================================================================================
    cout << fname << ". Differential correction procedure.              "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    pressEnter(refst.isFlagOn);

    int isPlotted   = 0;
    status = multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //status =multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //------------------------------------------------------------------------------------
    // If something went bad, an error is returned
    //------------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". Error during the differential correction procedure.";
        cerr << " ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }


    //------------------------------------------------------------------------------------
    // Final trajectory, on a grid
    //------------------------------------------------------------------------------------
    cout << fname << ". Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type, 0.0, 0.0,
                     coord_type, h2, comp_type, h3, 3, "Final guess");

    //====================================================================================
    // 6. Store the data, in sys units
    //====================================================================================
    cout << fname << ". Storing the refined solution in system units "  << endl;
    cout << " The solution will be stored in " << filenameCUM(OFTS_ORDER, TYPE_COMP_FOR_JPL, SEML.li_SEM) << endl;
    cout << "-----------------------------------------------------------"  << endl;
    writeCOMP_txt(t_traj_n, y_traj_n, final_index);


    //====================================================================================
    // 7. JPL refinement
    //====================================================================================
    if(refst.isJPL)
    {
        status = oojplrefft3d(coord_type, refst);

        //--------------------------------------------------------------------------------
        // If something went bad, an error is returned
        //--------------------------------------------------------------------------------
        if(status)
        {
            cerr << fname << ". Error during the differential correction procedure.";
            cerr << " ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }
    }



    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    pressEnter(refst.isFlagOn);
    gnuplot_close(h2);
    gnuplot_close(h3);


    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_eff);
    free_dvector(t_traj, 0, traj_grid_eff);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_eff);
    free_dvector(t_traj_comp, 0, traj_grid_eff);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_eff);
    free_dvector(t_traj_n, 0, traj_grid_eff);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v  , 0, mPlot*final_index);

    return FTC_SUCCESS;
}

//========================================================================================
//
//         Refinement of solutions: to JPL
//
//========================================================================================
/**
 *  \brief Refine a given output of oocomprefft3d into JPL ephemerides.
 **/
int oojplrefft3d(int coord_type, RefSt& refst)
{
    //====================================================================================
    // 1. Read the data, in sys units
    //====================================================================================
    string fname = "oojplrefft3d";
    cout << " oojplrefft3d. Read the data, in sys units..."  << endl;
    int status = 0;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    int final_index = getLengthCOMP_txt();
    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    readCOMP_txt(t_traj_n, y_traj_n, final_index);

    //====================================================================================
    // 2. Local parameters
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Complementary coordinate type
    //------------------------------------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //------------------------------------------------------------------------------------
    // State and time vectors
    //------------------------------------------------------------------------------------
    int mPlot = 10;
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);

    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    //int dcs  = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //------------------------------------------------------------------------------------
    //Notable points in SEM & EM systems
    //------------------------------------------------------------------------------------
    notablePoints(h2, h3, coord_type);

    //====================================================================================
    // 5. Initial trajectory, on a grid
    //====================================================================================
    cout << " oojplrefft3d. Initial trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type, 0.0, 0.0,
                     coord_type, h2, comp_type, h3, 3, "Initial trajectory");

    //------------------------------------------------------------------------------------
    //Free at this point, because final_index is gonna change
    //------------------------------------------------------------------------------------
    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v, 0, mPlot*final_index);

    //====================================================================================
    // 5. JPL refinement
    //====================================================================================
    pressEnter(refst.isFlagOn);

    //====================================================================================
    // Select between VECLI and J2000
    //====================================================================================
    int choice = 0;
    //If refst.djplcoord does not contain a valid choice (e.g. -1), we ask the user
    if(refst.djplcoord != 0 && refst.djplcoord != 1  && refst.djplcoord != 2)
    {
        cout << "Enter 0 for VECLI, 1 for J2000, 2 for NJ2000: ";
        cin >> choice;
    }
    else
    {
        choice = refst.djplcoord;
    }

    int coord_int = 0, fwrk_int = 0;
    switch(choice)
    {
    case 0:
        coord_int = VECLI;
        fwrk_int  = I_ECLI;
        break;
    case 1:
        coord_int = J2000;
        fwrk_int  = I_J2000;
        break;
    case 2:
        coord_int = NJ2000;
        fwrk_int  = I_NJ2000;
        changeDCS(SEML, fwrk);
        cout << "NJ2000 has been chosen: changeDCS is applied to match " << endl;
        cout << " the framework associated to coord_type.";
        break;
    default:
        cout << "wrong input. NJ2000 is used by default." << endl;
        cout << "NJ2000 has been chosen: changeDCS is applied to match " << endl;
        cout << " the framework associated to coord_type.";
        coord_int = NJ2000;
        fwrk_int  = I_NJ2000;
        changeDCS(SEML, fwrk);
        break;
    }

    //====================================================================================
    // Initialize SPICE kernerls & VF
    //====================================================================================
    gnuplot_ctrl* h4 = gnuplot_init();
    cout << fname << ". Initialize SPICE kernerls..." << endl;
    furnsh_c("spice/kernels/metakernel.furnsh");
    int shift = 0;


    //====================================================================================
    // Search for best fit in  JPL ephemerides
    //====================================================================================
    cout << fname << ". Search for best fit in JPL DE430..." << endl;

    //----------------------------------------------------------
    //Get the best fit at t = t_traj_n[shift]. et0 is in seconds
    //----------------------------------------------------------
    double et0, et;
    double tsys0 = t_traj_n[shift];
    qbcp2jpl(tsys0, &et0, coord_type);

    //----------------------------------------------------------
    //Display the value in "D" format
    //----------------------------------------------------------
    //ConstSpiceChar*  picture  = "Wkd Mon DD HH:MN:SC PDT YYYY ::UTC-7";
    SpiceChar output[50];
    et2utc_c (et0, "C", 6, 51, output);
    cout << fname << ". The best fit has been obtained for the following date:" << endl;
    printf ( "%s\n", output);
    cout << "tsys0 = " << tsys0 << endl;
    cout << "et0 = "   << et0   << endl;

    //------------------------------------------------------------------------------------
    // Complementary value ot tsys0
    //------------------------------------------------------------------------------------
    double tsys0_comp = tsys0;
    switch(eph_coord(comp_type))
    {
    case VEM:
        tsys0_comp /= SEML.us_em.ns;
        break;

    case VSEM:
        tsys0_comp *= SEML.us_em.ns;
        break;
    }

    //====================================================================================
    // Change of coordinates
    //====================================================================================
    cout << fname << ". Change of coordinates syn -> ecliptic..." << endl;

    //------------------------------------------------------------------------------------
    //Initialize the vectors that will be used throughout the computation
    //------------------------------------------------------------------------------------
    double** y_traj_jpl = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl  = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //    Computes a first guess in ephemerides coordinates, switching between the
    //    Earth-Moon and Sun-Earth plane of motion when the discrepancy between the two
    //    coordinates system is minimal.
    //------------------------------------------------------------------------------------
    //    final_index = oojplfg3d_switch(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index,
    //                                   coord_type, comp_type, coord_int,
    //                                   et0, tsys0, tsys0_comp);


    //------------------------------------------------------------------------------------
    //Same as oojplfg3d_switch, with a refinement of the position of the minimum.
    //------------------------------------------------------------------------------------
    //    int mRef = 50;
    //    final_index = oojplfg3d_super_switch(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index,
    //                                   coord_type, comp_type, coord_int, mRef,
    //                                   et0, tsys0, tsys0_comp);


    //------------------------------------------------------------------------------------
    //Same as oojplfg3d_switch, with an interpolation before & after the switching point.
    //------------------------------------------------------------------------------------
    int mRef = 5;
    final_index = oojplfg3d_interpolation(y_traj_n, t_traj_n, &y_traj_jpl, &t_traj_jpl, final_index,
                                          coord_type, comp_type, coord_int, mRef,
                                          et0, tsys0, tsys0_comp);
    cout << "final_index = " << final_index << endl;

    //------------------------------------------------------------------------------------
    //   Once the first guess is computed, we can safely initialize the other vectors and
    //   matrices, since final_index will not change anymore.
    //------------------------------------------------------------------------------------
    double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl_n  = dvector(0, final_index);

    double** y_jpl_temp   = dmatrix(0, 41, 0, final_index);
    double* t_jpl_temp    = dvector(0, final_index);

    double* et_traj_jpl   = dvector(0, final_index);

    double** y_earth_spice = dmatrix(0, 5, 0, final_index);
    double** y_moon_spice  = dmatrix(0, 5, 0, final_index);
    double** y_l2_spice    = dmatrix(0, 5, 0, final_index);

    ymc_v       = dmatrix(0, 5, 0, mPlot*final_index);
    tmc_v       = dvector(0, mPlot*final_index);
    ymc_comp_v  = dmatrix(0, 5, 0, mPlot*final_index);
    tmc_comp_v  = dvector(0, mPlot*final_index);


    //------------------------------------------------------------------------------------
    //Time in seconds: €€TODO: adapt the other cases, only NJ2000 for now!!
    //------------------------------------------------------------------------------------
    for(int p = 0; p <= final_index; p++)
    {
        et_traj_jpl[p] = t_traj_jpl[p]/SEML.ss.n;
    }

    //====================================================================================
    // Plotting Moon +  Earth + L2
    //====================================================================================
    cout << fname << ". Plotting Moon +  Earth + L2..." << endl;

    //------------------------------------------------------------------------------------
    //Position of the Moon +  Earth + L2
    //------------------------------------------------------------------------------------
    double lt, YV[6];
    for(int p = 0; p <= final_index; p++)
    {
        //Current time in seconds
        et = et_traj_jpl[p];
        //EARTH
        spkezr_c ("EARTH", et, DEFFRAME,  "NONE", DEFOBS, YV, &lt);
        for(int i = 0; i < 6; i++) y_earth_spice[i][p] = YV[i];
        //MOON
        spkezr_c ("MOON", et, DEFFRAME,  "NONE", DEFOBS, YV, &lt);
        for(int i = 0; i < 6; i++) y_moon_spice[i][p] = YV[i];
        //L2
        spkez_c (392, et, DEFFRAME,  "NONE", 0, YV, &lt);
        for(int i = 0; i < 6; i++) y_l2_spice[i][p] = YV[i];
    }

    //------------------------------------------------------------------------------------
    //Position of the Sun +  Earth + Initial guess in coord_type coordinates
    //------------------------------------------------------------------------------------
    ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    gnuplot_plot_X(h2, y_jpl_temp, final_index+1, (char*) "EARTH", "points", "3", "2", 8);

    ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    gnuplot_plot_X(h2, y_jpl_temp, final_index+1, (char*) "MOON", "points", "4", "2", 8);

    ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    gnuplot_plot_X(h2, y_jpl_temp, final_index+1, (char*) "SEMLi", "points", "5", "2", 8);

    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    pressEnter(refst.isFlagOn);
    cout << fname << ". Plotting Initial Guess..." << endl;

    //------------------------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    plottrajsegbyseg(y_traj_jpl, t_traj_jpl, final_index, mPlot, coord_int,
                     et0, tsys0, coord_type, h2, comp_type,  h3, 6,
                     "Initial guess in JPL ephemerides");

    //====================================================================================
    // 6.5 Differential correction
    //====================================================================================
    pressEnter(refst.isFlagOn);
    cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;
    status  = multiple_shooting_direct_variable_time(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int,  5e-10, true, h4);
    //status  = multiple_shooting_direct(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int, true, h4);
    //Plot
    gnuplot_plot_X(h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);

    //------------------------------------------------------------------------------------
    // If something went bad, an error is returned
    //------------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". Error during the differential correction procedure.";
        cerr << " ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }


    //====================================================================================
    //Final trajectory on lines, segment by segment
    //====================================================================================
    plottrajsegbyseg(y_traj_jpl_n, t_traj_jpl_n, final_index, mPlot, coord_int,
                     et0, tsys0, coord_type, h2, comp_type,  h3, 2,
                     "Final guess in JPL ephemerides");


    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    pressEnter(refst.isFlagOn);
    gnuplot_close(h2);
    gnuplot_close(h3);
    gnuplot_close(h4);

    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj_n, 0, 41, 0, final_index);
    free_dvector(t_traj_n, 0, final_index);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(y_jpl_temp, 0, 41, 0, final_index);
    free_dvector(t_jpl_temp, 0, final_index);
    free_dmatrix(y_traj_jpl_n, 0, 41, 0, final_index);
    free_dvector(t_traj_jpl_n, 0, final_index);
    free_dvector(et_traj_jpl, 0, final_index);

    return FTC_SUCCESS;
}


/**
 *  \brief Refine a given output of oocomprefft3d into Inertial Coordinates, then into
 *         JPL coordinates.
 **/
int oointojplrefft3d(int coord_type, RefSt& refst)
{
    //====================================================================================
    // 1. Read the data, in coord_type units
    //====================================================================================
    string fname = "ooinrefft3d";
    cout << fname << ". Read the data, in coord_type units..."  << endl;
    int status = 0;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    int final_index = getLengthCOMP_txt();
    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    readCOMP_txt(t_traj_n, y_traj_n, final_index);

    //------------------------------------------------------------------------------------
    //Arbitrarily reduce the number of points by a factor
    //------------------------------------------------------------------------------------
    cout << fname << ". Reduce the number of points..."  << endl;
    int factor = 3;
    final_index = final_index/factor;

    for(int k = 0; k <= final_index; k++)
    {
        for(int i = 0; i <6; i++) y_traj_n[i][k] = y_traj_n[i][factor*k];
        t_traj_n[k] = t_traj_n[factor*k];
    }


    //====================================================================================
    // 2. Local parameters
    //====================================================================================
    //Complementary coordinate type
    int comp_type = comp_coord_typ(coord_type);

    //------------------------------------------------------------------------------------
    // State and time vectors
    //------------------------------------------------------------------------------------
    int mPlot = 10;
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);

    double** y_jpl_temp   = dmatrix(0, 41, 0, final_index);
    double* t_jpl_temp    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //Get the default framework from the coord_type
    //------------------------------------------------------------------------------------
    //int dcs  = default_coordinate_system(coord_type);
    int fwrk = default_framework(coord_type);

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //Gnuplot window
    gnuplot_ctrl* h2, *h3, * h4;

    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");
    gnuplot_cmd(h4, "set title \"Continuation steps\" ");

    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");
    gnuplot_cmd(h4, "set grid");

    //Notable points in SEM & EM systems
    notablePoints(h2, h3, coord_type);

    //Color for plots
    int color = 1;

    //====================================================================================
    // 5. Initial trajectory, on a grid
    //====================================================================================
    cout << " oojplrefft3d. Initial trajectory, on a grid..."  << endl;
    //Initial trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type,
                     0.0, 0.0, coord_type, h2, comp_type,  h3,color++,
                     "Initial guess (NCSEM)");

    //====================================================================================
    // 5. INERTIAL refinement
    //====================================================================================
    //pressEnter(refst.isFlagOn);

    //====================================================================================
    // Select ECISEM
    //====================================================================================
    int coord_int = ECISEM;
    //int fwrk_int  = I_ECISEM;

    //====================================================================================
    // Change FOCUS
    //====================================================================================
    cout << "ECISEM has been chosen: changeDCS is applied to match the framework associated to coord_type.";
    changeDCS(SEML, fwrk);

    //====================================================================================
    // Change of coordinates
    //====================================================================================
    cout << fname << ". Change of coordinates syn -> ecliptic..." << endl;

    //------------------------------------------------------------------------------------
    //Initialize the vectors that will be used throughout the computation
    //------------------------------------------------------------------------------------
    double** y_traj_ecisem = dmatrix(0, 41, 0, final_index);
    double* t_traj_ecisem  = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //COC: coord_type -> coord_int
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_ecisem, t_traj_ecisem, final_index, coord_type, coord_int);


    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    //pressEnter(refst.isFlagOn);
    cout << fname << ". Plotting Initial Guess..." << endl;
    //------------------------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    plottrajsegbyseg(y_traj_ecisem, t_traj_ecisem, final_index, mPlot, coord_int,
                     0.0, 0.0, coord_type, h2, comp_type,  h3, color++,
                     "Initial guess (ECISEM)");


    //====================================================================================
    // 6.5 Differential correction, in ECISEM framework
    //====================================================================================
    //pressEnter(refst.isFlagOn);
    cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;
    //DiffCorr
    status  = multiple_shooting_direct(y_traj_ecisem, t_traj_ecisem, y_traj_ecisem, t_traj_ecisem, 42, final_index, coord_int, true, h4);
    //Plot
    gnuplot_plot_X(h4, y_traj_ecisem, final_index+1, (char*) "Final trajectory", "lines", "2", "2", color);

    //------------------------------------------------------------------------------------
    // If something went bad, an error is returned
    //------------------------------------------------------------------------------------
    if(status)
    {
        cerr << fname << ". Error during the differential correction procedure.";
        cerr << " ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }


    //====================================================================================
    //Final trajectory on lines, segment by segment, saved on txt file
    //====================================================================================
    plottrajsegbyseg(y_traj_ecisem, t_traj_ecisem, final_index, mPlot, coord_int,
                     0.0, 0.0, coord_type, h2, comp_type,  h3, color++,
                     "Final guess (ECISEM)");


    //====================================================================================
    // JPL: initialization and time selection
    //====================================================================================
    cout << fname << ". Initialize SPICE kernerls..." << endl;
    furnsh_c("spice/kernels/metakernel.furnsh");

    //------------------------------------------------------------------------------------
    //Get the best fit at tsys0 = 0.0. et0 is in seconds. After this step,
    //et0 contains the synchronized time defined by:
    //At epoch et0, the positions of the Sun, Earth and Moon in JPL ephemerides match
    //as much as possible the position of the Sun, Earth and Moon in the QBCP at tsys0 = 0
    //------------------------------------------------------------------------------------
    double et0;
    double tsys0 = 0.0;//t_traj_ecisem[0];
    //qbcp2jpl(tsys0, &et0, coord_int);
    qbcp2jpl_inertial(tsys0, &et0, coord_int);

    //------------------------------------------------------------------------------------
    //Display the value in "D" format
    //------------------------------------------------------------------------------------
    //ConstSpiceChar*  picture  = "Wkd Mon DD HH:MN:SC PDT YYYY ::UTC-7";
    SpiceChar output[50];
    et2utc_c (et0, "C", 6, 51, output);
    cout << fname << ". The best fit has been obtained for the following date:" << endl;
    printf ( "%s\n", output);

    //------------------------------------------------------------------------------------
    //We have to define tsys0 in complementary EM units,
    //------------------------------------------------------------------------------------
    //double tsys0_comp = tsys0/SEML.us_em.ns;

    //------------------------------------------------------------------------------------
    //Time shift stored in SEML: it is needed to synchronize the QBCP and the JPL ephemerides
    //------------------------------------------------------------------------------------
    SEML.ss.tshift = et0*SEML.ss.n - tsys0;

    //====================================================================================
    //Compute the equivalent of the state along the trajectory in NJ2000 coordinates, to
    //compare with ECI QBCP coordinates!
    //====================================================================================
    // ECISEM -> NJ2000
    coord2necistate_vec(y_traj_n, t_traj_n, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0,  tsys0, eph_coord(coord_int), SEML.ss);
    //Plot
    gnuplot_plot_X(h4, y_jpl_temp, final_index+1, (char*)"Inertial state (NJ2000)", "lines", "1", "1", color);

    //Moon's position (NJ2000)
    SpiceDouble lt, RS1[6], RS2[6], RE[6], REM[6];
    spkezr_c ("MOON",  et0, DEFFRAME,  "NONE",  DEFOBS, RS1, &lt);
    spkezr_c ("EARTH",  et0, DEFFRAME,  "NONE",  DEFOBS, RE, &lt);
    spkezr_c ("EARTH MOON BARYCENTER", et0, DEFFRAME,  "NONE",  DEFOBS, REM, &lt);
    ecl2neci(RS1, REM, RS2, SEML.ss);
    RS2[2] = 0.0;
    gnuplot_plot_xyz(h4, RS2, RS2+1, RS2+2, 1, (char*)"Moon (NJ2000, projected)", "points", "3", "5", color);

    //Moon's position (Custom)
    double me = SEML.us_sem.me;
    double mm = SEML.us_sem.mm;
    double n  = SEML.us_sem.n;
    double ni = SEML.us_sem.ni;
    double ai = SEML.us_sem.ai;
    double r1 = creal(evz(SEML.cs_sem.zt, tsys0, n, ni, ai));
    double r2 = cimag(evz(SEML.cs_sem.zt, tsys0, n, ni, ai));
    double Pm[3] = {- me/(mm + me)* r1, - me/(mm + me)* r2, 0.0};
    gnuplot_plot_xyz(h4, Pm, Pm+1, Pm+2, 1, (char*)"Moon (ECISEM)", "points", "4", "5", color++);


    //====================================================================================
    //Initial trajectory on lines, computed from NJ2000 first guess
    //====================================================================================
    plottrajsegbyseg(y_jpl_temp, t_jpl_temp, final_index, mPlot, NJ2000,
                     et0, tsys0, coord_type, h2, comp_type,  h3, color++,
                     "Initial guess (NJ2000)");

    //====================================================================================
    //Continutation procedure
    //====================================================================================
    cout << "Continuation procedure" << endl;
    pressEnter(refst.isFlagOn);


    //------------------------------------------------------------------------------------
    // Initialization of the containers
    //------------------------------------------------------------------------------------
    double** y_traj_ecisem_c = dmatrix(0, 47, 0, final_index);
    double* t_traj_ecisem_c  = dvector(0, final_index);
    double* t_traj_ecisem_s  = dvector(0, final_index);


    for(int k = 0; k <= final_index; k++)
    {
        //--------------------------------------------------------------------------------
        //Copy segment by segment : time and state
        //--------------------------------------------------------------------------------
        t_traj_ecisem_c[k] = t_traj_ecisem[k];
        for(int i = 0; i < 42; i++)
        {
            y_traj_ecisem_c[i][k] = y_traj_ecisem[i][k];
        }

        //Add Identity matrix
        y_traj_ecisem_c[6][k]  = 1.0;
        y_traj_ecisem_c[13][k] = 1.0;
        y_traj_ecisem_c[20][k] = 1.0;
        y_traj_ecisem_c[27][k] = 1.0;
        y_traj_ecisem_c[34][k] = 1.0;
        y_traj_ecisem_c[41][k] = 1.0;

        //--------------------------------------------------------------------------------
        //The last 6 components the variational equations for epsilon
        //--------------------------------------------------------------------------------
        for(int i = 42; i< 48; i++) y_traj_ecisem_c[i][k] = 0.0;
    }

    //Null vector
    int nfv = 6*(final_index+1)+1;
    double nullvector[nfv];
    double de = 5e-3;

    SEML.epsilon = 0.0;
    //------------------------------------------------------------------------------------
    // First differential correction
    //------------------------------------------------------------------------------------
    multiple_shooting_direct_deps(y_traj_ecisem_c, t_traj_ecisem_c, y_traj_ecisem_c, t_traj_ecisem_c,
                                  nullvector, true, 48, final_index, coord_int, true, h4);



    gnuplot_ctrl* h1;
    h1 = gnuplot_init();


    //------------------------------------------------------------------------------------
    // Solution with epsilon as a parameter + continuation
    //------------------------------------------------------------------------------------
    //Loop
    int Nmax = 3000;
    int iter = 0;
    do
    {
        //New de if necessary
        de = SEML.epsilon+de*nullvector[nfv-1] <= 1.0? de:(1.0-SEML.epsilon)/nullvector[nfv-1];

        //--------------------------------------------------------------------------------
        //  New state, with a correction from null vector
        //--------------------------------------------------------------------------------
        for(int k = 0; k <= final_index; k++)
        {
            //----------------------------------------------------------------------------
            //Copy segment by segment : time and state
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) y_traj_ecisem_c[i][k] += de*nullvector[6*k+i];
        }

        //--------------------------------------------------------------------------------
        //  New epsilon
        //--------------------------------------------------------------------------------
        cout << "Old epsilon = " << SEML.epsilon << endl;
        SEML.epsilon = SEML.epsilon+de*nullvector[nfv-1];
        cout << "New epsilon = " << SEML.epsilon << endl;

        //--------------------------------------------------------------------------------
        //  Correct the state with a differential correction procedure
        //--------------------------------------------------------------------------------
        multiple_shooting_direct_deps(y_traj_ecisem_c, t_traj_ecisem_c, y_traj_ecisem_c, t_traj_ecisem_c, nullvector, false, 48, final_index, coord_int, true, h4);

        //--------------------------------------------------------------------------------
        //The final time is shifted to be able to be plotted, using the right normalization
        //--------------------------------------------------------------------------------
        for(int k = 0; k <= final_index; k++) t_traj_ecisem_s[k] = t_traj_ecisem_c[k] + SEML.ss.tshift;

        //--------------------------------------------------------------------------------
        //Trajectory on lines, segment by segment, using SHIFTED time
        //--------------------------------------------------------------------------------
        if(SEML.epsilon > 0.90)
            plottrajsegbyseg(y_traj_ecisem_c, t_traj_ecisem_s, final_index, mPlot, NJ2000,
                             et0, tsys0, coord_type, h2, comp_type,  h3, color+1,
                             "");

        //--------------------------------------------------------------------------------
        // Plotting the direction of motion
        // in the (x, eps) space
        //--------------------------------------------------------------------------------
        gnuplot_plot_xy(h1, &y_traj_ecisem_c[0][0], &(SEML.epsilon), 1, "", "points", "7", "3", 1);


        //--------------------------------------------------------------------------------
        //Upload iter
        //--------------------------------------------------------------------------------
        iter++;

    }
    while(iter <= Nmax && SEML.epsilon < 1.0);


    //====================================================================================
    //Final trajectory on lines, segment by segment, using SHIFTED time
    //====================================================================================
    plottrajsegbyseg(y_traj_ecisem_c, t_traj_ecisem_s, final_index, mPlot, NJ2000,
                     et0, tsys0, coord_type, h2, comp_type,  h3, color++,
                     "Final guess (ecisem -> NJ2000)");

    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    pressEnter(refst.isFlagOn);
    gnuplot_close(h2);
    gnuplot_close(h3);
    gnuplot_close(h4);

    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj_n, 0, 41, 0, final_index);
    free_dvector(t_traj_n, 0, final_index);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(y_jpl_temp, 0, 41, 0, final_index);
    free_dvector(t_jpl_temp, 0, final_index);

    free_dmatrix(y_traj_ecisem_c, 0, 47, 0, final_index);
    free_dvector(t_traj_ecisem_c, 0, final_index);

    return FTC_SUCCESS;
}


//----------------------------------------------------------------------------------------
//         Refinement of solutions: to JPL - Subroutines
//----------------------------------------------------------------------------------------
/**
 * \brief Computes a first guess in ephemerides coordinates, switching between the
 *        Earth-Moon and Sun-Earth plane of motion when the discrepancy between the two
 *        coordinates system is minimal.
 **/
int oojplfg3d_switch(double** y_traj_n, double* t_traj_n,
                     double** y_traj_jpl, double* t_traj_jpl,
                     int final_index, int coord_type,
                     int comp_type, int coord_int,
                     double et0, double tsys0, double tsys0_comp)
{
    double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl_n  = dvector(0, final_index);

    //====================================================================================
    //Compute the change of coord: syn -> ecliptic, where we suppose that the current
    //plane of motion is the one of the primaries associated to coord_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), SEML.ss);
        break;
    }

    //====================================================================================
    //Half of the trajectory will be set in the other plane, i.e we suppose that the
    //current plane of motion is the one of the primaries associated to comp_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), SEML.ss);
        break;
    }


    //====================================================================================
    //Find the point that minimizes the discrepancy between the to sets of coordinates
    //The search is performed on the range [p0 p1]
    //====================================================================================
    int p0 = final_index/3, p1 = 2*final_index/3;
    int pmin = p0;
    int p = p0;
    double dmin = 0.0, dminmin = 0.0;
    do
    {
        dmin = 0.0;

        //The minimum is computed in 3 consecutive points: on p-1, p, p+1
        for(int i = 0; i < 6; i++) dmin += (y_traj_jpl[i][p] - y_traj_jpl_n[i][p])*(y_traj_jpl[i][p] - y_traj_jpl_n[i][p]);

        if(p == p0)
        {
            pmin = p;
            dminmin = dmin;
        }
        else
        {
            if(dmin < dminmin)
            {
                pmin = p;
                dminmin = dmin;
            }
        }
        p++;
    }
    while(p <= p1);
    cout << "oojplfg3d_switch. pmin = " << pmin << "/" << final_index << endl;

    //====================================================================================
    //Store the second half: at the beginning or at the end, depending on
    //the computation coordinates
    //====================================================================================
    int begind = 0, endind = 0;
    switch(eph_coord(coord_type))
    {
    case VEM:
        begind = pmin;
        endind = final_index;
        break;

    case VSEM:
        begind = 0;
        endind = pmin;
        break;
    }

    for(int p = begind; p <= endind; p++)
    {
        for(int i = 0; i <6; i++) y_traj_jpl[i][p] =  y_traj_jpl_n[i][p];
    }

    //====================================================================================
    //Free
    //====================================================================================
    free_dmatrix(y_traj_jpl_n, 0, 41, 0, final_index);
    free_dvector(t_traj_jpl_n, 0, final_index);


    return final_index;
}

/**
 * \brief Same as oojplfg3d_switch, with a refinement of the position of the minimum.
 **/
int oojplfg3d_super_switch(double** y_traj_n, double* t_traj_n,
                           double** y_traj_jpl, double* t_traj_jpl,
                           int final_index, int coord_type,
                           int comp_type, int coord_int,
                           int mRef,
                           double et0, double tsys0, double tsys0_comp)
{
    //====================================================================================
    //Initialize
    //====================================================================================
    double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl_n  = dvector(0, final_index);

    double** yms          = dmatrix(0, 5, 0, mRef);
    double* tms           = dvector(0, mRef);
    double** ymc          = dmatrix(0, 5, 0, mRef);
    double* tmc           = dvector(0, mRef);
    double** ymc_comp     = dmatrix(0, 5, 0, mRef);
    double* tmc_comp      = dvector(0, mRef);

    //Get the default coordinates system from the coord_type
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    //Compute the change of coord: syn -> ecliptic, where we suppose that the current
    //plane of motion is the one of the primaries associated to coord_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), SEML.ss);
        break;
    }

    //====================================================================================
    //Half of the trajectory will be set in the other plane, i.e we suppose that the
    //current plane of motion is the one of the primaries associated to comp_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), SEML.ss);
        break;
    }


    //====================================================================================
    //Find the point that minimizes the discrepancy between the to sets of coordinates
    //The search is performed on the range [p0 p1]
    //====================================================================================
    int p0 = final_index/3, p1 = 2*final_index/3;
    int pmin = p0;
    int p = p0;
    double dmin = 0.0, dminmin = 0.0;
    do
    {
        dmin = 0.0;
        for(int i = 0; i < 6; i++) dmin += (y_traj_jpl[i][p] - y_traj_jpl_n[i][p])*(y_traj_jpl[i][p] - y_traj_jpl_n[i][p]);

        if(p == p0)
        {
            pmin = p;
            dminmin = dmin;
        }
        else
        {
            if(dmin < dminmin)
            {
                pmin = p;
                dminmin = dmin;
            }
        }
        p++;
    }
    while(p <= p1);
    cout << "oojplfg3d_super_switch. pmin = " << pmin << "/" << final_index << endl;


    //====================================================================================
    //Once we have the point, we refine it by integration between pmin-1 and pmin+1
    //====================================================================================
    double tmin;
    double ymin[6], yv[6];
    int ode78coll;

    //------------------------------------------------------------------------------------
    // First subsearch range: pmin-1 to pmin
    //------------------------------------------------------------------------------------
    //Integration segment by segment, in coord_type coordinates
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][pmin-1];
    ode78(yms, tms, &ode78coll, t_traj_n[pmin-1], t_traj_n[pmin], yv, 6, mRef, dcs, coord_type, coord_type);

    //In ecliptic coordinates, in the plane associated to coord_type
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
        break;

    case J2000:
        coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
        break;

    case NJ2000:
        coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), SEML.ss);
        break;
    }

    //In ecliptic coordinates, in the plane associated to comp_type
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case J2000:
        coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case NJ2000:
        coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type), SEML.ss);
        break;
    }


    //Compute the minimum discrepancy
    for(p = 0; p <= mRef; p++)
    {
        dmin = 0.0;
        for(int i = 0; i < 3; i++) dmin += (ymc[i][p] - ymc_comp[i][p])*(ymc[i][p] - ymc_comp[i][p]);

        if(p == final_index/3)
        {
            dminmin = dmin;
            tmin = tmc[p];
            for(int i = 0; i <6; i++) ymin[i] = ymc[i][p];
        }
        else
        {
            if(dmin < dminmin)
            {
                dminmin = dmin;
                tmin = tmc[p];
                for(int i = 0; i <6; i++) ymin[i] = ymc[i][p];
            }
        }
    }

    //------------------------------------------------------------------------------------
    // Second subsearch range: pmin-1 to pmin
    //------------------------------------------------------------------------------------
    //Integration segment by segment, in coord_type coordinates
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][pmin];
    ode78(yms, tms, &ode78coll, t_traj_n[pmin], t_traj_n[pmin+1], yv, 6, mRef, dcs, coord_type, coord_type);

    //In ecliptic coordinates, in the plane associated to coord_type
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
        break;

    case J2000:
        coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
        break;

    case NJ2000:
        coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), SEML.ss);
        break;
    }

    //In ecliptic coordinates, in the plane associated to comp_type
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case J2000:
        coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case NJ2000:
        coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type), SEML.ss);
        break;
    }


    //Compute the minimum discrepancy
    tmin = tmc[0];
    for(p = 0; p <= mRef; p++)
    {
        dmin = 0.0;
        for(int i = 0; i < 3; i++) dmin += (ymc[i][p] - ymc_comp[i][p])*(ymc[i][p] - ymc_comp[i][p]);

        if(dmin < dminmin)
        {
            dminmin = dmin;
            tmin = tmc[p];
            for(int i = 0; i <6; i++) ymin[i] = ymc[i][p];
        }

    }

    cout << "-----------------------------------------------------------" << endl;
    cout << "oojplfg3d_super_switch. Switch part before refinment       " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = pmin-1; p <= pmin+1; p++)
    {
        cout << t_traj_jpl[p] << "  " << y_traj_jpl[0][p] << "   " << y_traj_jpl[1][p] <<  "   " << y_traj_jpl[2][p] << endl;
    }

    //The best fit replace the old point, both in position in time
    for(int i = 0; i <6; i++) y_traj_jpl_n[i][pmin] = ymin[i];
    t_traj_jpl[pmin] = tmin;

    cout << "-----------------------------------------------------------" << endl;
    cout << "oojplfg3d_super_switch. Switch part after refinment        " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = pmin-1; p <= pmin+1; p++)
    {
        cout << t_traj_jpl[p] << "  " << y_traj_jpl_n[0][p] << "   " << y_traj_jpl_n[1][p] <<  "   " << y_traj_jpl_n[2][p] << endl;
    }

    //====================================================================================
    //Store the second half: at the beginning or at the end, depending on
    //the computation coordinates
    //====================================================================================
    int begind = 0, endind = 0;
    switch(eph_coord(coord_type))
    {
    case VEM:
        begind = pmin;
        endind = final_index;
        break;

    case VSEM:
        begind = 0;
        endind = pmin;
        break;
    }

    for(int p = begind; p <= endind; p++)
    {
        for(int i = 0; i <6; i++) y_traj_jpl[i][p] =  y_traj_jpl_n[i][p];
    }


    //====================================================================================
    //Free
    //====================================================================================
    free_dmatrix(y_traj_jpl_n, 0, 41, 0, final_index);
    free_dvector(t_traj_jpl_n, 0, final_index);

    free_dmatrix(yms, 0, 5, 0, mRef);
    free_dvector(tms, 0, mRef);
    free_dmatrix(ymc, 0, 5, 0, mRef);
    free_dvector(tmc, 0, mRef);
    free_dmatrix(ymc_comp, 0, 5, 0, mRef);
    free_dvector(tmc_comp, 0, mRef);

    return final_index;
}


/**
 * \brief Same as oojplfg3d_switch, with an interpolation before and after the switching
 *        point.
 **/
int oojplfg3d_interpolation(double** y_traj_n, double* t_traj_n,
                            double** * y_traj_jpl, double** t_traj_jpl,
                            int final_index, int coord_type,
                            int comp_type, int coord_int,
                            int mRef,
                            double et0, double tsys0, double tsys0_comp)
{
    //====================================================================================
    //Initialize
    //====================================================================================
    double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl_n  = dvector(0, final_index);

    double** yms          = dmatrix(0, 5, 0, mRef);
    double* tms           = dvector(0, mRef);
    double** ymc          = dmatrix(0, 5, 0, mRef);
    double* tmc           = dvector(0, mRef);
    double** ymc_comp     = dmatrix(0, 5, 0, mRef);
    double* tmc_comp      = dvector(0, mRef);

    //Get the default coordinates system from the coord_type
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    //Compute the change of coord: syn -> ecliptic, where we suppose that the current
    //plane of motion is the one of the primaries associated to coord_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, *y_traj_jpl, *t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, *y_traj_jpl, *t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, *y_traj_jpl, *t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), SEML.ss);
        break;
    }

    //====================================================================================
    //Half of the trajectory will be set in the other plane, i.e we suppose that the
    //current plane of motion is the one of the primaries associated to comp_type
    //====================================================================================
    switch(coord_int)
    {
    case VECLI:
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case J2000:
        coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
        break;

    case NJ2000:
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), SEML.ss);
        break;
    }


    //====================================================================================
    //Find the point that minimizes the discrepancy between the to sets of coordinates
    //The search is performed on the range [p0 p1]
    //====================================================================================
    //gnuplot_ctrl *h4 = gnuplot_init();
    int p0 = 0, p1 = final_index;
    int pmin = p0;
    int p = p0;
    double dmin = 0.0, dminmin = 0.0;
    bool start = true, dte = false;
    double distanceToEarth = 0.0;
    double pe[3];
    do
    {
        //--------------------------------------------------------------------------------
        //Distance to Earth in NCSEM coordinates
        //--------------------------------------------------------------------------------
        distanceToEarth = 0.0;
        evaluateCoef(pe, t_traj_n[p], SEML.us_sem.n, SEML.nf, SEML.cs_sem.pe, 3);
        for(int i = 0; i < 3; i++) distanceToEarth += (y_traj_n[i][p] - pe[i])*(y_traj_n[i][p] - pe[i]);
        distanceToEarth = sqrt(distanceToEarth);
        dte = distanceToEarth > 0.4;

        //--------------------------------------------------------------------------------
        //Find best fit for a distance to Earth sufficiently large (dte == true)
        //--------------------------------------------------------------------------------
        dmin = 0.0;
        for(int i = 0; i < 6; i++) dmin += ((*y_traj_jpl)[i][p] - y_traj_jpl_n[i][p])*((*y_traj_jpl)[i][p] - y_traj_jpl_n[i][p]);
        if(start && dte)
        {
            pmin = p;
            dminmin = dmin;
            start = false;
        }
        else
        {
            //            if(dmin < dminmin && dte)
            //            {
            //                //cout << "Min! y_traj_n[0][p] = " << y_traj_n[0][p] << endl;
            //                pmin = p;
            //                dminmin = dmin;
            //            }
        }
        p++;

    }
    while(p <= p1);

    cout << "oojplfg3d_interpolation. pmin = " << pmin << "/" << final_index << endl;

    //    cout << "-----------------------------------------------------------" << endl;
    //    cout << "oojplfg3d_interpolation. Switch part before refinment      " << endl;
    //    cout << "-----------------------------------------------------------" << endl;
    //    for(int p = pmin-2; p <= pmin+1; p++)
    //    {
    //        cout << (*t_traj_jpl)[p] << "  " << (*y_traj_jpl)[0][p] << "   " << (*y_traj_jpl)[1][p] <<  "   " << (*y_traj_jpl)[2][p] << endl;
    //    }

    //pressEnter(true);

    //====================================================================================
    //Once we have the point, we interpolate it by integration between pmin-nsh and pmin+nsh
    //====================================================================================
    int nsh0 = 1;
    int nsh1 = 3;
    int nshs = nsh0+nsh1;

    double** y_traj_jpl_x   = dmatrix(0, 41, 0, final_index+nshs*mRef);
    double* t_traj_jpl_x    = dvector(0, final_index+nshs*mRef);
    double yv[6];
    int ode78coll;

    for(int nn = -nsh0; nn < nsh1; nn++)
    {
        //--------------------------------------------------------------------------------
        // First interpolation range: pmin-1 to pmin
        //--------------------------------------------------------------------------------
        //Integration segment by segment, in coord_type coordinates
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][pmin+nn];
        ode78(yms, tms, &ode78coll, t_traj_n[pmin+nn], t_traj_n[pmin+nn+1], yv, 6, mRef, dcs, coord_type, coord_type);

        //In ecliptic coordinates, in the plane associated to coord_type
        switch(coord_int)
        {
        case VECLI:
            coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
            break;

        case J2000:
            coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type));
            break;

        case NJ2000:
            coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), SEML.ss);
            break;
        }

        //In ecliptic coordinates, in the plane associated to comp_type
        switch(coord_int)
        {
        case VECLI:
            coord2eclstate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
            break;

        case J2000:
            coord2ecistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
            break;

        case NJ2000:
            coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type), SEML.ss);
            break;
        }

        //Interpolation from pmin-1 to pmin-1+mRef-1
        //cout << "- Storing from " << pmin-nsh0+(nsh0+nn)*mRef << "to " << pmin-nsh0+(nsh0+nn)*mRef+mRef-1 << endl;
        double epsilon = 0.0;
        for(p = 0; p <= mRef-1; p++)
        {
            epsilon = (1.0*p+(nsh0+nn)*mRef)/(nshs*(mRef)-1);
            //cout << "epsilon = " << epsilon << endl;
            for(int i = 0; i <6; i++) y_traj_jpl_x[i][p+pmin-nsh0+(nsh0+nn)*mRef] = epsilon*ymc[i][p] + (1-epsilon)*ymc_comp[i][p];
            t_traj_jpl_x[p+pmin-nsh0+(nsh0+nn)*mRef] = tmc[p];
        }
    }


    //====================================================================================
    //Store the parts before and after the interpolation part.
    //====================================================================================
    switch(eph_coord(coord_type))
    {
    case VEM:
        //--------------------------------------------------------------------------------
        // If we compute in VEM coordinates, the beginning of the trajectory is already
        // good, that last part should be set in the SE plane.
        //--------------------------------------------------------------------------------
        //Store the first part of the trajectory
        //cout << "- Storing from " << 0 << "to " << pmin-nsh0-1 << endl;
        for(p = 0; p <= pmin-nsh0-1; p++)
        {
            for(int i = 0; i <6; i++) y_traj_jpl_x[i][p] = (*y_traj_jpl)[i][p];
            t_traj_jpl_x[p] = (*t_traj_jpl)[p];
        }


        //Store the last part of the trajectory
        //cout << "- Storing from " << pmin-nsh0+nshs*mRef << "to " << final_index+nshs*mRef-nshs<< endl;
        for(p = pmin+nsh1; p <= final_index; p++)
        {
            for(int i = 0; i < 6; i++) y_traj_jpl_x[i][p+nshs*mRef-nshs] = y_traj_jpl_n[i][p];
            t_traj_jpl_x[p+nshs*mRef-nshs] = t_traj_jpl_n[p];
        }
        break;

    case VSEM:

        //--------------------------------------------------------------------------------
        // If we compute in VSEM coordinates, the end of the trajectory is already
        // good, that first part should be set in the EM plane.
        //--------------------------------------------------------------------------------
        //Store the first part of the trajectory
        //cout << "- Storing from " << 0 << "to " << pmin-nsh0-1 << endl;
        for(p = 0; p <= pmin-nsh0-1; p++)
        {
            for(int i = 0; i <6; i++) y_traj_jpl_x[i][p] = y_traj_jpl_n[i][p];
            t_traj_jpl_x[p] = t_traj_jpl_n[p];
        }


        //Store the last part of the trajectory
        //cout << "- Storing from " << pmin-nsh0+nshs*mRef << "to " << final_index+nshs*mRef-nshs<< endl;
        for(p = pmin+nsh1; p <= final_index; p++)
        {
            for(int i = 0; i < 6; i++) y_traj_jpl_x[i][p+nshs*mRef-nshs] = (*y_traj_jpl)[i][p];
            t_traj_jpl_x[p+nshs*mRef-nshs] = (*t_traj_jpl)[p];
        }
        break;
    }

    //====================================================================================
    //Free, before final_index is changed
    //====================================================================================
    free_dmatrix(y_traj_jpl_n, 0, 41, 0, final_index);
    free_dvector(t_traj_jpl_n, 0, final_index);

    free_dmatrix(yms, 0, 5, 0, mRef);
    free_dvector(tms, 0, mRef);
    free_dmatrix(ymc, 0, 5, 0, mRef);
    free_dvector(tmc, 0, mRef);
    free_dmatrix(ymc_comp, 0, 5, 0, mRef);
    free_dvector(tmc_comp, 0, mRef);

    //====================================================================================
    //New init. Once the previous computation has been done, we need to define again:
    //final_index, y_traj_jpl, and t_traj_jpl, so that the final trajectory incorporates
    //the interpolated part.
    //====================================================================================
    //Free from previous initialization
    free_dmatrix(*y_traj_jpl, 0, 41, 0, final_index);
    free_dvector(*t_traj_jpl, 0, final_index);

    //Redefinition of final_index
    int final_index_0 = final_index +nshs*(mRef-1);

    //Init again
    *y_traj_jpl  = dmatrix(0, 41, 0, final_index_0);
    *t_traj_jpl  = dvector(0, final_index_0);

    //Store
    for(p = 0; p <= final_index_0; p++)
    {
        for(int i = 0; i < 6; i++) (*y_traj_jpl)[i][p] = y_traj_jpl_x[i][p];
        (*t_traj_jpl)[p] = t_traj_jpl_x[p];
    }


    //    cout << "-----------------------------------------------------------" << endl;
    //    cout << "oojplfg3d_interpolation. Switch part after refinement      " << endl;
    //    cout << "-----------------------------------------------------------" << endl;
    //    for(int p = pmin-2; p <= pmin+2*mRef+1; p++)
    //    {
    //        cout << (*t_traj_jpl)[p] << "  " << (*y_traj_jpl)[0][p] << "   " << (*y_traj_jpl)[1][p] <<  "   " << (*y_traj_jpl)[2][p] << endl;
    //    }

    //====================================================================================
    //Return final_index_0
    //====================================================================================
    return final_index_0;
}

//----------------------------------------------------------------------------------------
//         Refinement of solutions: to JPL - Tests
//----------------------------------------------------------------------------------------
/**
 *  \brief Computes only a SEMLi orbit and test a JPL refinement.
 **/
int oocomprefft3d_test_seml_synjpl(int man_grid_size_t,
                                   int coord_type,
                                   Orbit& orbit_SEM,
                                   RefSt refst)

{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "oocomprefft3d_test_seml_synjpl";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_seml_SEM = refst.tspan_SEM;   //TOF on SEMLi orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");

    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type);

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    int max_grid    = 3000;
    int man_grid_2  = (refst.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;
    int traj_grid_2 = (refst.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;

    //------------------------------------------------------------------------------------
    //To store final data
    //------------------------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj  = dvector(0, traj_grid_2);

    double** y_traj_comp  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj_comp  = dvector(0, traj_grid_2);

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, traj_grid_2);
    double* t_traj_n    = dvector(0, traj_grid_2);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================

    //====================================================================================
    // 4.1 Initialize local variables
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Local variables to store the manifold leg
    //------------------------------------------------------------------------------------
    double** y_man_NCEM   = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_coord  = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_comp   = dmatrix(0, 5, 0, man_grid_2);
    double* t_man_EM      = dvector(0, man_grid_2);
    double* t_man_SEM     = dvector(0, man_grid_2);
    double* t_man_coord   = dvector(0, man_grid_2);
    double* t_man_comp    = dvector(0, man_grid_2);

    //====================================================================================
    // 6.4 Compute the final SEMLi orbit
    //====================================================================================
    cout << fname << ". Compute the final SEMLi orbit..."  << endl;
    //------------------------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //------------------------------------------------------------------------------------
    orbit_SEM.setSi(0.0, 4);
    orbit_SEM.update_ic(orbit_SEM.getSi());


    //------------------------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //------------------------------------------------------------------------------------
    int sem_index = man_grid_2;

    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    case REF_GIVEN_GRID:
    {
        int output = orbit_SEM.traj_int_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, true);
        //--------------------------------------------------------------------------------
        //If output is strictly greater than 0, then the projection procedure inside
        //the integration went wrong, and the new indix is output.
        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //is no point stored in the data, except the first one, and the new indix is zero
        //--------------------------------------------------------------------------------
        if(output == ORBIT_EPROJ) sem_index = 0;
        else if(output > 0) sem_index = output;
        break;
    }

    case REF_VAR_GRID:
    {
        sem_index = orbit_SEM.traj_int_var_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //To coord_type coordinates
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, sem_index, NCSEM, coord_type);
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_comp, t_man_comp, sem_index, NCSEM, comp_type);


    //------------------------------------------------------------------------------------
    // Save ALL points
    //------------------------------------------------------------------------------------
    for(int kman = 0; kman <= sem_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }

    //------------------------------------------------------------------------------------
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(h2, y_man_coord, sem_index+1, (char*)"SEMLi", "points", "1", "3", 6);
    gnuplot_plot_X(h3, y_man_comp, sem_index+1, (char*)"", "points", "1", "3", 6);

    //====================================================================================
    // Entire size is: no more than 3*man_grid_size_t points or so.
    //====================================================================================
    int final_index = 0;
    if(refst.grid == REF_VAR_GRID)
    {
        final_index = sem_index;
        int freq = final_index/(3*man_grid_size_t);
        cout << "----------------------------------"<< endl;
        cout << "Number of MSD points:        "     << endl;
        cout << "At SEML:          " << sem_index   << endl;
        cout << "Total:            " << final_index << endl;
        cout << "Selecting 1/" << freq << " points" << endl;
        cout << "Subtotal:         " << (int) floor(final_index/freq) << endl;
        cout << "----------------------------------"<< endl;

        //------------------------------------------------------------------------------------
        // Subset of points in y_traj_comp/t_traj_comp
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= floor(final_index/freq); kman++)
        {
            for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
            t_traj_comp[kman] = t_traj[freq*kman];
        }

        //New final index
        final_index = floor(final_index/freq);

        //------------------------------------------------------------------------------------
        // Subset of points back in y_traj/t_traj
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= final_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
            t_traj[kman] = t_traj_comp[kman];
        }
    }
    else
    {
        //If the grid size is fixed, we still return the variable indix, because
        //some integration procedure could have gone wrong so that the returned indices
        //are smaller than their optimal desired value
        final_index = sem_index;
    }

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    cout << fname << ". Initial trajectory, on a grid..."  << endl;
    int mPlot = 10;
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);
    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    double yv[6];
    int ode78coll;

    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, &ode78coll, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

        //Store
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }
    //To complementary coordinates
    qbcp_coc_vec(ymc_v, tmc_v, ymc_comp_v, tmc_comp_v, mPlot*final_index, coord_type, comp_type);

    //Plot on h2: lines
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);
    //Plot on h2: points
    gnuplot_plot_X(h2, y_traj, final_index+1, (char*)"", "points", "1", "2", 7);
    //Plot on h3: lines
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);


    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    pressEnter(refst.isFlagOn);

    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //------------------------------------------------------------------------------------
    // Final trajectory, on a grid
    //------------------------------------------------------------------------------------
    cout << fname << ". Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Store
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }

    //To complementary coordinates
    qbcp_coc_vec(ymc_v, tmc_v, ymc_comp_v, tmc_comp_v, mPlot*final_index, coord_type, comp_type);
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

    //====================================================================================
    // JPL
    //====================================================================================
    if(refst.isJPL)
    {
        //----------------------------------------------------------
        //Go on
        //----------------------------------------------------------
        pressEnter(refst.isFlagOn);

        //====================================================================================
        // Initialize SPICE kernerls
        //====================================================================================
        gnuplot_ctrl* h4 = gnuplot_init();
        cout << fname << ". Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //====================================================================================
        // Initialize VF
        //====================================================================================
        cout << fname << ". Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        int coll;
        OdeParams odeParams(&coll, &SEML);
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

        //====================================================================================
        // Search for best fit in  JPL ephemerides
        //====================================================================================
        cout << fname << ". Search for best fit in JPL DE430..." << endl;

        //----------------------------------------------------------
        //Get the best fit at t = t_traj_n[shift]. et0 is in seconds
        //----------------------------------------------------------
        double et0;
        double tsys0 = t_traj_n[shift];
        double tsys0_comp = tsys0;
        switch(eph_coord(comp_type))
        {
        case VEM:
            tsys0_comp /= SEML.us_em.ns;
            break;

        case VSEM:
            tsys0_comp *= SEML.us_em.ns;
            break;
        }
        //Best fit
        qbcp2jpl(tsys0, &et0, coord_type);

        //----------------------------------------------------------
        //Update the SEML
        //----------------------------------------------------------
        SEML.ss.et0 = et0;
        SEML.ss.t0  = tsys0;

        cout << "best fit (s) = " << et0 << endl;

        //====================================================================================
        // Change of coordinates
        //====================================================================================
        cout << fname << ". Change of coordinates syn -> ecliptic..." << endl;
        //------------------------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //------------------------------------------------------------------------------------
        double** y_traj_syn   = dmatrix(0, 41, 0, final_index);
        double* t_traj_syn    = dvector(0, final_index);
        double* et_traj_jpl   = dvector(0, final_index);

        double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double* t_traj_jpl_n  = dvector(0, final_index);

        //------------------------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //------------------------------------------------------------------------------------
        qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_syn, t_traj_syn, final_index, coord_type, eph_coord(coord_type));
        for(int p = 0; p <= final_index; p++)
        {
            et_traj_jpl[p] = et0 + (t_traj_n[p]-tsys0)/mean_motion(eph_coord(coord_type));
        }

        //====================================================================================
        // Plotting Initial Guess
        //====================================================================================
        pressEnter(refst.isFlagOn);
        cout << fname << ". Plotting Initial Guess..." << endl;
        //------------------------------------------------------------------------------------
        //Initial trajectory on lines, segment by segment
        //------------------------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
            ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);
        }

        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        pressEnter(refst.isFlagOn);
        cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
        //Plot
        gnuplot_plot_X(h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //--------------------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //--------------------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
        }

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        pressEnter(refst.isFlagOn);
        gnuplot_close(h4);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    pressEnter(refst.isFlagOn);
    gnuplot_close(h2);
    gnuplot_close(h3);


    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj, 0, traj_grid_2);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_comp, 0, traj_grid_2);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_n, 0, traj_grid_2);

    free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_2);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_2);
    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_2);
    free_dvector(t_man_EM, 0, man_grid_2);
    free_dvector(t_man_SEM, 0, man_grid_2);
    free_dvector(t_man_coord, 0, man_grid_2);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v  , 0, mPlot*final_index);

    return FTC_SUCCESS;
}


/**
 *  \brief Computes only an EML2 orbit and test a JPL refinement.
 **/
int oocomprefft3d_test_eml_synjpl(int man_grid_size_t,
                                  int coord_type,
                                  Orbit& orbit_EM,
                                  RefSt refst)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "oocomprefft3d_test_eml_synjpl";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_eml_EM = refst.tspan_EM;   //TOF on EML2 orbit
    //    double tof_seml_SEM = refst.tspan_SEM;   //TOF on SEMLi orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");

    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type);

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    int max_grid    = 3000;
    int man_grid_2  = (refst.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;
    int traj_grid_2 = (refst.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;

    //------------------------------------------------------------------------------------
    //To store final data
    //------------------------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj  = dvector(0, traj_grid_2);

    double** y_traj_comp  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj_comp  = dvector(0, traj_grid_2);

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, traj_grid_2);
    double* t_traj_n    = dvector(0, traj_grid_2);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================

    //====================================================================================
    // 4.1 Initialize local variables
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Local variables to store the manifold leg
    //------------------------------------------------------------------------------------
    double** y_man_NCEM   = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_coord  = dmatrix(0, 5, 0, man_grid_2);
    double** y_man_comp   = dmatrix(0, 5, 0, man_grid_2);
    double* t_man_EM      = dvector(0, man_grid_2);
    double* t_man_SEM     = dvector(0, man_grid_2);
    double* t_man_coord   = dvector(0, man_grid_2);
    double* t_man_comp    = dvector(0, man_grid_2);


    //====================================================================================
    // 4.2 Compute the initial orbit
    //====================================================================================
    cout << fname << ". Compute the initial orbit..."  << endl;
    //------------------------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //------------------------------------------------------------------------------------
    orbit_EM.setSi(0.0, 4);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //------------------------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //------------------------------------------------------------------------------------
    int em_index = man_grid_2;

    switch(refst.grid)
    {
    case REF_FIXED_GRID:
    case REF_GIVEN_GRID:
    {
        int output = orbit_EM.traj_int_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, true);
        //--------------------------------------------------------------------------------
        //If output is strictly greater than 0, then the projection procedure inside
        //the integration went wrong, and the new indix is output.
        //It output == ORBIT_EPROJ, the projection procedure failed so bad that there
        //is no point stored in the data, except the first one, and the new indix is zero
        //--------------------------------------------------------------------------------
        if(output == ORBIT_EPROJ) em_index = 0;
        else if(output > 0) em_index = output;
        break;
    }

    case REF_VAR_GRID:
    {
        em_index = orbit_EM.traj_int_var_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, true);
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //To coord_type coordinates
    //------------------------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, em_index, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp,  t_man_comp,  em_index, NCEM, comp_type);

    //------------------------------------------------------------------------------------
    // Save All points
    //------------------------------------------------------------------------------------
    for(int kman = 0; kman <= em_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][em_index-kman];
        t_traj[kman] = t_man_coord[em_index-kman];
    }

    //------------------------------------------------------------------------------------
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(h2, y_man_coord, em_index+1, (char*)"", "lines", "5", "3", 5);
    gnuplot_plot_X(h3, y_man_comp, em_index+1, (char*)"", "lines", "5", "3", 5);


    //====================================================================================
    // Entire size is: no more than 3*man_grid_size_t points or so.
    //====================================================================================
    int final_index = 0;
    if(refst.grid == REF_VAR_GRID)
    {
        final_index = em_index;
        int freq = final_index/(3*man_grid_size_t);
        cout << "----------------------------------"<< endl;
        cout << "Number of MSD points:        "     << endl;
        cout << "At EML:          " << em_index   << endl;
        cout << "Total:            " << final_index << endl;
        cout << "Selecting 1/" << freq << " points" << endl;
        cout << "Subtotal:         " << (int) floor(final_index/freq) << endl;
        cout << "----------------------------------"<< endl;

        //------------------------------------------------------------------------------------
        // Subset of points in y_traj_comp/t_traj_comp
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= floor(final_index/freq); kman++)
        {
            for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
            t_traj_comp[kman] = t_traj[freq*kman];
        }

        //New final index
        final_index = floor(final_index/freq);

        //------------------------------------------------------------------------------------
        // Subset of points back in y_traj/t_traj
        //------------------------------------------------------------------------------------
        for(int kman = 0; kman <= final_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
            t_traj[kman] = t_traj_comp[kman];
        }
    }
    else
    {
        //If the grid size is fixed, we still return the variable indix, because
        //some integration procedure could have gone wrong so that the returned indices
        //are smaller than their optimal desired value
        final_index = em_index;
    }

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    cout << fname << ". Initial trajectory, on a grid..."  << endl;
    int mPlot = 10;
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);
    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    double yv[6];
    int ode78coll;
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, &ode78coll, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

        //Store
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }
    //To complementary coordinates
    qbcp_coc_vec(ymc_v, tmc_v, ymc_comp_v, tmc_comp_v, mPlot*final_index, coord_type, comp_type);

    //Plot on h2: lines
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);
    //Plot on h2: points
    gnuplot_plot_X(h2, y_traj, final_index+1, (char*)"", "points", "1", "2", 7);
    //Plot on h3: lines
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);


    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //------------------------------------------------------------------------------------
    // Final trajectory, on a grid
    //------------------------------------------------------------------------------------
    cout << fname << ". Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Store
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }

    //To complementary coordinates
    qbcp_coc_vec(ymc_v, tmc_v, ymc_comp_v, tmc_comp_v, mPlot*final_index, coord_type, comp_type);
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

    //====================================================================================
    // JPL
    //====================================================================================
    if(refst.isJPL)
    {
        //----------------------------------------------------------
        //Go on
        //----------------------------------------------------------
        printf("Press ENTER to go on with the JPL ref");
        scanf("%c",&ch);


        //====================================================================================
        // Initialize SPICE kernerls
        //====================================================================================
        gnuplot_ctrl* h4 = gnuplot_init();
        cout << fname << ". Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //====================================================================================
        // Initialize VF
        //====================================================================================
        cout << fname << ". Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        int coll;
        OdeParams odeParams(&coll, &SEML);
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

        //====================================================================================
        // Search for best fit in  JPL ephemerides
        //====================================================================================
        cout << fname << ". Search for best fit in JPL DE430..." << endl;

        //----------------------------------------------------------
        //Get the best fit at t = t_traj_n[shift]. et0 is in seconds
        //----------------------------------------------------------
        double et0;
        double tsys0 = t_traj_n[shift];
        double tsys0_comp = tsys0;
        switch(eph_coord(comp_type))
        {
        case VEM:
            tsys0_comp /= SEML.us_em.ns;
            break;

        case VSEM:
            tsys0_comp *= SEML.us_em.ns;
            break;
        }
        //Best fit
        qbcp2jpl(tsys0, &et0, coord_type);

        //----------------------------------------------------------
        //Update the SEML
        //----------------------------------------------------------
        SEML.ss.et0 = et0;
        SEML.ss.t0  = tsys0;

        //====================================================================================
        // Change of coordinates
        //====================================================================================
        cout << fname << ". Change of coordinates syn -> ecliptic..." << endl;
        //------------------------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //------------------------------------------------------------------------------------
        double** y_traj_syn   = dmatrix(0, 41, 0, final_index);
        double* t_traj_syn    = dvector(0, final_index);
        double* et_traj_jpl   = dvector(0, final_index);

        double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double* t_traj_jpl_n  = dvector(0, final_index);

        //------------------------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //------------------------------------------------------------------------------------
        qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_syn, t_traj_syn, final_index, coord_type, eph_coord(coord_type));
        for(int p = 0; p <= final_index; p++)
        {
            et_traj_jpl[p] = et0 + (t_traj_n[p]-tsys0)/mean_motion(eph_coord(coord_type));
        }

        //================================================================================
        // Plotting Initial Guess
        //================================================================================
        printf("Press ENTER to plot the Initial Guess\n");
        scanf("%c",&ch);
        cout << fname << ". Plotting Initial Guess..." << endl;
        //--------------------------------------------------------------------------------
        //Initial trajectory on lines, segment by segment
        //--------------------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
            ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);

        }

        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
        //Plot
        gnuplot_plot_X(h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //--------------------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //--------------------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
        }

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
        gnuplot_close(h4);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);
    gnuplot_close(h3);


    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj, 0, traj_grid_2);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_comp, 0, traj_grid_2);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_n, 0, traj_grid_2);

    free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_2);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_2);
    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_2);
    free_dvector(t_man_EM, 0, man_grid_2);
    free_dvector(t_man_SEM, 0, man_grid_2);
    free_dvector(t_man_coord, 0, man_grid_2);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v  , 0, mPlot*final_index);

    return FTC_SUCCESS;
}


/**
 *  \brief Computes only a EML2-SEMLi connection and test a JPL refinement, in synodical coordinates
 **/
int oocomprefft3d_test_eml2seml_synjpl(int coord_type)
{
    //====================================================================================
    // 1. Read the data, in sys units
    //====================================================================================
    string fname = "oocomprefft3d_test_eml2seml_synjpl";
    cout << fname << ". Read the data, in sys units..."  << endl;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    int final_index = getLengthCOMP_txt();
    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    readCOMP_txt(t_traj_n, y_traj_n, final_index);

    //====================================================================================
    // 2. Local parameters
    //====================================================================================
    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //----------------------------------------------------------
    // State and time vectors
    //----------------------------------------------------------
    int mPlot = 100;

    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double* tmc_comp    = dvector(0, mPlot);

    double** ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_v       = dvector(0, mPlot*final_index);
    double** ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double* tmc_comp_v  = dvector(0, mPlot*final_index);

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type);

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    double yv[6];
    int ode78coll;

    cout << fname << ". Initial trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, &ode78coll, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Store
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }

    //To complementary coordinates
    qbcp_coc_vec(ymc_v, tmc_v, ymc_comp_v, tmc_comp_v, mPlot*final_index, coord_type, comp_type);
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

    //====================================================================================
    // JPL
    //====================================================================================
    //----------------------------------------------------------
    //Go on
    //----------------------------------------------------------
    char ch;
    printf("Press ENTER to go on with the JPL ref");
    scanf("%c",&ch);


    //====================================================================================
    // Initialize SPICE kernerls
    //====================================================================================
    gnuplot_ctrl* h4 = gnuplot_init();
    cout << fname << ". Initialize SPICE kernerls..." << endl;
    furnsh_c("spice/kernels/metakernel.furnsh");

    //====================================================================================
    // Initialize VF
    //====================================================================================
    cout << fname << ". Initialize associated VF..." << endl;
    int shift = 0;
    OdeStruct driver_JPL;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    int coll;
    OdeParams odeParams(&coll, &SEML);
    //Init ode structure
    init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

    //====================================================================================
    // Search for best fit in  JPL ephemerides
    //====================================================================================
    cout << fname << ". Search for best fit in JPL DE430..." << endl;

    //----------------------------------------------------------
    //Get the best fit at t = t_traj_n[shift]. et0 is in seconds
    //----------------------------------------------------------
    double et0;
    double tsys0      = t_traj_n[shift];
    double tsys0_comp = tsys0;
    switch(eph_coord(comp_type))
    {
    case VEM:
        tsys0_comp /= SEML.us_em.ns;
        break;

    case VSEM:
        tsys0_comp *= SEML.us_em.ns;
        break;
    }
    qbcp2jpl(tsys0, &et0, coord_type);

    //----------------------------------------------------------
    //Update the SEML
    //----------------------------------------------------------
    SEML.ss.et0 = et0;
    SEML.ss.t0  = tsys0;

    //====================================================================================
    // Change of coordinates
    //====================================================================================
    cout << fname << ". Change of coordinates syn -> ecliptic..." << endl;
    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    double** y_traj_syn   = dmatrix(0, 41, 0, final_index);
    double* t_traj_syn    = dvector(0, final_index);
    //double* et_traj_jpl   = dvector(0, final_index);
    double** y_traj_jpl   = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl    = dvector(0, final_index);
    double** y_traj_syn_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_syn_n  = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Compute the change of coord: syn -> ecliptic
    //------------------------------------------------------------------------------------
    //qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_syn, t_traj_syn, final_index, coord_type, eph_coord(coord_type));
    coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));

    cout << "-----------------------------------------------------------" << endl;
    cout << "Trajectory  before refinment                               " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = 0; p <= final_index; p++)
    {
        cout << y_traj_jpl[0][p] << "   " << y_traj_jpl[1][p] <<  "   " << y_traj_jpl[2][p] << endl;
    }

    //====================================================================================
    //Half of the trajectory is set in the other plane
    //====================================================================================
    double** y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_jpl_n  = dvector(0, final_index);
    coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0, tsys0_comp, eph_coord(comp_type));

    //----------------------------------------------------------
    //Find the point that minimizes the discrepancy
    //----------------------------------------------------------
    int pmin = final_index/3;
    int p = final_index/3;
    double dmin = 0.0, dminmin = 0.0;
    do
    {
        dmin = 0.0;
        for(int i = 0; i < 6; i++) dmin += (y_traj_jpl[i][p] - y_traj_jpl_n[i][p])*(y_traj_jpl[i][p] - y_traj_jpl_n[i][p]);

        if(p == final_index/3)
        {
            pmin = p;
            dminmin = dmin;
        }
        else
        {
            if(dmin < dminmin)
            {
                pmin = p;
                dminmin = dmin;
            }
        }
        p++;
    }
    while(p <= 2*final_index/3);

    cout << "pmin = " << pmin << "/" << final_index << endl;

    //----------------------------------------------------------
    //Store the second half: at the beginning or at the end, depending on the computation coordinates
    //----------------------------------------------------------
    int begind = 0, endind = 0;
    switch(eph_coord(coord_type))
    {
    case VEM:
        begind = pmin;
        endind = final_index;
        break;

    case VSEM:
        begind = 0;
        endind = pmin;
        break;
    }

    for(int p = begind; p <= endind; p++)
    {
        for(int i = 0; i <6; i++) y_traj_jpl[i][p] =  y_traj_jpl_n[i][p];
    }


    cout << "-----------------------------------------------------------" << endl;
    cout << "Trajectory  before refinment                               " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = 0; p <= final_index; p++)
    {
        cout << y_traj_jpl[0][p] << "   " << y_traj_jpl[1][p] <<  "   " << y_traj_jpl[2][p] << endl;
    }

    //====================================================================================
    //Whole trajectory back in syn coordinates
    //====================================================================================
    ecl2synstate_vec(y_traj_jpl, t_traj_jpl, y_traj_syn, t_traj_syn, final_index, et0, tsys0, eph_coord(coord_type));

    cout << "-----------------------------------------------------------" << endl;
    cout << "Trajectory (y_traj_syn)  before refinment                               " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = 0; p <= final_index; p++)
    {
        cout << y_traj_syn[0][p] << "   " << y_traj_syn[1][p] <<  "   " << y_traj_syn[2][p] << endl;
    }

    final_index = pmin;
    //------------------------------------------------------------------------------------
    //Time in seconds
    //------------------------------------------------------------------------------------
    //for(int p = 0; p <= final_index; p++)
    //{
    //        et_traj_jpl[p] = et0 + (t_traj_n[p]-tsys0)/mean_motion(eph_coord(coord_type));
    //}

    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    printf("Press ENTER to plot the Initial Guess\n");
    scanf("%c",&ch);
    cout << fname << ". Plotting Initial Guess..." << endl;
    //------------------------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
        ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

        //To ecl, then to NCSEM
        qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
        else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

        if(k == 0 || k == final_index-1)
        {
            for(int p = 0; p <= mPlot; p++)
            {
                cout << ymc[0][p] << "   " << ymc[1][p] <<  "   " << ymc[2][p] << endl;
            }
            cout << "--------------------------------------" << endl;
        }

        //Back to comp_type coordinates
        // @todo: put it NCEM  -> ECLIPTIC -> NCSEM, using ecl2coordstate_vec and its counterpart
        syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
        ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

        //Plot on h3
        if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
        else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);
    }

    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    printf("Press ENTER to refine\n");
    scanf("%c",&ch);
    cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

    int isPlotted   = 1;
    multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_syn_n, t_traj_syn_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
    //multiple_shooting_direct_variable_time(y_traj_syn, t_traj_syn, y_traj_syn_n, t_traj_syn_n,  42, final_index, eph_coord(coord_type), 5e-9, isPlotted, h4);
    //Plot
    gnuplot_plot_X(h4, y_traj_syn_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


    //------------------------------------------------------------------------------------
    //Final trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_syn_n[i][k];
        ode78_jpl(ymc_comp, tmc_comp, &ode78coll, t_traj_syn_n[k], t_traj_syn_n[k+1], yv, 6, mPlot,
                  eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));


        //Back to coord_type coordinates
        qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
        else gnuplot_plot_X(h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

        //Back to comp_type coordinates
        //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
        syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
        ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

        //Plot on h3
        if(k == 0) gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
        else gnuplot_plot_X(h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h4);


    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);
    gnuplot_close(h3);


    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(y_traj_n, 0, 41, 0, final_index);
    free_dvector(t_traj_n, 0, final_index);

    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    free_dmatrix(ymc_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_v, 0, mPlot*final_index);
    free_dmatrix(ymc_comp_v, 0, 5, 0, mPlot*final_index);
    free_dvector(tmc_comp_v  , 0, mPlot*final_index);

    return FTC_SUCCESS;
}


//========================================================================================
//
// Text format, read
//
//========================================================================================
/**
 * \brief Reads a given \c Ofsc  object within a \c Oftsc, in txt format.
 **/
void toCelestiaFormat(string filename)
{
    //Init
    ifstream readStream;
    string ct;
    //Reading
    readStream.open((filename).c_str(), ios::in);

    //Check that the opening went well
    if (!readStream.is_open())
    {
        cout << "toCelestiaFormat. Cannot open file " << filename << endl;
        cout << "Check the text data exist." << endl;
        return;
    }
    else
    {
        //--------------------------------------------------------------------------------
        // Count the number of lines
        //--------------------------------------------------------------------------------
        int count0 = -2;
        while(!readStream.eof())
        {
            getline(readStream, ct);
            count0 ++;
        }

        //--------------------------------------------------------------------------------
        // Go back to the beginning (not the best way, but enough for now)
        //--------------------------------------------------------------------------------
        //readStream.seekg (0, readStream.beg); --> should be better, but not working
        readStream.close();
        readStream.open((filename).c_str(), ios::in);

        //--------------------------------------------------------------------------------
        // Read the data
        //--------------------------------------------------------------------------------
        double** y1  = dmatrix(0, 5, 0, count0);
        double* t1   = dvector(0, count0);
        for(int i = 0; i<= count0; i++)
        {
            //Time
            readStream >> t1[i];
            //Position
            readStream >> y1[0][i];
            readStream >> y1[1][i];
            readStream >> y1[2][i];
        }
        readStream.close();

        //--------------------------------------------------------------------------------
        // COC towards VSEM + store data
        //--------------------------------------------------------------------------------
        double** y2  = dmatrix(0, 5, 0, count0);
        double* t2   = dvector(0, count0);

        furnsh_c("spice/kernels/metakernel.furnsh");
        eci2syndpos_vec(y1, t1, y2, t2, count0, VSEM);

        // Rewrite the data
        gnuplot_fplot_txyzv(t2, y2, count0+1, "jpltraj_SEM.xyz", "w");


        //--------------------------------------------------------------------------------
        // COC towards VEM + store data
        //--------------------------------------------------------------------------------
        eci2syndpos_vec(y1, t1, y2, t2, count0, VEM);
        // Rewrite the data
        gnuplot_fplot_txyzv(t2, y2, count0+1, "jpltraj_EM.xyz", "w");


        //--------------------------------------------------------------------------------
        // Position of the SEML1 point: may be wrong!
        //--------------------------------------------------------------------------------
        SpiceDouble lt;
        SpiceDouble RS[6], RB[6], RE[6];
        for(int i = 0; i<= count0; i++)
        {

            //Position
            spkez_c (m1name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RS, &lt);
            spkez_c (m2name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RB, &lt);
            spkez_c (EARTH, t1[i], DEFFRAME, "NONE", SSB, RE, &lt);

            //Position of SEML1
            y1[0][i] = (RS[0] - RB[0])*(SEML.cs_sem_l1.mu - 1 + SEML.cs_sem_l1.gamma) - RE[0];
            y1[1][i] = (RS[1] - RB[1])*(SEML.cs_sem_l1.mu - 1 + SEML.cs_sem_l1.gamma) - RE[1];
            y1[2][i] = (RS[2] - RB[2])*(SEML.cs_sem_l1.mu - 1 + SEML.cs_sem_l1.gamma) - RE[2];
        }

        //Towards VSEM coordinates
        eci2syndpos_vec(y1, t1, y2, t2, count0, VSEM);
        gnuplot_fplot_txyzv(t2, y2, count0+1, "SEML1_SEM.xyz", "w");


        //--------------------------------------------------------------------------------
        // Position of the SEMLi point: may be wrong!
        //--------------------------------------------------------------------------------
        for(int i = 0; i<= count0; i++)
        {

            //Position
            spkez_c (m1name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RS, &lt);
            spkez_c (m2name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RB, &lt);
            spkez_c (EARTH, t1[i], DEFFRAME, "NONE", SSB, RE, &lt);

            //Position of SEMLi
            y1[0][i] = (RS[0] - RB[0])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[0];
            y1[1][i] = (RS[1] - RB[1])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[1];
            y1[2][i] = (RS[2] - RB[2])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[2];
        }

        //Towards VSEM coordinates
        eci2syndpos_vec(y1, t1, y2, t2, count0, VSEM);
        gnuplot_fplot_txyzv(t2, y2, count0+1, "SEML2_SEM.xyz", "w");


        //--------------------------------------------------------------------------------
        // Position of the EML2 point
        //--------------------------------------------------------------------------------
        double k = 0.0;
        for(int i = 0; i<= count0; i++)
        {

            //Position
            spkez_c (m1name(VEM), t1[i], DEFFRAME, "NONE", SSB, RS, &lt);
            spkez_c (m2name(VEM), t1[i], DEFFRAME, "NONE", SSB, RB, &lt);
            spkez_c (EARTH, t1[i], DEFFRAME, "NONE", SSB, RE, &lt);

            k = sqrt((RS[0] - RB[0])*(RS[0] - RB[0]) + (RS[1] - RB[1])*(RS[1] - RB[1])+ (RS[2] - RB[2])*(RS[2] - RB[2]));
            //Position of EML2
            y2[0][i] = -(SEML.us_em.mu_EM - 1.0 - SEML.cs_em_l2.gamma)*k;
            y2[1][i] = 0.0;
            y2[2][i] = 0.0;
        }

        //Towards VEM coordinates
        gnuplot_fplot_txyzv(t2, y2, count0+1, "EML2_EM.xyz", "w");


        //--------------------------------------------------------------------------------
        // Position of the EML1 point
        //--------------------------------------------------------------------------------
        for(int i = 0; i<= count0; i++)
        {

            //Position
            spkez_c (m1name(VEM), t1[i], DEFFRAME, "NONE", SSB, RS, &lt);
            spkez_c (m2name(VEM), t1[i], DEFFRAME, "NONE", SSB, RB, &lt);
            spkez_c (EARTH, t1[i], DEFFRAME, "NONE", SSB, RE, &lt);

            k = sqrt((RS[0] - RB[0])*(RS[0] - RB[0]) + (RS[1] - RB[1])*(RS[1] - RB[1])+ (RS[2] - RB[2])*(RS[2] - RB[2]));
            //Position of EML1
            y2[0][i] = -(SEML.us_em.mu_EM - 1.0 + SEML.cs_em_l1.gamma)*k;
            y2[1][i] = 0.0;
            y2[2][i] = 0.0;
        }

        //Towards VEM coordinates
        gnuplot_fplot_txyzv(t2, y2, count0+1, "EML1_EM.xyz", "w");

    }

    //--------------------------------------------------------------------------------
    // Do the same for the Barycenters? The barycenter of the Earth-Moon system is easily retrieved.
    //--------------------------------------------------------------------------------
}
