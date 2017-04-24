#include "oolencon.h"


//========================================================================================
//
//          Computation of the CMU about SEMLi
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about SEMLi,
 *         in the QBCP model. The initial conditions (IC) are computed in a 5-dimensional
 *         box: one dimension for the starting time, four dimensions for the
 *         parameterization of the Center Manifold (s1 to s4 coordinates). The RCM
 *         coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:      the value in RCM coordinates on the unstable direction s5.
 *  \param projSt.TLIM:     the min/max starting time (in SEM units) in the IC box.
 *  \param projSt.TSIZE:    the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI:  the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                          in the IC box.
 *  \param projSt.GSIZE_SI: the number of points on the  s1, s2, s3, s4 values  grids
 *                          in the IC box.
 *  \param projSt.ISPAR;    if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/SEM/L2/cu_3d_order_16.bin"
 **/
int compute_grid_CMU_SEM_3D(double dist_to_cm, ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Splash screen
    //====================================================================================
    // Update the filename of type TYPE_CU_3D
    string filename_output = projSt.get_and_update_filename(projSt.FILE_CU, TYPE_CU_3D, ios::out);

    cout << resetiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;
    cout << "       Computation of center-unstable IC at SEMLi                  " << endl;
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
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us->T;
    cout << ", " << projSt.TLIM[1]/SEML.us->T << "] x T" << endl;
    cout << " The data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Building the working grids
    //------------------------------------------------------------------------------------
    double** grid_si_CMU_RCM = (double**) calloc(4, sizeof(double*));
    for(int i = 0; i <4; i++)
    {
        grid_si_CMU_RCM[i] = (double*) calloc(projSt.GSIZE_SI[i]+1, sizeof(double));
        init_grid(grid_si_CMU_RCM[i], projSt.GLIM_SI[i][0], projSt.GLIM_SI[i][1], projSt.GSIZE_SI[i]);
    }

    cout << " - The detailed values of s2 are:                                   " << endl;
    for(int i = 0; i < projSt.GSIZE_SI[1]; i++) cout << grid_si_CMU_RCM[1][i] << ", ";
    cout << grid_si_CMU_RCM[1][projSt.GSIZE_SI[1]]  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);


    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_SEM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_SEM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_SEM[i]/SEML.us->T << ", ";
    cout << grid_t_SEM[projSt.TSIZE]/SEML.us->T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.GSIZE_SI[1])*(1+projSt.GSIZE_SI[2])*(1+projSt.GSIZE_SI[3])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------------------------------------------------
    // Data structures
    //------------------------------------------------------------------------------------
    double** init_state_CMU_NCSEM = dmatrix(0, 5, 0, projSt.GSIZE_SI[2]);
    double** init_state_CMU_RCM   = dmatrix(0, 4, 0, projSt.GSIZE_SI[2]);

    //====================================================================================
    // Get the invariant manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //====================================================================================
    // Check that the invman is an unstable-manifold
    //====================================================================================
    if(invman.getManType() != MAN_CENTER_U)
    {
        cout << "compute_grid_CMU_SEM_3D. The invariant manifold must be of center-unstable type. return." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Reset the data file
    //------------------------------------------------------------------------------------
    initCU_bin_3D(projSt.GSIZE_SI, projSt.TSIZE, filename_output);

    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= projSt.TSIZE; kt++)
    {
        //--------------------------------------------------------------------------------
        //Append the time in data file
        //--------------------------------------------------------------------------------
        appTimeCU_bin_3D(grid_t_SEM, kt, filename_output);


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


                        //----------------------------------------------------------------
                        // Initialization on the center-unstable manifold
                        //----------------------------------------------------------------
                        //Init sti
                        sti[0] = grid_si_CMU_RCM[0][ks1];
                        sti[1] = grid_si_CMU_RCM[1][ks2];
                        sti[2] = grid_si_CMU_RCM[2][ks3];
                        sti[3] = grid_si_CMU_RCM[3][ks4];
                        sti[4] = dist_to_cm;


                        //----------------------------------------------------------------
                        //Equivalent state
                        //----------------------------------------------------------------
                        invman.evalRCMtoNC(sti, grid_t_SEM[kt], yvu, OFTS_ORDER, OFS_ORDER);

                        //----------------------------------------------------------------
                        //Save
                        //----------------------------------------------------------------
                        #pragma omp critical
                        {

                            for(int i = 0; i < 6; i++) init_state_CMU_NCSEM[i][ks3] = yvu[i];
                            for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][ks3] = sti[i];
                            //Display
                            displayCompletion("compute_grid_CMU_SEM", (double) iter++/noe*100);
                        }

                        free_dvector(yvu, 0, 5);
                        free_dvector(sti, 0, 4);

                    }

                    //--------------------------------------------------------------------
                    //Store values
                    //--------------------------------------------------------------------
                    writeCU_bin_3D(init_state_CMU_NCSEM, init_state_CMU_RCM,
                                   projSt.GSIZE_SI, filename_output);
                }
            }
        }
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_dmatrix(init_state_CMU_NCSEM, 0, 5, 0, projSt.GSIZE_SI[2]);
    free_dmatrix(init_state_CMU_RCM,   0, 4, 0, projSt.GSIZE_SI[2]);

    return FTC_SUCCESS;
}

//========================================================================================
//
//         Projection on the CM/CMS/CMU of EMLi
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM. Then, each point on the
 *         integration grid is projected on the Center Manifold CM_EM_NC about EMLi.
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param projSt.TM:          the maximum integration time on each leg, in SEM units.
 *
 *  \param projSt.MSIZE:       the number of points on each manifold leg.
 *
 *  \param projSt.NOD:         the number of dimensions on which the distance of
 *                             projection is computed (usually either 3 (the physical
 *                             distance) or 6 (the whole phase space)).
 *
 *  \param projSt.ISPAR:       if TRUE, the computation is parallelized.
 *
 *  \param projSt.YNMAX:       the maximum norm in NCEM coordinates for which a given
 *                             state on the integration grid is projected on CM_EM_NC
 *                             More precisely: for a given state y along the manifold leg,
 *                             if norm(y, 3) < projSt.YNMAX, the state is projected.
 *                             Otherwise, it is considered too far away from EMLi to be
 *                             a good candidate for projection.
 *
 *  \param projSt.SNMAX:       the maximum norm in RCM EM coordinates for which a given
 *                             projection state on the CM of EMLi (CM_SEM_NC) is
 *                             computed back in NCSEM coordinates. More precisely, for a
 *                             given state y in NCSEM coordinates, the result of the
 *                             projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                             coordinates.
 *                             if norm(sproj, 4) < projSt.SNMAX, the computation
 *                             yproj = CM_EM_NC(sproj, t) is performed. Otherwise, the
 *                             state sproj is considered too far away from the RCM origin
 *                             to be a good candidate - it is out of the domain of
 *                             practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form:
 *          "plot/QBCP/SEM/L2/projcu_3d_order_16.bin".
 **/
int int_proj_CMU_SEM_on_CM_EM_3D(ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Filenames
    //====================================================================================
    // Get the filename from the latest data of type TYPE_CU_3D
    string filename_input = projSt.get_filename(TYPE_CU_3D, ios::in);
    // Get the filename for the output
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ_3D, ios::out);


    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //Read data size
    int t_grid_size, si_grid_size[4], offset;
    offset = getLenghtCU_bin_3D(si_grid_size, &t_grid_size, filename_input);

    if(offset < 0)
    {
        cout << "int_proj_CMU_SEM_on_CM_EM_3D. Impossible to get length of the data file." << endl;
        return FTC_FAILURE;
    }


    //====================================================================================
    // Splash screen
    //====================================================================================
    cout << resetiosflags(ios::scientific) << setprecision(15);

    //------------------------------------------------------------------------------------
    //To store all data
    //------------------------------------------------------------------------------------
    double** init_state_CMU_NCSEM = dmatrix(0, 5, 0, si_grid_size[2]);
    double** init_state_CMU_RCM   = dmatrix(0, 4, 0, si_grid_size[2]);
    double* init_time_grid_SEM    = dvector(0, t_grid_size);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+si_grid_size[0])*(1+si_grid_size[1])*(1+si_grid_size[2])*(1+si_grid_size[3])*(1+t_grid_size);


    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @EML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_EM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_SEM_on_CM_EM_3D. The invariant manifold at EMLj ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_EM(OFTS_ORDER, OFS_ORDER, *SEML_EM.cs);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
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
        offset = readTCU_bin_3D(offset, init_time_grid_SEM, kt, filename_input);

        if(offset < 0)
        {
            cout << "int_proj_CMU_SEM_on_CM_EM_3D. Impossible to read the data file." << endl;
            return FTC_FAILURE;
        }

        for(int ks2 = 0; ks2 <= si_grid_size[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= si_grid_size[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= si_grid_size[0]; ks1++)
                {
                    //--------------------------------------------------------------------
                    //Read data from file
                    //--------------------------------------------------------------------
                    offset = readCU_bin_3D(offset, init_state_CMU_NCSEM,
                                           init_state_CMU_RCM, si_grid_size, filename_input);
                    if(offset < 0)
                    {
                        cout << "int_proj_CMU_SEM_on_CM_EM_3D. Impossible to read the data file." << endl;
                        return FTC_FAILURE;
                    }

                    //--------------------------------------------------------------------
                    //Most inner loop is parallelized
                    //--------------------------------------------------------------------
                    #pragma omp parallel for if(isPar)  shared(index)
                    for(int ks3 = 0; ks3 <= si_grid_size[2]; ks3++)
                    {

                        //----------------------------------------------------------------
                        //Inputs
                        //----------------------------------------------------------------
                        ProjResSt projResSt(NCSEM);

                        //----------------------------------------------------------------
                        //Initial conditions
                        //----------------------------------------------------------------
                        // Init the time in SEM coordinates
                        projResSt.init_time  = init_time_grid_SEM[kt];
                        // Init the state in NCSEM coordinates
                        for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = init_state_CMU_NCSEM[i][ks3];
                        // Init the state in RCM coordinates
                        for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = init_state_CMU_RCM[i][ks3];
                        //Label
                        projResSt.label = 0;

                        //----------------------------------------------------------------
                        //Projection on center manifold at EMLj
                        //----------------------------------------------------------------
                        proj_subroutine(projResSt, invman_EM, projSt);

                        //----------------------------------------------------------------
                        // Save outputs
                        //----------------------------------------------------------------
                        if(projResSt.min_proj_dist_SEM_o < ePdef)
                        {
                            //Save
                            #pragma omp critical
                            {
                                writeIntProjCU_bin(filename_output, projResSt);
                            }
                        }

                        //----------------------------------------------------------------
                        //Display completion
                        //----------------------------------------------------------------
                        #pragma omp critical
                        {
                            displayCompletion("int_proj_CMU_SEM_on_CM_EM_3D", 100.0*index++/noe);
                        }
                    }
                }
            }
        }
    }


    //------------------------------------------------------------------------------------
    //Free all data
    //------------------------------------------------------------------------------------
    free_dmatrix(init_state_CMU_NCSEM, 0, 5, 0, si_grid_size[2]);
    free_dmatrix(init_state_CMU_RCM, 0, 4, 0, si_grid_size[2]);
    free_dvector(init_time_grid_SEM, 0, t_grid_size);


    return FTC_SUCCESS;
}

//========================================================================================
//
//          Computation of the CMU about EMLi
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
int compute_grid_CMU_EM_3D(double dist_to_cm, ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Splash screen
    //====================================================================================
    // Update the filename of type TYPE_CU_3D
    string filename_output = projSt.get_and_update_filename(projSt.FILE_CU, TYPE_CU_3D, ios::out);

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
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us->T;
    cout << ", " << projSt.TLIM[1]/SEML.us->T << "] x T" << endl;
    cout << " The data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Building the working grids
    //------------------------------------------------------------------------------------
    double** grid_si_CMU_RCM = (double**) calloc(4, sizeof(double*));
    for(int i = 0; i <4; i++)
    {
        grid_si_CMU_RCM[i] = (double*) calloc(projSt.GSIZE_SI[i]+1, sizeof(double));
        init_grid(grid_si_CMU_RCM[i], projSt.GLIM_SI[i][0], projSt.GLIM_SI[i][1], projSt.GSIZE_SI[i]);
    }

    cout << " - The detailed values of s2 are:                                   " << endl;
    for(int i = 0; i < projSt.GSIZE_SI[1]; i++) cout << grid_si_CMU_RCM[1][i] << ", ";
    cout << grid_si_CMU_RCM[1][projSt.GSIZE_SI[1]]  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);


    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_EM[i]/SEML.us->T << ", ";
    cout << grid_t_EM[projSt.TSIZE]/SEML.us->T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.GSIZE_SI[1])*(1+projSt.GSIZE_SI[2])*(1+projSt.GSIZE_SI[3])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------------------------------------------------
    // Data structures
    //------------------------------------------------------------------------------------
    double** init_state_CMU_NCEM = dmatrix(0, 5, 0, projSt.GSIZE_SI[2]);
    double** init_state_CMU_RCM  = dmatrix(0, 4, 0, projSt.GSIZE_SI[2]);

    //====================================================================================
    // Get the invariant manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //====================================================================================
    // Check that the invman is an unstable-manifold
    //====================================================================================
    if(invman.getManType() != MAN_CENTER_U)
    {
        cout << "compute_grid_CMU_EM_3D. The invariant manifold must be of center-unstable type. return." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    // Reset the data file
    //------------------------------------------------------------------------------------
    initCU_bin_3D(projSt.GSIZE_SI, projSt.TSIZE, filename_output);

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
        appTimeCU_bin_3D(grid_t_EM, kt, filename_output);


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


                        //----------------------------------------------------------------
                        // Initialization on the center-unstable manifold
                        //----------------------------------------------------------------
                        //Init sti
                        sti[0] = grid_si_CMU_RCM[0][ks1];
                        sti[1] = grid_si_CMU_RCM[1][ks2];
                        sti[2] = grid_si_CMU_RCM[2][ks3];
                        sti[3] = grid_si_CMU_RCM[3][ks4];
                        sti[4] = dist_to_cm;


                        //----------------------------------------------------------------
                        //Equivalent state
                        //----------------------------------------------------------------
                        invman.evalRCMtoNC(sti, grid_t_EM[kt], yvu, OFTS_ORDER, OFS_ORDER);

                        //----------------------------------------------------------------
                        //Save
                        //----------------------------------------------------------------
                        #pragma omp critical
                        {

                            for(int i = 0; i < 6; i++) init_state_CMU_NCEM[i][ks3] = yvu[i];
                            for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][ks3] = sti[i];
                            //Display
                            displayCompletion("compute_grid_CMU_EM", (double) iter++/noe*100);
                        }

                        free_dvector(yvu, 0, 5);
                        free_dvector(sti, 0, 4);

                    }

                    //--------------------------------------------------------------------
                    //Store values
                    //--------------------------------------------------------------------
                    writeCU_bin_3D(init_state_CMU_NCEM, init_state_CMU_RCM,
                                   projSt.GSIZE_SI, filename_output);
                }
            }
        }
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
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
int compute_grid_CMU_EM(double dist_to_cm, ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Splash screen
    //====================================================================================
    // Update the filename of type TYPE_CU
    string filename_output = projSt.get_and_update_filename(projSt.FILE_CU, TYPE_CU, ios::out);

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
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us_em.T;
    cout << ", " << projSt.TLIM[1]/SEML.us_em.T << "] x T" << endl;
    cout << " The data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    if(projSt.dHd > 0)
    {
        cout << "Moreover, a fixed energy has been selected:                        " << endl;
        cout << "Hd - H0 = " << projSt.dHd                                            << endl;
        cout << "===================================================================" << endl;
    }

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Building the working grids
    //------------------------------------------------------------------------------------
    double* grid_s1_CMU_RCM = dvector(0,  projSt.GSIZE_SI[0]);
    double* grid_s3_CMU_RCM = dvector(0,  projSt.GSIZE_SI[2]);
    init_grid(grid_s1_CMU_RCM, projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1], projSt.GSIZE_SI[0]);
    init_grid(grid_s3_CMU_RCM, projSt.GLIM_SI[2][0], projSt.GLIM_SI[2][1], projSt.GSIZE_SI[2]);

    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_EM[i]/SEML.us_em.T << ", ";
    cout << grid_t_EM[projSt.TSIZE]/SEML.us_em.T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.GSIZE_SI[2])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------------------------------------------------
    // Data structures
    //------------------------------------------------------------------------------------
    double**** init_state_CMU_NCEM = d4tensor(0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    double**** init_state_CMU_RCM  = d4tensor(0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    double** *  init_dH_valid       = d3tensor(0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);

    //====================================================================================
    // Get the invariant center-unstable manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs_em);

    //====================================================================================
    // Get the invariant center manifold at EML2
    //====================================================================================
    CSYS cs_em_center;
    init_CSYS(&cs_em_center, &SEML, &SEM, F_EM, SEML.li_EM, SEML.numberOfCoefs, false, PMS_GRAPH, MAN_CENTER);
    Invman invman_center(OFTS_ORDER, OFS_ORDER, cs_em_center);

    //====================================================================================
    // Energy at the origin, at t = 0.0
    //====================================================================================
    double* st0 = dvector(0,4);
    for(int i = 0; i < 5; i++) st0[i] = 0.0;


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

                //------------------------------------------------------------------------
                // Initialization on the center-unstable manifold
                //------------------------------------------------------------------------
                //Init sti
                sti[0] = grid_s1_CMU_RCM[ks1];
                sti[1] = 0.0;
                sti[2] = grid_s3_CMU_RCM[ks3];
                sti[3] = 0.0;
                sti[4] = dist_to_cm;

                //------------------------------------------------------------------------
                // Energy: we impose H = Hd
                //------------------------------------------------------------------------
                if(projSt.dHd > 0)
                {
                    //Energy at the origin (at t = 0.0 for now)
                    double H0 = invman_center.H_SYS(st0, 0.0);

                    //Desired energy
                    double Hd = H0 + projSt.dHd;

                    //Refinement to get H = Hd
                    RefineH refineH(&invman_center, 2, sti, grid_t_EM[kt], Hd);
                    int status = init_s0_energy(&refineH, sti, grid_t_EM[kt]);

                    if(status == GSL_SUCCESS)
                    {
                        //IC are valid
                        init_dH_valid[kt][ks1][ks3] = 1;

                        //Display
                        //cout << "Convergence on dHd:" << invman_center.H_SYS(sti, grid_t_EM[kt]) - H0 << endl;
                        //cout << "sti = " << endl;
                        //vector_printf_prec(sti, 5);

                    }
                    else
                    {
                        //IC are  non valid
                        init_dH_valid[kt][ks1][ks3] = 0;

                        //Display
                        //cout << "No convergence on dHd." << endl;
                        //cout << "sti = " << endl;
                        //vector_printf_prec(sti, 5);
                    }
                }
                else
                {
                    //IC are always valid
                    init_dH_valid[kt][ks1][ks3] = 1;
                }

                //------------------------------------------------------------------------
                //Equivalent state
                //------------------------------------------------------------------------
                invman.evalRCMtoNC(sti, grid_t_EM[kt], yvu, OFTS_ORDER, OFS_ORDER);

                //------------------------------------------------------------------------
                //Save
                //------------------------------------------------------------------------
                #pragma omp critical
                {

                    for(int i = 0; i < 6; i++) init_state_CMU_NCEM[i][kt][ks1][ks3] = yvu[i];
                    for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][kt][ks1][ks3] = sti[i];
                    //Display
                    displayCompletion("compute_grid_CMU_EM", (double) iter++/noe*100);
                }

                free_dvector(yvu, 0, 5);
                free_dvector(sti, 0, 4);
            }
        }
    }

    //------------------------------------------------------------------------------------
    //Store values
    //------------------------------------------------------------------------------------
    writeCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, init_dH_valid, grid_t_EM,
                projSt.GSIZE_SI[0], projSt.GSIZE_SI[2], projSt.TSIZE, filename_output);

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    free_d4tensor(init_state_CMU_RCM,  0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    free_d3tensor(init_dH_valid,  0, projSt.TSIZE, 0, projSt.GSIZE_SI[0], 0, projSt.GSIZE_SI[2]);
    free_dvector(st0, 0, 4);

    return FTC_SUCCESS;
}

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         two-dimensional box: one dimension for the starting time (t0), one dimension
 *         for the parameterization of the Center Manifold (s1 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *          Two possibilities exist:
 *          (i) if projSt.dHd > 0, then a fixed value of energy
 *          at departure is desired by the user. Thereore, the value s3 = f(s1, t0)
 *          is computed and the IC are of the form (s1, 0, f(s1, t0), 0, dist_to_cm)
 *
 *          (ii) if projSt.dHd <= 0, the IC are of the form (s1, 0, 0 0, dist_to_cm)
 *
 *
 *  \param dist_to_cm:           the value in RCM coordinates on the unst. direction s5.
 *  \param projSt.TLIM[0]:       the minimum starting time (in EM units) in the IC box.
 *  \param projSt.TLIM[1]:       the maximum starting time (in EM units) in the IC box.
 *  \param projSt.TSIZE:         the number of points on the time grid in the IC box.
 *  \param projSt.GLIM_SI[0][0]: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GLIM_SI[0][1]: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param projSt.GSIZE_SI[0]:   the number of points on the s1 grid in the IC box.
 *  \param projSt.ISPAR:         if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16_dH_0.005.bin"
 **/
int compute_grid_CMU_EM_dH(double dist_to_cm, ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Splash screen
    //====================================================================================
    // Update the filename of type TYPE_CU
    string filename_output = projSt.get_and_update_filename(projSt.FILE_CU, TYPE_CU, ios::out);

    cout << resetiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;
    cout << "       Computation of center-unstable planar IC at EML2            " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << projSt.GSIZE_SI[0]+1  << " value(s) of s1 in [" << projSt.GLIM_SI[0][0];
    cout << ", " << projSt.GLIM_SI[0][1] << "]" << endl;
    cout << "  - " << projSt.TSIZE+1  << " value(s) of t in [" << projSt.TLIM[0]/SEML.us_em.T;
    cout << ", " << projSt.TLIM[1]/SEML.us_em.T << "] x T" << endl;
    cout << " The data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    if(projSt.dHd > 0)
    {
        cout << "Moreover, a fixed energy has been selected:                        " << endl;
        cout << "Hd - H0 = " << projSt.dHd                                            << endl;
        cout << "===================================================================" << endl;
    }

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Building the working grids
    //------------------------------------------------------------------------------------
    double* grid_s1_CMU_RCM = dvector(0,  projSt.GSIZE_SI[0]);
    init_grid(grid_s1_CMU_RCM, projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1], projSt.GSIZE_SI[0]);

    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_EM[i]/SEML.us_em.T << ", ";
    cout << grid_t_EM[projSt.TSIZE]/SEML.us_em.T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+projSt.GSIZE_SI[0])*(1+projSt.TSIZE);
    int iter = 1;

    //------------------------------------------------------------------------------------
    // Data structures
    //------------------------------------------------------------------------------------
    double** * init_state_CMU_NCEM = d3tensor(0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);
    double** * init_state_CMU_RCM  = d3tensor(0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);
    double**   init_dH_valid      = dmatrix(0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);

    //====================================================================================
    // Get the invariant center-unstable manifold at EML2
    //====================================================================================
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs_em);

    //====================================================================================
    // Get the invariant center manifold at EML2
    //====================================================================================
    CSYS cs_em_center;
    init_CSYS(&cs_em_center, &SEML, &SEM, F_EM, SEML.li_EM, SEML.numberOfCoefs, false, PMS_GRAPH, MAN_CENTER);
    Invman invman_center(OFTS_ORDER, OFS_ORDER, cs_em_center);

    //====================================================================================
    // Energy at the origin, at t = 0.0
    //====================================================================================
    double* st0 = dvector(0,4);
    for(int i = 0; i < 5; i++) st0[i] = 0.0;


    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= projSt.TSIZE; kt++)
    {
        #pragma omp parallel for if(isPar)  shared(iter)
        for(int ks1 = 0; ks1 <= projSt.GSIZE_SI[0]; ks1++)
        {
            //----------------------------------------------------------------------------
            // Inner variables
            //----------------------------------------------------------------------------
            Ofsc ofs(OFS_ORDER);
            double* yvu = dvector(0,5);
            double* sti = dvector(0,4);

            //----------------------------------------------------------------------------
            // Initialization on the center-unstable manifold
            //----------------------------------------------------------------------------
            //Init sti
            sti[0] = grid_s1_CMU_RCM[ks1];
            sti[1] = 0.0;
            sti[2] = 0.0;
            sti[3] = 0.0;
            sti[4] = dist_to_cm;

            //----------------------------------------------------------------------------
            // Energy: we impose H = Hd
            //----------------------------------------------------------------------------
            if(projSt.dHd > 0)
            {
                //Energy at the origin (at t = 0.0 for now)
                double H0 = invman_center.H_SYS(st0, grid_t_EM[kt]);

                //Desired energy
                double Hd = H0 + projSt.dHd;

                //Refinement to get H = Hd
                RefineH refineH(&invman_center, 3, sti, grid_t_EM[kt], Hd);
                int status = init_s0_energy(&refineH, sti, grid_t_EM[kt]);

                if(status == GSL_SUCCESS)
                {
                    //IC are valid
                    init_dH_valid[kt][ks1] = 1;

                    //Display
                    //cout << "Convergence on dHd:" << invman_center.H_SYS(sti, grid_t_EM[kt]) - H0 << endl;
                    //cout << "sti = " << endl;
                    //vector_printf_prec(sti, 5);

                }
                else
                {
                    //IC are  non valid
                    init_dH_valid[kt][ks1] = 0;

                    //Display
                    //cout << "No convergence on dHd." << endl;
                    //cout << "sti = " << endl;
                    //vector_printf_prec(sti, 5);
                }
            }
            else
            {
                //IC are always valid
                init_dH_valid[kt][ks1] = 1;
            }

            //----------------------------------------------------------------------------
            //Equivalent state
            //----------------------------------------------------------------------------
            invman.evalRCMtoNC(sti, grid_t_EM[kt], yvu, OFTS_ORDER, OFS_ORDER);

            //----------------------------------------------------------------------------
            //Save
            //----------------------------------------------------------------------------
            #pragma omp critical
            {

                for(int i = 0; i < 6; i++) init_state_CMU_NCEM[i][kt][ks1] = yvu[i];
                for(int i = 0; i < 5; i++) init_state_CMU_RCM[i][kt][ks1] = sti[i];
                //Display
                displayCompletion("compute_grid_CMU_EM", (double) iter++/noe*100);
            }

            free_dvector(yvu, 0, 5);
            free_dvector(sti, 0, 4);
        }
    }

    //------------------------------------------------------------------------------------
    //Store values
    //------------------------------------------------------------------------------------
    writeCU_bin_dH(init_state_CMU_NCEM, init_state_CMU_RCM, init_dH_valid, grid_t_EM,
                   projSt.GSIZE_SI[0], projSt.TSIZE, filename_output);

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_d3tensor(init_state_CMU_NCEM, 0, 5, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);
    free_d3tensor(init_state_CMU_RCM,  0, 4, 0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);
    free_dmatrix(init_dH_valid,  0, projSt.TSIZE, 0, projSt.GSIZE_SI[0]);
    free_dvector(st0, 0, 4);

    return FTC_SUCCESS;
}

//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEMLi
//
//========================================================================================
/**
 *  \brief Projection subroutine, used in all the subsequent routines in this section.
 *         All outputs in projResSt are updated
 *         (see the declaration of this structure for details).
 *
 *         This routine performs the following steps:
 *
 *          - It integrates the initial state at emli, stored in projResSt
 *          - It compute the minimum and the argminimum distance of projection along the
 *            corresponding trajectory, targeting the center manifold invman_target
 *          - The information relative to this min/argmin are stored in projResSt.
 *
 *         The coordinates within which the results are saved are encoded in ProjResSt.
 **/
int proj_subroutine(ProjResSt& projResSt, Invman& invman_target, ProjSt& projSt)
{
    //====================================================================================
    // 1. Local variables
    //====================================================================================
    //------------------------------------------------------------------------------------
    // State and time vectors, at initialization, in IC_COORD
    //------------------------------------------------------------------------------------
    double** y_man_PR = dmatrix(0, 5, 0, projSt.MSIZE);
    double* t_man_PR  = dvector(0, projSt.MSIZE);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //------------------------------------------------------------------------------------
    // Inner variables
    //------------------------------------------------------------------------------------
    // Variables in IC_COORD coordinates
    double yv_IC[6], tv_IC, tf_IC;

    // Variables in PR_COORD (projection) coordinates
    double yvproj_PR[6], sproj[4], y_man_norm_PR = 0.0;
    double yv_PR[6], tv_PR;

    // Variables in FC_COORD & FV_COORD coordinates
    double yv_FC[6], yvproj_FC[6], yv_FV[6], yvproj_FV[6];
    double dv_at_projection_FV = 0.0;

    // Variables for the distance of projection, in DP_COORD
    double yv_DP[6], yvproj_DP[6];
    double proj_dist_DP, min_proj_dist_DP = ePdef;

    // Variables with their own, fixed, coordinate system
    double crossings_NCSEM = 0.0, crossings_lunar = 0.0;
    int collision_NCEM = 0;

    //Temp variables
    int kmin = 0;
    Ofsc ofs(OFS_ORDER);

    //====================================================================================
    // 2. Integration of the manifold leg
    //====================================================================================
    //Initial conditions
    tv_IC  = projResSt.init_time;
    for(int i = 0; i < 6; i++) yv_IC[i] = projResSt.init_state_CMU_NC[i];

    // Integration on projSt.MSIZE+1 fixed grid, from IC_COORD to PR_COORD coordinates
    // Note that the integration framework is always I_NCEM, hence collision_NCEM always
    // makes sense
    OdeEvent odeEvent(false, false);
    int status = ode78(y_man_PR, t_man_PR, &odeEvent, tv_IC, tv_IC+projSt.TM, yv_IC, 6, projSt.MSIZE, I_NCEM, projResSt.IC_COORD, projResSt.PR_COORD);

    //====================================================================================
    // 3. Projection on the center manifold of EMLi.
    // No use of SEML_SEM or SEML after this point!
    // We need first to check that the integration went well
    //====================================================================================
    //------------------------------------------------------------------------------------
    // We need first to check that the integration went well
    //------------------------------------------------------------------------------------
    if(status)
    {
        // If status != 0, something went wrong during the integration
        // in ode78 and we use the maximum value ePdef (the solution
        // is basically discarded)
        proj_dist_DP = ePdef;
    }
    else
    {

        //--------------------------------------------------------------------------------
        //Loop on trajectory
        //--------------------------------------------------------------------------------
        for(int kman = 0; kman <= projSt.MSIZE; kman++)
        {
            //Current state
            for(int i = 0; i < 6; i++) yv_PR[i] = y_man_PR[i][kman];
            tv_PR = t_man_PR[kman];

            //Current distance from EMLi in PR_COORD coordinates
            y_man_norm_PR = 0.0;
            for(int i = 0; i < 2; i++) y_man_norm_PR += yv_PR[i]*yv_PR[i];
            y_man_norm_PR = sqrt(y_man_norm_PR);

            //cout << "y_man_norm_PR = " << y_man_norm_PR << endl;

            //----------------------------------------------------------------------------
            //Check n°1: the current state is close enough to EMLi
            //----------------------------------------------------------------------------
            if(y_man_norm_PR < projSt.YNMAX)
            {
                // Projection on the center manifold
                invman_target.NCprojCCMtoCM(yv_PR, tv_PR, sproj);

                //cout << "s_man_norm_RCM = " << ENorm(sproj, 4) << endl;
                //------------------------------------------------------------------------
                //Check n°2: the projection is close enough to EMLi
                //------------------------------------------------------------------------
                if(ENorm(sproj, 4)< projSt.SNMAX)
                {
                    //yvproj_PR = W(sproj, tv)
                    invman_target.evalRCMtoNC(sproj, tv_PR, yvproj_PR, OFTS_ORDER, OFS_ORDER);

                    //Distance of projection in DP_COORD coordinates
                    qbcp_coc(tv_PR, yv_PR,     yv_DP,     projResSt.PR_COORD, projResSt.DP_COORD);
                    qbcp_coc(tv_PR, yvproj_PR, yvproj_DP, projResSt.PR_COORD, projResSt.DP_COORD);

                    proj_dist_DP = 0.0;
                    for(int i = 0; i < projSt.NOD; i++) proj_dist_DP += (yvproj_DP[i] - yv_DP[i])*(yvproj_DP[i] - yv_DP[i]);
                    proj_dist_DP = sqrt(proj_dist_DP);

                    //cout << "proj_dist_DP = " << proj_dist_DP << endl;
                }
                else proj_dist_DP = ePdef;
            }
            else proj_dist_DP = ePdef;

            //----------------------------------------------------------------------------
            //Update distance min if necessary
            //----------------------------------------------------------------------------
            if(proj_dist_DP < min_proj_dist_DP)
            {
                //Distance of projection
                min_proj_dist_DP = proj_dist_DP;

                //Indix (argmin)
                kmin = kman;

                //In FC coordinates
                qbcp_coc(tv_PR, yv_PR,     yv_FC,     projResSt.PR_COORD, projResSt.FC_COORD);
                qbcp_coc(tv_PR, yvproj_PR, yvproj_FC, projResSt.PR_COORD, projResSt.FC_COORD);

                //In FC coordinates, with velocity
                qbcp_coc(tv_PR, yv_PR,     yv_FV,     projResSt.PR_COORD, projResSt.FV_COORD);
                qbcp_coc(tv_PR, yvproj_PR, yvproj_FV, projResSt.PR_COORD, projResSt.FV_COORD);

                // Saving states
                for(int i = 0; i < 6; i++) projResSt.final_state_CMU_FC_o[i]     = yv_FC[i];
                for(int i = 0; i < 6; i++) projResSt.projected_state_CMU_FC_o[i] = yvproj_FC[i];
                for(int i = 0; i < 4; i++) projResSt.projected_state_CMU_RCM_o[i]  = sproj[i];

                //Associated DV
                dv_at_projection_FV = 0.0;
                for(int i = 3; i < 6; i++) dv_at_projection_FV += (yv_FV[i] - yvproj_FV[i])*(yv_FV[i] - yvproj_FV[i]);
                dv_at_projection_FV = sqrt(dv_at_projection_FV);

            }else if(proj_dist_DP > min_proj_dist_DP && projSt.PRIMARY)
            {
                // If we only look for the primary family, we can stop at the first
                // minimum without any risk of losing good solutions. If
                // proj_dist_DP > min_proj_dist_DP, we know that we have passed the
                // first minimum, and we can break the loop.
                break;
            }

        }

        //--------------------------------------------------------------------------------
        // We check if we have the primary family, if necessary
        // To do so, we check that we have made two clockwise/counterclockwise
        //  crossings of x = -1.
        //
        // To do a clockwise turn:
        //  1. x1 > 0 -> x2 < 0 && y1 > 0
        //  2. x1 < 0 -> x2 > 0 && y1 < 0
        //
        // To do a counterclockwise turn:
        //  1. x1 > 0 -> x2 < 0 && y1 < 0
        //  2. x1 < 0 -> x2 > 0 && y1 > 0
        //--------------------------------------------------------------------------------
        if(min_proj_dist_DP < ePdef)
        {
            // Init temporary objects
            OdeEvent odeEvent_purge(true, false);
            double** y_purge = dmatrix(0, 5, 0, 2);
            double* t_purge  = dvector(0, 2);

            // Initial condition in IC_COORD coordinates
            tv_IC  = projResSt.init_time;
            for(int i = 0; i < 6; i++) yv_IC[i] = projResSt.init_state_CMU_NC[i];

            // Compute the event using ode78 on [tv_IC, tf_IC]
            // Note that the integration framework is always I_NCSEM, so that
            // crossings_NCSEM always makes sense
            qbcp_coc_time(&t_man_PR[kmin], &tf_IC, 0, projResSt.PR_COORD, projResSt.IC_COORD);
            ode78(y_purge, t_purge, &odeEvent_purge, tv_IC, tf_IC, yv_IC, 6, 2, I_NCSEM, projResSt.IC_COORD, NCSEM);

            // Save value in crossings_NCSEM
            crossings_NCSEM = odeEvent_purge.crossings;

            // Save value in crossings_lunar
            crossings_lunar = odeEvent_purge.crossings_moon;

            // Free temporary objects
            free_dmatrix(y_purge, 0, 5, 0, 2);
            free_dvector(t_purge, 0, 2);
        }

        //--------------------------------------------------------------------------------
        //We check for collisions
        //--------------------------------------------------------------------------------
        if(odeEvent.coll)
        {
            collision_NCEM = odeEvent.coll;
        }

        //--------------------------------------------------------------------------------
        //If we save ONLY the primary family
        //--------------------------------------------------------------------------------
        if(projSt.PRIMARY)
        {
            if(collision_NCEM != 0 || crossings_NCSEM != 2.0) min_proj_dist_DP = ePdef;
        }

        //--------------------------------------------------------------------------------
        //If we save ONLY the primary family
        //--------------------------------------------------------------------------------
        if(crossings_lunar > 0)
        {
            //cout << "There has been " << crossings_lunar << "crossings of the lunar orbits" << endl;
        }
    }


    //====================================================================================
    // 3.3. Save outputs
    //====================================================================================
    //Initial position in FC_COORD coordinates
    for(int i = 0; i < 6; i++) yv_PR[i] = y_man_PR[i][0];
    qbcp_coc(t_man_PR[0], yv_PR, yv_FC, projResSt.PR_COORD, projResSt.FC_COORD);
    for(int i = 0; i < 6; i++) projResSt.init_state_CMU_FC_o[i] = yv_FC[i];

    //Minimum projection distance
    projResSt.min_proj_dist_SEM_o    = min_proj_dist_DP;

    //DV at projection
    projResSt.dv_at_projection_FC_o  = dv_at_projection_FV;

    //Crossings
    projResSt.crossings_NCSEM_o      = crossings_NCSEM;

    //Collisions
    projResSt.collision_NCEM_o       = collision_NCEM;

    //Final time
    projResSt.final_time_FC_o        = t_man_PR[kmin];


    //====================================================================================
    //Free
    //====================================================================================
    free_dmatrix(y_man_PR, 0, 5, 0, projSt.MSIZE);
    free_dvector(t_man_PR, 0, projSt.MSIZE);


    return GSL_SUCCESS;
}

/**
 *  \brief Projection subroutine, used in all the subsequent routines in this section.
 *         All outputs in projResSt are updated
 *         (see the declaration of this structure for details).
 *
 *         This routine performs the following steps:
 *
 *          - It integrates the initial state at emli, stored in projResSt
 *          - It compute the minimum and the argminimum distance of projection along the
 *            corresponding trajectory, targeting the center manifold invman_SEM
 *          - The information relative to this min/argmin are stored in projResSt.
 *
 *
 *         Note: this routine is the old version from the EML to SEML case.
 **/
int proj_subroutine_old_from_EML(ProjResSt& projResSt, Invman& invman_SEM, ProjSt& projSt)
{
    //====================================================================================
    // 1. Local variables
    //====================================================================================
    //------------------------------------------------------------------------------------
    // State and time vectors
    //------------------------------------------------------------------------------------
    double** y_man_NCSEM = dmatrix(0, 5, 0, projSt.MSIZE);
    double* t_man_SEM    = dvector(0, projSt.MSIZE);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //------------------------------------------------------------------------------------
    // Misc
    //------------------------------------------------------------------------------------
    double yvproj_NCSEM[6], sproj[4], yv_SEM[6], yvproj_SEM[6], yv_VSEM[6], yvproj_VSEM[6];
    double proj_dist_SEM, min_proj_dist_SEM = ePdef;
    double dv_at_projection_SEM = 0.0, y_man_norm_NCSEM = 0.0;
    double crossings_NCSEM = 0.0, crossings_lunar = 0.0;
    int collision_NCEM = 0;
    int kmin = 0;
    Ofsc ofs(OFS_ORDER);

    //====================================================================================
    // 2. Integration of the manifold leg
    //====================================================================================
    //Initial conditions
    double yv[6], tv;
    tv  = projResSt.init_time;
    for(int i = 0; i < 6; i++) yv[i] = projResSt.init_state_CMU_NC[i];

    //Integration on projSt.MSIZE+1 fixed grid
    OdeEvent odeEvent(false, false);
    int status = ode78(y_man_NCSEM, t_man_SEM, &odeEvent, tv, tv+projSt.TM, yv, 6, projSt.MSIZE, I_NCEM, NCEM, NCSEM);

    //====================================================================================
    // 3. Projection on the center manifold of SEMLi.
    // No use of SEML_EM or SEML after this point!
    // We need first to check that the integration went well
    //====================================================================================
    //------------------------------------------------------------------------------------
    // We need first to check that the integration went well
    //------------------------------------------------------------------------------------
    if(status)
    {
        // If status != 0, something went wrong during the integration
        // in ode78 and we use the maximum value ePdef (the solution
        // is basically discarded)
        proj_dist_SEM = ePdef;
    }
    else
    {

        //--------------------------------------------------------------------------------
        //Loop on trajectory
        //--------------------------------------------------------------------------------
        for(int kman = 0; kman <= projSt.MSIZE; kman++)
        {
            //Current state
            for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
            tv = t_man_SEM[kman];

            //Current distance from SEMLi in NCSEM units
            y_man_norm_NCSEM = 0.0;
            for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
            y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

            //----------------------------------------------------------------------------
            //Check n°1: the current state is close enough to SEMLi
            //----------------------------------------------------------------------------
            if(y_man_norm_NCSEM < projSt.YNMAX)
            {
                // Projection on the center manifold
                invman_SEM.NCprojCCMtoCM(yv, tv, sproj);

                //------------------------------------------------------------------------
                //Check n°2: the projection is close enough to SEMLi
                //------------------------------------------------------------------------
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

            //----------------------------------------------------------------------------
            //Update distance min if necessary
            //----------------------------------------------------------------------------
            if(proj_dist_SEM < min_proj_dist_SEM)
            {
                //Distance of projection
                min_proj_dist_SEM = proj_dist_SEM;

                //Indix (argmin)
                kmin = kman;

                // Misc states
                for(int i = 0; i < 6; i++) projResSt.final_state_CMU_FC_o[i]     = yv_SEM[i];
                for(int i = 0; i < 6; i++) projResSt.projected_state_CMU_FC_o[i] = yvproj_SEM[i];
                for(int i = 0; i < 4; i++) projResSt.projected_state_CMU_RCM_o[i]  = sproj[i];

                //Associated DV
                dv_at_projection_SEM = 0.0;
                for(int i = 3; i < 6; i++) dv_at_projection_SEM += (yvproj_VSEM[i] - yv_VSEM[i])*(yvproj_VSEM[i] - yv_VSEM[i]);
                dv_at_projection_SEM = sqrt(dv_at_projection_SEM);
            }else if(proj_dist_SEM > min_proj_dist_SEM && projSt.PRIMARY)
            {
                // If we only look for the primary family, we can stop at the first
                // minimum without any risk of losing good solutions. If
                // proj_dist_SEM > min_proj_dist_SEM, we know that we have passed the
                // first minimum, and we can break the loop.
                break;
            }

        }

        //--------------------------------------------------------------------------------
        // We check if we have the primary family, if necessary
        // To do so, we check that we have made two clockwise/counterclockwise
        //  crossings of x = -1.
        //
        // To do a clockwise turn:
        //  1. x1 > 0 -> x2 < 0 && y1 > 0
        //  2. x1 < 0 -> x2 > 0 && y1 < 0
        //
        // To do a counterclockwise turn:
        //  1. x1 > 0 -> x2 < 0 && y1 < 0
        //  2. x1 < 0 -> x2 > 0 && y1 > 0
        //--------------------------------------------------------------------------------
        if(min_proj_dist_SEM < ePdef)
        {
            // Init temporary objects
            OdeEvent odeEvent_purge(true, false);
            double** y_purge = dmatrix(0, 5, 0, 2);
            double* t_purge  = dvector(0, 2);

            // Initial condition in NCEM coordinates
            tv  = projResSt.init_time;
            for(int i = 0; i < 6; i++) yv[i] = projResSt.init_state_CMU_NC[i];

            // Compute the event using ode78 on [tv, t_man_EM[kmin]]
            ode78(y_purge, t_purge, &odeEvent_purge, tv, t_man_SEM[kmin]/SEML.us_em.ns, yv, 6, 2, I_NCSEM, NCEM, NCSEM);

            // Save value in crossings_NCSEM
            crossings_NCSEM = odeEvent_purge.crossings;

            // Save value in crossings_lunar
            crossings_lunar = odeEvent_purge.crossings_moon;

            // Free temporary objects
            free_dmatrix(y_purge, 0, 5, 0, 2);
            free_dvector(t_purge, 0, 2);
        }

        //--------------------------------------------------------------------------------
        //We check for collisions
        //--------------------------------------------------------------------------------
        if(odeEvent.coll)
        {
            collision_NCEM = odeEvent.coll;
        }

        //--------------------------------------------------------------------------------
        //If we save ONLY the primary family
        //--------------------------------------------------------------------------------
        if(projSt.PRIMARY)
        {
            if(collision_NCEM != 0 || crossings_NCSEM != 2.0) min_proj_dist_SEM = ePdef;
        }

        //--------------------------------------------------------------------------------
        //If we save ONLY the primary family
        //--------------------------------------------------------------------------------
        if(crossings_lunar > 0)
        {
            //cout << "There has been " << crossings_lunar << "crossings of the lunar orbits" << endl;
        }
    }


    //====================================================================================
    // 3.3. Save outputs
    //====================================================================================
    //Initial position in SEM coordinates
    for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][0];
    NCSEMmtoSEMm(t_man_SEM[0], yv, yv_SEM, &SEML_SEM);
    for(int i = 0; i < 6; i++) projResSt.init_state_CMU_FC_o[i] = yv_SEM[i];

    //Minimum projection distance
    projResSt.min_proj_dist_SEM_o     = min_proj_dist_SEM;

    //DV at projection
    projResSt.dv_at_projection_FC_o = dv_at_projection_SEM;

    //Crossings
    projResSt.crossings_NCSEM_o       = crossings_NCSEM;

    //Collisions
    projResSt.collision_NCEM_o        = collision_NCEM;

    //Final time
    projResSt.final_time_FC_o       = t_man_SEM[kmin];


    //====================================================================================
    //Free
    //====================================================================================
    free_dmatrix(y_man_NCSEM, 0, 5, 0, projSt.MSIZE);
    free_dvector(t_man_SEM, 0, projSt.MSIZE);


    return GSL_SUCCESS;
}

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
int int_proj_CMU_EM_on_CM_SEM_3D(ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Filenames
    //====================================================================================
    // Get the filename from the latest data of type TYPE_CU_3D
    string filename_input = projSt.get_filename(TYPE_CU_3D, ios::in);
    // Get the filename for the output
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ_3D, ios::out);

    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //Read data size
    int t_grid_size, si_grid_size[4], offset;
    offset = getLenghtCU_bin_3D(si_grid_size, &t_grid_size, filename_input);

    if(offset < 0)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM_3D. Impossible to get length of the data file." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Init
    //====================================================================================
    cout << resetiosflags(ios::scientific) << setprecision(15);

    //------------------------------------------------------------------------------------
    //To store all data
    //------------------------------------------------------------------------------------
    double** init_state_CMU_NCEM = dmatrix(0, 5, 0, si_grid_size[2]);
    double** init_state_CMU_RCM  = dmatrix(0, 4, 0, si_grid_size[2]);
    double* init_time_grid_EM    = dvector(0, t_grid_size);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+si_grid_size[0])*(1+si_grid_size[1])*(1+si_grid_size[2])*(1+si_grid_size[3])*(1+t_grid_size);


    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM_3D. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
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
        offset = readTCU_bin_3D(offset, init_time_grid_EM, kt, filename_input);

        if(offset < 0)
        {
            cout << "int_proj_CMU_EM_on_CM_SEM_3D. Impossible to read the data file." << endl;
            return FTC_FAILURE;
        }

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
                                           init_state_CMU_RCM, si_grid_size,
                                           filename_input);
                    if(offset < 0)
                    {
                        cout << "int_proj_CMU_EM_on_CM_SEM_3D. Impossible to read the data file." << endl;
                        return FTC_FAILURE;
                    }

                    //--------------------------------------------------------------------
                    //Most inner loop is parallelized
                    //--------------------------------------------------------------------
                    #pragma omp parallel for if(isPar)  shared(index)
                    for(int ks3 = 0; ks3 <= si_grid_size[2]; ks3++)
                    {

                        //----------------------------------------------------------------
                        //Inputs
                        //----------------------------------------------------------------
                        ProjResSt projResSt(NCEM);

                        //----------------------------------------------------------------
                        //Initial conditions
                        //----------------------------------------------------------------
                        // Init the time in EM coordinates
                        projResSt.init_time  = init_time_grid_EM[kt];
                        // Init the state in NCEM coordinates
                        for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = init_state_CMU_NCEM[i][ks3];
                        // Init the state in RCM coordinates
                        for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = init_state_CMU_RCM[i][ks3];
                        //Label
                        projResSt.label = 0;

                        //----------------------------------------------------------------
                        //Projection on center manifold at SEML2
                        //----------------------------------------------------------------
                        proj_subroutine(projResSt, invman_SEM, projSt);

                        //----------------------------------------------------------------
                        // Save outputs
                        //----------------------------------------------------------------
                        if(projResSt.min_proj_dist_SEM_o < ePdef)
                        {
                            //Save
                            #pragma omp critical
                            {
                                writeIntProjCU_bin(filename_output, projResSt);
                            }
                        }

                        //----------------------------------------------------------------
                        //Display completion
                        //----------------------------------------------------------------
                        #pragma omp critical
                        {
                            displayCompletion("int_proj_CMU_EM_on_CM_SEM_3D", 100.0*index++/noe);
                        }

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
 *      "plot/QBCP/EM/L2/projcu_order_16.bin"
 **/
int int_proj_CMU_EM_on_CM_SEM(ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Filenames
    //====================================================================================
    // Get the filename from the latest data of type TYPE_CU
    string filename_input = projSt.get_filename(TYPE_CU, ios::in);
    // Get the filename for the output
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ, ios::out);

    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Read data size
    //------------------------------------------------------------------------------------
    int t_grid_size, s1_grid_size, s3_grid_size;
    int status = getLenghtCU_bin(&s1_grid_size, &s3_grid_size, &t_grid_size, filename_input);

    if(status != FTC_SUCCESS)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. Impossible to get length of the data file." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //To store all data
    //------------------------------------------------------------------------------------
    double**** init_state_CMU_NCEM = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double**** init_state_CMU_RCM  = d4tensor(0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double** *  init_dH_valid      = d3tensor(0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double* init_time_grid_EM      = dvector(0, t_grid_size);

    //------------------------------------------------------------------------------------
    //Read data from file
    //------------------------------------------------------------------------------------
    status = readCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, init_dH_valid,
                        init_time_grid_EM, s1_grid_size, s3_grid_size, t_grid_size,
                        filename_input);

    if(status != FTC_SUCCESS)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. Impossible to read data file." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Splash screen
    //====================================================================================
    cout << "===================================================================" << endl;
    cout << "              Computation of the connections between:              " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << s1_grid_size+1  << " value(s) of s1" << endl;
    cout << "  - " << s3_grid_size+1  << " value(s) of s3" << endl;
    cout << "  - " << t_grid_size+1   << " value(s) of t"  << endl;
    cout << " The CMU data are taken from " << filename_input << endl;
    cout << " The projection data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;


    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
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
                // We check that the energy condition is valid
                //========================================================================
                if(init_dH_valid[kt][ks1][ks3])
                {
                    //--------------------------------------------------------------------
                    //Inputs
                    //--------------------------------------------------------------------
                    ProjResSt projResSt(NCEM);

                    //--------------------------------------------------------------------
                    //Initial conditions
                    //--------------------------------------------------------------------
                    // Init the time in EM coordinates
                    projResSt.init_time  = init_time_grid_EM[kt];
                    // Init the state in NCEM coordinates
                    for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = init_state_CMU_NCEM[i][kt][ks1][ks3];
                    // Init the state in RCM coordinates
                    for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = init_state_CMU_RCM[i][kt][ks1][ks3];
                    //Label
                    projResSt.label = 0;

                    //--------------------------------------------------------------------
                    //Projection on center manifold at SEML2
                    //--------------------------------------------------------------------
                    proj_subroutine(projResSt, invman_SEM, projSt);

                    //--------------------------------------------------------------------
                    // Save outputs
                    //--------------------------------------------------------------------
                    if(projResSt.min_proj_dist_SEM_o < ePdef)
                    {
                        //Save
                        #pragma omp critical
                        {
                            writeIntProjCU_bin(filename_output, projResSt);
                        }
                    }
                }

                //------------------------------------------------------------------------
                //Display completion
                //------------------------------------------------------------------------
                #pragma omp critical
                {
                    displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/((1+s1_grid_size)*(1+s3_grid_size)*(1+t_grid_size)));
                }
            }
        }
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_RCM, 0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d3tensor(init_dH_valid,  0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_dvector(init_time_grid_EM, 0, t_grid_size);


    return FTC_SUCCESS;
}

/**
 *  \brief Same as int_proj_CMU_EM_on_CM_SEM but for the outputs from the routine
 *         compute_grid_CMU_EM_dH.
 *
 *         Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM_dH.
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
 *      "plot/QBCP/EM/L2/projcu_order_16_dH_0.005.bin"
 **/
int int_proj_CMU_EM_on_CM_SEM_dH(ProjSt& projSt)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;

    //====================================================================================
    // Filenames
    //====================================================================================
    // Get the filename from the latest data of type TYPE_CU
    string filename_input = projSt.get_filename(TYPE_CU, ios::in);
    // Get the filename for the output
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ, ios::out);

    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Read data size
    //------------------------------------------------------------------------------------
    int t_grid_size, s1_grid_size;
    int status = getLenghtCU_bin_dH(&s1_grid_size, &t_grid_size, filename_input);

    if(status != FTC_SUCCESS)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. Impossible to get length of the data file." << endl;
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //To store all data
    //------------------------------------------------------------------------------------
    double** * init_state_CMU_NCEM = d3tensor(0, 5, 0, t_grid_size, 0, s1_grid_size);
    double** * init_state_CMU_RCM  = d3tensor(0, 4, 0, t_grid_size, 0, s1_grid_size);
    double**   init_dH_valid       = dmatrix(0, t_grid_size, 0, s1_grid_size);
    double* init_time_grid_EM      = dvector(0, t_grid_size);

    //------------------------------------------------------------------------------------
    //Read data from file
    //------------------------------------------------------------------------------------
    status = readCU_bin_dH(init_state_CMU_NCEM, init_state_CMU_RCM, init_dH_valid,
                           init_time_grid_EM, s1_grid_size, t_grid_size, filename_input);

    if(status != FTC_SUCCESS)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. Impossible to read data file." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Splash screen
    //====================================================================================
    cout << "===================================================================" << endl;
    cout << "              Computation of the connections between:              " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER << endl;
    cout << "  - " << s1_grid_size+1  << " value(s) of s1" << endl;
    cout << "  - " << t_grid_size+1   << " value(s) of t"  << endl;
    cout << " The data will be stored in " << filename_output << endl;
    cout << setiosflags(ios::scientific) << setprecision(15);
    cout << "===================================================================" << endl;

    //====================================================================================
    // 2.2. Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
    filestream.close();

    //====================================================================================
    // 3. Loop: only the inner loop is parallelized, since
    //    open_mp doest not allow nested // loops by default
    //====================================================================================
    COMPLETION = 0;
    int index  = 0;
    for(int kt = 0; kt <= t_grid_size; kt++)
    {
        #pragma omp parallel for if(isPar) shared(index)
        for(int ks1 = 0; ks1 <= s1_grid_size; ks1++)
        {
            //============================================================================
            // We check that the energy condition is valid
            //============================================================================
            if(init_dH_valid[kt][ks1])
            {
                //------------------------------------------------------------------------
                //Inputs
                //------------------------------------------------------------------------
                ProjResSt projResSt(NCEM);

                //------------------------------------------------------------------------
                //Initial conditions
                //------------------------------------------------------------------------
                // Init the time in EM coordinates
                projResSt.init_time  = init_time_grid_EM[kt];
                // Init the state in NCEM coordinates
                for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = init_state_CMU_NCEM[i][kt][ks1];
                // Init the state in RCM coordinates
                for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = init_state_CMU_RCM[i][kt][ks1];
                //Label
                projResSt.label = 0;

                //------------------------------------------------------------------------
                //Projection on center manifold at SEML2
                //------------------------------------------------------------------------
                proj_subroutine(projResSt, invman_SEM, projSt);

                //------------------------------------------------------------------------
                // Save outputs
                //------------------------------------------------------------------------
                if(projResSt.min_proj_dist_SEM_o < ePdef)
                {
                    //Save
                    #pragma omp critical
                    {
                        writeIntProjCU_bin(filename_output, projResSt);
                    }
                }
            }

            //----------------------------------------------------------------------------
            //Display completion
            //----------------------------------------------------------------------------
            #pragma omp critical
            {
                displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/((1+s1_grid_size)*(1+t_grid_size)));
            }
        }
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_d3tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size);
    free_d3tensor(init_state_CMU_RCM, 0, 4, 0, t_grid_size, 0, s1_grid_size);
    free_dmatrix(init_dH_valid,  0, t_grid_size, 0, s1_grid_size);
    free_dvector(init_time_grid_EM, 0, t_grid_size);


    return FTC_SUCCESS;
}

/**
 *  \brief Detection of the possible connections between a set of orbits about EMLi
 *         and the center manifold at SEMLj.
 *         The routine works with a system of seeds in RCM coordinates.
 *
 *         For all projSt.GSIZE_SI[0] values of s1 in
 *         [projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1]], we define the following initial
 *         conditions:
 *                       s0 = (s1, projSt.GLIM_SI[1][0], -s1, projSt.GLIM_SI[3][0])
 *
 *         Note that s3 = -s1 guarantees that the IC are on the y = 0 plane in NC coord.
 *         (at least in the planar case).
 *
 *
 *
 *         Then, for all projSt.TSIZE value of initial time t0 in
 *         [projSt.TLIM[0], projSt.TLIM[1] ], we cen compute z0 = W(s0, t0). We can then
 *         project z0 on the center manifold, in order to get
 *                  si = (s1, s3, s3, s4) = W^{-1}(z0, t0),
 *
 *         Roughly speaking, si provides IC for a similar orbit but a different starting
 *         phase for the Sun-Earth-Moon system. These IC define a given orbit
 *         whose period is estimated via cmu_orbit_estimate_period. We can then integrate
 *         this orbit on one of its period and look for connections on a fine grid along
 *         the resulting trajectory
 *
 *  The output data are saved in a binary file of the form
 *      "plot/QBCP/EM/L2/projcu_order_16_Orbit.bin"
 **/
int int_proj_ORBIT_EM_on_CM_SEM(ProjSt& projSt, int Nperiods)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;
    double dt  = projSt.dt;

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Invariant manifold
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);
    int ncs = invman.getNCS();
    int dcs = default_coordinate_system(ncs);

    //------------------------------------------------------------------------------------
    // ODE
    //------------------------------------------------------------------------------------
    //Driver
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------------------------------------------------------------------------
    //Orbit
    //------------------------------------------------------------------------------------
    Orbit orbit(&invman, &SEML, &odestruct, OFTS_ORDER, OFS_ORDER, projSt.TLIM[0], projSt.TLIM[0]+10);


    //====================================================================================
    // Building the initial conditions on the orbit
    //====================================================================================
    double t0;
    double* st0 = dvector(0,4);

    // Initial conditions @ t0
    t0 =  projSt.TLIM[0];
    for(int i = 0; i < 4; i++) st0[i] = projSt.GLIM_SI[i][0];
    st0[4] = 0.0;

    //====================================================================================
    // Filename
    //====================================================================================
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ_ORBIT, ios::out);

    //====================================================================================
    // Splash screen
    //====================================================================================
    cout << resetiosflags(ios::scientific) << setprecision(4);
    cout << "===================================================================" << endl;
    cout << "              Computation of the connections of a single orbit:    " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER                                      << endl;
    cout << " - " << projSt.TSIZE+1  << " values of times"                                 ;
    cout << " in ["  << projSt.TLIM[0]/SEML.us_em.T                                      ;
    cout << ", " << projSt.TLIM[1]/SEML.us_em.T  << "]"                           << endl;
    cout << " - " << projSt.GSIZE_SI[0]+1  << " values of s1/s3"                           ;
    cout << " in ["  << projSt.GLIM_SI[0][0]                                             ;
    cout << ", " << projSt.GLIM_SI[0][1]  << "]"                                  << endl;
    cout << " The data will be stored in " << filename_output                     << endl;
    cout << "===================================================================" << endl;

    //------------------------------------------------------------------------------------
    //Building the working grids
    //------------------------------------------------------------------------------------
    double* grid_s1_CMU_RCM = dvector(0,  projSt.GSIZE_SI[0]);
    //double* grid_s3_CMU_RCM = dvector(0,  projSt.GSIZE_SI[2]);
    init_grid(grid_s1_CMU_RCM, projSt.GLIM_SI[0][0], projSt.GLIM_SI[0][1], projSt.GSIZE_SI[0]);
    //init_grid(grid_s3_CMU_RCM, projSt.GLIM_SI[2][0], projSt.GLIM_SI[2][1], projSt.GSIZE_SI[2]);

    cout << " - The detailed values of s1 are:                                   " << endl;
    for(int i = 0; i < projSt.GSIZE_SI[0]; i++) cout << grid_s1_CMU_RCM[i] << ", ";
    cout << grid_s1_CMU_RCM[projSt.GSIZE_SI[0]]  << endl;
    cout << "===================================================================" << endl;

    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_EM[i]/SEML.us->T << ", ";
    cout << grid_t_EM[projSt.TSIZE]/SEML.us->T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);
    cout << setiosflags(ios::scientific) << setprecision(15);

    //====================================================================================
    // Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);


    //====================================================================================
    // Misc initialization. Require better presentation?
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    //double ePdef = 1e5;

    //====================================================================================
    // Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
    filestream.close();

    //====================================================================================
    // Loop: only the inner loop is parallelized, since
    //    open_mp doest not allow nested // loops by default
    //====================================================================================
    int index  = 0;
    int label  = 0;
    COMPLETION = 0;
    double stp[5];

    //s0 loop
    for(int ks1 = 0; ks1 <= projSt.GSIZE_SI[0]; ks1++)
    {

        //--------------------------------------------------------------------------------
        // Changing the initial conditions on the orbit
        //--------------------------------------------------------------------------------
        st0[0] =  grid_s1_CMU_RCM[ks1];
        st0[1] =  projSt.GLIM_SI[1][0];
        st0[2] = -grid_s1_CMU_RCM[ks1];
        st0[3] =  projSt.GLIM_SI[3][0];
        st0[4] =  0.0;

        t0 = grid_t_EM[0]; //just for the specific purpose of estimating the period, right after

        //--------------------------------------------------------------------------------
        // Estimate the period and the number of points on the time grid to match a
        // frequency of dt over this period
        //--------------------------------------------------------------------------------
        double Tp = 0.0;
        int N = 0.0;
        cmu_orbit_estimate_period(st0, t0, &Tp, &N, dt, orbit);

        // We multiply Tp and N by the number of desired periods
        Tp *= Nperiods;
        N  *= Nperiods;

        //Save the initial position
        double z0[6];
        state_memcpy(z0, orbit.getZ0());

        //--------------------------------------------------------------------------------
        // Storing initial conditions along the orbit
        //--------------------------------------------------------------------------------
        //To store data
        double** yNCE = dmatrix(0, 5, 0, N);
        double** sRCM = dmatrix(0, 4, 0, N);
        double* tNCE  = dvector(0, N);

        //--------------------------------------------------------------------------------
        //Number of elements
        //--------------------------------------------------------------------------------
        int noe = (1+N)*(1+projSt.TSIZE)*(1+projSt.GSIZE_SI[0]);


        // Time loop
        for(int kt = 0; kt <= projSt.TSIZE; kt++)
        {
            //----------------------------------------------------------------------------
            //Strob map & unstable directions
            //----------------------------------------------------------------------------
            t0 = grid_t_EM[kt];

            //Projection
            orbit.NCprojCCMtoCM(z0, t0, stp);
            stp[4] =  0.0;

            cmu_grid_orbit_on_one_period(orbit, tNCE, yNCE, sRCM, stp, t0, Tp, N, isPar, projSt.hyp_epsilon_eml2);

            //----------------------------------------------------------------------------
            //Once the unstable directions are obtained, we propagate & project
            //----------------------------------------------------------------------------
            #pragma omp parallel for if(isPar) shared(index)
            for(int ks = 0; ks <= N; ks++)
            {
                //------------------------------------------------------------------------
                //Inputs
                //------------------------------------------------------------------------
                ProjResSt projResSt(NCEM);

                //------------------------------------------------------------------------
                //Seeds
                //------------------------------------------------------------------------
                //Init time (global)
                projResSt.seed_time = t0;
                //RCM state (seed)
                for(int i = 0; i < 4; i++) projResSt.seed_state_CMU_RCM[i] = st0[i];

                //------------------------------------------------------------------------
                //Initial conditions
                //------------------------------------------------------------------------
                // Init the time in EM coordinates
                projResSt.init_time  = tNCE[ks];
                // Init the state in NCEM coordinates
                for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = yNCE[i][ks];
                // Init the state in RCM coordinates
                for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = sRCM[i][ks];
                //Label
                projResSt.label = label;

                //------------------------------------------------------------------------
                //Projection on center manifold at SEML2
                //------------------------------------------------------------------------
                proj_subroutine(projResSt, invman_SEM, projSt);

                //------------------------------------------------------------------------
                // Save outputs
                //------------------------------------------------------------------------
                //if(projResSt.min_proj_dist_SEM_o < ePdef)
                //{
                //Save
                #pragma omp critical
                {
                    writeIntProjCUSeed_bin(filename_output, projResSt);
                }
                //}


                //------------------------------------------------------------------------
                //Display completion
                //------------------------------------------------------------------------
                #pragma omp critical
                {
                    displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/noe);
                }
            }

            label++;
        }

        //--------------------------------------------------------------------------------
        //Free
        //--------------------------------------------------------------------------------
        free_dmatrix(yNCE, 0, 5, 0, N);
        free_dmatrix(sRCM, 0, 4, 0, N);
        free_dvector(tNCE, 0, N);
    }

    return FTC_SUCCESS;
}

/**
 *  \brief Detection of the possible connections between a set of orbits about EMLi
 *         and the center manifold at SEMLj.
 *         We define the following initial
 *         conditions:
 *         s0 = (projSt.GLIM_SI[0][0], projSt.GLIM_SI[1][0], projSt.GLIM_SI[2][0], projSt.GLIM_SI[3][0])
 *
 *         Then, for all projSt.TSIZE value of initial time t0 in
 *         [projSt.TLIM[0], projSt.TLIM[1] ], we cen compute z0 = W(s0, t0). We can then
 *         project z0 on the center manifold, in order to get
 *                  si = (s1, s3, s3, s4) = W^{-1}(z0, t0),
 *
 *         Roughly speaking, si provides IC for a similar orbit but a different starting
 *         phase for the Sun-Earth-Moon system. These IC define a given orbit
 *         whose period is estimated via cmu_orbit_estimate_period. We can then integrate
 *         this orbit on one of its period and look for connections on a fine grid along
 *         the resulting trajectory
 *
 *  The output data are saved in a binary file of the form
 *      "plot/QBCP/EM/L2/projcu_order_16_Orbit.bin"
 **/
int int_proj_SINGLE_ORBIT_EM_on_CM_SEM(ProjSt& projSt, int Nperiods)
{
    //====================================================================================
    // Retrieve the parameters in the projection structure
    //====================================================================================
    bool isPar = projSt.ISPAR;
    double dt  = projSt.dt;

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Invariant manifold
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);
    int ncs = invman.getNCS();
    int dcs = default_coordinate_system(ncs);

    //------------------------------------------------------------------------------------
    // ODE
    //------------------------------------------------------------------------------------
    //Driver
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------------------------------------------------------------------------
    //Orbit
    //------------------------------------------------------------------------------------
    Orbit orbit(&invman, &SEML, &odestruct, OFTS_ORDER, OFS_ORDER, projSt.TLIM[0], projSt.TLIM[0]+10);


    //====================================================================================
    // Building the initial conditions on the orbit
    //====================================================================================
    double t0;
    double* st0 = dvector(0,4);

    // Initial conditions @ t0
    t0 =  projSt.TLIM[0];
    for(int i = 0; i < 4; i++) st0[i] = projSt.GLIM_SI[i][0];
    st0[4] = 0.0;

    //====================================================================================
    // Filename
    //====================================================================================
    string filename_output = projSt.get_and_update_filename(projSt.FILE_PCU, TYPE_MAN_PROJ_ORBIT, ios::out);

    //====================================================================================
    // Splash screen
    //====================================================================================
    cout << resetiosflags(ios::scientific) << setprecision(4);
    cout << "===================================================================" << endl;
    cout << "              Computation of the connections of a single orbit:    " << endl;
    cout << "===================================================================" << endl;
    cout << " The computation domain is the following:                          " << endl;
    cout << " - OFTS_ORDER = " << OFTS_ORDER                                      << endl;
    cout << " - " << projSt.TSIZE+1  << " values of times"                               ;
    cout << " in ["  << projSt.TLIM[0]/SEML.us_em.T                                      ;
    cout << ", " << projSt.TLIM[1]/SEML.us_em.T  << "]"                           << endl;
    cout << " - The initial condition in RCM coordinates are " << endl;
    cout << " (" << st0[0]<<", "<<st0[1] << ", " <<st0[2]<< ", " <<st0[3]<< ")"   << endl;
    cout << " The data will be stored in " << filename_output                     << endl;
    cout << "===================================================================" << endl;

    //------------------------------------------------------------------------------------
    //Building the time grid
    //------------------------------------------------------------------------------------
    double* grid_t_EM = dvector(0,  projSt.TSIZE);
    init_grid(grid_t_EM, projSt.TLIM[0], projSt.TLIM[1], projSt.TSIZE);

    cout << " - The detailed values of t are:                                   " << endl;
    for(int i = 0; i < projSt.TSIZE; i++) cout << grid_t_EM[i]/SEML.us->T << ", ";
    cout << grid_t_EM[projSt.TSIZE]/SEML.us->T  << endl;
    cout << "===================================================================" << endl;
    pressEnter(true);
    cout << setiosflags(ios::scientific) << setprecision(15);

    //====================================================================================
    // Initialize tools for the projection phase. Namely the center manifold @SEML
    //====================================================================================
    //First, check that the type of manifold provided is good:
    if(SEML_SEM.cs->manType != MAN_CENTER)
    {
        cout << "int_proj_CMU_EM_on_CM_SEM. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return FTC_FAILURE;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);


    //====================================================================================
    // Misc initialization. Require better presentation?
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Notable points in SEM system
    //------------------------------------------------------------------------------------
    double** semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //------------------------------------------------------------------------------------
    // projection by default is arbitrary big
    //------------------------------------------------------------------------------------
    //double ePdef = 1e5;

    //====================================================================================
    // Reset the data file (projection)
    //====================================================================================
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename_output.c_str(), ios::binary | ios::out);
    filestream.close();

    //====================================================================================
    // Loop: only the inner loop is parallelized, since
    //    open_mp doest not allow nested // loops by default
    //====================================================================================
    int index  = 0;
    int label  = 0;
    COMPLETION = 0;
    double stp[5];

    //------------------------------------------------------------------------------------
    //just for the specific purpose of estimating the period, right after
    //------------------------------------------------------------------------------------
    t0 = grid_t_EM[0];

    //------------------------------------------------------------------------------------
    // Estimate the period and the number of points on the time grid to match a
    // frequency of dt over this period
    //------------------------------------------------------------------------------------
    double Tp = 0.0;
    int N = 0.0;
    cmu_orbit_estimate_period(st0, t0, &Tp, &N, dt, orbit);

    // We multiply Tp and N by the number of desired periods
    Tp *= Nperiods;
    N  *= Nperiods;

    //Save the initial position
    double z0[6];
    state_memcpy(z0, orbit.getZ0());

    //------------------------------------------------------------------------------------
    // Storing initial conditions along the orbit
    //------------------------------------------------------------------------------------
    //To store data
    double** yNCE = dmatrix(0, 5, 0, N);
    double** sRCM = dmatrix(0, 4, 0, N);
    double* tNCE  = dvector(0, N);

    //------------------------------------------------------------------------------------
    //Number of elements
    //------------------------------------------------------------------------------------
    int noe = (1+N)*(1+projSt.TSIZE)*(1+projSt.GSIZE_SI[0]);


    // Time loop
    for(int kt = 0; kt <= projSt.TSIZE; kt++)
    {
        //--------------------------------------------------------------------------------
        //Strob map & unstable directions
        //--------------------------------------------------------------------------------
        t0 = grid_t_EM[kt];

        //Projection
        orbit.NCprojCCMtoCM(z0, t0, stp);
        stp[4] =  0.0;

        cmu_grid_orbit_on_one_period(orbit, tNCE, yNCE, sRCM, stp, t0, Tp, N, isPar, projSt.hyp_epsilon_eml2);

        //--------------------------------------------------------------------------------
        //Once the unstable directions are obtained, we propagate & project
        //--------------------------------------------------------------------------------
        #pragma omp parallel for if(isPar) shared(index)
        for(int ks = 0; ks <= N; ks++)
        {
            //----------------------------------------------------------------------------
            //Inputs
            //----------------------------------------------------------------------------
            ProjResSt projResSt(NCEM);

            //----------------------------------------------------------------------------
            //Seeds
            //----------------------------------------------------------------------------
            //Init time (global)
            projResSt.seed_time = t0;
            //RCM state (seed)
            for(int i = 0; i < 4; i++) projResSt.seed_state_CMU_RCM[i] = st0[i];

            //----------------------------------------------------------------------------
            //Initial conditions
            //----------------------------------------------------------------------------
            // Init the time in EM coordinates
            projResSt.init_time  = tNCE[ks];
            // Init the state in NCEM coordinates
            for(int i = 0; i < 6; i++) projResSt.init_state_CMU_NC[i] = yNCE[i][ks];
            // Init the state in RCM coordinates
            for(int i = 0; i < 5; i++) projResSt.init_state_CMU_RCM[i] = sRCM[i][ks];
            //Label
            projResSt.label = label;

            //----------------------------------------------------------------------------
            //Projection on center manifold at SEML2
            //----------------------------------------------------------------------------
            proj_subroutine(projResSt, invman_SEM, projSt);

            //----------------------------------------------------------------------------
            // Save outputs
            //----------------------------------------------------------------------------
            //if(projResSt.min_proj_dist_SEM_o < ePdef)
            //{
            //Save
            #pragma omp critical
            {
                writeIntProjCUSeed_bin(filename_output, projResSt);
            }
            //}


            //----------------------------------------------------------------------------
            //Display completion
            //----------------------------------------------------------------------------
            #pragma omp critical
            {
                displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/noe);
            }
        }

        label++;
    }

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(sRCM, 0, 4, 0, N);
    free_dvector(tNCE, 0, N);


    return FTC_SUCCESS;
}

//========================================================================================
//
//         Initial conditions for projection of single orbits
//
//========================================================================================
/**
 *  \brief Computes the stroboscopic map on N+1 points along a trajectory in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] = t0 + k*SEML.us->T, for all k in [0, N],
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_strob(double* tNCE, double** yNCE, double** sRCM, double st0[], double t0, int N, int isPar, double hyp_epsilon)
{
    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Invariant manifold
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);
    int ncs = invman.getNCS();
    int dcs = default_coordinate_system(ncs);

    //------------------------------------------------------------------------------------
    // ODE
    //------------------------------------------------------------------------------------
    //Hard precision
    Config::configManager().C_PREC_HARD();

    //Driver
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------------------------------------------------------------------------
    //Orbit
    //------------------------------------------------------------------------------------
    double tf = t0 + N*SEML.us->T;
    Orbit orbit(&invman, &SEML, &odestruct, OFTS_ORDER, OFS_ORDER, t0, tf);

    //Orbit IC
    orbit.update_ic(st0, t0);

    //====================================================================================
    //Integration
    //====================================================================================
    int status = orbit.traj_int_grid(tf, yNCE, tNCE, N, true);

    //====================================================================================
    //Projection
    //====================================================================================
    orbit.proj_traj_grid(sRCM, yNCE, tNCE, N);

    //====================================================================================
    // Building the initial conditions on the unstable manifold of the orbit
    //====================================================================================
    int iter = 1;
    COMPLETION = 0;
    #pragma omp parallel for if(isPar)  shared(iter)
    for(int k = 0; k <= N; k++)
    {
        Ofsc ofs(OFS_ORDER);
        double* yvu = dvector(0,5);
        double* sti = dvector(0,4);


        //--------------------------------------------------------------------------------
        // Initialization on the center-unstable manifold
        //--------------------------------------------------------------------------------
        //Init sti
        sti[0] = sRCM[0][k];
        sti[1] = sRCM[1][k];
        sti[2] = sRCM[2][k];
        sti[3] = sRCM[3][k];
        sti[4] = hyp_epsilon;


        //--------------------------------------------------------------------------------
        //Equivalent state
        //--------------------------------------------------------------------------------
        invman.evalRCMtoNC(sti, tNCE[k], yvu, OFTS_ORDER, OFS_ORDER);

        //--------------------------------------------------------------------------------
        //Save
        //--------------------------------------------------------------------------------
        #pragma omp critical
        {

            for(int i = 0; i < 6; i++) yNCE[i][k] = yvu[i];
            for(int i = 0; i < 5; i++) sRCM[i][k] = sti[i];
            //Display
            displayCompletion("compute_grid_CMU_EM", (double) iter++/(N+1)*100);
        }

        free_dvector(yvu, 0, 5);
        free_dvector(sti, 0, 4);

    }

    //====================================================================================
    //Back to original precision
    //====================================================================================
    Config::configManager().C_PREC_BACK();


    return status;
}

/**
 *  \brief Computes N+1 points along a trajectory in the QBCP, on NPeriods periods T
 *         of the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] is the time, in [0 NPeriods*T]
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_orbit(double* tNCE, double** yNCE, double** sRCM, double st0[], double t0, int N, int NPeriods,  int isPar, double hyp_epsilon)
{
    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Invariant manifold
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);
    int ncs = invman.getNCS();
    int dcs = default_coordinate_system(ncs);

    //------------------------------------------------------------------------------------
    // ODE
    //------------------------------------------------------------------------------------
    //Hard precision
    Config::configManager().C_PREC_HARD();

    //Driver
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------------------------------------------------------------------------
    //Orbit
    //------------------------------------------------------------------------------------
    double tf = t0 + NPeriods*SEML.us->T;
    Orbit orbit(&invman, &SEML, &odestruct, OFTS_ORDER, OFS_ORDER, t0, tf);

    //Orbit IC
    orbit.update_ic(st0, t0);

    //====================================================================================
    //Integration
    //====================================================================================
    int status = orbit.traj_int_grid(tf, yNCE, tNCE, N, true);

    //====================================================================================
    //Projection
    //====================================================================================
    orbit.proj_traj_grid(sRCM, yNCE, tNCE, N);

    //====================================================================================
    // Building the initial conditions on the unstable manifold of the orbit
    //====================================================================================
    int iter = 1;
    COMPLETION = 0;
    #pragma omp parallel for if(isPar)  shared(iter)
    for(int k = 0; k <= N; k++)
    {
        Ofsc ofs(OFS_ORDER);
        double* yvu = dvector(0,5);
        double* sti = dvector(0,4);


        //--------------------------------------------------------------------------------
        // Initialization on the center-unstable manifold
        //--------------------------------------------------------------------------------
        //Init sti
        sti[0] = sRCM[0][k];
        sti[1] = sRCM[1][k];
        sti[2] = sRCM[2][k];
        sti[3] = sRCM[3][k];
        sti[4] = hyp_epsilon;


        //--------------------------------------------------------------------------------
        //Equivalent state
        //--------------------------------------------------------------------------------
        invman.evalRCMtoNC(sti, tNCE[k], yvu, OFTS_ORDER, OFS_ORDER);

        //--------------------------------------------------------------------------------
        //Save
        //--------------------------------------------------------------------------------
        #pragma omp critical
        {

            for(int i = 0; i < 6; i++) yNCE[i][k] = yvu[i];
            for(int i = 0; i < 5; i++) sRCM[i][k] = sti[i];
            //Display
            displayCompletion("compute_grid_CMU_EM", (double) iter++/(N+1)*100);
        }

        free_dvector(yvu, 0, 5);
        free_dvector(sti, 0, 4);

    }

    //====================================================================================
    //Back to original precision
    //====================================================================================
    Config::configManager().C_PREC_BACK();


    return status;
}

/**
 *  \brief Estimates the period T (in adimensionalized units) of a given orbit Orbit
 *         in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *         The period is estimated as follows:
 *         The initial condition are z0v = (x0, y0, z0, px0, py0, pz0) = W(s0, t0).
 *         The polar argument of (x0, y0) is computed.
 *         And the integration is stopped when a similar argument is obtained.
 *         The time of the event is the desired time T.
 *
 *         Moreover, the number of points N is computed as the desired number of points to
 *         achieve a frequency of dt in [0 T].
 *         Hence, N = floor(T/dt);
 **/
int cmu_orbit_estimate_period(const double st0[], double t0, double* T, int* N, double dt, Orbit& orbit)
{
    //====================================================================================
    //Orbit IC
    //====================================================================================
    int ncs = orbit.getInvman()->getNCS();
    int dcs = default_coordinate_system(ncs);

    orbit.update_ic(st0, t0);

    //====================================================================================
    // Estimate the period of the orbit via event function
    //====================================================================================
    //Hard precision
    Config::configManager().C_PREC_HARD();
    //------------------------------------------------------------------------------------
    //Event structure
    //------------------------------------------------------------------------------------
    double angle0 = atan(orbit.getZ0()[1]/orbit.getZ0()[0]);
    double center[3] = {0,0,0}; //center is the origin
    struct value_params val_par;
    val_par.max_events = 3;
    val_par.direction  = 0;
    val_par.dim        = 0;
    val_par.value      = angle0;
    val_par.center     = center;
    val_par.type       = 'A';

    double** ye_NCSEM = dmatrix(0, 5, 0, val_par.max_events);
    double* te_NCSEM  = dvector(0, val_par.max_events);

    //------------------------------------------------------------------------------------
    // Event detection
    //------------------------------------------------------------------------------------
    int ode78coll;
    double yv[6];
    for(int i = 0; i < 6; i++) yv[i] = orbit.getZ0()[i];
    ode78_qbcp_event(ye_NCSEM, te_NCSEM, &ode78coll, t0, t0+10, yv, 6, dcs, ncs, ncs, &val_par);

    //------------------------------------------------------------------------------------
    //Then, we can redefine tf and N
    //------------------------------------------------------------------------------------
    //Find the closest point to the initial one
    int kargmin = 1;
    double dmin = fabs(orbit.getZ0()[0] - ye_NCSEM[0][1]), dm = 0.0;
    for(int k = min(val_par.max_events, 2); k <  val_par.max_events; k++)
    {
        dm = fabs(orbit.getZ0()[0] - ye_NCSEM[0][k]);
        if(dm < dmin)
        {
            dmin = dm;
            kargmin = k;
        }
    }

    // Estimated period, as a ratio
    double r0 = (te_NCSEM[kargmin] - t0)/SEML.us_em.T;
    cout << "Estimated period is : " << r0 << " xT " << endl;

    // Approximation as a multiple of dt
    double rEst = ( dt*floor(r0/dt) + dt );
    cout << "We use the following approximation : " << rEst << endl;

    // Period, and number of points to match a frequency of "dt" on a period
    *T = rEst * SEML.us_em.T;
    *N = floor(rEst/dt);

    //====================================================================================
    //Back to original precision
    //====================================================================================
    Config::configManager().C_PREC_BACK();

    return GSL_SUCCESS;
}

/**
 *  \brief Computes N+1 points along a trajectory in the QBCP, on one period T
 *         of the orbit Orbit, in the QBCP.
 *
 *         The initial conditions (IC) are (s0, t0), in RCM coordinates.
 *
 *         The routine stores N+1 points corresponding to the following positions:
 *              - tNCE[k] is the time, in [0 T]
 *              - yNCE[k] is the NC state at tNCE[k],
 *                (NC coordinates given by SEML.cs (either NCSEM or NCEM)
 *              - sRCM[k] is the RCM state at tNCE[k].
 *
 *         If isPar == true, the computation and storage of sRCM is made parallel
 *         (but careful with that, may not work with higher parallelized loops).
 **/
int cmu_grid_orbit_on_one_period(Orbit& orbit, double* tNCE, double** yNCE, double** sRCM, const double st0[], double t0, double T, int N, int isPar, double hyp_epsilon)
{
    //====================================================================================
    //Orbit IC
    //====================================================================================
    orbit.update_ic(st0, t0);

    //====================================================================================
    //Integration with Hard precision
    //====================================================================================
    //Config::configManager().C_PREC_HARD();
    int status = orbit.traj_int_grid(t0+T, yNCE, tNCE, N, true);

    //====================================================================================
    //Projection
    //====================================================================================
    orbit.proj_traj_grid(sRCM, yNCE, tNCE, N);

    //====================================================================================
    // Building the initial conditions on the unstable manifold of the orbit
    //====================================================================================
    //int iter = 1;
    //COMPLETION = 0;
    #pragma omp parallel for if(isPar) //shared(iter)
    for(int k = 0; k <= N; k++)
    {
        Ofsc ofs(OFS_ORDER);
        double* yvu = dvector(0,5);
        double* sti = dvector(0,4);


        //--------------------------------------------------------------------------------
        // Initialization on the center-unstable manifold
        //--------------------------------------------------------------------------------
        //Init sti
        sti[0] = sRCM[0][k];
        sti[1] = sRCM[1][k];
        sti[2] = sRCM[2][k];
        sti[3] = sRCM[3][k];
        sti[4] = hyp_epsilon;


        //--------------------------------------------------------------------------------
        //Equivalent state
        //--------------------------------------------------------------------------------
        orbit.getInvman()->evalRCMtoNC(sti, tNCE[k], yvu, OFTS_ORDER, OFS_ORDER);

        //--------------------------------------------------------------------------------
        //Save
        //--------------------------------------------------------------------------------
        #pragma omp critical
        {

            for(int i = 0; i < 6; i++) yNCE[i][k] = yvu[i];
            for(int i = 0; i < 5; i++) sRCM[i][k] = sti[i];
            //Display
            //displayCompletion("cmu_grid_orbit_on_one_period", (double) iter++/(N+1)*100);
        }

        free_dvector(yvu, 0, 5);
        free_dvector(sti, 0, 4);

    }

    //====================================================================================
    //Back to original precision
    //====================================================================================
    //Config::configManager().C_PREC_BACK();


    return status;
}


//========================================================================================
//
//         Refinement of solutions: CMU to CMS - general routine
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EMLi.
 *         The final conditions vary in the paramerization of the CMS of SEMLi.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,
 *         the refSt.coord_type must be NCSEM. However, the user can put
 *         other coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines
 *         are not used.
 **/
int refemlisemli(RefSt& refSt)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    string fname = "refemlisemli";
    cout << "===================================================================" << endl;
    cout << "   refemlisemli. Refinement of EMLi-SEMLj arc                      " << endl;
    cout << "===================================================================" << endl;

    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int coord_type    = refSt.coord_type;
    int man_grid_size = refSt.gridSize;
    int dcs  = default_coordinate_system(coord_type);
    int status = 0;

    //====================================================================================
    // 2. Structures to compute the invariant manifolds
    //====================================================================================
    Invman invman_EM( OFTS_ORDER, OFS_ORDER, SEML.cs_em);
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML.cs_sem);

    //====================================================================================
    // 3. User-defined number of continuation steps, if necessary
    //====================================================================================
    int cont_steps_MAX = 0;
    if(refSt.type == REF_CONT || refSt.type == REF_CONT_D ||
            refSt.type == REF_CONT_D_HARD_CASE ||  refSt.type == REF_COMP)
    {
        if(refSt.isLimUD)
        {
            cout << "refemlisemli. User-defined number of continuation steps" << endl;
            cout << "Enter a value for cont_steps_MAX: ";
            cin >> cont_steps_MAX;
            refSt.cont_step_max    = cont_steps_MAX;

            cout << "refemlisemli. User-defined number of continuation steps" << endl;
            cout << "Enter a value for cont_steps_MAX_vt: ";
            cin >> cont_steps_MAX;
            refSt.cont_step_max_vt = cont_steps_MAX;
        }
    }
    else
    {
        cout << "refemlisemli. No continuation will be performed." << endl;
        refSt.cont_step_max    = 0;
        refSt.cont_step_max_vt = 0;
    }

    //====================================================================================
    // 4. Select the good IC for EML2-to-SEMLi connections in data files
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    status = selectemlisemli(refSt, st_EM, st_SEM, t_EM, &t0_SEM, &pmin_dist_SEM_out);

    if(status == FTC_FAILURE)
    {
        cout << "No solution satisfying all the constraints has been found." << endl;
        cout << "End of computation." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Display the selected first guess
    //====================================================================================
    coutmp();
    cout << "===================================================================" << endl;
    cout << " refemlisemli. A solution has been selected,                       " << endl;
    cout << " with the following characteristics:                               " << endl;
    cout << " 1. Estimated error at patch point (km):                           " << endl;
    cout << " ep = " << pmin_dist_SEM_out* SEML.cs_sem.cr3bp.L  << " km, ";
    cout <<  pmin_dist_SEM_out  << " in SEMSU coord."                             << endl;
    cout << " 2. Initial conditions at EML2 (RCM):                              " << endl;
    cout << " s0 = (" << st_EM[0] << ", " << st_EM[1] << ", ";
    cout              << st_EM[2] << ", " << st_EM[3] << ", " << st_EM[4] << ")"  << endl;
    cout << " 3. Final conditions at SEML2 (RCM):                               " << endl;
    cout << " qf = (" << st_SEM[0] << ", " << st_SEM[1] << ", ";
    cout              << st_SEM[2] << ", " << st_SEM[3] << ")"                    << endl;
    cout << " 4. Estimated TOF (x T):                                           " << endl;
    cout << " TOF = " << (t_EM[1] - t_EM[0])/SEML.us_em.T                         << endl;
    cout << "===================================================================" << endl;
    coutlp();
    pressEnter(refSt.isFlagOn);

    //====================================================================================
    // 5.1 Initialize local variables: EM
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Initialisation of the orbit structure
    //------------------------------------------------------------------------------------
    OdeStruct odestruct_EM;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams_EM(&SEML_EM, I_NCEM);
    //Init ode structure
    init_ode_structure(&odestruct_EM, T, T_root, 6, qbcp_vfn, &odeParams_EM);
    //Init routine
    Orbit orbit_EM(&invman_EM, &SEML_EM, &odestruct_EM, OFTS_ORDER, OFS_ORDER, t_EM[0], t_EM[1]);

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
    OdeStruct odestruct_SEM;
    //Parameters
    OdeParams odeParams_SEM(&SEML_SEM, I_NCSEM);
    //Init ode structure
    init_ode_structure(&odestruct_SEM, T, T_root, 6, qbcp_vfn, &odeParams_SEM);
    //Init routine
    Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &odestruct_SEM, OFTS_ORDER, OFS_ORDER,
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
    switch(refSt.type)
    {
    case REF_CONT:
    case REF_SINGLE:
    {
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init(refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type, refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //Multiple shooting procedure
        //--------------------------------------------------------------------------------
        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        gnuplot_close(refSt.isPlotted, h2);
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
        h2 = gnuplot_init(refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type, refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //First step: variable tn
        //--------------------------------------------------------------------------------
        refSt.time = REF_VAR_TN;

        //There is no need to save at this point,
        //since everything is gonna be erased by the second step
        int isSaved = refSt.isSaved;
        refSt.isSaved = 0;

        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }


        //--------------------------------------------------------------------------------
        //Second step: fixed time + saved, if desired by the user
        //--------------------------------------------------------------------------------
        refSt.time    = REF_FIXED_TIME;
        refSt.isSaved = isSaved;

        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with fixed time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        gnuplot_close(refSt.isPlotted, h2);
        break;
    }

    case REF_CONT_D_HARD_CASE:
    {
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init(refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type, refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //First step: variable tn
        //--------------------------------------------------------------------------------
        refSt.time = REF_VAR_TN;

        //There is no need to save at this point,
        //since everything is gonna be erased by the second step
        int isSaved = refSt.isSaved;
        refSt.isSaved = 0;

        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
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
        refSt.time    = REF_FIXED_TIME;
        refSt.grid    = REF_GIVEN_GRID;
        refSt.isSaved = isSaved;

        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with fixed time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        gnuplot_close(refSt.isPlotted, h2);
        break;
    }


    case REF_COMP:
    {
        cout << "===============================================================" << endl;
        cout << "refemlisemli.                                                 " << endl;
        cout << "The computation of the entire trajectory has been selected.    " << endl;
        cout << "===============================================================" << endl;
        cout << "refemlisemli. First part: "                                     << endl;
        cout << "Refinement of the connection leg in the parameterization space." << endl;
        //--------------------------------------------------------------------------------
        //First part: REF_CONT with variable tn
        //--------------------------------------------------------------------------------
        refSt.type = REF_CONT;
        refSt.time = REF_VAR_TN;

        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init(refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type, refSt.isPlotted);

        //Multiple shooting procedure
        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        //Second part: REF_SINGLE
        //--------------------------------------------------------------------------------
        refSt.type = REF_SINGLE;
        refSt.time = REF_FIXED_TIME;


        //Multiple shooting procedure
        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        gnuplot_close(refSt.isPlotted, h2);

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
        refSt.type = REF_COMP;
        int sf_eml2, sf_man, sf_seml2;
        cout << "===============================================================" << endl;
        cout << "refemlisemli. Second part: "                                    << endl;
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
        status = comprefemlisemli3d(sampfreq, coord_type, orbit_EM, orbit_SEM, refSt, 0, 1);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure on the whole trajectory. ref_errno = " << ref_strerror(status) << endl;
            return FTC_FAILURE;
        }

        //--------------------------------------------------------------------------------
        //CAREFUL!! the rest of the routines are still using msd_grid_size
        //(a number of points) and NOT the sampling frequencies!! To be adapted...
        //--------------------------------------------------------------------------------
        //compref3d_test_seml_synjpl(msd_grid_size, coord_type, orbit_SEM, refSt);
        //compref3d_test_eml_synjpl(msd_grid_size, coord_type, orbit_EM, refSt);
        //comprefft3d_test_eml2seml_synjpl(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refSt);
        //comprefft3d_test_eml2seml_insem(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refSt);
        break;
    }
    }

    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    pressEnter(refSt.isFlagOn);


    return FTC_SUCCESS;
}


/**
 *  \brief Computes the best trajectories from int_proj_ORBIT_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of i.
 *         The final conditions vary in the paramerization of the CMS of SEMLi.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 *
 *         Because of the possible use of msvt3d and msvtplan inside this routine,
 *         the refSt.coord_type must be NCSEM. However, the user can put
 *         other coordinate systems (VNCSEM, PSEM, PEM...), as long as these two routines
 *         are not used.
 *
 *         so stands for single orbit.
 **/
int sorefemlisemli(RefSt& refSt)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    string fname = "sorefemlisemli";
    cout << "===================================================================" << endl;
    cout << "   sorefemlisemli. Refinement of EMLi-SEMLi arc                    " << endl;
    cout << "===================================================================" << endl;

    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int coord_type    = refSt.coord_type;
    int man_grid_size = refSt.gridSize;
    int dcs           = default_coordinate_system(coord_type);
    int status        = 0;

    //====================================================================================
    // 2. Structures to compute the invariant manifolds
    //====================================================================================
    Invman invman_EM( OFTS_ORDER, OFS_ORDER, SEML.cs_em);
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML.cs_sem);


    //====================================================================================
    // 3. Select the good IC for EML2-to-SEMLi connections in data files
    //====================================================================================
    ProjResClass projRes;
    status = soselectemlisemli(refSt, projRes);

    if(status == FTC_FAILURE)
    {
        cout << "No solution satisfying all the constraints has been found." << endl;
        cout << "End of computation." << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // Display
    //====================================================================================
    coutmp();
    cout << "===================================================================" << endl;
    cout << " refemlisemli. " << projRes.size() << " solutions have been found: " << endl;
    //Display
    coutmp();
    cout << "--------------------------------------" << endl;
    cout << "The first entry is:" << endl;
    projRes.displayFirstEntry();

    cout << "--------------------------------------" << endl;
    cout << "The last entry is:" << endl;
    projRes.displayLastEntry();
    cout << "===================================================================" << endl;
    coutlp();
    pressEnter(refSt.isFlagOn);


    //====================================================================================
    // 5. We initialize the local variables with the first instance in projRes
    //====================================================================================
    double st_EM[5], st_SEM[5], t_EM[2], t0_SEM = 0.0, pmin_dist_SEM_out = 0.0;
    double st_EM_seed[4], t_EM_seed = 0.0;
    int label = 0;
    projRes.update_ic(st_EM, st_SEM, t_EM, st_EM_seed, &t_EM_seed, &t0_SEM, &pmin_dist_SEM_out, &label, 0);

    //====================================================================================
    // 5.1 Initialize local variables: EM
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Initialisation of the orbit structure
    //------------------------------------------------------------------------------------
    OdeStruct odestruct_EM;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams_EM(&SEML_EM, I_NCEM);
    //Init ode structure
    init_ode_structure(&odestruct_EM, T, T_root, 6, qbcp_vfn, &odeParams_EM);
    //Init routine
    Orbit orbit_EM(&invman_EM, &SEML_EM, &odestruct_EM, OFTS_ORDER, OFS_ORDER, t_EM[0], t_EM[1]);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ifc(st_EM,  t_EM[0],  t_EM[1]);

    //====================================================================================
    // 5.2 Initialize local variables: SEM
    //====================================================================================

    //------------------------------------------------------------------------------------
    // Initialisation of the orbit structure
    //------------------------------------------------------------------------------------
    OdeStruct odestruct_SEM;
    //Parameters
    OdeParams odeParams_SEM(&SEML_SEM, I_NCSEM);
    //Init ode structure
    init_ode_structure(&odestruct_SEM, T, T_root, 6, qbcp_vfn, &odeParams_SEM);
    //Init routine
    Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &odestruct_SEM, OFTS_ORDER, OFS_ORDER,
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
    // 6 Multiple shooting procedure, for all values in projRes
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h2;
    h2 = gnuplot_init(refSt.isPlotted);

    string filename      = refSt.get_and_update_filename(refSt.FILE_CONT, TYPE_CONT_ATF, ios::out);           //filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, TYPE_CONT_ATF, SEML.li);
    string filename_traj = refSt.get_and_update_filename(refSt.FILE_CONT_TRAJ, TYPE_CONT_ATF_TRAJ, ios::out); //filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, TYPE_CONT_ATF_TRAJ, SEML.li);

    //------------------------------------------------------------------------------------
    // Storage of the type
    //------------------------------------------------------------------------------------
    int refst_type = refSt.type;

    //------------------------------------------------------------------------------------
    // Loop on all first guesses
    //------------------------------------------------------------------------------------
    int isfirst = 1;
    for(int k = 0; k < projRes.size(); k++)
    {
        //--------------------------------------------------------------------------------
        // Update initial conditions
        //--------------------------------------------------------------------------------
        projRes.update_ic(st_EM, st_SEM, t_EM, st_EM_seed, &t_EM_seed, &t0_SEM, &pmin_dist_SEM_out, &label, k);

        orbit_EM.update_ifc(st_EM,  t_EM[0],  t_EM[1]);
        orbit_SEM.update_ic(st_SEM, t0_SEM);

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        notablePoints_sem(h2, coord_type, refSt.isPlotted);

        //--------------------------------------------------------------------------------
        //First step: variable tn, no continuation
        //--------------------------------------------------------------------------------
        refSt.type    = REF_ORBIT;
        refSt.time    = REF_FIXED_TIME;
        refSt.isSaved = 0;

        status = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);

        //--------------------------------------------------------------------------------
        //Second step: variable tn & continuation, if necessary
        //--------------------------------------------------------------------------------
        if(status == FTC_SUCCESS  && refst_type == REF_CONT_ORBIT)
        {
            refSt.type  = REF_CONT_ORBIT;
            refSt.time  = REF_VAR_TN;
            status      = subrefemlisemli(orbit_EM, orbit_SEM, y_traj, t_traj, dcs, coord_type, &man_grid_size, refSt, h2);
        }

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure with variable time. ref_errno = " << ref_strerror(status) << endl;
            //return FTC_FAILURE;
        }
        else if(refSt.last_error < refSt.inner_prec)
        {
            //----------------------------------------------------------------------------
            // Find the intersection with a certain Poincaré Section (PS)
            //----------------------------------------------------------------------------
            double ye[6], te = 0.0;
            xpkemlisemli(ye, &te, t_traj, y_traj, man_grid_size, refSt);

            //----------------------------------------------------------------------------
            // Saving
            //----------------------------------------------------------------------------
            // Only the txt results
            writeCONT_txt(filename, orbit_EM, orbit_SEM, te, ye, projRes, isfirst,  k);
            // Entire trajectory
            writeCONT_bin(refSt, filename_traj, y_traj, t_traj, man_grid_size,
                          orbit_EM, orbit_SEM, isfirst, 0, 0, projRes, k);
            isfirst = 0;
        }
        else
        {
            cout << "PRECISION NOT GOOD ENOUGH" << endl;
        }

        //--------------------------------------------------------------------------------
        // press Enter
        //--------------------------------------------------------------------------------
        pressEnter(refSt.isFlagOn);
    }

    gnuplot_close(refSt.isPlotted, h2);

    return FTC_SUCCESS;
}

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - subroutines
//
//========================================================================================
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
int subrefemlisemli(Orbit& orbit_EM, Orbit& orbit_SEM, double** y_traj, double* t_traj,
                    int dcs, int coord_type, int* man_grid_size_t,
                    RefSt& refSt, gnuplot_ctrl* h2)
{
    //====================================================================================
    // 0. Check on  coord_type (for now)
    //====================================================================================
    if(coord_type !=  NCSEM)
    {
        cout << "subrefemlisemli. WARNING: the coord_type inside refSt is different  " << endl;
        cout << "from NCSEM, which means that the variable-time capability cannot   " << endl;
        cout << "be used. If you still wish to continue, press enter.               " << endl;
        cout << "===================================================================" << endl;
        char ch;
        scanf("%c",&ch);
    }

    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "subrefemlisemli";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = 2;//floor((double)refSt.cont_step_max/floor(0.1*GP_MAX_TMP_FILES))+1;

    //------------------------------------------------------------------------------------
    //If the grid is fixed, then we use the user-provided size. Else, we use a
    //big value (max_grid_size) so that the arrays would not be saturated.
    //------------------------------------------------------------------------------------
    int max_grid_size = 1000;
    int man_grid_size = (refSt.grid == REF_VAR_GRID) ? max_grid_size:*man_grid_size_t;

    //------------------------------------------------------------------------------------
    //Local variables for plotting
    //------------------------------------------------------------------------------------
    int mPlot  = refSt.mplot;
    int isPlot = refSt.isPlotted;
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
    diffcorrptr  diffcorr  = ftc_select_diffcorr(refSt);
    predictorptr predictor = ftc_select_predictor(refSt);

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
    int man_index = icmanemlisemli(y_traj, t_traj, orbit_EM, orbit_SEM, dcs, coord_type, man_grid_size, refSt);

    // Plot the resulting trajectory
    gnuplot_plot_X(isPlot, h2, y_traj, man_index,   (char*)"", "lines", "1", "1", 4);
    gnuplot_plot_X(isPlot, h2, y_traj, man_index+1, (char*)"", "points", "1", "1", 4);


    //====================================================================================
    // 3.  Differential correction
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    pressEnter(refSt.isFlagOn);

    //------------------------------------------------------------------------------------
    //Free variables
    //------------------------------------------------------------------------------------
    int nfv = nfreevariables(refSt, man_index);
    double* nullvector = dvector(0, nfv-1);

    //------------------------------------------------------------------------------------
    //GNUPLOT:
    //  1. if there is a continuation procedure based on variable time,
    //     the component s5 at SEML is plotted.
    //  2. if there is a continuation procedure based on fixed time,
    //     several computations are produced.
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h3 = 0, *h4 = 0, *h5 = 0, *h6 = 0, *h7 = 0;
    if(refSt.isCont())
    {
        switch(refSt.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:

            h5 = gnuplot_init(refSt.isPlotted);
            gnuplot_cmd(refSt.isPlotted, h5,  "set title \"s5_SEM vs steps\" ");
            gnuplot_cmd(refSt.isPlotted, h5, "set grid");
            break;


        case REF_FIXED_TIME:

            h3 = gnuplot_init(refSt.isPlotted);
            gnuplot_cmd(refSt.isPlotted, h3,  "set title \"s3_EM vs s1_EM\" ");
            gnuplot_cmd(refSt.isPlotted, h3, "set grid");

            h4 = gnuplot_init(refSt.isPlotted);
            gnuplot_cmd(refSt.isPlotted, h4,  "set title \"s3_SEM vs s1_SEM\" ");
            gnuplot_cmd(refSt.isPlotted, h4, "set grid");

            h5 = gnuplot_init(refSt.isPlotted);
            gnuplot_cmd(refSt.isPlotted, h5,  "set title \"t0_SEM vs steps\" ");
            gnuplot_cmd(refSt.isPlotted, h5, "set grid");

            if(refSt.type == REF_3D)
            {
                h6 = gnuplot_init(refSt.isPlotted);
                gnuplot_cmd(refSt.isPlotted, h6,  "set title \"s4_EM vs s2_EM\" ");
                gnuplot_cmd(refSt.isPlotted, h6, "set grid");
            }

            if(refSt.is3D())
            {
                h7 = gnuplot_init(refSt.isPlotted);
                gnuplot_cmd(refSt.isPlotted, h7,  "set title \"s4_SEM vs s2_SEM\" ");
                gnuplot_cmd(refSt.isPlotted, h7, "set grid");
            }

            break;


        }
    }


    //====================================================================================
    //6.1. First step of the continuation procedure.
    //====================================================================================
    int status = 0;
    int niter = 1;

    status = diffcorr(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                      man_index, coord_type, true, orbit_EM, orbit_SEM,
                      h2, refSt, &niter);




    //====================================================================================
    // If it is not a success, we print the value, and we return
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Error during the first refinement procedure. ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }

    //====================================================================================
    // If it is a sucess, we go on.
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Some inner variables
    //------------------------------------------------------------------------------------
    double yv[6];
    int kn = 0;

    //====================================================================================
    // Find the intersection with a certain Poincaré Section (PS)
    //====================================================================================
    double ye[6], te = 0.0;
    //int newpos;
    xpkemlisemli(ye, &te, t_traj_n, y_traj_n, man_index, refSt);
    //xpkemlisemli(ye, &te, t_traj_n, y_traj_n, &newpos, man_index, refSt);

    //====================================================================================
    // Save first entry if the continuation process is on
    //====================================================================================
    string filename = "", filename_traj = "";

    if(refSt.isSaved && refSt.isCont())
    {
        filename      = refSt.get_and_update_filename(refSt.FILE_CONT, TYPE_CONT_ATF, ios::out);
        filename_traj = refSt.get_and_update_filename(refSt.FILE_CONT_TRAJ, TYPE_CONT_ATF_TRAJ, ios::out);

        cout << "-----------------------------------------------------------"  << endl;
        cout << fname << ". First refinement was a success.                 "  << endl;
        cout << " We can save the first entry in txt and binary files:      "  << endl;
        cout << "-" << filename                                                << endl;
        cout << "-" << filename_traj                                           << endl;
        cout << "-----------------------------------------------------------"  << endl;

        // Main parameters
        writeCONT_txt(kn, filename, orbit_EM, orbit_SEM, te, ye, true);

        // Entire trajectory
        writeCONT_bin(refSt, filename_traj, dcs, coord_type,
                      y_traj_n, t_traj_n, man_index, mPlot,
                      orbit_EM, orbit_SEM, kn++, true,
                      refSt.isSaved_EM, refSt.isSaved_SEM);
    }

    //====================================================================================
    //Continuation procedure
    //====================================================================================
    if(refSt.isCont())
    {
        cout << fname << ". Continuation procedure.                         "  << endl;
        cout << "-----------------------------------------------------------"  << endl;

        //================================================================================
        // Initialization
        //================================================================================

        //--------------------------------------------------------------------------------
        // Initialize the stepper & desired number of iterations in Newton's method
        //--------------------------------------------------------------------------------
        double ds0    = 5e-2;
        double niterd = 4;
        double dsmin  = refSt.dsmin, dsmax = refSt.dsmax;

        switch(refSt.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:
            ds0    = refSt.ds0_vt;
            niterd = refSt.nu0_vt;
            dsmin  = refSt.dsmin_vt;
            dsmax  = refSt.dsmax_vt;
            break;
        case REF_FIXED_TIME:
            ds0    = refSt.ds0;
            niterd = refSt.nu0;
            dsmin  = refSt.dsmin;
            dsmax  = refSt.dsmax;
            break;
        }

        refSt.dsc  = ds0;
        double ds  = ds0;
        double dkn = 0.0;

        //--------------------------------------------------------------------------------
        //An additional condition can be set:
        //
        // REF_COND_S5: if the time is not fixed, we want to decrease the
        //s5 component at SEMLi. Therefore, we add a stopping condition on the loop.
        //
        // REF_COND_T: we want to make enough "turns" about SEMLi
        //--------------------------------------------------------------------------------
        bool addCondition = true;
        double theta_acc  = 0.0, theta_old = 0.0, theta_max  = refSt.thetaMax;
        int mani_old = man_index;
        switch(refSt.termination)
        {
        case REF_COND_S5:
        {
            switch(refSt.time)
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

        //================================================================================
        // Continuation loop
        //================================================================================

        //--------------------------------------------------------------------------------
        // The maximum number of steps depends on the strategy
        //--------------------------------------------------------------------------------
        int cont_step_max_local = refSt.cont_step_max;
        switch(refSt.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:
            cont_step_max_local = refSt.cont_step_max_vt;
            break;


        case REF_FIXED_TIME:
            cont_step_max_local = refSt.cont_step_max;
            break;
        }

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        while(kn < cont_step_max_local && status == GSL_SUCCESS && addCondition)
        {
            //============================================================================
            // Updating the free variables via the predictor
            //============================================================================
            predictor(y_traj_n, t_traj_n, &ds, ds, nullvector,
                      orbit_EM, orbit_SEM, man_index, coord_type, refSt);

            //============================================================================
            // Diff correction
            //============================================================================
            status = diffcorr(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                              nullvector, 42, man_index, coord_type,
                              false, orbit_EM, orbit_SEM, h2, refSt, &niter);

            //============================================================================
            // Find the intersection with a certain Poincaré Section (PS)
            //============================================================================
            double ye[6], te = 0.0;
            if(status == GSL_SUCCESS) xpkemlisemli(ye, &te, t_traj_n, y_traj_n, man_index, refSt);

            //============================================================================
            // Save
            //============================================================================
            if(refSt.isSaved && status == GSL_SUCCESS)
            {
                //Update the filenames, if necessary
                if(filename.empty())      filename      = refSt.get_and_update_filename(refSt.FILE_CONT, TYPE_CONT_ATF, ios::out);
                if(filename_traj.empty()) filename_traj = refSt.get_and_update_filename(refSt.FILE_CONT_TRAJ, TYPE_CONT_ATF_TRAJ, ios::out);

                // Main parameters
                writeCONT_txt(kn, filename, orbit_EM, orbit_SEM, te, ye, false);

                // Entire trajectory
                writeCONT_bin(refSt, filename_traj, dcs, coord_type,
                              y_traj_n, t_traj_n, man_index, mPlot,
                              orbit_EM, orbit_SEM, kn, false,
                              refSt.isSaved_EM, refSt.isSaved_SEM);
            }


            //============================================================================
            // Display
            //============================================================================
            dkn = (double) kn;
            double t0 = t_traj_n[0]/SEML.us_sem.T;//t_traj_n[0]/SEML.us_sem.T;
            if(refSt.isCont() && status == GSL_SUCCESS && kn % plotfreq == 0)
            {
                switch(refSt.time)
                {
                case REF_VAR_TIME:
                case REF_VAR_TN:
                {
                    gnuplot_plotc_xy(isPlot, h5, &dkn, &orbit_SEM.getSi()[4], 1, (char*)"", "points", "1", "2", 0);
                    break;
                }

                case REF_FIXED_TIME:
                {
                    gnuplot_plotc_xy(isPlot, h3, &orbit_EM.getSi()[0],  &orbit_EM.getSi()[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plotc_xy(isPlot, h4, &orbit_SEM.getSi()[0],  &orbit_SEM.getSi()[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plotc_xy(isPlot, h5, &dkn, &t0, 1, (char*)"", "points", "1", "2", 0);

                    if(refSt.type == REF_3D)
                    {
                        gnuplot_plotc_xy(isPlot, h6, &orbit_EM.getSi()[1],  &orbit_EM.getSi()[3], 1, (char*)"", "points", "1", "2", 0);
                    }

                    if(refSt.is3D())
                    {
                        gnuplot_plotc_xy(isPlot, h7, &orbit_SEM.getSi()[1],  &orbit_SEM.getSi()[3], 1, (char*)"", "points", "1", "2", 0);
                    }

                    break;
                }

                }
            }

            //============================================================================
            // Update the additionnal condition, if necessary
            //============================================================================
            switch(refSt.termination)
            {
            case REF_COND_S5:
            {
                //------------------------------------------------------------------------
                // First type of condition: stop when we are close enough to the
                // center manifold (unstable component is small enough)
                // note: this condition is not used when all the times are fixed.
                //------------------------------------------------------------------------
                switch(refSt.time)
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
                //------------------------------------------------------------------------
                // Another possible condition: enough turns around SEMLi.
                //------------------------------------------------------------------------
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
                    //--------------------------------------------------------------------
                    // Each time that theta_old*theta < 0, it means that we have made
                    // about 180° around the last point. We then do 2 actions:
                    //  1. We update the accumulated angle (+90°)
                    //  2. We add a few patch points, to ensure numerical stability.
                    //--------------------------------------------------------------------
                    if(theta_old*theta < 0)
                    {
                        //  1. We update the accumulated angle (+90°)
                        theta_acc += 90;
                    }

                    theta_old = theta;
                }

                //------------------------------------------------------------------------
                // The condition is only applied when some times are left free to vary
                //------------------------------------------------------------------------
                switch(refSt.time)
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


            //============================================================================
            // Advance
            //============================================================================
            kn++;
            cout << "Step n°" << kn << "/" << cont_step_max_local << " complete"  << endl;
            cout << "----------------------------------------------------------"  << endl;


            //============================================================================
            // Update the stepper. For now, only in the REF_VAR_TIME case !
            //============================================================================
            cout << "Updating the stepper:                                     "  << endl;
            cout << "ds_old = " << ds                                             << endl;

            //Prior to updating ds, we ensure that, at least, niter = 1
            niter = max(1, niter);

            switch(refSt.time)
            {
            case REF_VAR_TIME:
            case REF_VAR_TN:
                ds = min(dsmax, max(dsmin, ds*niterd/niter));
                break;
            case REF_FIXED_TIME:
                ds = min(dsmax, max(dsmin, ds*niterd/niter));
                break;
            }

            refSt.dsc  = ds;
            cout << "niter  = " << niter << endl;
            cout << "ds_new = " << ds                                             << endl;
            cout << "----------------------------------------------------------"  << endl;


        }

        //================================================================================
        // If the trajectory has been extended... We add a few points.
        // Careful with this piece of code: if nnew is big, then a big number of points
        // is added to the state vectors, which can make the Newton procedure unstable.
        //
        // Consequently, this option is taken only if refSt.type = REF_CONT_D_HARD_CASE
        // i.e. basically if the user has tried refSt.type = REF_CONT_D before
        // but it did not work...
        //
        // IMPORTANT: after this step, the number of points man_grid_size_t is changed!
        // this means that any further use of this routine with man_grid_size_t will
        // take into account this change.
        //
        //================================================================================
        if(theta_acc >= theta_max && refSt.type == REF_CONT_D_HARD_CASE)
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
            gnuplot_plot_X(isPlot, h2, y_traj_n, man_index+1, (char*)"", "points", "1", "1", 8);
        }
    }

    //====================================================================================
    // Final correction with better precision
    // CAREFUL: it is NOT done for now because some elements depends on man_index for their size,
    // and man_index CAN move if REF_CONT_D_HARD_CASE is used !!
    // Need to make a check on who is using man_index...
    //====================================================================================
    /*
    //    inner_prec = 1e-10;                    //harder precision in diff corr procedures
    //    Config::configManager().C_PREC_HARD(); //harder precision in numerical integration
    //    switch(refSt.dim)
    //    {
    //    case REF_3D:
    //        switch(refSt.time)
    //        {
    //        case REF_FIXED_TIME:
    //            status = msft3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                            nullvector, 42,
    //                            man_index, coord_type, inner_prec, false,
    //                            orbit_EM, orbit_SEM,
    //                            h2, refSt, &niter);
    //            break;
    //
    //        case REF_VAR_TIME:
    //        case REF_VAR_TN:
    //            status = msvt3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                            nullvector, 42,
    //                            man_index, coord_type, inner_prec, false,
    //                            orbit_EM, orbit_SEM,
    //                            h2, refSt, &niter);
    //            break;
    //        }
    //        break;
    //    case REF_PLANAR:
    //        switch(refSt.time)
    //        {
    //        case REF_FIXED_TIME:
    //            status = msftplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                              nullvector, 42,
    //                              man_index, coord_type, inner_prec, false,
    //                              orbit_EM, orbit_SEM,
    //                              h2, refSt, &niter);
    //            break;
    //
    //        case REF_VAR_TIME:
    //            status = msvtplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                              nullvector, 42,
    //                              man_index, coord_type, inner_prec, false,
    //                              orbit_EM, orbit_SEM,
    //                              h2, refSt, &niter);
    //            break;
    //
    //        case REF_VAR_TN:
    //            status = msvltplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
    //                               nullvector, 42,
    //                               man_index, coord_type, inner_prec, false,
    //                               orbit_EM, orbit_SEM,
    //                               h2, refSt, &niter);
    //            break;
    //        }
    //        break;
    //    }
    //
    //    Config::configManager().C_PREC_BACK(); //back to initial precision
    */

    //====================================================================================
    // Save final solution, only if REF_CONT. The savings of the orbit is FORCED.
    //====================================================================================
    if(refSt.type == REF_CONT)
    {
        cout << "-----------------------------------------------------------"  << endl;
        cout << fname << ".Since REF_CONT is used, the final orbit is saved."  << endl;
        cout << " In txt and binary files:                                  "  << endl;

        //Update the filenames, if necessary
        if(filename.empty())      filename      = refSt.get_and_update_filename(refSt.FILE_CONT, TYPE_CONT_ATF, ios::out);
        if(filename_traj.empty()) filename_traj = refSt.get_and_update_filename(refSt.FILE_CONT_TRAJ, TYPE_CONT_ATF_TRAJ, ios::out);

        cout << "-" << filename                                                << endl;
        cout << "-" << filename_traj                                           << endl;
        cout << "-----------------------------------------------------------"  << endl;

        // If !refSt.isSaved, this is the first solution to be solved, hence label = 0
        double label = (refSt.isSaved)? kn+1:0;

        // Intersection with the Pk section
        xpkemlisemli(ye, &te, t_traj_n, y_traj_n, man_index, refSt);

        // Main parameters
        writeCONT_txt(label, filename, orbit_EM, orbit_SEM, te, ye, true);

        // Entire trajectory
        writeCONT_bin(refSt, filename_traj, dcs, coord_type,
                      y_traj_n, t_traj_n, man_index, mPlot,
                      orbit_EM, orbit_SEM, label, !refSt.isSaved, 1, 1);
    }


    //====================================================================================
    // Display final solution
    //====================================================================================
    if(status)
    {
        cerr << fname << ". Error during the continuation procedure. ref_errno = " << ref_strerror(status) << endl;
        return FTC_FAILURE;
    }
    else if(refSt.isPlotted)
    {
        //================================================================================
        // Final trajectory, on a grid
        //================================================================================
        cout << fname << ". Final trajectory, on a grid.                    "  << endl;
        cout << "-----------------------------------------------------------"  << endl;


        //--------------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------------
        double** ymc   = dmatrix(0, 5, 0, mPlot);
        double* tmc    = dvector(0, mPlot);
        double yv[6];
        int ode78coll;

        //--------------------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //--------------------------------------------------------------------------------
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
                gnuplot_plot_X(isPlot, h2, ymc, mPlot+1, (char*)"Final.", "lines", "1", "2", 0);
            else
                gnuplot_plot_X(isPlot, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 0);
        }

        //--------------------------------------------------------------------------------
        // Free
        //--------------------------------------------------------------------------------
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
    }


    //====================================================================================
    // 9. Update the orbits for next step
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
        //orbit_EM.setSi(PROJ_EPSILON, 4);
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

    cout << fname << ". "                                                  << endl;
    cout << "          End of computation.                              "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    pressEnter(refSt.isFlagOn);

    //====================================================================================
    // 10. Reset pk section
    //====================================================================================
    refSt.pkpos = 0;


    //====================================================================================
    // 11. Free
    //====================================================================================
    if(refSt.isCont())
    {
        switch(refSt.time)
        {
        case REF_VAR_TIME:
        case REF_VAR_TN:

            gnuplot_close(refSt.isPlotted, h5);
            break;

        case REF_FIXED_TIME:

            gnuplot_close(refSt.isPlotted, h3);
            gnuplot_close(refSt.isPlotted, h4);
            gnuplot_close(refSt.isPlotted, h5);

            if(refSt.type == REF_3D)
            {
                gnuplot_close(refSt.isPlotted, h6);
            }

            if(refSt.is3D())
            {
                gnuplot_close(refSt.isPlotted, h7);
            }

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
 *         through data files produced by int_proj_CMU_EM_on_CM_SEM_3D or
 *         int_proj_CMU_EM_on_CM_SEM.
 **/
int selectemlisemli(RefSt& refSt, double st_EM[5], double st_SEM[5], double t_EM[2],
                    double* t0_SEM, double* pmin_dist_SEM_out)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    cout << "-------------------------------------------------------------------" << endl;
    cout << "   selectemlisemli. Selects good IC for EML2-SEMLi arc"                << endl;
    cout << "-------------------------------------------------------------------" << endl;
    //====================================================================================
    // 1. Init the data containers
    //====================================================================================
    //------------------------------------------------------------------------------------
    //To store data from the data file
    //------------------------------------------------------------------------------------
    ProjResClass sortSt;

    //------------------------------------------------------------------------------------
    //To store data from a sub selection of the previous elements
    //------------------------------------------------------------------------------------
    ProjResClass subSt;

    //====================================================================================
    // 2. Getting back the data
    //====================================================================================
    string filename;
    if(refSt.isFromServer)
    {
        //Getting the desired t0 as a fraction of T
        double t0_des_mod = fmod(refSt.t0xT_des, 1.0);

        //Getting the name of the data file, depending on the parameters in refSt
        switch(refSt.dim)
        {

        case REF_PLANAR:
        {
            if(SEML.li_EM == 2)
            {
                if(SEML.li_SEM == 2) filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_dest_L2_tspan_";
                else filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_dest_L1_tspan_";

                if(refSt.crossings > 0)
                {
                    filename += "0T_T_crossings.bin";
                    //filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_dest_L2_t0_0065T.bin";
                }
                else
                {
                    if(t0_des_mod >= 0 && t0_des_mod < 0.25)
                    {
                        filename += "0T_025T_FINAL.bin";
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
                }

                //filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_t0_099T.bin";
            }
            else
            {
                if(SEML.li_SEM == 2) filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_dest_L2_t0_0.bin";
                else filename = SEML.cs->F_PLOT+"Serv/projcu_order_20_dest_L1.bin";
            }

            break;
        }

        case REF_3D:
        case REF_MIXED:
        {

            if(SEML.li_SEM == 2) filename = SEML.cs->F_PLOT+"Serv/projcu_3d_order_20_dest_L2_t0_0995T.bin";
            else filename = SEML.cs->F_PLOT+"Serv/projcu_3d_order_20_dest_L1_t0_0995T.bin";

            break;
        }

        }


        cout << "selectemlisemli. The data will be retrieved from the server-computed file named:" << endl;
        cout << filename << endl;

        //Read data file
        sortSt.readProjRes_t0(filename, refSt.t0_des, refSt.typeOfTimeSelection);
    }
    else
    {
        int type = TYPE_MAN_PROJ;
        switch(refSt.dim)
        {
        case REF_PLANAR:
            type = TYPE_MAN_PROJ;
            break;

        case REF_3D:
        case REF_MIXED:
            type = TYPE_MAN_PROJ_3D;
            break;
        }

        filename = refSt.get_and_update_filename(refSt.FILE_PCU, type, ios::in); //filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, type, SEML.li_SEM);

        //Read data file
        sortSt.readProjRes_t0(filename, refSt.t0_des, refSt.typeOfTimeSelection);
    }



    //====================================================================================
    // Select the subselection
    //====================================================================================
    bool flag = subSt.push_back_conditional(sortSt, refSt);

    //====================================================================================
    // If there at least one solution that satisfies all the criteria
    //====================================================================================
    if(flag)
    {
        //--------------------------------------------------------------------------------
        //Sort the subselection by pmin distance
        //--------------------------------------------------------------------------------
        subSt.sort_pmin_dist_SEM();

        //================================================================================
        // 5. The best solution in the subselection will serve as the first guess.
        //================================================================================
        subSt.update_ic(st_EM, st_SEM, t_EM, t0_SEM, pmin_dist_SEM_out, 0);
        return FTC_SUCCESS;
    }
    else
    {
        //================================================================================
        // Else, we return a failure
        //================================================================================
        return FTC_FAILURE;

    }
}


/**
 *  \brief Selects good initial conditions for EML2-to-SEMLi connections, searching
 *         through data files produced by int_proj_ORBIT_EM_on_CM_SEM
 **/
int soselectemlisemli(RefSt& refSt, ProjResClass& subSt)
{
    //====================================================================================
    // 0. Splash screen
    //====================================================================================
    cout << "-------------------------------------------------------------------" << endl;
    cout << "   selectemlisemli. Selects good IC for EMLi-SEMLi arc"              << endl;
    cout << "-------------------------------------------------------------------" << endl;

    //====================================================================================
    // 1. Init the data containers
    //====================================================================================
    //------------------------------------------------------------------------------------
    //To store data from the data file
    //------------------------------------------------------------------------------------
    ProjResClass readSt;

    //====================================================================================
    // 2. Getting back the data
    //====================================================================================
    string filename;
    if(refSt.isFromServer)
    {
        /// @TODO: remains to be adapted when the time comes

        int type = TYPE_MAN_PROJ_ORBIT;
        switch(refSt.dim)
        {
        case REF_PLANAR:
            type = TYPE_MAN_PROJ_ORBIT;
            break;

        case REF_3D:
        case REF_MIXED:
            type = TYPE_MAN_PROJ_3D;
            break;
        }

        filename = refSt.get_and_update_filename(refSt.FILE_PCU, type, ios::in); //filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, type, SEML.li_SEM);

        //Read data file
        readSt.readProjRes(filename);
    }
    else
    {
        int type = TYPE_MAN_PROJ_ORBIT;
        switch(refSt.dim)
        {
        case REF_PLANAR:
            type = TYPE_MAN_PROJ_ORBIT;
            break;

        case REF_3D:
        case REF_MIXED:
            type = TYPE_MAN_PROJ_3D;
            break;
        }

        filename = refSt.get_and_update_filename(refSt.FILE_PCU, type, ios::in); //filenameCUM(SEML.cs->F_PLOT, OFTS_ORDER, type, SEML.li_SEM);

        //Read data file
        readSt.readProjRes(filename);
    }

    //------------------------------------------------------------------------------------
    // Subselection
    //------------------------------------------------------------------------------------
    bool flag = subSt.push_back_conditional(readSt, refSt);

    if(!flag)
    {
        cout << "soselectemlisemli. Warning: empty subselection." << endl;
        pressEnter(true);
        return FTC_FAILURE;
    }


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
int icmanemlisemli(double** y_traj, double* t_traj,
                   Orbit& orbit_EM, Orbit& orbit_SEM,
                   int dcs, int coord_type, int man_grid_size,
                   RefSt& refSt)
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
    cout << " icmanemlisemli. Compute the manifold leg..."  << endl;
    //------------------------------------------------------------------------------------
    // Integration: from orbit_EM.t0 to orbit_EM.tf
    // Note that the collisionner is not used at this step. It is used later, after
    // the refinement procedures
    //------------------------------------------------------------------------------------
    int man_index = man_grid_size;
    int ode78coll;

    switch(refSt.grid)
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
    cout << " icmanemlisemli. Compute the first point of the final SEMLi orbit..."<< endl;

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
 *         is computed. The grid size is returned.
 **/
int iccompemlisemli(double** y_traj, double* t_traj,
                    double** y_traj_comp, double* t_traj_comp,
                    Orbit& orbit_EM, Orbit& orbit_SEM,
                    int dcs, int coord_type, int grid_points_des[3],
                    int grid_points_eff[3], int max_grid,
                    RefSt& refSt, gnuplot_ctrl* h2, gnuplot_ctrl* h3)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "iccompemlisemli";

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
    double tof_eml_EM   = refSt.tspan_EM;    //TOF on EML2 orbit
    double tof_seml_SEM = refSt.tspan_SEM;   //TOF on SEMLi orbit

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 2. Compute the initial orbit
    //====================================================================================
    cout << fname << ". Compute the initial orbit..."  << endl;
    //------------------------------------------------------------------------------------
    //For the computation of the initial orbit, we save & kill the unstable part.
    //------------------------------------------------------------------------------------
    double hyp_epsilon = orbit_EM.getSi()[4];
    orbit_EM.setSi(0.0, 4);

    //------------------------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //------------------------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //------------------------------------------------------------------------------------
    //Integration on des_grid_size+1 fixed grid
    //------------------------------------------------------------------------------------
    int em_index = grid_points_eff[0];
    switch(refSt.grid)
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
    gnuplot_plot_X(refSt.isPlotted, h2, y_man_coord, em_index+1, (char*)"", "points", "5", "3", 5);
    gnuplot_plot_X(refSt.isPlotted, h3, y_man_comp, em_index+1, (char*)"", "points", "5", "3", 5);

    //====================================================================================
    // 4.3 Compute the manifold leg
    //====================================================================================
    cout << fname << ". Compute the manifold leg..."  << endl;

    double tf_EM = orbit_SEM.getT0()/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //------------------------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //------------------------------------------------------------------------------------
    orbit_EM.setSi(hyp_epsilon, 4);

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

    switch(refSt.grid)
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
    gnuplot_plot_X(refSt.isPlotted, h2, y_man_coord, man_index+1, (char*)"Man", "points", "1", "3", 4);
    gnuplot_plot_X(refSt.isPlotted, h3, y_man_comp, man_index+1, (char*)"", "points", "1", "3", 4);

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

    switch(refSt.grid)
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
    gnuplot_plot_X(refSt.isPlotted, h2, y_man_coord, sem_index+1, (char*)"SEMLi", "points", "1", "3", 6);
    gnuplot_plot_X(refSt.isPlotted, h3, y_man_comp, sem_index+1, (char*)"", "points", "1", "3", 6);

    //====================================================================================
    // Entire size is: no more than sum(grid_points_des) points or so.
    //====================================================================================
    int final_index = 0;
    if(refSt.grid == REF_VAR_GRID)
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

//----------------------------------------------------------------------------------------
//         Brick C: Find the intersection of a EML2-SEMLi connection
//                     with a certain Pk section x = cst
//----------------------------------------------------------------------------------------
/**
 *  \brief Find the intersection of a EML2-SEMLi connection contained in
 *         y_traj_n/t_traj_n with a certain Pk section x = cst defined by refSt.
 **/
int xpkemlisemli(double ye[6], double* te, double* t_traj_n, double** y_traj_n,
                 int man_index, RefSt& refSt)
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
        test = (y_traj_n[0][k2] > refSt.xps && y_traj_n[0][k1] < refSt.xps) ||
               (y_traj_n[0][k2] < refSt.xps && y_traj_n[0][k1] > refSt.xps);
    }
    while(!test && k1 >= 0);


    //====================================================================================
    //2. If a solution has been found:
    //====================================================================================
    if(test && k1 >= 0)
    {
        //cout << "xpkemlisemli. Range accross the PS found:" << endl;
        //cout << y_traj_n[0][k1] << " < " << refSt.xps << " < " << y_traj_n[0][k2] << endl;

        //================================================================================
        //2. Integrate until x = xps
        //================================================================================
        double center[3];
        struct value_params val_par;
        val_par.max_events = 1;
        val_par.direction  = 0;
        val_par.dim        = 0;
        val_par.value      = refSt.xps;
        val_par.center     = center;
        val_par.type       = 'X';

        double** ye_NCSEM = dmatrix(0, 5, 0, 1);
        double* te_NCSEM  = dvector(0, 1);

        //================================================================================
        //3. After this step: ye_NCSEM[*][0] & te_NCSEM[0] contains the intersection
        //   Note: the (possible) collisions with the primaries are not taken into account here:
        //   We suppose that for this particular application (intersection with a Pk section
        //   far away from any primary), such collision are very unlikely!
        //================================================================================
        int ode78coll;
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k1];
        ode78_qbcp_event(ye_NCSEM, te_NCSEM, &ode78coll, t_traj_n[k1], t_traj_n[k2], yv, 6, dcs,
                         NCSEM, NCSEM, &val_par);

        //================================================================================
        //4. Store
        //================================================================================
        *te = te_NCSEM[0];
        for(int i = 0; i < 6; i++) ye[i] = ye_NCSEM[i][0];

        //================================================================================
        //5. Free
        //================================================================================
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
 *  \brief Find the intersection of a EML2-SEMLi connection contained in
 *         y_traj_n/t_traj_n with a certain Pk section x = cst defined by refSt.
 *         Once the intersection is found, it is incorporated in the sequence of patch points,
 *         in place of a given point, at position newpos
 **/
int xpkemlisemli(double ye[6], double* te, double* t_traj_n, double** y_traj_n,
                 int* newpos, int man_index, RefSt& refSt)
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
        test = (y_traj_n[0][k2] > refSt.xps && y_traj_n[0][k1] < refSt.xps) ||
               (y_traj_n[0][k2] < refSt.xps && y_traj_n[0][k1] > refSt.xps);
    }
    while(!test && k1 >= 0);


    //====================================================================================
    //2. If a solution has been found:
    //====================================================================================
    if(test && k1 >= 0)
    {
        //cout << "xpkemlisemli. Range accross the PS found:" << endl;
        //cout << y_traj_n[0][k1] << " < " << refSt.xps << " < " << y_traj_n[0][k2] << endl;


        //================================================================================
        //2. Integrate until x = xps
        //================================================================================
        double center[3];
        struct value_params val_par;
        val_par.max_events = 1;
        val_par.direction  = 0;
        val_par.dim        = 0;
        val_par.value      = refSt.xps;
        val_par.center     = center;
        val_par.type       = 'X';

        double** ye_NCSEM = dmatrix(0, 5, 0, 1);
        double* te_NCSEM  = dvector(0, 1);

        //================================================================================
        //3. After this step: ye_NCSEM[*][0] & te_NCSEM[0] contains the intersection
        //   Note: the (possible) collisions with the primaries are not taken into account here:
        //   We suppose that for this particular application (intersection with a Pk section
        //   far away from any primary), such collision are very unlikely!
        //================================================================================
        int ode78coll;
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k1];
        ode78_qbcp_event(ye_NCSEM, te_NCSEM, &ode78coll, t_traj_n[k1], t_traj_n[k2], yv, 6, dcs,
                         NCSEM, NCSEM, &val_par);

        //================================================================================
        //4. Store
        //================================================================================
        *te = te_NCSEM[0];
        for(int i = 0; i < 6; i++) ye[i] = ye_NCSEM[i][0];

        //================================================================================
        //5. Switch in the sequence of patch points
        //================================================================================
        *newpos = k2;
        refSt.pkpos = k2;
        for(int i = 0; i < 6; i++) y_traj_n[i][*newpos] = ye[i];
        t_traj_n[*newpos] = *te;

        //================================================================================
        //6. Free
        //================================================================================
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
//         Refining trajectories from continuation procedures
//
//========================================================================================
/**
 *  \brief Refine complete EMLi-SEMLi connections from semi-analyticalcontinuation results
 *         The orbits at both ends are included in the refinement process. Each trajectory
 *         is then refined in a higher-fidelity model (JPL DE430).
 **/
int reffromcontemlisemli(RefSt& refSt)
{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int coord_type = refSt.coord_type;
    int status = 0;

    //====================================================================================
    //         Initialization of the containers
    //====================================================================================
    //Name
    string fname = "reffromcontemlisemli";

    // Get size of the data file
    string filename = refSt.get_and_update_filename(refSt.FILE_CONT, TYPE_CONT_ATF, ios::in);
    int fsize = getLengthCONT_txt(filename);

    // Init the vectors
    double* t0_CMU_EM     = dvector(0, fsize-1);
    double* tf_CMU_EM     = dvector(0, fsize-1);
    double** si_CMU_EM    = dmatrix(0, 4, 0, fsize-1);
    double** si_CMS_SEM   = dmatrix(0, 4, 0, fsize-1);
    double** z0_CMU_NCEM  = dmatrix(0, 5, 0, fsize-1);
    double** z0_CMS_NCSEM = dmatrix(0, 5, 0, fsize-1);
    double* thetae        = dvector(0, fsize-1);
    double** ze_NCSEM     = dmatrix(0, 5, 0, fsize-1);
    double* H0_NCEM       = dvector(0, fsize-1);
    double* He_NCEM       = dvector(0, fsize-1);
    double* H0_NCSEM      = dvector(0, fsize-1);
    double* He_NCSEM      = dvector(0, fsize-1);

    //====================================================================================
    //         Read
    //====================================================================================
    status = readCONT_txt(t0_CMU_EM, tf_CMU_EM, si_CMU_EM, si_CMS_SEM,
                          z0_CMU_NCEM, z0_CMS_NCSEM, thetae, ze_NCSEM, H0_NCEM, He_NCEM,
                          H0_NCSEM, He_NCSEM, fsize, filename);

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
    Invman invman_EM(OFTS_ORDER, OFS_ORDER, *SEML.cs);
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, *SEML_SEM.cs);

    //====================================================================================
    // Select the parameters
    //====================================================================================
    int isFirst = 1;
    for(int k = 0; k < fsize; k++)
    {
        cout << "##################################" << endl;
        cout << "k = " << k  << "/" << fsize         << endl;
        cout << "##################################" << endl;
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
        OdeStruct odestruct_EM;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        OdeParams odeParams_EM(&SEML_EM, I_NCEM);
        //Init ode structure
        init_ode_structure(&odestruct_EM, T, T_root, 6, qbcp_vfn, &odeParams_EM);
        //Init routine
        Orbit orbit_EM(&invman_EM, &SEML_EM, &odestruct_EM, OFTS_ORDER, OFS_ORDER, t_EM[0], t_EM[1]);

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
        OdeStruct odestruct_SEM;
        //Parameters
        OdeParams odeParams_SEM(&SEML_SEM, I_NCSEM);
        //Init ode structure
        init_ode_structure(&odestruct_SEM, T, T_root, 6, qbcp_vfn, &odeParams_SEM);
        //Init routine
        Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &odestruct_SEM, OFTS_ORDER, OFS_ORDER,
                        t0_SEM, t0_SEM+5*SEML.us_sem.T);

        //--------------------------------------------------------------------------------
        // Update the initial state in the orbit_SEM, with the RCM coordinates
        //--------------------------------------------------------------------------------
        orbit_SEM.update_ic(st_SEM, t0_SEM);

        //--------------------------------------------------------------------------------
        // The max projection error is greatly widen, in order to accept bigger orbits
        //--------------------------------------------------------------------------------
        orbit_SEM.setEPmaxx(1e-1);

        //--------------------------------------------------------------------------------
        // We can advance to REF_COMP for the rest of the computation
        //--------------------------------------------------------------------------------
        refSt.type = REF_COMP;
        int sampfreq[3] = {refSt.sf_eml2, refSt.sf_man, refSt.sf_seml2};
        status = comprefemlisemli3d(sampfreq, coord_type, orbit_EM, orbit_SEM, refSt, k, isFirst);

        if(status)
        {
            cerr << fname << ". Error during the continuation procedure on the whole trajectory. ref_errno = " << ref_strerror(status) << endl;
            //return FTC_FAILURE;
        }

        if(isFirst) isFirst = 0;

    }

    return FTC_SUCCESS;
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
int comprefemlisemli3d(int grid_freq_days[3], int coord_type,
                       Orbit& orbit_EM, Orbit& orbit_SEM,
                       RefSt& refSt, int label, int isFirst)

{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "comprefemlisemli3d";
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
    h2 = gnuplot_init(refSt.isPlotted);
    h3 = gnuplot_init(refSt.isPlotted);

    gnuplot_cmd(refSt.isPlotted, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(refSt.isPlotted, h2, "set grid");
    gnuplot_cmd(refSt.isPlotted, h3, "set grid");

    //------------------------------------------------------------------------------------
    //Notable points in SEM & EM systems
    //------------------------------------------------------------------------------------
    notablePoints(h2, h3, coord_type, refSt.isPlotted);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Initialize the number of points on the grids. If the number of points is fixed,
    // the solution is straightforward. If the number of points is desired variable,
    // the default number are set to an arbitrarily high number of points.
    //------------------------------------------------------------------------------------
    cout << " comprefemlisemli3d. Initialize the number of points on the grids..."  << endl;
    int max_grid    = 50000;
    int grid_points_des[3];

    //Desired number of points
    grid_points_des[0] = refSt.tspan_EM/(2*M_PI*86400*grid_freq_days[0]/SEML.cs_em.cr3bp.T);
    grid_points_des[1] = (orbit_SEM.getT0()/SEML.us_em.ns - orbit_EM.getT0())/(2*M_PI*86400*grid_freq_days[1]/SEML.cs_em.cr3bp.T);
    grid_points_des[2] = refSt.tspan_SEM/(2*M_PI*86400*grid_freq_days[2]/SEML.cs_sem.cr3bp.T);

    cout << "Desired frequency on leg 1: "  << grid_freq_days[0]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[0] << endl    << endl;
    cout << "Desired frequency on leg 2: "  << grid_freq_days[1]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[1] << endl    << endl;
    cout << "Desired frequency on leg 3: "  << grid_freq_days[2]  << "days " << endl;
    cout << "==> the number of points  = "  << grid_points_des[2] << endl    << endl;


    //Effective number of points for the following computation: equal to grid_points_des
    //if the number of points if fixed, equal to max_grid otherwise
    int grid_points_eff[3];
    grid_points_eff[0]  = (refSt.grid == REF_VAR_GRID) ? max_grid:grid_points_des[0];
    grid_points_eff[1]  = (refSt.grid == REF_VAR_GRID) ? max_grid:grid_points_des[1];
    grid_points_eff[2]  = (refSt.grid == REF_VAR_GRID) ? max_grid:grid_points_des[2];

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
    int final_index = iccompemlisemli(y_traj, t_traj, y_traj_comp, t_traj_comp,
                                      orbit_EM, orbit_SEM, dcs, coord_type,
                                      grid_points_des,
                                      grid_points_eff,
                                      max_grid,
                                      refSt, h2, h3);


    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    cout << " comprefemlisemli3d. Initial trajectory, on a grid..."  << endl;
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
    pressEnter(refSt.isFlagOn);

    int isPlotted   = 0;
    status = multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, PREC_GSM, isPlotted, h2, 0);

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
    if(refSt.isPlotted)
    {
        cout << fname << ". Final trajectory, on a grid..."  << endl;
        //Final trajectory on lines, segment by segment
        plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type, 0.0, 0.0,
                         coord_type, h2, comp_type, h3, 3, "Final guess");
    }

    //====================================================================================
    // 6. Store the data, in sys units
    //====================================================================================
    string filename = refSt.get_and_update_filename(refSt.FILE_JPL_TXT, TYPE_COMP_FOR_JPL, ios::out);
    cout << fname << ". Storing the refined solution in system units "     << endl;
    cout << " The solution will be stored in " << filename                 << endl;
    cout << "-----------------------------------------------------------"  << endl;
    writeCOMP_txt(t_traj_n, y_traj_n, final_index, filename);


    //====================================================================================
    // 7. JPL refinement
    //====================================================================================
    if(refSt.isJPL)
    {
        status = jplref3d(coord_type, refSt, label, isFirst, filename);

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
    pressEnter(refSt.isFlagOn);
    gnuplot_close(refSt.isPlotted, h2);
    gnuplot_close(refSt.isPlotted, h3);


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
 *  \brief Refine a given output of comprefemlisemli3d into JPL ephemerides.
 **/
int jplref3d(int coord_type, RefSt& refSt, int label, int isFirst, string filename_in)
{
    //====================================================================================
    // 1. Read the data, in sys units
    //====================================================================================
    string fname = "jplref3d";
    cout << " jplref3d. Read the data, in sys units..."  << endl;
    int status = 0;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    int final_index = getLengthCOMP_txt(filename_in);

    if(final_index == FTC_ENOENT)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    status = readCOMP_txt(t_traj_n, y_traj_n, final_index, filename_in);

    if(status != FTC_SUCCESS)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

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
    int mPlot = refSt.mplot;
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

    //Filename for saving
    string filename = refSt.get_and_update_filename(refSt.FILE_JPL_BIN, TYPE_CONT_JPL_TRAJ, ios::out);

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init(refSt.isPlotted);
    h3 = gnuplot_init(refSt.isPlotted);

    gnuplot_cmd(refSt.isPlotted, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(refSt.isPlotted, h2, "set view equal xyz");
    //gnuplot_cmd(refSt.isPlotted, h3, "set view equal xyz");
    gnuplot_cmd(refSt.isPlotted, h2, "set grid");
    gnuplot_cmd(refSt.isPlotted, h3, "set grid");

    //------------------------------------------------------------------------------------
    //Notable points in SEM & EM systems
    //------------------------------------------------------------------------------------
    notablePoints(h2, h3, coord_type, refSt.isPlotted);

    //====================================================================================
    // 5. Initial trajectory, on a grid
    //====================================================================================
    cout << " jplref3d. Initial trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type, 0.0, 0.0,
                     coord_type, h2, comp_type, h3, 3, "Initial trajectory");


    if(refSt.isSaved)
    {
        cout << "The initial trajectory is saved in " << filename << endl;
        savetrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type, 0.0, 0.0, coord_type, comp_type, filename, label, isFirst);
    }

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
    pressEnter(refSt.isFlagOn);

    //====================================================================================
    // Select between VECLI, J2000, and NJ2000 (best choice is NJ2000)
    //====================================================================================
    int choice = 0;
    //If refSt.djplcoord does not contain a valid choice (e.g. -1), we ask the user
    if(refSt.djplcoord != VECLI && refSt.djplcoord != J2000  && refSt.djplcoord != NJ2000)
    {
        cout << "Enter " <<  VECLI << "for VECLI, ";
        cout << J2000  << "for J2000, " << NJ2000 << "for NJ2000: ";
        cin >> choice;
    }
    else
    {
        choice = refSt.djplcoord;
    }

    //------------------------------------------------------------------------------------
    // Select the coordinate system and the frawemork for future integration
    //------------------------------------------------------------------------------------
    int coord_int = 0, fwrk0 = SEML.fwrk;// fwrk_int = 0;
    switch(choice)
    {
    case VECLI:
        coord_int = VECLI;
        //fwrk_int  = I_ECLI;
        break;
    case J2000:
        coord_int = J2000;
        //fwrk_int  = I_J2000;
        break;
    case NJ2000:
        coord_int = NJ2000;
        //fwrk_int  = I_NJ2000;

        //Change focus, if necessary
        if(fwrk0 != fwrk) changeDCS(SEML, fwrk);
        cout << "NJ2000 has been chosen: changeDCS is applied to match " << endl;
        cout << " the framework associated to coord_type.";
        break;
    default:
        cout << "wrong input. NJ2000 is used by default." << endl;
        cout << "NJ2000 has been chosen: changeDCS is applied to match " << endl;
        cout << " the framework associated to coord_type.";
        coord_int = NJ2000;
        //fwrk_int  = I_NJ2000;
        if(fwrk0 != fwrk) changeDCS(SEML, fwrk);
        break;
    }

    //====================================================================================
    // Initialize SPICE kernerls & VF
    //====================================================================================
    gnuplot_ctrl* h4 = gnuplot_init(refSt.isPlotted);
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
    //    final_index = jplfg3d_switch(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index,
    //                                   coord_type, comp_type, coord_int,
    //                                   et0, tsys0, tsys0_comp);


    //------------------------------------------------------------------------------------
    //Same as jplfg3d_switch, with a refinement of the position of the minimum.
    //------------------------------------------------------------------------------------
    //    int mRef = 50;
    //    final_index = jplfg3d_super_switch(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index,
    //                                   coord_type, comp_type, coord_int, mRef,
    //                                   et0, tsys0, tsys0_comp);


    //------------------------------------------------------------------------------------
    //Same as jplfg3d_switch, with an interpolation before & after the switching point.
    //------------------------------------------------------------------------------------
    int mRef = 5;
    final_index = jplfg3d_interpolation(y_traj_n, t_traj_n, &y_traj_jpl, &t_traj_jpl, final_index,
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
    //double** y_moon_spice  = dmatrix(0, 5, 0, final_index);
    double** y_l2_spice    = dmatrix(0, 5, 0, final_index);

    ymc_v       = dmatrix(0, 5, 0, mPlot*final_index);
    tmc_v       = dvector(0, mPlot*final_index);
    ymc_comp_v  = dmatrix(0, 5, 0, mPlot*final_index);
    tmc_comp_v  = dvector(0, mPlot*final_index);


    //------------------------------------------------------------------------------------
    //Time in seconds: @TODO: adapt the other cases, only NJ2000 for now!!
    //------------------------------------------------------------------------------------
    for(int p = 0; p <= final_index; p++)
    {
        et_traj_jpl[p] = t_traj_jpl[p]/SEML.ss->n;
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

        //L2
        spkez_c (392, et, DEFFRAME,  "NONE", DEFOBSINT, YV, &lt);
        for(int i = 0; i < 6; i++) y_l2_spice[i][p] = YV[i];

        //MOON - The Moon does NOT work for now, for an unknown reason... Linked to SPICE kernels?
        //spkezr_c ("MOON", et, DEFFRAME,  "NONE", DEFOBS, YV, &lt);
        //spkez_c (301, et, DEFFRAME,  "NONE", DEFOBSINT, YV, &lt);
        //for(int i = 0; i < 6; i++) y_moon_spice[i][p] = YV[i];
    }

    //------------------------------------------------------------------------------------
    //Position of the Sun +  Earth + Initial guess in coord_type coordinates
    //------------------------------------------------------------------------------------
    ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    gnuplot_plot_X(refSt.isPlotted, h2, y_jpl_temp, final_index+1, (char*) "EARTH", "points", "3", "2", 8);

    //ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    //gnuplot_plot_X(refSt.isPlotted, h2, y_jpl_temp, final_index+1, (char*) "MOON", "points", "4", "2", 8);

    ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
    gnuplot_plot_X(refSt.isPlotted, h2, y_jpl_temp, final_index+1, (char*) "SEMLi", "points", "5", "2", 8);

    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    //    pressEnter(refSt.isFlagOn);
    //    cout << fname << ". Plotting Initial Guess..." << endl;
    //
    //    //------------------------------------------------------------------------------
    //    //Initial trajectory on lines, segment by segment
    //    //------------------------------------------------------------------------------
    //    plottrajsegbyseg(y_traj_jpl, t_traj_jpl, final_index, mPlot, coord_int,
    //                     et0, tsys0, coord_type, h2, comp_type,  h3, 6,
    //                     "Initial guess in JPL ephemerides");


    //====================================================================================
    // 6.5 Differential correction
    //====================================================================================
    pressEnter(refSt.isFlagOn);
    cout << fname << ". Refine trajectory in JPL ephemerides, with fixed times..." << endl;
    status  = multiple_shooting_direct(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int, PREC_GSM, true, h4, 1);

    //------------------------------------------------------------------------------------
    // If something went bad, we try the same but with variable times
    //------------------------------------------------------------------------------------
    if(status)
    {
        cout << fname << ". Refine trajectory in JPL ephemerides, with variables times..." << endl;
        status  = multiple_shooting_direct_variable_time(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int,  PREC_GSM, true, h4, 1);

    }

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
    //Plot
    //------------------------------------------------------------------------------------
    gnuplot_plot_X(refSt.isPlotted, h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


    //====================================================================================
    //Final trajectory on lines, segment by segment
    //====================================================================================
    plottrajsegbyseg(y_traj_jpl_n, t_traj_jpl_n, final_index, mPlot, coord_int,
                     et0, tsys0, coord_type, h2, comp_type,  h3, 2,
                     "Final guess in JPL ephemerides");

    if(refSt.isSaved)
    {
        cout << "The final trajectory is saved in " << filename << endl;
        savetrajsegbyseg(y_traj_jpl_n, t_traj_jpl_n, final_index, mPlot, coord_int, et0, tsys0, coord_type, comp_type, filename, label, 0);
    }


    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    pressEnter(refSt.isFlagOn);
    gnuplot_close(refSt.isPlotted, h2);
    gnuplot_close(refSt.isPlotted, h3);
    gnuplot_close(refSt.isPlotted, h4);

    //====================================================================================
    // Reset the focus in SEML, if necessary
    //====================================================================================
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);

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
 *  \brief Refine a given output of comprefemlisemli3d into Inertial Coordinates, then into
 *         JPL coordinates.
 **/
int comptojplref3d(int coord_type, RefSt& refSt, string filename_in)
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
    int final_index = getLengthCOMP_txt(filename_in);

    if(final_index == FTC_ENOENT)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    status = readCOMP_txt(t_traj_n, y_traj_n, final_index, filename_in);

    if(status != FTC_SUCCESS)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

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

    h2 = gnuplot_init(refSt.isPlotted);
    h3 = gnuplot_init(refSt.isPlotted);
    h4 = gnuplot_init(refSt.isPlotted);

    gnuplot_cmd(refSt.isPlotted, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h3, "set title \"Complementary coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h4, "set title \"Continuation steps\" ");

    gnuplot_cmd(refSt.isPlotted, h2, "set grid");
    gnuplot_cmd(refSt.isPlotted, h3, "set grid");
    gnuplot_cmd(refSt.isPlotted, h4, "set grid");

    //Notable points in SEM & EM systems
    notablePoints(h2, h3, coord_type, refSt.isPlotted);

    //Color for plots
    int color = 1;

    //====================================================================================
    // 5. Initial trajectory, on a grid
    //====================================================================================
    cout << " jplref3d. Initial trajectory, on a grid..."  << endl;
    //Initial trajectory on lines, segment by segment
    plottrajsegbyseg(y_traj_n, t_traj_n, final_index, mPlot, coord_type,
                     0.0, 0.0, coord_type, h2, comp_type,  h3,color++,
                     "Initial guess (NCSEM)");

    //====================================================================================
    // 5. INERTIAL refinement
    //====================================================================================
    //pressEnter(refSt.isFlagOn);

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
    //pressEnter(refSt.isFlagOn);
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
    //pressEnter(refSt.isFlagOn);
    cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;
    //DiffCorr
    status  = multiple_shooting_direct(y_traj_ecisem, t_traj_ecisem, y_traj_ecisem, t_traj_ecisem, 42, final_index, coord_int, PREC_GSM, true, h4, 0);
    //Plot
    gnuplot_plot_X(refSt.isPlotted, h4, y_traj_ecisem, final_index+1, (char*) "Final trajectory", "lines", "2", "2", color);

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
    SEML.ss->tshift = et0*SEML.ss->n - tsys0;

    //====================================================================================
    //Compute the equivalent of the state along the trajectory in NJ2000 coordinates, to
    //compare with ECI QBCP coordinates!
    //====================================================================================
    // ECISEM -> NJ2000
    coord2necistate_vec(y_traj_n, t_traj_n, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0,  tsys0, eph_coord(coord_int), *SEML.ss);
    //Plot
    gnuplot_plot_X(refSt.isPlotted, h4, y_jpl_temp, final_index+1, (char*)"Inertial state (NJ2000)", "lines", "1", "1", color);

    //Moon's position (NJ2000)
    SpiceDouble lt, RS1[6], RS2[6], RE[6], REM[6];
    spkezr_c ("MOON",  et0, DEFFRAME,  "NONE",  DEFOBS, RS1, &lt);
    spkezr_c ("EARTH",  et0, DEFFRAME,  "NONE",  DEFOBS, RE, &lt);
    spkezr_c ("EARTH MOON BARYCENTER", et0, DEFFRAME,  "NONE",  DEFOBS, REM, &lt);
    ecl2neci(RS1, REM, RS2, *SEML.ss);
    RS2[2] = 0.0;
    gnuplot_plot_xyz(refSt.isPlotted, h4, RS2, RS2+1, RS2+2, 1, (char*)"Moon (NJ2000, projected)", "points", "3", "5", color);

    //Moon's position (Custom)
    double me = SEML.us_sem.me;
    double mm = SEML.us_sem.mm;
    double n  = SEML.us_sem.n;
    double ni = SEML.us_sem.ni;
    double ai = SEML.us_sem.ai;
    double r1 = creal(evz(SEML.cs_sem.zt, tsys0, n, ni, ai));
    double r2 = cimag(evz(SEML.cs_sem.zt, tsys0, n, ni, ai));
    double Pm[3] = {- me/(mm + me)* r1, - me/(mm + me)* r2, 0.0};
    gnuplot_plot_xyz(refSt.isPlotted, h4, Pm, Pm+1, Pm+2, 1, (char*)"Moon (ECISEM)", "points", "4", "5", color++);


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
    pressEnter(refSt.isFlagOn);


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
    h1 = gnuplot_init(refSt.isPlotted);


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
        for(int k = 0; k <= final_index; k++) t_traj_ecisem_s[k] = t_traj_ecisem_c[k] + SEML.ss->tshift;

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
        gnuplot_plot_xy(refSt.isPlotted, h1, &y_traj_ecisem_c[0][0], &(SEML.epsilon), 1, "", "points", "7", "3", 1);


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
    pressEnter(refSt.isFlagOn);
    gnuplot_close(refSt.isPlotted, h2);
    gnuplot_close(refSt.isPlotted, h3);
    gnuplot_close(refSt.isPlotted, h4);

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
int jplfg3d_switch(double** y_traj_n, double* t_traj_n,
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
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), *SEML.ss);
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
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), *SEML.ss);
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
    cout << "jplfg3d_switch. pmin = " << pmin << "/" << final_index << endl;

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
 * \brief Same as jplfg3d_switch, with a refinement of the position of the minimum.
 **/
int jplfg3d_super_switch(double** y_traj_n, double* t_traj_n,
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
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), *SEML.ss);
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
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), *SEML.ss);
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
    cout << "jplfg3d_super_switch. pmin = " << pmin << "/" << final_index << endl;


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
        coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), *SEML.ss);
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
        coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type), *SEML.ss);
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
        coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), *SEML.ss);
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
        coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type), *SEML.ss);
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
    cout << "jplfg3d_super_switch. Switch part before refinment       " << endl;
    cout << "-----------------------------------------------------------" << endl;
    for(int p = pmin-1; p <= pmin+1; p++)
    {
        cout << t_traj_jpl[p] << "  " << y_traj_jpl[0][p] << "   " << y_traj_jpl[1][p] <<  "   " << y_traj_jpl[2][p] << endl;
    }

    //The best fit replace the old point, both in position in time
    for(int i = 0; i <6; i++) y_traj_jpl_n[i][pmin] = ymin[i];
    t_traj_jpl[pmin] = tmin;

    cout << "-----------------------------------------------------------" << endl;
    cout << "jplfg3d_super_switch. Switch part after refinment        " << endl;
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
 * \brief Same as jplfg3d_switch, with an interpolation before and after the switching
 *        point.
 **/
int jplfg3d_interpolation(double** y_traj_n, double* t_traj_n,
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
        coord2necistate_vec(y_traj_n, t_traj_n, *y_traj_jpl, *t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type), *SEML.ss);
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
        coord2necistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type), *SEML.ss);
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
    double dmin = 0.0;//, dminmin = 0.0;
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
            //dminmin = dmin;
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

    cout << "jplfg3d_interpolation. pmin = " << pmin << "/" << final_index << endl;

    //    cout << "-----------------------------------------------------------" << endl;
    //    cout << "jplfg3d_interpolation. Switch part before refinment      " << endl;
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
            coord2necistate_vec(yms, tms, ymc, tmc, mRef, coord_type, et0, tsys0, eph_coord(coord_type), *SEML.ss);
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
            coord2necistate_vec(yms, tms, ymc_comp, tmc_comp, mRef, coord_type, et0, tsys0_comp, eph_coord(comp_type),*SEML.ss);
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
    //    cout << "jplfg3d_interpolation. Switch part after refinement      " << endl;
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
int compref3d_test_seml_synjpl(int man_grid_size_t,
                               int coord_type,
                               Orbit& orbit_SEM,
                               RefSt refSt)

{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "compref3d_test_seml_synjpl";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_seml_SEM = refSt.tspan_SEM;   //TOF on SEMLi orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init(refSt.isPlotted);
    h3 = gnuplot_init(refSt.isPlotted);

    gnuplot_cmd(refSt.isPlotted, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(refSt.isPlotted, h2, "set view equal xyz");
    //gnuplot_cmd(refSt.isPlotted, h3, "set view equal xyz");
    gnuplot_cmd(refSt.isPlotted, h2, "set grid");
    gnuplot_cmd(refSt.isPlotted, h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type, refSt.isPlotted);

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    int max_grid    = 3000;
    int man_grid_2  = (refSt.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;
    int traj_grid_2 = (refSt.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;

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

    switch(refSt.grid)
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
    gnuplot_plot_X(refSt.isPlotted, h2, y_man_coord, sem_index+1, (char*)"SEMLi", "points", "1", "3", 6);
    gnuplot_plot_X(refSt.isPlotted, h3, y_man_comp, sem_index+1, (char*)"", "points", "1", "3", 6);

    //====================================================================================
    // Entire size is: no more than 3*man_grid_size_t points or so.
    //====================================================================================
    int final_index = 0;
    if(refSt.grid == REF_VAR_GRID)
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
    gnuplot_plot_X(refSt.isPlotted, h2, ymc_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);
    //Plot on h2: points
    gnuplot_plot_X(refSt.isPlotted, h2, y_traj, final_index+1, (char*)"", "points", "1", "2", 7);
    //Plot on h3: lines
    gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);


    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    pressEnter(refSt.isFlagOn);

    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, PREC_GSM, isPlotted, h2, 0);

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
    gnuplot_plot_X(refSt.isPlotted, h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(refSt.isPlotted, h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

    //====================================================================================
    // JPL
    //====================================================================================
    if(refSt.isJPL)
    {
        //----------------------------------------------------------
        //Go on
        //----------------------------------------------------------
        pressEnter(refSt.isFlagOn);

        //====================================================================================
        // Initialize SPICE kernerls
        //====================================================================================
        gnuplot_ctrl* h4 = gnuplot_init(refSt.isPlotted);
        cout << fname << ". Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //====================================================================================
        // Initialize VF
        //====================================================================================
        cout << fname << ". Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct odestruct_JPL;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        OdeParams odeParams(&SEML, dcs);
        //Init ode structure
        init_ode_structure(&odestruct_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

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
        SEML.ss->et0 = et0;
        SEML.ss->t0  = tsys0;

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
        pressEnter(refSt.isFlagOn);
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
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);
        }

        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        pressEnter(refSt.isFlagOn);
        cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), PREC_GSM, isPlotted, h4, 0);
        //Plot
        gnuplot_plot_X(refSt.isPlotted, h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


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
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
        }

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        pressEnter(refSt.isFlagOn);
        gnuplot_close(refSt.isPlotted, h4);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    pressEnter(refSt.isFlagOn);
    gnuplot_close(refSt.isPlotted, h2);
    gnuplot_close(refSt.isPlotted, h3);


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
int compref3d_test_eml_synjpl(int man_grid_size_t,
                              int coord_type,
                              Orbit& orbit_EM,
                              RefSt refSt)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    string fname = "compref3d_test_eml_synjpl";
    cout << fname << ". Initialization of the local variables..."  << endl;

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_eml_EM = refSt.tspan_EM;   //TOF on EML2 orbit
    //    double tof_seml_SEM = refSt.tspan_SEM;   //TOF on SEMLi orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h2, *h3;
    h2 = gnuplot_init(refSt.isPlotted);
    h3 = gnuplot_init(refSt.isPlotted);

    gnuplot_cmd(refSt.isPlotted, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(refSt.isPlotted, h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(refSt.isPlotted, h2, "set view equal xyz");
    //gnuplot_cmd(refSt.isPlotted, h3, "set view equal xyz");
    gnuplot_cmd(refSt.isPlotted, h2, "set grid");
    gnuplot_cmd(refSt.isPlotted, h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type, refSt.isPlotted);

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
    int comp_type = comp_coord_typ(coord_type);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    int max_grid    = 3000;
    int man_grid_2  = (refSt.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;
    int traj_grid_2 = (refSt.grid == REF_VAR_GRID) ? max_grid:man_grid_size_t;

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

    switch(refSt.grid)
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
    gnuplot_plot_X(refSt.isPlotted, h2, y_man_coord, em_index+1, (char*)"", "lines", "5", "3", 5);
    gnuplot_plot_X(refSt.isPlotted, h3, y_man_comp, em_index+1, (char*)"", "lines", "5", "3", 5);


    //====================================================================================
    // Entire size is: no more than 3*man_grid_size_t points or so.
    //====================================================================================
    int final_index = 0;
    if(refSt.grid == REF_VAR_GRID)
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
    gnuplot_plot_X(refSt.isPlotted, h2, ymc_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);
    //Plot on h2: points
    gnuplot_plot_X(refSt.isPlotted, h2, y_traj, final_index+1, (char*)"", "points", "1", "2", 7);
    //Plot on h3: lines
    gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp_v, mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);


    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    cout << fname << ". Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, PREC_GSM, h2, 0);

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
    gnuplot_plot_X(refSt.isPlotted, h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(refSt.isPlotted, h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

    //====================================================================================
    // JPL
    //====================================================================================
    if(refSt.isJPL)
    {
        //----------------------------------------------------------
        //Go on
        //----------------------------------------------------------
        printf("Press ENTER to go on with the JPL ref");
        scanf("%c",&ch);


        //====================================================================================
        // Initialize SPICE kernerls
        //====================================================================================
        gnuplot_ctrl* h4 = gnuplot_init(refSt.isPlotted);
        cout << fname << ". Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //====================================================================================
        // Initialize VF
        //====================================================================================
        cout << fname << ". Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct odestruct_JPL;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Parameters
        OdeParams odeParams(&SEML, dcs);
        //Init ode structure
        init_ode_structure(&odestruct_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

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
        SEML.ss->et0 = et0;
        SEML.ss->t0  = tsys0;

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
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);

        }

        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), PREC_GSM, isPlotted, h4, 0);
        //Plot
        gnuplot_plot_X(refSt.isPlotted, h4, y_traj_jpl_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


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
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(refSt.isPlotted, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_X(refSt.isPlotted, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
        }

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
        gnuplot_close(refSt.isPlotted, h4);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(refSt.isPlotted, h2);
    gnuplot_close(refSt.isPlotted, h3);


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
int compref3d_test_eml2seml_synjpl(int coord_type, string filename_in)
{
    //====================================================================================
    // 1. Read the data, in sys units
    //====================================================================================
    string fname = "compref3d_test_eml2seml_synjpl";
    cout << fname << ". Read the data, in sys units..."  << endl;

    //------------------------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //------------------------------------------------------------------------------------
    int final_index = getLengthCOMP_txt(filename_in);

    if(final_index == FTC_ENOENT)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

    double** y_traj_n   = dmatrix(0, 41, 0, final_index);
    double* t_traj_n    = dvector(0, final_index);

    //------------------------------------------------------------------------------------
    //Read the data stored in a data file
    //------------------------------------------------------------------------------------
    int status = readCOMP_txt(t_traj_n, y_traj_n, final_index, filename_in);

    if(status != FTC_SUCCESS)
    {
        cout << fname << ". Impossible to read the data file." << endl;
        return FTC_FAILURE;
    }

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
    h2 = gnuplot_init(true);
    h3 = gnuplot_init(true);
    gnuplot_cmd(true, h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(true, h3, "set title \"Complementary coordinates\" ");
    //gnuplot_cmd(true, h2, "set view equal xyz");
    //gnuplot_cmd(true, h3, "set view equal xyz");
    gnuplot_cmd(true, h2, "set grid");
    gnuplot_cmd(true, h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM & EM systems
    //----------------------------------------------------------
    notablePoints(h2, h3, coord_type, true);

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
    gnuplot_plot_X(true, h2, ymc_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_X(true, h2, y_traj_n, final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_X(true, h3, ymc_comp_v, mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

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
    gnuplot_ctrl* h4 = gnuplot_init(true);
    cout << fname << ". Initialize SPICE kernerls..." << endl;
    furnsh_c("spice/kernels/metakernel.furnsh");

    //====================================================================================
    // Initialize VF
    //====================================================================================
    cout << fname << ". Initialize associated VF..." << endl;
    int shift = 0;
    OdeStruct odestruct_JPL;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct_JPL, T, T_root, 6, jpl_vf_syn, &odeParams);

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
    SEML.ss->et0 = et0;
    SEML.ss->t0  = tsys0;

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
        if(k == 0) gnuplot_plot_X(true, h2, ymc, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
        else gnuplot_plot_X(true, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 1);

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
        if(k == 0) gnuplot_plot_X(true, h3, ymc_comp, mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
        else gnuplot_plot_X(true, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 1);
    }

    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    printf("Press ENTER to refine\n");
    scanf("%c",&ch);
    cout << fname << ". Refine trajectory in JPL ephemerides..." << endl;

    int isPlotted   = 1;
    multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_syn_n, t_traj_syn_n, 42, final_index, eph_coord(coord_type), PREC_GSM, isPlotted, h4, 0);
    //Plot
    gnuplot_plot_X(true, h4, y_traj_syn_n, final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


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
        if(k == 0) gnuplot_plot_X(true, h2, ymc, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
        else gnuplot_plot_X(true, h2, ymc, mPlot+1, (char*)"", "lines", "1", "2", 2);

        //Back to comp_type coordinates
        //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
        syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
        ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

        //Plot on h3
        if(k == 0) gnuplot_plot_X(true, h3, ymc_comp, mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
        else gnuplot_plot_X(true, h3, ymc_comp, mPlot+1, (char*)"", "lines", "1", "2", 2);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(true, h4);


    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(true, h2);
    gnuplot_close(true, h3);


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
//         PLOTTING TRAJECTORIES
//
//========================================================================================
/**
 *  \brief Plotting the notable points in the SEM system, in coord_type coordinates
 **/
int notablePoints_sem(gnuplot_ctrl* h2, int coord_type, int isPlot)
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
        semPlot(h2, semP_coord, isPlot);
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord, isPlot);
        break;
    case PEM:
    case VEM:
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord, isPlot);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord, isPlot);
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
int notablePoints(gnuplot_ctrl* h2, gnuplot_ctrl* h3, int coord_type, int isPlot)
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
        semPlot(h2, semP_coord, isPlot);
        //Comp
        emPoints(0.0, semP_comp);
        emPlot(h3, semP_comp, isPlot);

        //h2 displays the system in SEM coordinates
        gnuplot_set_xlabel(isPlot, h2, (char*) "Xsem");
        gnuplot_set_ylabel(isPlot, h2, (char*) "Ysem");
        gnuplot_set_zlabel(isPlot, h2, (char*) "Zsem");

        //h3 displays the system in EM coordinates
        gnuplot_set_xlabel(isPlot, h3, (char*) "Xem");
        gnuplot_set_ylabel(isPlot, h3, (char*) "Yem");
        gnuplot_set_zlabel(isPlot, h3, (char*) "Zem");
        break;
    }
    case NCSEM:
    case VNCSEM:
    {
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord, isPlot);
        //Comp
        emNCPoints(0.0, semP_comp);
        emPlot(h3, semP_comp, isPlot);

        //h2 displays the system in NCSEM coordinates
        gnuplot_set_xlabel(isPlot, h2, (char*) "xsem");
        gnuplot_set_ylabel(isPlot, h2, (char*) "ysem");
        gnuplot_set_zlabel(isPlot, h2, (char*) "zsem");

        //h3 displays the system in NCEM coordinates
        gnuplot_set_xlabel(isPlot, h3, (char*) "xem");
        gnuplot_set_ylabel(isPlot, h3, (char*) "yem");
        gnuplot_set_zlabel(isPlot, h3, (char*) "zem");
        break;
    }
    case PEM:
    case VEM:
    {
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord, isPlot);
        //Comp
        semPoints(0.0, semP_comp);
        semPlot(h3, semP_comp, isPlot);

        //h3 displays the system in SEM coordinates
        gnuplot_set_xlabel(isPlot, h3, (char*) "Xsem");
        gnuplot_set_ylabel(isPlot, h3, (char*) "Ysem");
        gnuplot_set_zlabel(isPlot, h3, (char*) "Zsem");

        //h2 displays the system in EM coordinates
        gnuplot_set_xlabel(isPlot, h2, (char*) "Xem");
        gnuplot_set_ylabel(isPlot, h2, (char*) "Yem");
        gnuplot_set_zlabel(isPlot, h2, (char*) "Zem");
        break;
    }
    case NCEM:
    case VNCEM:
    {
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord, isPlot);
        //Comp
        semNCPoints(0.0, semP_comp);
        semPlot(h3, semP_comp, isPlot);

        //h3 displays the system in NCSEM coordinates
        gnuplot_set_xlabel(isPlot, h3, (char*) "xsem");
        gnuplot_set_ylabel(isPlot, h3, (char*) "ysem");
        gnuplot_set_zlabel(isPlot, h3, (char*) "zsem");

        //h2 displays the system in NCEM coordinates
        gnuplot_set_xlabel(isPlot, h2, (char*) "xem");
        gnuplot_set_ylabel(isPlot, h2, (char*) "yem");
        gnuplot_set_zlabel(isPlot, h2, (char*) "zem");
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
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys1, et0, tsys0, eph_coord(coordsys1), *SEML.ss);
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
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys2, et0, tsys0_comp, eph_coord(coordsys2), *SEML.ss);
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
    gnuplot_plot_X(true, h2, ymc_v, mPlot*final_index+1, (char*)title.c_str(), "lines", "1", "2", color);

    //Plot on h3
    gnuplot_plot_X(true, h3, ymc_comp_v, mPlot*final_index+1, (char*)title.c_str(), "lines", "1", "2", color);

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


/**
 *  \brief Save a trajectory, segment by segment, in to complementary coordinate systems.
 **/
int savetrajsegbyseg(double** y_traj, double* t_traj,
                     int final_index, int mPlot,
                     int coord_int,
                     double et0,      double tsys0,
                     int coordsys1,   int coordsys2,
                     string filename, int label, bool isFirst)
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

    double* tmc_original_v  = dvector(0, mPlot*final_index);

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
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys1, et0, tsys0, eph_coord(coordsys1), *SEML.ss);
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
            tmc_original_v[k*mPlot + p] = tmc_comp[p];
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
            neci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coordsys2, et0, tsys0_comp, eph_coord(coordsys2), *SEML.ss);
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
    //Save mPlot*final_index+1 points
    //------------------------------------------------------------------------------------
    fstream filestream;
    double res;
    if(isFirst) filestream.open (filename.c_str(), ios::out | ios::binary);
    else filestream.open (filename.c_str(), ios::out | ios::binary | ios::app);

    for(int k = 0; k < mPlot*final_index+1; k++)
    {
        //1. Label of the solution
        res = label;
        filestream.write((char*) &res, sizeof(double));

        //2. coord_int
        res = coord_int;
        filestream.write((char*) &res, sizeof(double));

        //3. Time in original units
        res = tmc_original_v[k];
        filestream.write((char*) &res, sizeof(double));

        //3. Time in new units (e.g. NCSEM)
        res = tmc_v[k];
        filestream.write((char*) &res, sizeof(double));

        //4. Current state in new coordinates (e.g. NCSEM)
        for(int i = 0; i < 6; i++)
        {
            res  = ymc_v[i][k];
            filestream.write((char*) &res, sizeof(double));
        }

        //5. Time in comp units (e.g. NCEM)
        res = tmc_comp_v[k];
        filestream.write((char*) &res, sizeof(double));

        //6. Current state in comp coordinates (e.g. NCEM)
        for(int i = 0; i < 6; i++)
        {
            res  = ymc_comp_v[i][k];
            filestream.write((char*) &res, sizeof(double));
        }

    }
    filestream.close();

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
//          I/O (misc)
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
