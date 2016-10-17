#include "oolencon.h"


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
 *  \param dist_to_cm:     the value in RCM coordinates on the unstable direction s5.
 *  \param tlim_CMU_EM:    the min/max starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param si_LIM_CMU_RCM: the min/max of s1, s2, s3, s4 values (in RCM coordinates)
 *                         in the IC box.
 *  \param si_grid_size:   the number of points on the  s1, s2, s3, s4 values  grids
 *                         in the IC box.
 *  \param invman:         the center-unstable manifold, if the type of manifold provided
 *                         is not center-unstable, a warning message is displayed and
 *                         nothing is done.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *                   "plot/QBCP/EM/L2/cu_3d_order_16.bin"
 **/
int oo_compute_grid_CMU_EM_3D(double dist_to_cm, double *tlim_CMU_EM,
                              int t_grid_size, double si_LIM_CMU_RCM[4][2],
                              int *si_grid_size, Invman &invman, bool isPar)
{
    //====================================================================================
    // Check that the invman is an unstable-manifold
    //====================================================================================
    if(invman.getManType() != MAN_CENTER_U)
    {
        cout << "oo_compute_grid_CMU_EM_3D. The invariant manifold must be of center-unstable type. return." << endl;
        return -1;
    }

    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------
    //Building the working grids
    //------------------------------------------
    double **grid_si_CMU_RCM = (double**) calloc(4, sizeof(double*));
    for(int i = 0; i <4; i++)
    {
        grid_si_CMU_RCM[i] = (double*) calloc(si_grid_size[i]+1, sizeof(double));
        init_grid(grid_si_CMU_RCM[i], si_LIM_CMU_RCM[i][0], si_LIM_CMU_RCM[i][1], si_grid_size[i]);
    }


    //------------------------------------------
    //Building the time grid
    //------------------------------------------
    double *grid_t_EM = dvector(0,  t_grid_size);
    init_grid(grid_t_EM, tlim_CMU_EM[0], tlim_CMU_EM[1], t_grid_size);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = (1+si_grid_size[0])*(1+si_grid_size[1])*(1+si_grid_size[2])*(1+si_grid_size[3])*(1+t_grid_size);
    int iter = 1;

    //------------------------------------------
    // Data structures
    //------------------------------------------
    double **init_state_CMU_NCEM = dmatrix(0, 5, 0, si_grid_size[2]);
    double **init_state_CMU_RCM  = dmatrix(0, 4, 0, si_grid_size[2]);

    //------------------------------------------
    // Reset the data file
    //------------------------------------------
    initCU_bin_3D(si_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU_3D);

    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= t_grid_size; kt++)
    {
        //----------------------
        //Append the time in data file
        //----------------------
        appTimeCU_bin_3D(grid_t_EM, kt, OFTS_ORDER, TYPE_CU_3D);


        for(int ks2 = 0; ks2 <= si_grid_size[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= si_grid_size[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= si_grid_size[0]; ks1++)
                {
                    #pragma omp parallel for if(isPar)  shared(iter)
                    for(int ks3 = 0; ks3 <= si_grid_size[2]; ks3++)
                    {
                        Ofsc ofs(OFS_ORDER);
                        double *yvu = dvector(0,5);
                        double *sti = dvector(0,4);

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
                    writeCU_bin_3D(init_state_CMU_NCEM, init_state_CMU_RCM, si_grid_size, OFTS_ORDER, TYPE_CU_3D);
                }
            }
        }
    }

    //------------------------------------------
    //Free
    //------------------------------------------
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, si_grid_size[2]);
    free_dmatrix(init_state_CMU_RCM,  0, 4, 0, si_grid_size[2]);

    return 1;
}

/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2,
 *         in the QBCP model. The initial conditions (IC) are computed in a
 *         three-dimensional box: one dimension for the starting time, two dimensions for
 *         the parameterization of the Center Manifold (s1 and s3 coordinates).
 *         The RCM coordinate s5 along the unstable direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates on the unstable direction s5.
 *  \param tmin_CMU_EM:    the minimum starting time (in EM units) in the IC box.
 *  \param tmax_CMU_EM:    the maximum starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param s1_MIN_CMU_RCM: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param s1_MAX_CMU_RCM: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param s3_MIN_CMU_RCM: the minimum s3 value (in RCM coordinates) in the IC box.
 *  \param s3_MAX_CMU_RCM: the maximum s3 value (in RCM coordinates) in the IC box.
 *  \param s1_grid_size:   the number of points on the s1 grid in the IC box.
 *  \param s3_grid_size:   the number of points on the s3 grid in the IC box.
 *  \param invman:         the center-unstable manifold, if the type of manifold provided
 *                         is not center-unstable, a warning message is displayed and
 *                         nothing is done.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form:
 *               "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int oo_compute_grid_CMU_EM(double dist_to_cm, double tmin_CMU_EM, double tmax_CMU_EM,
                           int t_grid_size, double s1_MIN_CMU_RCM, double s1_MAX_CMU_RCM,
                           double s3_MIN_CMU_RCM, double s3_MAX_CMU_RCM,
                           int s1_grid_size, int s3_grid_size, Invman &invman, bool isPar)

{
    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------
    //Building the working grids
    //------------------------------------------
    double *grid_s1_CMU_RCM = dvector(0,  s1_grid_size);
    double *grid_s3_CMU_RCM = dvector(0,  s3_grid_size);
    init_grid(grid_s1_CMU_RCM, s1_MIN_CMU_RCM, s1_MAX_CMU_RCM, s1_grid_size);
    init_grid(grid_s3_CMU_RCM, s3_MIN_CMU_RCM, s3_MAX_CMU_RCM, s3_grid_size);

    //------------------------------------------
    //Building the time grid
    //------------------------------------------
    double *grid_t_EM = dvector(0,  t_grid_size);
    init_grid(grid_t_EM, tmin_CMU_EM, tmax_CMU_EM, t_grid_size);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = (1+s1_grid_size)*(1+s3_grid_size)*(1+t_grid_size);
    int iter = 1;

    //------------------------------------------
    // Data structures
    //------------------------------------------
    double ****init_state_CMU_NCEM = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double ****init_state_CMU_RCM  = d4tensor(0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);

    //====================================================================================
    // Loop on all elements.
    //
    // Note that openMP is only used on the inner loop, since
    // it is useless to use on nested loops.
    //====================================================================================
    COMPLETION = 0;
    for(int kt = 0; kt <= t_grid_size; kt++)
    {
        for(int ks1 = 0; ks1 <= s1_grid_size; ks1++)
        {
            #pragma omp parallel for if(isPar)  shared(iter)
            for(int ks3 = 0; ks3 <= s3_grid_size; ks3++)
            {
                Ofsc ofs(OFS_ORDER);
                double *yvu = dvector(0,5);
                double *sti = dvector(0,4);

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
    writeCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, grid_t_EM, s1_grid_size, s3_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_RCM,  0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);

    return 1;
}


//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEML2
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM. Then, each point on the
 *         integration grid is projected on the Center Manifold CM_SEM_NC about SEML2.
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param tmax_on_manifold_EM: the maximum integration time on each leg, in EM units.
 *
 *  \param man_grid_size:       the number of points on each manifold leg.
 *
 *  \param nod:                 the number of dimensions on which the distance of
 *                              projection is computed (usually either 3 (the physical
 *                              distance) or 6 (the whole phase space)).
 *
 *  \param isPar:               if TRUE, the computation is parallelized.
 *
 *  \param ynormMax:            the maximum norm in NCSEM coordinates for which a given
 *                              state on the integration grid is projected on CM_SEM_NC
 *                              More precisely: for a given state y along the manifold leg,
 *                              if norm(y, 3) < ynormMax, the state is projected.
 *                              Otherwise, it is considered too far away from SEML2 to be
 *                              a good candidate for projection.
 *
 *  \param snormMax:            the maximum norm in RCM SEM coordinates for which a given
 *                              projection state on the CM of SEML2 (CM_SEM_NC) is
 *                              computed back in NCSEM coordinates. More precisely, for a
 *                              given state y in NCSEM coordinates, the result of the
 *                              projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                              coordinates.
 *                              if norm(sproj, 4) < snormMax, the computation
 *                              yproj = CM_SEM_NC(sproj, t) is performed. Otherwise, the
 *                              state sproj is considered too far away from the RCM origin
 *                              to be a good candidate - it is out of the domain of
 *                              practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form:
 *          "plot/QBCP/EM/L2/projcu_3d_order_16.bin".
 **/
int oo_int_proj_CMU_EM_on_CM_SEM_3D(double tmax_on_manifold_EM,
                                    int man_grid_size, int nod,
                                    int isPar, double ynormMax,
                                    double snormMax)
{
    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, si_grid_size[4], offset;
    offset = getLenghtCU_bin_3D(si_grid_size, &t_grid_size, OFTS_ORDER, TYPE_CU_3D);

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double **init_state_CMU_NCEM = dmatrix(0, 5, 0, si_grid_size[2]);
    double **init_state_CMU_RCM  = dmatrix(0, 4, 0, si_grid_size[2]);
    double *init_time_grid_EM    = dvector(0, t_grid_size);

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
        return -1;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);


    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    double **final_state_CMU_SEM      = dmatrix(0, 5, 0, si_grid_size[2]);
    double **projected_state_CMU_SEM  = dmatrix(0, 5, 0, si_grid_size[2]);
    double **projected_state_CMU_RCM  = dmatrix(0, 3, 0, si_grid_size[2]);
    double **init_state_CMU_SEM       = dmatrix(0, 5, 0, si_grid_size[2]);
    double **init_state_CMU_NCEM_0    = dmatrix(0, 5, 0, si_grid_size[2]);
    double *min_proj_dist_tens_SEM    = dvector(0, si_grid_size[2]);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================
    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //----------------------------------------------------------
    // projection by default is arbitrary big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Filename
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ_3D);
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
        //----------------------------------------------------------
        //Read time from file
        //----------------------------------------------------------
        offset = readTCU_bin_3D(offset, init_time_grid_EM, kt, OFTS_ORDER, TYPE_CU_3D);

        for(int ks2 = 0; ks2 <= si_grid_size[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= si_grid_size[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= si_grid_size[0]; ks1++)
                {
                    //----------------------------------------------------------
                    //Read data from file
                    //----------------------------------------------------------
                    offset = readCU_bin_3D(offset, init_state_CMU_NCEM, init_state_CMU_RCM, si_grid_size, OFTS_ORDER, TYPE_CU_3D);

                    //----------------------------------------------------------
                    //Most inner loop is parallelized
                    //----------------------------------------------------------
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
                        double **y_man_NCSEM = dmatrix(0, 5, 0, man_grid_size);
                        double *t_man_SEM    = dvector(0, man_grid_size);


                        //----------------------------------------------------------------
                        //Integration on man_grid_size+1 fixed grid
                        // PB: when Release + F_NCSEM!!
                        //----------------------------------------------------------------
                        int status = ode78_qbcp(y_man_NCSEM, t_man_SEM, tv, tv+tmax_on_manifold_EM, yv, 6, man_grid_size, F_NCEM, NCEM, NCSEM);

                        //================================================================
                        // 3.2. Projection on the center manifold of SEML2.
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
                        if(status != -1)
                        {
                            //------------------------------------------------------------
                            //Loop on trajectory
                            //------------------------------------------------------------
                            for(int kman = 0; kman <= man_grid_size; kman++)
                            {
                                //Current state
                                for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
                                tv = t_man_SEM[kman];

                                //Current distance from SEMLi in NCSEM units
                                y_man_norm_NCSEM = 0.0;
                                for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
                                y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

                                //--------------------------------------------------------
                                //Check n°1: the current state is close enough to SEMLi
                                //--------------------------------------------------------
                                if(y_man_norm_NCSEM < ynormMax)
                                {
                                    // Projection on the center manifold
                                    invman_SEM.NCprojCCMtoCM(yv, tv, sproj);

                                    //----------------------------------------------------
                                    //Check n°2: the projection is close enough to SEMLi
                                    //----------------------------------------------------
                                    if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
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
                                        for(int i = 0; i < nod; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
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

                        }
                        else proj_dist_SEM = ePdef;


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
                        //free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
                        free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
                        //free_dvector(t_man_EM, 0, man_grid_size);
                        free_dvector(t_man_SEM, 0, man_grid_size);


                    }
                }
            }
        }
    }
    return 1;
}


/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions
 *         obtained using the routine compute_grid_CMU_EM.
 *         Then, each point on the integration grid is projected on the center Manifold
 *         about SEML2 (denoted here CM_SEM_NC).
 *         The best solution (minimum distance of projection) is stored.
 *
 *  \param tmax_on_manifold_EM: the maximum integration time on each leg, in EM units.
 *  \param t_grid_size_x:.......the number of points on the time grid in the IC box.
 *                              If -1, the value used in compute_grid_CMU_EM is used.
 *  \param s1_grid_size_x:......the number of points on the s1 grid in the IC box.
 *                              If -1, the value used in compute_grid_CMU_EM is used.
 *  \param s3_grid_size_x:......the number of points on the s3 grid in the IC box.
 *                              If -1, the value used in compute_grid_CMU_EM is used.
 *  \param man_grid_size:.......the number of points on each manifold leg.
 *  \param NsortMin:............the number of best solutions that are kept
 *  \param nod:.................the number of dimensions on which the distance of
 *                              projection is computed - usually either 3
 *                              (the physical distance) or 6 (the whole phase space).
 *  \param isPar:...............if TRUE, the computation is parallelized.
 *  \param ynormMax:............the maximum norm in NCSEM coordinates for which a given
 *                              state on the integration grid is projected on CM_SEM_NC.
 *                              More precisely: for a given state y along the manifold leg,
 *                              if norm(y, 3) < ynormMax, the state is projected.
 *                              Otherwise, it is considered too far away from SEML2 to
 *                              be a good candidate for projection.
 *  \param snormMax:............the maximum norm in RCM SEM coordinates for which a given
 *                              projection state on the CM of SEML2 (CM_SEM_NC) is
 *                              computed back in NCSEM coordinates. More precisely, for a
 *                              given state y in NCSEM coordinates, the result of the
 *                              projection on CM_SEM_NC gives a state sproj in RCM SEM
 *                              coordinates. If norm(sproj, 4) < snormMax, the computation
 *                              yproj = CM_SEM_NC(sproj, t) is performed. Otherwise,
 *                              the state sproj is considered too far away from the RCM
 *                              origin to be a good candidate - it is out of the domain of
 *                              practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form "plot/QBCP/EM/L2/projcu_order_16.bin", and "plot/QBCP/EM/L2/sortprojcu_order_16.bin" for the NsortMin best solutions.
 **/
int oo_int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM, int t_grid_size_x,
                                 int s1_grid_size_x, int s3_grid_size_x,
                                 int man_grid_size, int NsortMin, int nod, int isPar,
                                 double ynormMax, double snormMax)
{
    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, s1_grid_size, s3_grid_size;
    getLenghtCU_bin(&s1_grid_size, &s3_grid_size, &t_grid_size, OFTS_ORDER, TYPE_CU);

    //----------------------------------------------------------
    //Initialize the actual sizes comparing t_grid_size/t_grid_size_x
    //----------------------------------------------------------
    int t_grid_size_t  = t_grid_size_x  < 0 ? t_grid_size:t_grid_size_x;
    int s1_grid_size_t = s1_grid_size_x < 0 ? s1_grid_size:s1_grid_size_x;
    int s3_grid_size_t = s3_grid_size_x < 0 ? s3_grid_size:s3_grid_size_x;

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double ****init_state_CMU_NCEM = d4tensor(0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double ****init_state_CMU_RCM  = d4tensor(0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    double *init_time_grid_EM      = dvector(0, t_grid_size);

    //----------------------------------------------------------
    //Read data from file
    //----------------------------------------------------------
    readCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM, init_time_grid_EM,
               s1_grid_size, s3_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU);

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
        cout << "oo_int_proj_CMU_EM_on_CM_SEM_3D. The invariant manifold at SEMLi ";
        cout << "must be of center type. return." << endl;
        return -1;
    }
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);


    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    double ****final_state_CMU_SEM      = d4tensor(0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    double ****projected_state_CMU_SEM  = d4tensor(0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    double ****projected_state_CMU_RCM  = d4tensor(0, 3, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    double ****init_state_CMU_SEM       = d4tensor(0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    double ****init_state_CMU_NCEM_0    = d4tensor(0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    double ***min_proj_dist_tens_SEM    = d3tensor(0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);


    //====================================================================================
    // 2.3. Misc initialization. Require better presentation?
    //====================================================================================

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //----------------------------------------------------------
    // projection by default is arbitrary big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Filename
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ);
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
    for(int kt = 0; kt <= t_grid_size_t; kt++)
    {
        for(int ks1 = 0; ks1 <= s1_grid_size_t; ks1++)
        {
            #pragma omp parallel for if(isPar) shared(index)
            for(int ks3 = 0; ks3 <= s3_grid_size_t; ks3++)
            {
                //====================================================================================
                // 3.1. Integration of the manifold leg
                //====================================================================================
                //---------------------------------------------------------------------
                //Initialize the initial conditions (both NC and RCM coordinates)
                //---------------------------------------------------------------------
                double yv[6], tv;
                tv  = init_time_grid_EM[kt];
                for(int i = 0; i < 6; i++) yv[i] = init_state_CMU_NCEM[i][kt][ks1][ks3];

                //---------------------------------------------------------------------
                //Local variables to store the manifold leg
                //---------------------------------------------------------------------
                double **y_man_NCSEM = dmatrix(0, 5, 0, man_grid_size);
                double *t_man_SEM    = dvector(0, man_grid_size);


                //---------------------------------------------------------------------
                //Integration on man_grid_size+1 fixed grid
                // PB: when Release + F_NCSEM!!
                //---------------------------------------------------------------------
                int status = ode78_qbcp(y_man_NCSEM, t_man_SEM, tv, tv+tmax_on_manifold_EM, yv, 6, man_grid_size, F_NCEM, NCEM, NCSEM);

                //====================================================================================
                // 3.2. Projection on the center manifold of SEML2. No use of SEML_EM or SEML after this point!
                // We need first to check that the integration went well
                //====================================================================================
                //----------------------------------------------------------
                //Temp variables
                //----------------------------------------------------------
                Ofsc ofs(OFS_ORDER);
                double yvproj_NCSEM[6], sproj[4], yv_SEM[6], yvproj_SEM[6], yvEM[6], yv_VSEM[6], yvproj_VSEM[6];
                double proj_dist_SEM, min_proj_dist_SEM = ePdef, dv_at_projection_SEM = 0.0, y_man_norm_NCSEM = 0.0;
                int ksort = (s1_grid_size_t+1)*(s3_grid_size_t+1)*kt + (s3_grid_size_t+1)*ks1 + ks3;
                int kmin = 0;

                //----------------------------------------------------------
                // We need first to check that the integration went well
                //----------------------------------------------------------
                if(status != -1)
                {
                    //----------------------------------------------------------
                    //Loop on trajectory
                    //----------------------------------------------------------
                    for(int kman = 0; kman <= man_grid_size; kman++)
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
                        if(y_man_norm_NCSEM < ynormMax)
                        {
                            // Projection on the center manifold
                            invman_SEM.NCprojCCMtoCM(yv, tv, sproj);

                            //------------------------------------------------------------
                            //Check n°2: the projection is close enough to SEMLi
                            //------------------------------------------------------------
                            if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
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
                                for(int i = 0; i < nod; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
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

                }
                else proj_dist_SEM = ePdef;


                //====================================================================================
                // 3.3. Save outputs
                //====================================================================================
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
                        distMin[ksort]  = min_proj_dist_SEM;

                        //----------------------------------------------------------
                        //Open datafile
                        //----------------------------------------------------------
                        writeIntProjCU_bin(filename, init_time_grid_EM, init_state_CMU_NCEM,
                        init_state_CMU_SEM, init_state_CMU_RCM, final_state_CMU_SEM, projected_state_CMU_SEM,
                        projected_state_CMU_RCM, min_proj_dist_SEM, dv_at_projection_SEM, t_man_SEM, kmin, ks1, ks3, kt);
                    }
                }

                //----------------------------------------------------------
                //Display completion
                //----------------------------------------------------------
                #pragma omp critical
                {
                    displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/((1+s1_grid_size_t)*(1+s3_grid_size_t)*(1+t_grid_size_t)));
                }

                //---------------------------------------------------------------------
                //Free
                //---------------------------------------------------------------------
                //free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
                free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
                //free_dvector(t_man_EM, 0, man_grid_size);
                free_dvector(t_man_SEM, 0, man_grid_size);
            }
        }
    }

    //====================================================================================
    // 4. Sorting the best solutions
    //====================================================================================
    vector<size_t> sortId = sort_indexes(distMin);

    //----------------------------------------------------------
    //Saving the NsortMin best results or all results if less than 50 have been computed
    //----------------------------------------------------------
    int number_of_sol  = min(NsortMin, Nsort-2);


    //---------------------
    //Filename
    //---------------------
    filename = filenameCUM(OFTS_ORDER, TYPE_MAN_SORT);

    //---------------------
    //Write sorted solutions
    //---------------------
    writeIntProjSortCU_bin(filename, init_state_CMU_NCEM, init_state_CMU_RCM,
                           final_state_CMU_SEM,  projected_state_CMU_SEM, projected_state_CMU_RCM,
                           min_proj_dist_tens_SEM, sortId, ktMin, ks1Min, ks3Min, t0_min_EM, tf_min_EM, distMin, number_of_sol);


    //----------------------------------------------------------
    //Free
    //----------------------------------------------------------
    free_d4tensor(final_state_CMU_SEM, 0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    free_d4tensor(projected_state_CMU_SEM, 0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    free_d4tensor(projected_state_CMU_RCM, 0, 3, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    free_d4tensor(init_state_CMU_SEM, 0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    free_d4tensor(init_state_CMU_NCEM_0, 0, 5, 0, t_grid_size_t, 0, s1_grid_size_t, 0, s3_grid_size_t);
    free_d3tensor(min_proj_dist_tens_SEM, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);


    //---------------------------------------------------------------------
    //Free
    //---------------------------------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_RCM, 0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_dvector(init_time_grid_EM, 0, t_grid_size);


    return 1;
}


//========================================================================================
//
//         Refinement of solutions: CMU to CMS - general routines
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2.
 *         The final conditions vary in the paramerization of the CMS of SEML2.
 *         The time at each point except the first one is allowed to vary.
 *         A continuation procedure can be performed to get more than one refined solution.
 **/
int refeml2seml(int man_grid_size,
                int coord_type,
                Invman &invman,
                RefSt &refst)
{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 2. Structures to compute the final orbit about SEML2
    //====================================================================================
    Invman invman_SEM(OFTS_ORDER, OFS_ORDER, SEML_SEM.cs);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    //---------------------------------------------------------------------
    //To store data from the data file
    //---------------------------------------------------------------------
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

    //---------------------------------------------------------------------
    //To store data from a sub selection of the previous elements
    //---------------------------------------------------------------------
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
    // 4. Getting back the data
    //====================================================================================
    int type = TYPE_MAN_PROJ;
    switch(refst.dim)
    {
        case REF_PLANAR:
            type = refst.isFromServer? TYPE_MAN_PROJ_FROM_SERVER:TYPE_MAN_PROJ;
        break;

        case REF_3D:
            type = refst.isFromServer? TYPE_MAN_PROJ_3D_FROM_SERVER:TYPE_MAN_PROJ_3D;
        break;
    }
    string filename = filenameCUM(OFTS_ORDER, type);

    readIntProjCU_bin(filename, t0_EM, tf_EM, s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM,
            s5_CMU_EM, pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM, sortId);

    //====================================================================================
    // 5. Select given intervals for the inputs
    //====================================================================================
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;

    if(refst.isLimUD)
    {
        cout << "-----------------------------------------------------" << endl;
        cout << "   refeml2seml. User-defined interval of research    " << endl;
        cout << "-----------------------------------------------------" << endl;
        cout << "Enter a value for s1_CMU_EM_MIN: ";
        cin >> s1_CMU_EM_MIN;
        cout << "Enter a value for s1_CMU_EM_MAX: ";
        cin >> s1_CMU_EM_MAX;
        cout << "Enter a value for s3_CMU_EM_MIN: ";
        cin >> s3_CMU_EM_MIN;
        cout << "Enter a value for s3_CMU_EM_MAX: ";
        cin >> s3_CMU_EM_MAX;
    }else
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
    // 7. User-defined number of continuation steps, if necessary
    //====================================================================================
    int cont_steps_MAX = 0;
    if(refst.type == REF_CONT || refst.type == REF_CONT_D)
    {
        if(refst.isLimUD)
        {
            cout << "-----------------------------------------------------" << endl;
            cout << "refeml2seml. User-defined number of continuation steps" << endl;
            cout << "-----------------------------------------------------" << endl;
            cout << "Enter a value for cont_steps_MAX: ";
            cin >> cont_steps_MAX;
            refst.cont_step_max = cont_steps_MAX;
        }
    }else
    {
        cout << "refeml2seml. No continuation will be performed." << endl;
        refst.cont_step_max = 0;
    }

    //====================================================================================
    // 7. The best solution in the subselection will serve as the first guess.
    //====================================================================================
    int kpos = sortId_R[0];

    cout << "-----------------------------------------------------" << endl;
    cout << "   refeml2seml. Refinement of EML2-SEML2 arc         " << endl;
    cout << "-----------------------------------------------------" << endl;
    cout << "Estimated error at patch point (km):       "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]*SEML.cs_sem.cr3bp.L << endl;
    cout << "Estimated error at patch point (SEMSU):    "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]                           << endl;
    cout << "-----------------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // 7.1 Initialize local variables: EM
    //====================================================================================
    //---------------------------------------------------------------------
    //Initialize the initial conditions (both NC and RCM coordinates)
    //---------------------------------------------------------------------
    double st_EM[5], st_SEM[5];
    st_EM[0] = s1_CMU_EM_R[kpos];
    st_EM[1] = s2_CMU_EM_R[kpos];
    st_EM[2] = s3_CMU_EM_R[kpos];
    st_EM[3] = s4_CMU_EM_R[kpos];
    st_EM[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    OdeStruct driver_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbfbp_vfn_novar, &SEML_EM);
    //Init routine
    Orbit orbit_EM(&invman, &SEML_EM, &driver_EM, OFTS_ORDER, OFS_ORDER, t0_EM_R[kpos], tf_EM_R[kpos]);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_EM.update_ic(st_EM,  t0_EM_R[kpos]);

    //====================================================================================
    // 7.2 Initialize local variables: SEM
    //====================================================================================
    //---------------------------------------------------------------------
    //Initialize the initial conditions (both NC and RCM coordinates)
    //---------------------------------------------------------------------
    st_SEM[0] = s1_CM_SEM_R[kpos];
    st_SEM[1] = s2_CM_SEM_R[kpos];
    st_SEM[2] = s3_CM_SEM_R[kpos];
    st_SEM[3] = s4_CM_SEM_R[kpos];
    st_SEM[4] = 0.0;

    //Initial time in SEM units
    double t0_SEM = tf_EM_R[kpos]*SEML.us_em.ns;

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    OdeStruct driver_SEM;
    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbfbp_vfn_novar,&SEML_SEM);
    //Init routine
    Orbit orbit_SEM(&invman_SEM, &SEML_SEM, &driver_SEM, OFTS_ORDER, OFS_ORDER, t0_SEM, t0_SEM+5*SEML.us_sem.T);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_SEM.update_ic(st_SEM, t0_SEM);

    //====================================================================================
    // 8 Multiple shooting procedure
    //====================================================================================
    double **semP_coord = dmatrix(0, 6, 0, 2);
    gnuplot_ctrl *h2;
    switch(refst.type)
    {
    case REF_CONT:
    case REF_SINGLE:
        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        switch(coord_type)
        {
        case PSEM:
            semPoints(0.0, semP_coord);
            semPlot(h2, semP_coord);
            break;
        case NCSEM:
        case VNCSEM:
            semNCPoints(0.0, semP_coord);
            semPlot(h2, semP_coord);
            break;

        case PEM:
            emPoints(0.0, semP_coord);
            emPlot(h2, semP_coord);
            break;
        case NCEM:
        case VNCEM:
            emNCPoints(0.0, semP_coord);
            emPlot(h2, semP_coord);
            break;
        }

        //--------------------------------------------------------------------------------
        //Multiple shooting procedure
        //--------------------------------------------------------------------------------
        srefeml2seml(orbit_EM, orbit_SEM, dcs, coord_type, man_grid_size, refst, h2);
        gnuplot_close(h2);
        break;

    case REF_CONT_D:
        //--------------------------------------------------------------------------------
        //REF_CONT for the rest of the computation
        //--------------------------------------------------------------------------------
        refst.type = REF_CONT;

        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        h2 = gnuplot_init();

        //--------------------------------------------------------------------------------
        //Notable points in SEM system
        //--------------------------------------------------------------------------------
        switch(coord_type)
        {
        case PSEM:
            semPoints(0.0, semP_coord);
            semPlot(h2, semP_coord);
            break;
        case NCSEM:
        case VNCSEM:
            semNCPoints(0.0, semP_coord);
            semPlot(h2, semP_coord);
            break;

        case PEM:
            emPoints(0.0, semP_coord);
            emPlot(h2, semP_coord);
            break;
        case NCEM:
        case VNCEM:
            emNCPoints(0.0, semP_coord);
            emPlot(h2, semP_coord);
            break;
        }

        //--------------------------------------------------------------------------------
        //First step: variable time
        //--------------------------------------------------------------------------------
        refst.time = REF_VAR_TIME;
        srefeml2seml(orbit_EM, orbit_SEM, dcs, coord_type, man_grid_size, refst, h2);


        //--------------------------------------------------------------------------------
        //Second step: fixed time
        //--------------------------------------------------------------------------------
        refst.time = REF_FIXED_TIME;
        srefeml2seml(orbit_EM, orbit_SEM, dcs, coord_type, man_grid_size, refst, h2);

        gnuplot_close(h2);
        break;

    case REF_COMP:
        cout << "-------------------------------------------------------" << endl;
        cout << "refeml2seml. User-defined number of points in REF_COMP " << endl;
        cout << "-------------------------------------------------------" << endl;
        cout << "Enter a value for msd_grid_size: ";
        int msd_grid_size;
        cin >> msd_grid_size;
        char ch;
        scanf("%c",&ch);


        comprefft3d(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);

        //comprefft3d_test_eml2seml_synjpl(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        //comprefft3d_test_seml_synjpl(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        //comprefft3d_test_seml_synjpl(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        //comprefft3d_test_eml2seml_insem(msd_grid_size, coord_type, orbit_EM, orbit_SEM, refst);
        break;
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);


    return 0;
}

//========================================================================================
//
//         Refinement of solutions: CMU to CMS - subroutines
//
//========================================================================================
/**
 *  \brief Continuation of a single of EML2-to-SEML2 connection, between orbit_EM and orbit_SEM.
 *         The Jacobian of the parameterization of the manifolds are contained in  DCM_EM_TFC and DCMS_SEM_TFC.
**/
void srefeml2seml(Orbit &orbit_EM,
                  Orbit &orbit_SEM,
                  int dcs,
                  int coord_type,
                  int man_grid_size_t,
                  RefSt &refst,
                  gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    cout << " srefeml2seml. Initialization of the local variables..."  << endl;
    //Local plotting condition
    int isPlotted   = refst.isPlotted;
    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = floor((double)refst.cont_step_max/floor(0.1*GP_MAX_TMP_FILES))+1;

    //Integration on the final orbit
    double tof_seml_SEM = 30*SEML.us_sem.T;   //TOF on SEML2 orbit

    //---------------------------------------------------------------------
    //If the grid is fixed, then we use the user-provided size. Else, we use a
    //big value (max_grid_size) so that the arrays would not be saturated.
    //---------------------------------------------------------------------
    int max_grid_size = 1000;
    int man_grid_size = refst.grid == REF_FIXED_GRID? man_grid_size_t:max_grid_size;

    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double  *t_man_coord  = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 500;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    double **y_traj       = dmatrix(0, 41, 0, man_grid_size);
    double  *t_traj       = dvector(0, man_grid_size);
    double **y_traj_comp  = dmatrix(0, 41, 0, man_grid_size);
    double  *t_traj_comp  = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, man_grid_size);
    double *t_traj_n    = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    cout << " srefeml2seml. Compute the manifold leg..."  << endl;
    //---------------------------------------------------------------------
    // Integration: from orbit_EM.t0 to orbit_EM.tf
    //---------------------------------------------------------------------
    int man_index = man_grid_size;
    if(refst.grid == REF_FIXED_GRID) //fixed time grid
        ode78_qbcp(y_man_coord, t_man_coord, orbit_EM.getT0(), orbit_EM.getTf(), orbit_EM.getZ0(), 6, man_grid_size, dcs, NCEM, coord_type);
    else //variable time grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.getT0(), orbit_EM.getTf(), orbit_EM.getZ0(), 6, man_grid_size, dcs, NCEM, coord_type, man_grid_size_t);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    // Plot the resulting trajectory
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], man_index, (char*)"", "lines", "1", "1", 4);

    //====================================================================================
    // 5. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    cout << " srefeml2seml. Compute the first point of the final SEML2 orbit..."  << endl;
    //---------------------------------------------------------------------
    // Save the last point at position man_index
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.getT0(), orbit_SEM.getZ0(), yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_index] = yout[i];
    t_traj[man_index] = tout;

    //---------------------------------------------------------------------
    // Plot the resulting trajectory
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], man_index+1, (char*)"", "points", "1", "1", 4);


    //====================================================================================
    // 6.  Differential correction
    //====================================================================================
    cout << " srefeml2seml. Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //======================================================================
    //Free variables
    //======================================================================
    int nfv = nfreevariables(refst, man_index);
    double *nullvector = dvector(0, nfv-1);

    //======================================================================
    //GNUPLOT:
    //  1. if there is a continuation procedure based on variable time,
    //     the component s5 at SEML is plotted.
    //  2. if there is a continuation procedure based on fixed time,
    //     several computations are produced.
    //======================================================================
    gnuplot_ctrl *h3 = 0, *h4 = 0, *h5 = 0;
    if(refst.type == REF_CONT)
    {
        if(refst.time == REF_VAR_TIME)
        {
            h5 = gnuplot_init();
            gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
            gnuplot_cmd(h5, "set grid");
        }
        if(refst.time == REF_FIXED_TIME)
        {
            h3 = gnuplot_init();
            gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");
            gnuplot_cmd(h3, "set grid");

            h4 = gnuplot_init();
            gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");
            gnuplot_cmd(h4, "set grid");

            h5 = gnuplot_init();
            gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
            gnuplot_cmd(h5, "set grid");
        }
    }

    cout << "NewNew. y_traj[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][0] << endl;

    cout << "NewNew. y_traj[end-1] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][man_grid_size-1] << endl;

    cout << "NewNew. y_traj[end] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][man_grid_size] << endl;


    //======================================================================
    //6.1. First step of the continuation procedure.
    //======================================================================
    int status = 0;
    double inner_prec = 5e-8;
    switch(refst.dim)
    {
    case REF_3D:
        switch(refst.time)
        {
        case REF_FIXED_TIME:
            status = msft3d(y_traj, t_traj, y_traj_n, t_traj_n,
                            nullvector, 42,
                            man_index, coord_type, inner_prec, true,
                            orbit_EM, orbit_SEM,
                            h2, refst);

            break;

        case REF_VAR_TIME:
            status = msvt3d(y_traj, t_traj, y_traj_n, t_traj_n,
                            nullvector, 42,
                            man_index, coord_type, inner_prec, true,
                            orbit_EM, orbit_SEM,
                            h2, refst);

            break;
        default:
            perror("srefeml2seml. Unknown refst.time.");
            break;
        }
        break;
    case REF_PLANAR:
        switch(refst.time)
        {
        case REF_FIXED_TIME:
            status = msftplan(y_traj, t_traj, y_traj_n, t_traj_n,
                              nullvector, 42,
                              man_index, coord_type, inner_prec, true,
                              orbit_EM, orbit_SEM,
                              h2, refst);

            break;

        case REF_VAR_TIME:
            status = msvtplan(y_traj, t_traj, y_traj_n, t_traj_n,
                              nullvector, 42,
                              man_index, coord_type, inner_prec, true,
                              orbit_EM, orbit_SEM,
                              h2, refst);

            break;
        default:
            perror("srefeml2seml. Unknown refst.time.");
            break;
        }
         break;
    default:
        perror("srefeml2seml. Unknown refst.dim.");
        break;
    }

    cout << "NewNew. y_traj_n[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj_n[i][0] << endl;

    //====================================================================================
    //6.2. If it is a sucess, we go on.
    //====================================================================================
    if(status == GSL_SUCCESS)
    {

        //================================================================================
        // Save first entry
        //================================================================================
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);
        double yv[42];
        string filename      = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, orbit_EM.getT0()/SEML.us_em.T);
        string filename_traj = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF_TRAJ, orbit_EM.getT0()/SEML.us_em.T);
        fstream filestream;
        if(refst.isSaved)
        {
            cout << " srefeml2seml. Save first entry in " << filename_traj         << endl;
            cout << "-----------------------------------------------------------"  << endl;
            //============================================================================
            // Main parameters
            //============================================================================
            filestream.open (filename.c_str(), ios::out);
            //Title
            filestream << "t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
            filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
            filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
            filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
            filestream << endl;
            //Data
            filestream << setprecision(15) <<  setiosflags(ios::scientific);
            filestream << orbit_EM.getT0() << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_EM.getSi()[i]  << "  ";
            filestream << orbit_EM.getTf() << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_SEM.getSi()[i] << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_EM.getZ0()[i]  << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_SEM.getZ0()[i] << "  ";
            filestream << endl;
            filestream.close();


            //============================================================================
            // Entire trajectory
            //============================================================================
            filestream.open (filename_traj.c_str(), ios::out | ios::binary);
            double res;
            //Final trajectory on lines, segment by segment
            for(int k = 0; k < man_index; k++)
            {
                //Integration segment by segment
                for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                ode78_qbcp(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
            //            cout << "-------------------------------------------"  << endl;
            //            cout << "Compute the initial SEML2 orbit            "  << endl;
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
            //            //---------------------------------------------------------------------
            //            //For dot(s) = fh(s)
            //            //---------------------------------------------------------------------
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
            //            //---------------------------------------------------------------------
            //            // Temp variables
            //            //---------------------------------------------------------------------
            //            double t0_SEM = t_traj_n[man_index];
            //            double t1_SEM = t0_SEM+tof_seml_SEM;
            //            double z[6];
            //            double t2 = t0_SEM;
            //            int k  = 0;
            //            double  s1ccm8[2*rvf.reduced_nv]; //CCM8
            //
            //            //---------------------------------------------------------------------
            //            // Initial state in CCM8 form
            //            //---------------------------------------------------------------------
            //            RCMtoCCM8(orbit_SEM.getSi(), s1ccm8, 5);
            //
            //            //---------------------------------------------------------------------
            //            // Reopen the file
            //            //---------------------------------------------------------------------
            //            filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
            //
            //            //---------------------------------------------------------------------
            //            // Loop
            //            //---------------------------------------------------------------------
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
            //                gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);
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
            //            filestream.close();

        }


        cout << " srefeml2seml. Continuation procedure.                     "  << endl;
        cout << "-----------------------------------------------------------"  << endl;
        //======================================================================
        //Continuation procedure
        //======================================================================
        double ds0 = 1e-1;
        switch(refst.time )
        {
        case REF_VAR_TIME:
            ds0 = 1e-1;
            break;
        case REF_FIXED_TIME:
            ds0 = 5e-1;
            break;
        }
        double ds  = ds0;
        double dkn = 0.0;
        int kn = 0;

        //An additional condition can be set: if the time is not fixed, we want to decrease the
        //s5 component at SEMLi. Therefore, we add a stopping condition on the loop.
        bool addCondition = refst.time == REF_VAR_TIME? fabs(orbit_SEM.getSi()[4]) > ORBIT_SEM_UNSTABLE_MIN: true;

        //======================================================================
        // 5.2.1. Updating the free variables
        //======================================================================
        while(kn < refst.cont_step_max && status == GSL_SUCCESS && addCondition)
        {
            switch(refst.dim)
            {
            case REF_3D:
                switch(refst.time)
                {
                case REF_FIXED_TIME:
                //                    ufvarft3d(y_traj_n, t_traj_n, ds, nullvector,
                //                              orbit_EM, orbit_SEM,
                //                              man_grid_size, coord_type);
                    break;

                case REF_VAR_TIME:
                //                    ufvarvt3d(y_traj_n, t_traj_n, ds, nullvector,
                //                              orbit_EM, orbit_SEM,
                //                              man_grid_size, coord_type);
                    break;
                default:
                    perror("srefeml2seml. Unknown refst.time.");
                    break;
                }
                break;
            case REF_PLANAR:
                switch(refst.time)
                {
                case REF_FIXED_TIME:
                        ufvarftplan(y_traj_n, t_traj_n, ds, nullvector,
                                    orbit_EM, orbit_SEM,
                                    man_grid_size, coord_type);
                    break;

                case REF_VAR_TIME:
                    ufvarvtplan(y_traj_n, t_traj_n, &ds, ds0, nullvector,
                                orbit_EM, orbit_SEM,
                                man_grid_size, coord_type);
                    break;
                default:
                    perror("srefeml2seml. Unknown refst.time.");
                    break;
                }
                break;
            default:
                perror("srefeml2seml. Unknown refst.dim.");

                break;
            }

            //======================================================================
            // 5.2.2. Diff correction
            //======================================================================
            isPlotted = (kn % plotfreq == 0) && refst.isPlotted ? 1:0;
            switch(refst.dim)
            {
            case REF_3D:
                switch(refst.time)
                {
                case REF_FIXED_TIME:
                    status = msft3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                    nullvector, 42,
                                    man_index, coord_type, inner_prec, false,
                                    orbit_EM, orbit_SEM,
                                    h2, refst);

                    break;

                case REF_VAR_TIME:
                    status = msvt3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                    nullvector, 42,
                                    man_index, coord_type, inner_prec, false,
                                    orbit_EM, orbit_SEM,
                                    h2, refst);

                    break;
                default:
                    perror("srefeml2seml. Unknown refst.time.");
                    break;
                }
                break;
            case REF_PLANAR:
                switch(refst.time)
                {
                case REF_FIXED_TIME:
                    status = msftplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                      nullvector, 42,
                                      man_index, coord_type, inner_prec, false,
                                      orbit_EM, orbit_SEM,
                                      h2, refst);

                    break;

                case REF_VAR_TIME:
                    status = msvtplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                      nullvector, 42,
                                      man_index, coord_type, inner_prec, false,
                                      orbit_EM, orbit_SEM,
                                      h2, refst);

                    break;
                default:
                    perror("srefeml2seml. Unknown refst.time.");
                    break;
                }
                break;
            default:
                perror("srefeml2seml. Unknown refst.dim.");
                break;

            }


            //======================================================================
            //Save
            //======================================================================
            if(refst.isSaved && status == GSL_SUCCESS)
            {
                filestream.open (filename.c_str(), ios::out | ios::app);
                //Data
                filestream << setprecision(15) <<  setiosflags(ios::scientific);
                filestream << orbit_EM.getT0() << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_EM.getSi()[i]  << "  ";
                filestream << orbit_EM.getTf() << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_SEM.getSi()[i] << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_EM.getZ0()[i]  << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_SEM.getZ0()[i] << "  ";
                filestream << endl;
                filestream.close();

                //=================================
                // Entire trajectory
                //=================================
                filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
                double res;
                //Final trajectory on lines, segment by segment
                for(int k = 0; k < man_index; k++)
                {
                    //Integration segment by segment
                    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                    ode78_qbcp(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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


            //======================================================================
            // 5.2.3. Display
            //======================================================================
            dkn = (double) kn;
            if(kn % plotfreq == 0)
            if(refst.type == REF_CONT)
            {
                if(refst.time == REF_VAR_TIME)
                gnuplot_plotc_xy(h5, &dkn, &orbit_SEM.getSi()[4], 1, (char*)"", "points", "1", "2", 0);

                if(refst.time == REF_FIXED_TIME)
                {
                    gnuplot_plotc_xy(h3, &orbit_EM.getSi()[0],  &orbit_EM.getSi()[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plotc_xy(h4, &orbit_SEM.getSi()[0], &orbit_SEM.getSi()[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plotc_xy(h5, &dkn, &orbit_SEM.getSi()[4], 1, (char*)"", "points", "1", "2", 0);
                }
            }


            //======================================================================
            // 5.2.4. Advance
            //======================================================================
            kn++;
            cout << "Step n°" << kn << "/" << refst.cont_step_max << " completed.               "  << endl;
            cout << "---------------------------------------------"  << endl;
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);

            //======================================================================
            // 5.2.5. update the additionnal condition, if necessary
            //======================================================================
            addCondition = refst.time == REF_VAR_TIME? fabs(orbit_SEM.getSi()[4]) > ORBIT_SEM_UNSTABLE_MIN: true;

        }

    }

    //====================================================================================
    // 7. Display final solution
    //====================================================================================
    if(status == GSL_SUCCESS)
    {
        //=============================================================
        // 7.1. Compute the initial orbit
        //=============================================================
        //TOF on EML2 orbit
        double tof_eml_EM   = 5*SEML.us.T;
        //TOF on SEML2 orbit
        //double tof_seml_SEM = 30*SEML.us_sem.T;
        cout << "---------------------------------------------------------"  << endl;
        cout << " srefeml2seml. Computation of the initial EML2 orbit     "  << endl;
        cout << "---------------------------------------------------------"  << endl;

        //---------------------------------------------------------------------
        // Reset the unstable direction
        //---------------------------------------------------------------------
        orbit_EM.setSi(0, 4);
        orbit_EM.setTf(orbit_EM.getT0()-tof_eml_EM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        orbit_EM.traj_int_grid(orbit_EM.getTf(), y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
        //gnuplot_plot_xyz(h6, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"EML2 orbit", "lines", "1", "2", 2);

        //================================================================================
        // 7.2. Final trajectory, on a grid
        //================================================================================
        cout << "----------------------------------------------"  << endl;
        cout << " srefeml2seml Final trajectory, on a grid.    "  << endl;
        cout << "----------------------------------------------"  << endl;
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);
        double yv[42];
        //Final trajectory on lines, segment by segment
        for(int k = 0; k < man_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78_qbcp(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
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
    if(status == GSL_SUCCESS)
    {
        //====================================
        // At EML2. The first 4 RCM components are good, as well as the initial time.
        // Hence, we need to update:
        // 1. The last RCM component (unstable part),
        // 2. The final time.
        //====================================
        orbit_EM.setSi(PROJ_EPSILON, 4);
        orbit_EM.setTf(t_traj_n[man_index]/SEML.us_em.ns);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_EM.update_ic(orbit_EM.getSi(), orbit_EM.getT0());

        //====================================
        // At SEML2. The first 5 RCM components are good.
        // Hence, we need to update:
        // 1. The initial time.
        //====================================
        orbit_SEM.setT0(t_traj_n[man_index]);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_SEM.update_ic(orbit_SEM.getSi(), orbit_SEM.getT0());
    }

    cout << "---------------------------------------------"  << endl;
    cout << " crefvtplan.                                 "  << endl;
    cout << "          End of computation.                "  << endl;
    cout << "---------------------------------------------"  << endl;
    printf("Press ENTER to close all\n");
    scanf("%c",&ch);
    scanf("%c",&ch);

    //====================================================================================
    // 9. Free
    //====================================================================================
    if(refst.type == REF_CONT)
    {
        if(refst.time == REF_VAR_TIME)
            gnuplot_close(h5);

        if(refst.time == REF_FIXED_TIME)
        {
            gnuplot_close(h3);
            gnuplot_close(h4);
            gnuplot_close(h5);
        }
    }
    //gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, man_grid_size);
    free_dvector(t_traj, 0, man_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, man_grid_size);
    free_dvector(t_traj_comp, 0, man_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, man_grid_size);
    free_dvector(t_traj_n, 0, man_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}


//========================================================================================
//
//         Refinement of solutions: Complete trajectory
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 *         A multiple_shooting_direct is applied on the WHOLE trajectory
 *         (EML2 orbit + manifold leg + SEML2 orbit).
 **/
int comprefft3d(int man_grid_size_t,
                int coord_type,
                Orbit &orbit_EM,
                Orbit &orbit_SEM,
                RefSt refst)

{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    cout << " comprefft3d. Initialization of the local variables..."  << endl;

    //----------------------------------------------------------
    //Get the default coordinates system from the coord_type
    //----------------------------------------------------------
    int dcs  = default_coordinate_system(coord_type);

    //----------------------------------------------------------
    // Define the time of flight on each orbit
    //----------------------------------------------------------
    double tof_eml_EM   = refst.tspan_EM;    //TOF on EML2 orbit
    double tof_seml_SEM = refst.tspan_SEM;   //TOF on SEML2 orbit


    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h2, "set title \"Computation coordinates\" ");
    gnuplot_cmd(h3, "set title \"Complementary coordinates\" ");

    //gnuplot_cmd(h2, "set view equal xyz");
    //gnuplot_cmd(h3, "set view equal xyz");
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP_coord = dmatrix(0, 6, 0, 2);
    double **semP_comp  = dmatrix(0, 6, 0, 2);
    switch(coord_type)
    {
    case PSEM:
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
    case NCSEM:
    case VNCSEM:
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

    case PEM:
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
    case NCEM:
    case VNCEM:
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

    //----------------------------------------------------------
    //Complementary coordinate type
    //----------------------------------------------------------
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

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    int max_grid    = 50000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? 3*man_grid_size_t:max_grid;

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    double **y_traj  = dmatrix(0, 41, 0, traj_grid_2);
    double  *t_traj  = dvector(0, traj_grid_2);

    double **y_traj_comp  = dmatrix(0, 41, 0, traj_grid_2);
    double  *t_traj_comp  = dvector(0, traj_grid_2);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_2);
    double *t_traj_n    = dvector(0, traj_grid_2);
    //TBC
    double **yma        = dmatrix(0, 41, 0, traj_grid_2);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================

    //====================================================================================
    // 4.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_NCEM   = dmatrix(0, 5, 0, man_grid_2);
    double **y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_2);
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_2);
    double **y_man_comp   = dmatrix(0, 5, 0, man_grid_2);
    double *t_man_EM      = dvector(0, man_grid_2);
    double *t_man_SEM     = dvector(0, man_grid_2);
    double *t_man_coord   = dvector(0, man_grid_2);
    double *t_man_comp    = dvector(0, man_grid_2);


    //====================================================================================
    // 4.2 Compute the initial orbit
    //====================================================================================
    cout << " comprefft3d. Compute the initial orbit..."  << endl;

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.setSi(0.0, 4);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int em_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        orbit_EM.traj_int_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, true);
    else //variable grid
        em_index = orbit_EM.traj_int_var_grid(orbit_EM.getT0()-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, true);

    //---------------------------------------------------------------------
    //To coord_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, em_index, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp,  t_man_comp,  em_index, NCEM, comp_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < em_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][em_index-kman];
        t_traj[kman] = t_man_coord[em_index-kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], em_index+1, (char*)"", "points", "5", "3", 5);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], em_index+1, (char*)"", "points", "5", "3", 5);

    //====================================================================================
    // 4.3 Compute the manifold leg
    //====================================================================================
    double yv[6];
    cout << " comprefft3d. Compute the manifold leg..."  << endl;

    double tf_EM = orbit_SEM.getT0()/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.setSi(PROJ_EPSILON, 4);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_EM.update_ic(orbit_EM.getSi());

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    int man_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        ode78_qbcp(y_man_coord, t_man_coord, orbit_EM.getT0(), tf_EM, orbit_EM.getZ0(), 6, man_grid_2, dcs, NCEM, coord_type);
    else //variable grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.getT0(), tf_EM, orbit_EM.getZ0(), 6, man_grid_2, dcs, NCEM, coord_type, -1);

    //---------------------------------------------------------------------
    //To comp_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_coord, t_man_coord, y_man_comp, t_man_comp, man_index, coord_type, comp_type);


    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][em_index + kman] = y_man_coord[i][kman];
        t_traj[em_index + kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_index+1, (char*)"Man", "points", "1", "3", 4);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], man_index+1, (char*)"", "points", "1", "3", 4);

    //====================================================================================
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    cout << " comprefft3d. Compute the final SEML2 orbit..."  << endl;
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.setSi(0.0, 4);
    orbit_SEM.update_ic(orbit_SEM.getSi());


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        orbit_SEM.traj_int_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, true);
    else //variable grid
        sem_index = orbit_SEM.traj_int_var_grid(orbit_SEM.getT0()+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

    //---------------------------------------------------------------------
    //To coord_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, sem_index, NCSEM, coord_type);
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_comp, t_man_comp, sem_index, NCSEM, comp_type);


    //---------------------------------------------------------------------
    // Save ALL points
    //---------------------------------------------------------------------
    for(int kman = 0; kman <= sem_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][em_index + man_index + kman] = y_man_coord[i][kman];
        t_traj[em_index + man_index + kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"SEML2", "points", "1", "3", 6);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], sem_index+1, (char*)"", "points", "1", "3", 6);

    //====================================================================================
    // Entire size is: no more than 3*man_grid_size_t points or so.
    //====================================================================================
    int final_index = 0;
    if(refst.grid == REF_VAR_GRID)
    {
        final_index = em_index + man_index + sem_index;
        int freq = final_index/(3*man_grid_size_t);
        cout << "----------------------------------"<< endl;
        cout << "Number of MSD points:        "     << endl;
        cout << "At EML:           " << em_index    << endl;
        cout << "During coast arc: " << man_index   << endl;
        cout << "At SEML:          " << sem_index   << endl;
        cout << "Total:            " << final_index << endl;
        cout << "Selecting 1/" << freq << " points" << endl;
        cout << "Subtotal:         " << (int) floor(final_index/freq) << endl;
        cout << "----------------------------------"<< endl;

        //---------------------------------------------------------------------
        // Subset of points in y_traj_comp/t_traj_comp
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= floor(final_index/freq); kman++)
        {
            for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
            t_traj_comp[kman] = t_traj[freq*kman];
        }

        //New final index
        final_index = floor(final_index/freq);

        //---------------------------------------------------------------------
        // Subset of points back in y_traj/t_traj
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= final_index; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
            t_traj[kman] = t_traj_comp[kman];
        }
    }
    else
    {
        final_index = traj_grid_2;
    }

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    cout << " comprefft3d. Initial trajectory, on a grid..."  << endl;
    int mPlot = 100;

    double **ymc        = dmatrix(0, 5, 0, mPlot);
    double *tmc         = dvector(0, mPlot);
    double **ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp    = dvector(0, mPlot);

    double **ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_v       = dvector(0, mPlot*final_index);
    double **ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78_qbcp(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    gnuplot_plot_xyz(h2, ymc_v[0], ymc_v[1],  ymc_v[2], mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);
    //Plot on h2: points
    gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], final_index+1, (char*)"", "points", "1", "2", 7);
    //Plot on h3: lines
    gnuplot_plot_xyz(h3, ymc_comp_v[0], ymc_comp_v[1],  ymc_comp_v[2], mPlot*final_index+1, (char*)"Initial guess", "lines", "1", "2", 7);


    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    cout << " comprefft3d. Differential correction procedure.          "  << endl;
    cout << "-----------------------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78_qbcp(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
    gnuplot_plot_xyz(h2, ymc_v[0], ymc_v[1],  ymc_v[2], mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);
    gnuplot_plot_xyz(h2, y_traj_n[0], y_traj_n[1],  y_traj_n[2], final_index+1, (char*)"", "points", "1", "2", 3);
    //Plot on h3
    gnuplot_plot_xyz(h3, ymc_comp_v[0], ymc_comp_v[1],  ymc_comp_v[2], mPlot*final_index+1, (char*)"Final guess", "lines", "1", "2", 3);

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

        //================================================================================
        // Select between VECLI and J2000
        //================================================================================
        int choice = 0;
        cout << "Enter 0 for VECLI, 1 for J2000: ";
        cin >> choice;
        int coord_int = choice == 0? VECLI:J2000;
        int fwrk_int  = choice == 0? F_ECLI:F_J2000;

        //Focus on SEM for normalized units
        //if(choice == 1) changeDCS(SEML, F_SEM);

        //================================================================================
        // Initialize SPICE kernerls & VF
        //================================================================================
        gnuplot_ctrl *h4 = gnuplot_init();
        cout << " comprefft3d. Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");
        int shift = 0;


        //================================================================================
        // Search for best fit in  JPL ephemerides
        //================================================================================
        cout << " comprefft3d. Search for best fit in JPL DE430..." << endl;

        //----------------------------------------------------------
        //Get the best fit at t = t_traj_n[shift]. et0 is in seconds
        //----------------------------------------------------------
        double et0, et;
        double tsys0 = t_traj_n[shift];
        qbcp2jpl(tsys0, &et0, coord_type);

        //----------------------------------------------------------
        //Display the value in "D" format
        //----------------------------------------------------------
        ConstSpiceChar  *pictur  = "Wkd Mon DD HH:MN:SC PDT YYYY ::UTC-7";
        SpiceChar output[50];
        et2utc_c (et0, "C", 6, 51, output);
        cout << "comprefft3d. The best fit has been obtained for the following date:" << endl;
        printf ( "%s\n", output);

        cout << "tsys0 = " << tsys0 << endl;
        cout << "et0 = "   << et0   << endl;

        //----------------------------------------------------------
        // Complementary value ot tsys0
        //----------------------------------------------------------
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

        //================================================================================
        // Change of coordinates
        //================================================================================
        cout << " comprefft3d. Change of coordinates syn -> ecliptic..." << endl;
        //---------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //---------------------------------------------------------------------
        double **y_traj_jpl   = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl   = dvector(0, final_index);
        double *et_traj_jpl   = dvector(0, final_index);
        double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl_n  = dvector(0, final_index);
        double **y_jpl_temp   = dmatrix(0, 41, 0, final_index);
        double *t_jpl_temp    = dvector(0, final_index);

        //---------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //---------------------------------------------------------------------
        switch(coord_int)
        {
            case VECLI:
            coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));
            break;

            case J2000:
            coord2eclnstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));//, SEML.ss);
            break;
        }

        //---------------------------------------------------------------------
        //Time in seconds
        //---------------------------------------------------------------------
        for(int p = 0; p <= final_index; p++)
        {
            et_traj_jpl[p] = et0 + (t_traj_n[p]-tsys0)/mean_motion(eph_coord(coord_type));
        }

        //================================================================================
        //Half of the trajectory is set in the other plane
        //================================================================================
        switch(coord_int)
        {
            case VECLI:
            coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));
            break;

            case J2000:
            coord2eclnstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));//, SEML.ss);
            break;
        }

        //----------------------------------------------------------
        //Store the second half: at the beginning or at the end, depending on the computation coordinates
        //----------------------------------------------------------
        int begind = 0, endind = 0;
        switch(eph_coord(coord_type))
        {
                case VEM:
                begind = floor(final_index/2)-6; //-6 for 16, 3D
                endind = final_index;
                break;

                case VSEM:
                begind = 0;
                endind = floor(final_index/2);
                break;
        }

        for(int p = begind; p <= endind; p++)
        {
            for(int i = 0; i <6; i++) y_traj_jpl[i][p] =  y_traj_jpl_n[i][p];
        }

        //================================================================================
        // Plotting Moon +  Earth + L2
        //================================================================================
        cout << " comprefft3d. Plotting Moon +  Earth + L2..." << endl;

        //---------------------------------------------------------------------
        //Position of the Moon +  Earth + L2
        //---------------------------------------------------------------------
        double **y_earth_spice = dmatrix(0, 5, 0, final_index);
        double **y_moon_spice  = dmatrix(0, 5, 0, final_index);
        double **y_l2_spice    = dmatrix(0, 5, 0, final_index);
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

        //---------------------------------------------------------------------
        //Position of the Sun +  Earth + Initial guess in coord_type coordinates
        //---------------------------------------------------------------------
        ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
        gnuplot_plot_xyz(h2, y_jpl_temp[0], y_jpl_temp[1],  y_jpl_temp[2], final_index+1, (char*) "EARTH", "points", "3", "2", 8);

        ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
        gnuplot_plot_xyz(h2, y_jpl_temp[0], y_jpl_temp[1],  y_jpl_temp[2], final_index+1, (char*) "MOON", "points", "4", "2", 8);

        ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, tsys0, eph_coord(coord_type));
        gnuplot_plot_xyz(h2, y_jpl_temp[0], y_jpl_temp[1], y_jpl_temp[2], final_index+1, (char*) "SEML2", "points", "5", "2", 8);


        //================================================================================
        // Plotting Initial Guess
        //================================================================================
        printf("Press ENTER to plot the Initial Guess\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Plotting Initial Guess..." << endl;
        //---------------------------------------------------------------------
        //Initial trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl[k], t_traj_jpl[k+1], yv, 6, mPlot, fwrk_int, coord_int, coord_int);

            //Back to coord_type coordinates
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));
                break;

                case J2000:
                ecln2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));//, SEML.ss);
                break;
            }

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);


            //----------------------------------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------------------------------
            //ECLIPTIC -> comp_type
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));
                break;

                case J2000:
                ecln2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));//, SEML.ss);
                break;
            }


            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

        }


        //================================================================================
        // 6.5 Differential correction
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct_variable_time(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, yma, 42, final_index, coord_int, isPlotted, h4);
        //multiple_shooting_direct(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, yma, 42, final_index, coord_int, isPlotted, h4);
        //Plot
        gnuplot_plot_xyz(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], y_traj_jpl_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //================================================================================
        //Final trajectory on lines, segment by segment, saved on txt file
        //================================================================================
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, fwrk_int, coord_int, coord_int);

            //            //----------------------------------------------------------------------------
            //            //To coord_type coordinates, non normalized
            //            //----------------------------------------------------------------------------
            //            switch(coord_int)
            //            {
            //                case VECLI:
            //                break;
            //
            //                case J2000:
            //                ecln2syndpos_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, et0, tsys0, eph_coord(coord_type));
            //                break;
            //            }


            //----------------------------------------------------------------------------
            //Store, in km, km/s and julian date
            //----------------------------------------------------------------------------
            for(int p = 0; p <= mPlot; p++)
            {
                tmc_v[k*mPlot + p] = tmc_comp[p]; //unitim_c(tmc_comp[p], "TDB", "JDTDB");
                //position in km, velocity not taken into account for now.
                ymc_v[0][k*mPlot + p] = ymc_comp[0][p];
                ymc_v[1][k*mPlot + p] = ymc_comp[1][p];
                ymc_v[2][k*mPlot + p] = ymc_comp[2][p];
                //for(int i = 3; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p]*SEML.ss.a*SEML.ss.n; //velocity in km/s
            }


            //----------------------------------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------------------------------
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));
                break;

                case J2000:
                ecln2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));//, SEML.ss);
                break;
            }


            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);


            ////----------------------------------------------------------------------------
            ////To comp_type coordinates, non normalized
            ////----------------------------------------------------------------------------
            //switch(coord_int)
            //            {
            //                case VECLI:
            //                break;
            //
            //                case J2000:
            //                ecln2syndpos_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, et0, tsys0, eph_coord(comp_type));
            //                break;
            //            }
            //
            //            //----------------------------------------------------------------------------
            //            //Store, in km, km/s and julian date
            //            //----------------------------------------------------------------------------
            //            for(int p = 0; p <= mPlot; p++)
            //            {
            //                tmc_comp_v[k*mPlot + p] = unitim_c(tmc_comp[p], "TDB", "JDTDB");                                    //seconds in TDB to Julian Date
            //                //position in km, velocity not taken into account for now.
            //                ymc_comp_v[0][k*mPlot + p] = -ymc[0][p];
            //                ymc_comp_v[1][k*mPlot + p] = -ymc[1][p];
            //                ymc_comp_v[2][k*mPlot + p] = +ymc[2][p];
            //                //for(int i = 3; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p]*SEML.ss.a*SEML.ss.n; //velocity in km/s
            //            }


            //----------------------------------------------------------------------------
            //To comp_type coordinates
            //----------------------------------------------------------------------------
            //ECLIPTIC -> comp_type
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));
                break;

                case J2000:
                ecln2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));//, SEML.ss);
                break;
            }

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

        }

        //--------------------------------------------------------------------------------
        //Store in files
        //--------------------------------------------------------------------------------
        //string fileEM  = "jpltraj_EM.xyz";
        string fileSEM = "jpltraj.xyz";
        //Save
        //gnuplot_fplot_txyzv(tmc_comp_v, ymc_comp_v, mPlot*final_index+1, fileEM.c_str(), "w");
        gnuplot_fplot_txyzv(tmc_v, ymc_v, mPlot*final_index+1, fileSEM.c_str(), "w");

        //--------------------------------------------------------------------------------
        //Gnuplot window
        //--------------------------------------------------------------------------------
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
        gnuplot_close(h4);
    }

    //------------------------------------------------------------------------------------
    //Gnuplot window
    //------------------------------------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);
    gnuplot_close(h3);


    //====================================================================================
    // Free the data containers
    //====================================================================================
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj, 0, traj_grid_2);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_comp, 0, traj_grid_2);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_2);
    free_dvector(t_traj_n, 0, traj_grid_2);
    free_dmatrix(yma, 0, 41, 0, traj_grid_2);

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

    return 0;
}


//----------------------------------------------------------------------------------------
// Text format, read
//----------------------------------------------------------------------------------------
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
        double **y1  = dmatrix(0, 5, 0, count0);
        double *t1   = dvector(0, count0);
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
        double **y2  = dmatrix(0, 5, 0, count0);
        double *t2   = dvector(0, count0);

        furnsh_c("spice/kernels/metakernel.furnsh");
        ecln2syndpos_vec(y1, t1, y2, t2, count0, VSEM);

        // Rewrite the data
        gnuplot_fplot_txyzv(t2, y2, count0+1, "jpltraj_SEM.xyz", "w");


        //--------------------------------------------------------------------------------
        // COC towards VEM + store data
        //--------------------------------------------------------------------------------
        ecln2syndpos_vec(y1, t1, y2, t2, count0, VEM);
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
        ecln2syndpos_vec(y1, t1, y2, t2, count0, VSEM);
        gnuplot_fplot_txyzv(t2, y2, count0+1, "SEML1_SEM.xyz", "w");


        //--------------------------------------------------------------------------------
        // Position of the SEML2 point: may be wrong!
        //--------------------------------------------------------------------------------
        for(int i = 0; i<= count0; i++)
        {

            //Position
            spkez_c (m1name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RS, &lt);
            spkez_c (m2name(VSEM), t1[i], DEFFRAME, "NONE", SSB, RB, &lt);
            spkez_c (EARTH, t1[i], DEFFRAME, "NONE", SSB, RE, &lt);

            //Position of SEML2
            y1[0][i] = (RS[0] - RB[0])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[0];
            y1[1][i] = (RS[1] - RB[1])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[1];
            y1[2][i] = (RS[2] - RB[2])*(SEML.cs_sem_l2.mu - 1 - SEML.cs_sem_l2.gamma) - RE[2];
        }

        //Towards VSEM coordinates
        ecln2syndpos_vec(y1, t1, y2, t2, count0, VSEM);
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
