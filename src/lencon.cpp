#include "lencon.h"



//========================================================================================
//
//          Computation of the CMU about EML2
//
//========================================================================================
/**
 *  \brief Computes initial conditions in the Planar Center-Unstable Manifold about EML2, in the QBCP model.
 *         The initial conditions (IC) are computed in a three-dimensional box: one dimension for the starting time,
 *         two dimensions for the parameterization of the Center Manifold (s1 and s3 coordinates). The RCM coordinate s5 along the unstable
 *         direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates applied on the unstable direction s5.
 *  \param tmin_CMU_EM:    the minimum starting time (in EM units) in the IC box.
 *  \param tmax_CMU_EM:    the maximum starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param s1_MIN_CMU_RCM: the minimum s1 value (in RCM coordinates) in the IC box.
 *  \param s1_MAX_CMU_RCM: the maximum s1 value (in RCM coordinates) in the IC box.
 *  \param s3_MIN_CMU_RCM: the minimum s3 value (in RCM coordinates) in the IC box.
 *  \param s3_MAX_CMU_RCM: the maximum s3 value (in RCM coordinates) in the IC box.
 *  \param s1_grid_size:   the number of points on the s1 grid in the IC box.
 *  \param s3_grid_size:   the number of points on the s3 grid in the IC box.
 *  \param CM_TFC:         the Fourier-Taylor representation of the Center-Unstable Manifold about EML2, in TFC coordinates
 *  \param Mcoc:           the Matrix that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param MIcoc:          the invarse of Mcoc.
 *  \param Vcoc:           the vector that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form "plot/QBCP/EM/L2/cu_order_16.bin"
 **/
int compute_grid_CMU_EM(double dist_to_cm,
                        double tmin_CMU_EM,
                        double tmax_CMU_EM,
                        int t_grid_size,
                        double s1_MIN_CMU_RCM,
                        double s1_MAX_CMU_RCM,
                        double s3_MIN_CMU_RCM,
                        double s3_MAX_CMU_RCM,
                        int s1_grid_size,
                        int s3_grid_size,
                        vector<Oftsc> &CM_TFC,
                        matrix<Ofsc>  &Mcoc,
                        vector<Ofsc>  &Vcoc,
                        bool isPar)

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
                RCMtoNCbyTFC(sti, grid_t_EM[kt], SEML.us.n, OFTS_ORDER, OFS_ORDER, 5, CM_TFC, Mcoc, Vcoc, yvu, 1);

                //Save
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

    //------------------------------------------
    //Store values
    //------------------------------------------
    writeCU_bin(init_state_CMU_NCEM, init_state_CMU_RCM,
                grid_t_EM, s1_grid_size, s3_grid_size,
                t_grid_size, OFTS_ORDER, TYPE_CU, SEML.li_SEM);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_d4tensor(init_state_CMU_NCEM, 0, 5, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);
    free_d4tensor(init_state_CMU_RCM,  0, 4, 0, t_grid_size, 0, s1_grid_size, 0, s3_grid_size);

    return 1;
}

/**
 *  \brief Computes initial conditions in the 3D Center-Unstable Manifold about EML2, in the QBCP model.
 *         The initial conditions (IC) are computed in a 5-dimensional box: one dimension for the starting time,
 *         four dimensions for the parameterization of the Center Manifold (s1 to s4 coordinates). The RCM coordinate s5 along the unstable
 *         direction is fixed to dist_to_cm.
 *
 *  \param dist_to_cm:     the value in RCM coordinates applied on the unstable direction s5.
 *  \param tlim_CMU_EM:    the min/max starting time (in EM units) in the IC box.
 *  \param t_grid_size:    the number of points on the time grid in the IC box.
 *  \param si_LIM_CMU_RCM: the min/max of s1, s2, s3, s4 values (in RCM coordinates) in the IC box.
 *  \param si_grid_size:   the number of points on the  s1, s2, s3, s4 values  grids in the IC box.
 *  \param CM_TFC:         the Fourier-Taylor representation of the Center-Unstable Manifold about EML2, in TFC coordinates
 *  \param Mcoc:           the Matrix that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param MIcoc:          the invarse of Mcoc.
 *  \param Vcoc:           the vector that appears in the TFC to NC change of coordinates: CM_NC = Mcoc*CM_TFC + Vcoc.
 *  \param isPar:          if TRUE, the computation is parallelized.
 *
 *
 * The output data are saved in a binary file of the form "plot/QBCP/EM/L2/cu_3d_order_16.bin"
 **/
int compute_grid_CMU_EM_3D(double dist_to_cm,
                           double *tlim_CMU_EM,
                           int t_grid_size,
                           double si_LIM_CMU_RCM[4][2],
                           int *si_grid_size,
                           vector<Oftsc> &CM_TFC,
                           matrix<Ofsc>  &Mcoc,
                           vector<Ofsc>  &Vcoc,
                           bool isPar)
{
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

    cout << "grid_t_EM[0] = " << grid_t_EM[0] << endl;

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
    initCU_bin_3D(si_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);

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
        appTimeCU_bin_3D(grid_t_EM, kt, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);


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
                        RCMtoNCbyTFC(sti, grid_t_EM[kt], SEML.us.n, OFTS_ORDER, OFS_ORDER, 5, CM_TFC, Mcoc, Vcoc, yvu, 1);

                        //Save
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

                    //------------------------------------------
                    //Store values
                    //------------------------------------------
                    writeCU_bin_3D(init_state_CMU_NCEM, init_state_CMU_RCM,
                                   si_grid_size, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);
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



//========================================================================================
//
//         Projection on the CM/CMS/CMU of SEML2
//
//========================================================================================
/**
 *  \brief Integrates the central-unstable legs from a discrete set of unstable directions obtained using the routine compute_grid_CMU_EM.
 *         Then, each point on the integration grid is projected on the Center Manifold CM_SEM_NC about SEML2. The best solution (minimum distance of projection)
 *         is stored.
 *
 *  \param tmax_on_manifold_EM: the maximum integration time on each leg, in EM units.
 *  \param t_grid_size_x:       the number of points on the time grid in the IC box. If -1, the value used in compute_grid_CMU_EM is used.
 *  \param s1_grid_size_x:      the number of points on the s1 grid in the IC box. If -1, the value used in compute_grid_CMU_EM is used.
 *  \param s3_grid_size_x:      the number of points on the s3 grid in the IC box. If -1, the value used in compute_grid_CMU_EM is used.
 *  \param man_grid_size:       the number of points on each manifold leg.
 *  \param NsortMin:            the number of best solutions that are kept
 *  \param nod:                 the number of dimensions on which the distance of projection is computed (usually either 3 (the physical distance) or 6 (the whole phase space)).
 *  \param isPar:               if TRUE, the computation is parallelized.
 *  \param ynormMax:            the maximum norm in NCSEM coordinates for which a given state on the integration grid is projected on CM_SEM_NC More precisely: for a given state
 *                              y along the manifold leg, iv norm(y, 3) < ynormMax, the state is projected. Otherwise, it is considered too far away from SEML2 to be a good candidate for projection.
 *  \param snormMax:            the maximum norm in RCM SEM coordinates for which a given projection state on the CM of SEML2 (CM_SEM_NC) is computed back in NCSEM coordinates.
 *                              More precisely, for a given state y in NCSEM coordinates, the result of the projection on CM_SEM_NC gives a state sproj in RCM SEM coordinates.
 *                              if norm(sproj, 4) < snormMax, the computation yproj = CM_SEM_NC(sproj, t) is performed. Otherwise, the state sproj is considered too far away from the RCM origin to be a
 *                              good candidate - it is out of the domain of practical convergence of CM_SEM_NC.
 *
 * The output data are saved in a binary file of the form "plot/QBCP/EM/L2/projcu_order_16.bin", and "plot/QBCP/EM/L2/sortprojcu_order_16.bin" for the NsortMin best solutions.
 **/
int int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM,
                              int t_grid_size_x,
                              int s1_grid_size_x,
                              int s3_grid_size_x,
                              int man_grid_size,
                              int NsortMin,
                              int nod,
                              int isPar,
                              double ynormMax,
                              double snormMax)
{
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
               s1_grid_size, s3_grid_size, t_grid_size, OFTS_ORDER, TYPE_CU, SEML.li_SEM);

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
    // 2.2. Initialize tools for the projection phase
    //====================================================================================

    //--------------------------------------
    // Center-manifold around SEML2
    //--------------------------------------
    vector<Oftsc>  CM_SEM_NC;     //center manifold in NC coordinates
    vector<Oftsc>  CM_SEM_TFC;    //center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);


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
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ, SEML.li_SEM);
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
                // PB: when Release + I_NCSEM!!
                //---------------------------------------------------------------------
                int status = ode78(y_man_NCSEM, t_man_SEM, tv, tv+tmax_on_manifold_EM, yv, 6, man_grid_size, I_NCEM, NCEM, NCSEM);

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



                        //If the current state is close enough to SEMLi
                        if(y_man_norm_NCSEM < ynormMax)
                        {
                            // Projection on the center manifold
                            NCprojCCMtoCM(yv, tv, SEML_SEM.us.n, sproj, CM_SEM_TFC, MIcoc_SEM, Vcoc_SEM);

                            //If the projection is close enough to SEMLi
                            if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
                            {
                                //yvproj_NCSEM = W(sproj, tv)
                                RCMtoNCbyTFC(sproj, tv, SEML_SEM.us.n, OFTS_ORDER, OFS_ORDER, 4, CM_SEM_TFC, Mcoc_SEM, Vcoc_SEM, yvproj_NCSEM, 1);

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
int int_proj_CMU_EM_on_CM_SEM_3D(double tmax_on_manifold_EM,
                                 int man_grid_size, int nod,
                                 int isPar, double ynormMax, double snormMax)
{
    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, si_grid_size[4], offset;
    offset = getLenghtCU_bin_3D(si_grid_size, &t_grid_size,
                                OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);

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
    // 2.2. Initialize tools for the projection phase
    //====================================================================================

    //--------------------------------------
    // Center-manifold around SEML2
    //--------------------------------------
    vector<Oftsc>  CM_SEM_NC;     //center manifold in NC coordinates
    vector<Oftsc>  CM_SEM_TFC;    //center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);


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
    string filename = filenameCUM(OFTS_ORDER, TYPE_MAN_PROJ_3D, SEML.li_SEM);
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
        offset = readTCU_bin_3D(offset, init_time_grid_EM,
                                kt, OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);
        cout << "init_time_grid_EM[kt] = " << init_time_grid_EM[kt] << endl;

        for(int ks2 = 0; ks2 <= si_grid_size[1]; ks2++)
        {
            for(int ks4 = 0; ks4 <= si_grid_size[3]; ks4++)
            {
                for(int ks1 = 0; ks1 <= si_grid_size[0]; ks1++)
                {
                    //----------------------------------------------------------
                    //Read data from file
                    //----------------------------------------------------------
                    offset = readCU_bin_3D(offset, init_state_CMU_NCEM,
                                           init_state_CMU_RCM, si_grid_size,
                                           OFTS_ORDER, TYPE_CU_3D, SEML.li_SEM);

                    //----------------------------------------------------------
                    //Most inner loop is parallelized
                    //----------------------------------------------------------
                    #pragma omp parallel for if(isPar)  shared(index)
                    for(int ks3 = 0; ks3 <= si_grid_size[2]; ks3++)
                    {

                        //====================================================================================
                        // 3.1. Integration of the manifold leg
                        //====================================================================================
                        //---------------------------------------------------------------------
                        //Initialize the initial conditions (both NC and RCM coordinates)
                        //---------------------------------------------------------------------
                        double yv[6], tv;
                        tv  = init_time_grid_EM[kt];
                        for(int i = 0; i < 6; i++) yv[i] = init_state_CMU_NCEM[i][ks3];

                        //---------------------------------------------------------------------
                        //Local variables to store the manifold leg
                        //---------------------------------------------------------------------
                        double **y_man_NCSEM = dmatrix(0, 5, 0, man_grid_size);
                        double *t_man_SEM    = dvector(0, man_grid_size);


                        //---------------------------------------------------------------------
                        //Integration on man_grid_size+1 fixed grid
                        // PB: when Release + I_NCSEM!!
                        //---------------------------------------------------------------------
                        int status = ode78(y_man_NCSEM, t_man_SEM, tv, tv+tmax_on_manifold_EM, yv, 6, man_grid_size, I_NCEM, NCEM, NCSEM);

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


                                //If the current state is close enough to SEMLi
                                if(y_man_norm_NCSEM < ynormMax)
                                {
                                    // Projection on the center manifold
                                    NCprojCCMtoCM(yv, tv, SEML_SEM.us.n, sproj, CM_SEM_TFC, MIcoc_SEM, Vcoc_SEM);

                                    //If the projection is close enough to SEMLi
                                    if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
                                    {
                                        //yvproj_NCSEM = W(sproj, tv)
                                        RCMtoNCbyTFC(sproj, tv, SEML_SEM.us.n, OFTS_ORDER, OFS_ORDER, 4, CM_SEM_TFC, Mcoc_SEM, Vcoc_SEM, yvproj_NCSEM, 1);

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
                            displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*index++/noe);
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
        }
    }
    return 1;
}



/**
 *  \brief Same int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM, but using the Center-Stable Manifold about SEML2 as a target.
 *         Note that, for now, the coordinates along the stable direction is fixed to a small arbitrary value.
 *         Since this routine does not yield substantially better results thant int_proj_CMU_EM_on_CM_SEM(double tmax_on_manifold_EM, it is not currently used.
 *
 **/
int int_proj_CMU_EM_on_CMS_SEM(double tmax_on_manifold_EM,
                               int ofts_order,
                               int t_grid_size_x,
                               int s1_grid_size_x,
                               int s3_grid_size_x,
                               int man_grid_size,
                               int NsortMin,
                               int isPar,
                               double ynormMax,
                               double snormMax)
{
    //====================================================================================
    // 1. Get initial condition in the center-unstable manifold from a data file
    //====================================================================================
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int t_grid_size, s1_grid_size, s3_grid_size;
    getLenghtCU_bin(&s1_grid_size, &s3_grid_size, &t_grid_size,
                    ofts_order, TYPE_CU, SEML.li_SEM);

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
               s1_grid_size, s3_grid_size, t_grid_size, ofts_order, TYPE_CU, SEML.li_SEM);

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
    // 2.2. Initialize tools for the projection phase
    //====================================================================================

    //--------------------------------------
    // Center-stable manifold around SEML2
    //--------------------------------------
    vector<Oftsc>  CMS_SEM_NC;     //center-stable manifold in NC coordinates
    vector<Oftsc>  CMS_SEM_TFC;    //center-stable manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);


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
    // projection by default is arbitrarily big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 2.4. Reset the data file (projection)
    //====================================================================================
    //Filename
    string filename = filenameCUM(ofts_order, TYPE_MAN_PROJ, SEML.li_SEM);
    //Open whitout append datafile and therefore erase its content.
    fstream filestream;
    filestream.open (filename.c_str(), ios::binary | ios::out);
    filestream.close();


    //====================================================================================
    // 3. Loop: only the inner loop is parallelized, since
    //    open_mp doest not allow nested // loops by default
    //====================================================================================
    int index  = 0;
    COMPLETION = 0;
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
                //---------------------------------------------------------------------
                int status = ode78(y_man_NCSEM, t_man_SEM, tv, tv+tmax_on_manifold_EM, yv, 6, man_grid_size, I_NCEM, NCEM, NCSEM);

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

                        //If the current state is close enough to SEMLi
                        if(y_man_norm_NCSEM < ynormMax)
                        {
                            // Projection on the center manifold
                            NCprojCCMtoCUS(yv, tv, SEML_SEM.us.n, sproj, CMS_SEM_TFC, PROJ_EPSILON_SEM, MIcoc_SEM, Vcoc_SEM);

                            //If the projection is close enough to SEMLi
                            if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
                            {
                                //yvproj_NCSEM = W(sproj, tv)
                                RCMtoNCbyTFC(sproj, tv, SEML_SEM.us.n, OFTS_ORDER, OFS_ORDER, 5, CMS_SEM_TFC, Mcoc_SEM, Vcoc_SEM, yvproj_NCSEM, 1);

                                //Back to SEM coordinates
                                NCSEMmtoSEMm(tv, yv, yv_SEM, &SEML_SEM);
                                NCSEMmtoSEMm(tv, yvproj_NCSEM, yvproj_SEM, &SEML_SEM);

                                //Back to SEM coordinates with velocity
                                SEMmtoSEMv(tv, yv_SEM, yv_VSEM, &SEML_SEM);
                                SEMmtoSEMv(tv, yvproj_SEM, yvproj_VSEM, &SEML_SEM);

                                //Distance of projection in SEM coordinates
                                proj_dist_SEM = 0.0;
                                for(int i = 0; i < 6; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
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
                // 3.3. Save a first batch of outputs
                //====================================================================================
                #pragma omp critical
                {
                    //Initial position in SEM coordinates
                    for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][0];
                    NCSEMmtoSEMm(t_man_SEM[0], yv, yv_SEM, &SEML_SEM);
                    for(int i = 0; i < 6; i++) init_state_CMU_SEM[i][kt][ks1][ks3] = yv_SEM[i];

                    //Same in NCEM coordinates
                    NCSEMmtoNCEMm(t_man_SEM[0], yv, yvEM, &SEML_SEM);
                    for(int i = 0; i < 6; i++) init_state_CMU_NCEM_0[i][kt][ks1][ks3] = yvEM[i];

                    //Initial time in EM coordinates
                    t0_min_EM[ksort]= t_man_SEM[0]/SEML.us_em.ns;        //in EM coordinates
                }


                //====================================================================================
                // 3.4. Refine the solutions, if desired. Not very convincing. Abandoned for now
                //====================================================================================

                //====================================================================================
                // 3.5. Save the rest of the outputs
                //====================================================================================
                #pragma omp critical
                {
                    //----------------------------------------------------------
                    //Save
                    //----------------------------------------------------------
                    //Minimum projection distance
                    min_proj_dist_tens_SEM[kt][ks1][ks3] = min_proj_dist_SEM;

                    //For sorting
                    indexMin[ksort] = (int) kmin;
                    ks1Min[ksort]   = (int) ks1;
                    ks3Min[ksort]   = (int) ks3;
                    ktMin[ksort]    = (int) kt;
                    distMin[ksort]  = min_proj_dist_SEM;
                    tf_min_EM[ksort]= t_man_SEM[kmin]/SEML.us_em.ns;    //in EM coordinates


                    //----------------------------------------------------------
                    //Open datafile
                    //----------------------------------------------------------
                    writeIntProjCU_bin(filename, init_time_grid_EM, init_state_CMU_NCEM,
                                       init_state_CMU_SEM, init_state_CMU_RCM, final_state_CMU_SEM, projected_state_CMU_SEM,
                                       projected_state_CMU_RCM, min_proj_dist_SEM, dv_at_projection_SEM, t_man_SEM, kmin, ks1, ks3, kt);

                    //----------------------------------------------------------
                    //Display completion
                    //----------------------------------------------------------
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
    filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

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
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE(int cont_grid_size,
        int coord_type,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        int isPlanar)

{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);
    //int coordsys = default_framework(coord_type);

    //====================================================================================
    // 2. Time on the orbits
    //====================================================================================
    double tof_seml_SEM = 5*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 3. Init the gnuplot windows
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP_coord = dmatrix(0, 6, 0, 2);
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


    //====================================================================================
    // 4. Init the data containers
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
    // 5. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCMS_SEM_TFC(6, 5, 5, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCMS_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 6. Getting back the data
    //====================================================================================
    int type = isPlanar? TYPE_MAN_PROJ:TYPE_MAN_PROJ_3D;
    string filename = filenameCUM(OFTS_ORDER, type, SEML.li_SEM);

    readIntProjCU_bin(filename, t0_EM, tf_EM, s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM, s5_CMU_EM,
                      pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM, sortId);

    //====================================================================================
    // 7. Select given intervals for the inputs
    //====================================================================================
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;

    cout << "----------------------------------------" << endl;
    cout << "   User-defined interval of research    " << endl;
    cout << "----------------------------------------" << endl;
    cout << "Enter a value for s1_CMU_EM_MIN: ";
    cin >> s1_CMU_EM_MIN;
    cout << "Enter a value for s1_CMU_EM_MAX: ";
    cin >> s1_CMU_EM_MAX;
    cout << "Enter a value for s3_CMU_EM_MIN: ";
    cin >> s3_CMU_EM_MIN;
    cout << "Enter a value for s3_CMU_EM_MAX: ";
    cin >> s3_CMU_EM_MAX;

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
    // 8. The best solution in the subselection will serve as the first guess.
    //====================================================================================
    int kpos = sortId_R[0];


    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;
    cout << "Estimated error at patch point (km):       "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]*SEML.cs_sem.cr3bp.L << endl;
    cout << "Estimated error at patch point (SEMSU):    "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]                         << endl;
    cout << "-------------------------------------------"  << endl;

    //====================================================================================
    // 8.1 Initialize local variables: EM
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
    Ofsc orbit_ofs_EM(OFS_ORDER);
    OdeStruct driver_EM;
    SingleOrbit orbit_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
    //Init routine
    init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
               &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
               1, t0_EM_R[kpos], tf_EM_R[kpos], SEML.us_em.T/5, &driver_EM, &SEML_EM);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, st_EM, t0_EM_R[kpos]);

    //====================================================================================
    // 8.2 Initialize local variables: SEM
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
    Ofsc orbit_ofs_SEM(OFS_ORDER);
    OdeStruct driver_SEM;
    SingleOrbit orbit_SEM;

    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
    //Init routine
    init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
               &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
               1, t0_SEM, t0_SEM+tof_seml_SEM, SEML.us_sem.T, &driver_SEM, &SEML_SEM);


    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_SEM, st_SEM, t0_SEM);


    //====================================================================================
    // 8.3 Display info
    //====================================================================================
    coutmp();
    cout << "The initial orbit is:                      "  << endl;
    cout << "- orbit_EM.si = (" << orbit_EM.si[0];
    cout << ", "<< orbit_EM.si[1] << ", " << orbit_EM.si[2];
    cout << ", "<< orbit_EM.si[3] << ", " << orbit_EM.si[4] << ")" << endl;
    cout << "- orbit_EM.t0 = " << orbit_EM.t0/SEML.us_em.T << " T" << endl;
    cout << "-------------------------------------------"  << endl;
    cout << "The final orbit is:                        "  << endl;
    cout << "- orbit_SEM.si = (" << orbit_SEM.si[0];
    cout << ", "<< orbit_SEM.si[1] << ", " << orbit_SEM.si[2];
    cout << ", "<< orbit_SEM.si[3] << ", " << orbit_SEM.si[4] << ")"  << endl;
    cout << "- orbit_SEM.t0 = " << orbit_SEM.t0/SEML.us_sem.T << " T" << endl;
    coutlp();
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //====================================================================================
    // 8.4 Multiple shooting procedure
    //====================================================================================
    if(isPlanar)
    {
        //====================================================================================
        // 8.4.1 First step: continuation with free final time, to decrease the stable component at SEML2
        //====================================================================================
        cout << "First step: continuation with free final time, to decrease the stable component at SEML2" << endl;
        crefvtplan(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, cont_grid_size, h2);

        //====================================================================================
        // 8.4.2 Display info
        //====================================================================================
        coutmp();
        cout << "The initial orbit is:                      "  << endl;
        cout << "- orbit_EM.si = (" << orbit_EM.si[0];
        cout << ", "<< orbit_EM.si[1] << ", " << orbit_EM.si[2];
        cout << ", "<< orbit_EM.si[3] << ", " << orbit_EM.si[4] << ")" << endl;
        cout << "- orbit_EM.t0 = " << orbit_EM.t0/SEML.us_em.T << " T" << endl;
        cout << "-------------------------------------------"  << endl;
        cout << "The final orbit is:                        "  << endl;
        cout << "- orbit_SEM.si = (" << orbit_SEM.si[0];
        cout << ", "<< orbit_SEM.si[1] << ", " << orbit_SEM.si[2];
        cout << ", "<< orbit_SEM.si[3] << ", " << orbit_SEM.si[4] << ")"  << endl;
        cout << "- orbit_SEM.t0 = " << orbit_SEM.t0/SEML.us_sem.T << " T" << endl;
        coutlp();
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //====================================================================================
        // 8.4.3 Second step: refinement of the complete trajectory
        //====================================================================================
        cout << "Second step: refinement of the complete trajectory" << endl;
        //ref_CMU_EM_to_CM_SEM_MSD_COMP(msd_grid_size, coord_type, orbit_EM, orbit_SEM);
        ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID(coord_type, orbit_EM, orbit_SEM);
        //ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID_TEST(coord_type, orbit_EM, orbit_SEM);
        //srefvtplan(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }
    else
    {
        //ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, cont_grid_size, h2);
        cout << "3D case is not handled." << endl;
    }



    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    return 0;
}


/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_2(int man_grid_size,
                                    int coord_type,
                                    vector<Oftsc> &CM_EM_NC,
                                    vector<Oftsc> &CM_EM_TFC,
                                    matrix<Oftsc> &DCM_EM_TFC,
                                    matrix<Ofsc>  &Mcoc_EM,
                                    matrix<Ofsc>  &MIcoc_EM,
                                    vector<Ofsc>  &Vcoc_EM,
                                    int isPlanar,
                                    int isSaved,
                                    int isUserDefined)
{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);
    //int coordsys = default_framework(coord_type);

    //====================================================================================
    // 2. Time on the orbits
    //====================================================================================
    double tof_seml_SEM = 5*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 3. Init the gnuplot windows
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP_coord = dmatrix(0, 6, 0, 2);
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


    //====================================================================================
    // 4. Init the data containers
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
    // 5. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCMS_SEM_TFC(6, 5, 5, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCMS_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 6. Getting back the data
    //====================================================================================
    int type = isPlanar? TYPE_MAN_PROJ:TYPE_MAN_PROJ_3D;
    string filename = filenameCUM(OFTS_ORDER, type, SEML.li_SEM);

    readIntProjCU_bin(filename, t0_EM, tf_EM, s1_CMU_EM, s2_CMU_EM, s3_CMU_EM, s4_CMU_EM, s5_CMU_EM,
                      pmin_dist_SEM, s1_CM_SEM, s2_CM_SEM, s3_CM_SEM, s4_CM_SEM, sortId);

    //====================================================================================
    // 7. Select given intervals for the inputs
    //====================================================================================
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;

    cout << "----------------------------------------" << endl;
    cout << "   User-defined interval of research    " << endl;
    cout << "----------------------------------------" << endl;
    cout << "Enter a value for s1_CMU_EM_MIN: ";
    cin >> s1_CMU_EM_MIN;
    cout << "Enter a value for s1_CMU_EM_MAX: ";
    cin >> s1_CMU_EM_MAX;
    cout << "Enter a value for s3_CMU_EM_MIN: ";
    cin >> s3_CMU_EM_MIN;
    cout << "Enter a value for s3_CMU_EM_MAX: ";
    cin >> s3_CMU_EM_MAX;

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

    //--------------------------------------
    // User-defined number of continuation steps
    //--------------------------------------
    int cont_steps_MAX = 0;
    cout << "-----------------------------------------" << endl;
    cout << "User-defined number of continuation steps" << endl;
    cout << "-----------------------------------------" << endl;
    cout << "Enter a value for cont_steps_MAX: ";
    cin >> cont_steps_MAX;

    //====================================================================================
    // 8. The best solution in the subselection will serve as the first guess.
    //====================================================================================
    int kpos = sortId_R[0];


    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;
    cout << "Estimated error at patch point (km):       "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]*SEML.cs_sem.cr3bp.L << endl;
    cout << "Estimated error at patch point (SEMSU):    "  << endl;
    cout <<  pmin_dist_SEM_R[kpos]                           << endl;
    cout << "-------------------------------------------"  << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // 8.1 Initialize local variables: EM
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
    Ofsc orbit_ofs_EM(OFS_ORDER);
    OdeStruct driver_EM;
    SingleOrbit orbit_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
    //Init routine
    init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
               &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
               1, t0_EM_R[kpos], tf_EM_R[kpos], SEML.us_em.T/5, &driver_EM, &SEML_EM);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, st_EM, t0_EM_R[kpos]);

    //====================================================================================
    // 8.2 Initialize local variables: SEM
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
    Ofsc orbit_ofs_SEM(OFS_ORDER);
    OdeStruct driver_SEM;
    SingleOrbit orbit_SEM;

    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
    //Init routine
    init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
               &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
               1, t0_SEM, t0_SEM+tof_seml_SEM, SEML.us_sem.T, &driver_SEM, &SEML_SEM);


    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_SEM, st_SEM, t0_SEM);


    //====================================================================================
    // 8.3 Multiple shooting procedure
    //====================================================================================
    if(isPlanar)
    {

        //====================================================================================
        // 8.3.1 First step: continuation with free final time, to decrease the stable component at SEML2
        //====================================================================================
        crefvtplan(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);

        //====================================================================================
        // 8.3.2 Second step: contination with fixed final time.
        //====================================================================================
        crefftplan(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC,
                dcs, coord_type, man_grid_size, cont_steps_MAX,
                isSaved, isUserDefined, h2);

        //====================================================================================
        // ANNEX: only one solution
        //====================================================================================
        //srefvtplan(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }
    else
    {
        //ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
        crefft3d(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC,
                dcs, coord_type, man_grid_size, cont_steps_MAX,
                isUserDefined, h2);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    return 0;
}


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
                vector<Oftsc> &CM_EM_NC,
                vector<Oftsc> &CM_EM_TFC,
                matrix<Oftsc> &DCM_EM_TFC,
                matrix<Ofsc>  &Mcoc_EM,
                matrix<Ofsc>  &MIcoc_EM,
                vector<Ofsc>  &Vcoc_EM,
                RefSt refst)
{
    //====================================================================================
    // 1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);


    //====================================================================================
    // 2. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCMS_SEM_TFC(6, 5, 5, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCMS_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


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
    int type = refst.dim == REF_PLANAR? TYPE_MAN_PROJ:TYPE_MAN_PROJ_3D;
    string filename = filenameCUM(OFTS_ORDER, type, SEML.li_SEM);

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
    if(refst.type == REF_CONT)
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
    Ofsc orbit_ofs_EM(OFS_ORDER);
    OdeStruct driver_EM;
    SingleOrbit orbit_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
    //Init routine
    init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC, &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
               1, t0_EM_R[kpos], tf_EM_R[kpos], SEML.us_em.T/5, &driver_EM, &SEML_EM);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, st_EM, t0_EM_R[kpos]);

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
    Ofsc orbit_ofs_SEM(OFS_ORDER);
    OdeStruct driver_SEM;
    SingleOrbit orbit_SEM;

    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn,&SEML_SEM);
    //Init routine
    init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
               &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
               1, t0_SEM, t0_SEM+5*SEML.us_sem.T, SEML.us_sem.T, &driver_SEM, &SEML_SEM);


    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_SEM, st_SEM, t0_SEM);



    //====================================================================================
    // 8 Multiple shooting procedure
    //====================================================================================
    double **semP_coord = dmatrix(0, 6, 0, 2);
    switch(refst.type)
    {
    case REF_CONT:
    case REF_SINGLE:
        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        gnuplot_ctrl *h2;
        h2 = gnuplot_init();

        //----------------------------------------------------------
        //Notable points in SEM system
        //----------------------------------------------------------
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

        //----------------------------------------------------------
        //Multiple shooting procedure
        //----------------------------------------------------------
        srefeml2seml(orbit_EM, orbit_SEM, DCM_EM_TFC, DCMS_SEM_TFC,
                     dcs, coord_type, man_grid_size, refst, h2);
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
void srefeml2seml(SingleOrbit &orbit_EM,
                  SingleOrbit &orbit_SEM,
                  matrix<Oftsc> &DCM_EM_TFC,
                  matrix<Oftsc> &DCMS_SEM_TFC,
                  int dcs,
                  int coord_type,
                  int man_grid_size_t,
                  RefSt refst,
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
        ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);
    else //variable time grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type, man_grid_size_t);


    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }

    cout << "initial point = " << endl;
    cout << y_traj[0][0] << endl;
    cout << y_traj[1][0] << endl;
    cout << y_traj[2][0] << endl;
    cout << y_traj[3][0] << endl;
    cout << y_traj[4][0] << endl;
    cout << y_traj[5][0] << endl;

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
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_index] = yout[i];
    t_traj[man_index] = tout;


    cout << "final point = " << endl;
    cout << y_traj[0][man_index] << endl;
    cout << y_traj[1][man_index] << endl;
    cout << y_traj[2][man_index] << endl;
    cout << y_traj[3][man_index] << endl;
    cout << y_traj[4][man_index] << endl;
    cout << y_traj[5][man_index] << endl;

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

    cout << "New. y_traj[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][0] << endl;

    cout << "New. y_traj[end-1] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][man_grid_size-1] << endl;

    cout << "New. y_traj[end] = " << endl;
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
                            DCM_EM_TFC, DCMS_SEM_TFC,
                            CCM_R_RCM_EM,
                            CCM_R_RCM_SEM,
                            orbit_EM, orbit_SEM,
                            h2, isPlotted, refst.isDirUD, false);

            break;

        case REF_VAR_TIME:
            status = msvt3d(y_traj, t_traj, y_traj_n, t_traj_n,
                            nullvector, 42,
                            man_index, coord_type, inner_prec, true,
                            DCM_EM_TFC, DCMS_SEM_TFC,
                            CCM_R_RCM_EM,
                            CCM_R_RCM_SEM,
                            orbit_EM, orbit_SEM,
                            h2, isPlotted, false);

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
                              DCM_EM_TFC, DCMS_SEM_TFC,
                              CCM_R_RCM_EM,
                              CCM_R_RCM_SEM,
                              orbit_EM, orbit_SEM,
                              h2, isPlotted, refst.isDirUD, false);

            break;

        case REF_VAR_TIME:
            status = msvtplan(y_traj, t_traj, y_traj_n, t_traj_n,
                              nullvector, 42,
                              man_index, coord_type, inner_prec, true,
                              DCM_EM_TFC, DCMS_SEM_TFC,
                              CCM_R_RCM_EM,
                              CCM_R_RCM_SEM,
                              orbit_EM, orbit_SEM,
                              h2, isPlotted, false);

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


    cout << "New. y_traj_n[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj_n[i][0] << endl;

    //======================================================================
    //6.2. If it is a sucess, we go on.
    //======================================================================
    if(status == GSL_SUCCESS)
    {
        //======================================================================
        // Save first entry
        //======================================================================
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);
        double yv[42];
        string filename      = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, orbit_EM.t0/SEML.us_em.T);
        string filename_traj = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF_TRAJ, SEML.li_SEM, orbit_EM.t0/SEML.us_em.T);
        fstream filestream;
        if(refst.isSaved)
        {
            cout << " srefeml2seml. Save first entry.                           "  << endl;
            cout << "-----------------------------------------------------------"  << endl;
            //=================================
            // Main parameters
            //=================================
            filestream.open (filename.c_str(), ios::out);
            //Title
            filestream << "t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
            filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
            filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
            filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
            filestream << endl;
            //Data
            filestream << setprecision(15) <<  setiosflags(ios::scientific);
            filestream << orbit_EM.t0 << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_EM.si[i]  << "  ";
            filestream << orbit_EM.tf << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_SEM.si[i] << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_EM.z0[i]  << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_SEM.z0[i] << "  ";
            filestream << endl;
            filestream.close();


            //=================================
            // Entire trajectory
            //=================================
            filestream.open (filename_traj.c_str(), ios::out | ios::binary);
            double res;
            //Final trajectory on lines, segment by segment
            for(int k = 0; k < man_index; k++)
            {
                //Integration segment by segment
                for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        }


        cout << " srefeml2seml. Continuation procedure.                     "  << endl;
        cout << "-----------------------------------------------------------"  << endl;
        //======================================================================
        //Continuation procedure
        //======================================================================
        double ds = 1e-1, ds0 = 1e-1;
        double dkn = 0.0;
        int kn = 0;

        //An additional condition can be set: if the time is not fixed, we want to decrease the
        //s5 component at SEMLi. Therefore, we add a stopping condition on the loop.
        bool addCondition = refst.time == REF_VAR_TIME? fabs(orbit_SEM.si[4]) > ORBIT_SEM_UNSTABLE_MIN: true;

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
                    ufvarft3d(y_traj_n, t_traj_n, ds, nullvector,
                              orbit_EM, orbit_SEM,
                              man_grid_size, coord_type);
                    break;

                case REF_VAR_TIME:
                    ufvarvt3d(y_traj_n, t_traj_n, ds, nullvector,
                              orbit_EM, orbit_SEM,
                              man_grid_size, coord_type);
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
                                    DCM_EM_TFC, DCMS_SEM_TFC,
                                    CCM_R_RCM_EM,
                                    CCM_R_RCM_SEM,
                                    orbit_EM, orbit_SEM,
                                    h2, isPlotted, refst.isDirUD, false);

                    break;

                case REF_VAR_TIME:
                    status = msvt3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                    nullvector, 42,
                                    man_index, coord_type, inner_prec, false,
                                    DCM_EM_TFC, DCMS_SEM_TFC,
                                    CCM_R_RCM_EM,
                                    CCM_R_RCM_SEM,
                                    orbit_EM, orbit_SEM,
                                    h2, isPlotted, false);

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
                                      DCM_EM_TFC, DCMS_SEM_TFC,
                                      CCM_R_RCM_EM,
                                      CCM_R_RCM_SEM,
                                      orbit_EM, orbit_SEM,
                                      h2, isPlotted, refst.isDirUD, false);

                    break;

                case REF_VAR_TIME:
                    status = msvtplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n,
                                      nullvector, 42,
                                      man_index, coord_type, inner_prec, false,
                                      DCM_EM_TFC, DCMS_SEM_TFC,
                                      CCM_R_RCM_EM,
                                      CCM_R_RCM_SEM,
                                      orbit_EM, orbit_SEM,
                                      h2, isPlotted, false);

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
                filestream << orbit_EM.t0 << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_EM.si[i]  << "  ";
                filestream << orbit_EM.tf << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_SEM.si[i] << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_EM.z0[i]  << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_SEM.z0[i] << "  ";
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
                    ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
                gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);

                if(refst.time == REF_FIXED_TIME)
                {
                    gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
                    gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);
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
            addCondition = refst.time == REF_VAR_TIME? fabs(orbit_SEM.si[4]) > ORBIT_SEM_UNSTABLE_MIN: true;

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
        orbit_EM.si[4] = 0.0;
        orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

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
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        orbit_EM.si[4] = PROJ_EPSILON;
        orbit_EM.tf    = t_traj_n[man_index]/SEML.us_em.ns;

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        //====================================
        // At SEML2. The first 5 RCM components are good.
        // Hence, we need to update:
        // 1. The initial time.
        //====================================
        orbit_SEM.t0    = t_traj_n[man_index];

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);
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

//====================================
// Continuation
//====================================
/**
 *  \brief Continuation of a single of EML2-to-SEML2 connection, between orbit_EM and orbit_SEM.
 *         The Jacobian of the parameterization of the manifolds are contained in  DCM_EM_TFC and DCMS_SEM_TFC.
 *         The patch and final times are free to vary. Hence, the continuation takes
 *         place along a given connection, rather than an entire family. The continuation is forced to make the stable component at SEML2 decrease.
 **/
void crefvtplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;       //TOF on EML2 orbit

    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //====================================================================================
    // 4. Plot the resulting trajectory
    //====================================================================================
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 5.  Differential correction
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "---------------------------------------------"  << endl;
    cout << " crefvtplan. "  << endl;
    cout << " Differential correction with variable times "  << endl;
    cout << "---------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //Free variables
    int nfv = 7*(man_grid_size+1)-4;
    double *nullvector = dvector(0, nfv-1);

    //Maximum number of steps allowed
    int maxStep = 50;

    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
    gnuplot_cmd(h5, "set grid");

    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = floor((double)maxStep/floor(0.1*GP_MAX_TMP_FILES))+1;
    cout << "plotfreq = " << plotfreq << endl;

    //======================================================================
    //5.1. First step of the continuation procedure.
    //======================================================================
    int status = 0;
    double inner_prec = 5e-8;
    status = msvtplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                      traj_grid_size, coord_type, inner_prec, true,
                      DCM_EM_TFC, DCMS_SEM_TFC,
                      CCM_R_RCM_EM,
                      CCM_R_RCM_SEM,
                      orbit_EM, orbit_SEM,
                      h2, isPlotted, false);

    //======================================================================
    //5.2. If it is a sucess, we go on.
    //======================================================================
    if(status == GSL_SUCCESS)
    {
        cout << "---------------------------------------------"  << endl;
        cout << " crefvtplan. "  << endl;
        cout << "    Continuation with variable times         "  << endl;
        cout << "---------------------------------------------"  << endl;

        //======================================================================
        //Continuation procedure
        //======================================================================
        double ds = 1e-1, ds0 = 1e-1;
        double yv[42], ye[42];
        double dkn = 0.0;
        int kn = 0;
        do
        {
            //======================================================================
            // 5.2.1. Updating the free variables
            //======================================================================
            //--------------------------------------------
            // The objective of this continuation is to bring orbit_SEM.si[4] to 0.0.
            // Prior to updating, we check that the continuation is not going "to far"
            // (orbit_SEM.si[4] may change sign).
            //--------------------------------------------
            dkn = orbit_SEM.si[4] + ds0*nullvector[5*man_grid_size-1];
            if(dkn * orbit_SEM.si[4] < 0) //if there is a change of sign, we reduce the stepsize
            {
                ds = -orbit_SEM.si[4]/nullvector[5*man_grid_size-1];
            }
            else ds = ds0;

            //--------------------------------------------
            // Then we can go on
            //--------------------------------------------
            //Updating CM_EM_RCM coordinates
            orbit_EM.si[0] += ds*nullvector[0];
            orbit_EM.si[2] += ds*nullvector[1];

            //Updating CM_EM_NCEM coordinates
            orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
            //To CM_EM_NCSEM coordinates
            for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
            qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
            for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


            //The middle (patch points) is classical cartesian coordinates at patch points
            for(int k = 1; k < man_grid_size; k++)
            {
                y_traj_n[0][k] += ds*nullvector[0 + 5*k-3];
                y_traj_n[1][k] += ds*nullvector[1 + 5*k-3];
                y_traj_n[3][k] += ds*nullvector[2 + 5*k-3];
                y_traj_n[4][k] += ds*nullvector[3 + 5*k-3];
                t_traj_n[k]    += ds*nullvector[5*k+1];
            }

            //Last time:
            t_traj_n[man_grid_size] += ds*nullvector[5*man_grid_size];

            //Last 4 correction variables is orbit.si
            //Updating CM_SEM_RCM coordinates
            orbit_SEM.si[0] += ds*nullvector[5*man_grid_size-3];
            orbit_SEM.si[2] += ds*nullvector[5*man_grid_size-2];
            orbit_SEM.si[4]  = max(0.0, orbit_SEM.si[4] + ds*nullvector[5*man_grid_size-1]);

            //Updating CM_SEM_NCSEM coordinates
            orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

            //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
            for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];

            //======================================================================
            // 5.2.2. Diff correction
            //======================================================================
            isPlotted = (kn % plotfreq == 0) ? 1:0;
            status = msvtplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n, nullvector, 42,
                              traj_grid_size, coord_type, inner_prec, false,
                              DCM_EM_TFC, DCMS_SEM_TFC,
                              CCM_R_RCM_EM,
                              CCM_R_RCM_SEM,
                              orbit_EM, orbit_SEM,
                              h2, isPlotted, false);

            //======================================================================
            // 5.2.3. Display
            //======================================================================
            dkn = (double) kn;
            if(kn % plotfreq == 0) gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);


            //======================================================================
            // 5.2.4. Advance
            //======================================================================
            kn++;
            cout << "Step n°" << kn << " completed.               "  << endl;
            cout << "---------------------------------------------"  << endl;
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);

        }
        while(kn < maxStep && status == GSL_SUCCESS && fabs(orbit_SEM.si[4]) > ORBIT_SEM_UNSTABLE_MIN);
    }

    //====================================================================================
    // 6. Display final solution
    //====================================================================================
    if(status == GSL_SUCCESS)
    {
        //=============================================================
        // 6.1. Compute the initial orbit
        //=============================================================
        cout << "---------------------------------------------"  << endl;
        cout << " crefvtplan. "  << endl;
        cout << "   Computation of the initial EML2 orbit     "  << endl;
        cout << "---------------------------------------------"  << endl;

        //---------------------------------------------------------------------
        // Reset the unstable direction
        //---------------------------------------------------------------------
        orbit_EM.si[4] = 0.0;
        orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

        cout << "orbit_EM.si = " << endl;
        vector_printf_prec(orbit_EM.si, 5);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);


        //====================================================================================
        // 6.2. Final trajectory, on a grid
        //====================================================================================
        cout << "---------------------------------------------"  << endl;
        cout << " crefvtplan. "  << endl;
        cout << "   Final trajectory, on a grid               "  << endl;
        cout << "---------------------------------------------"  << endl;

        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);
        double yv[42];
        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
        }

        //====================================================================================
        // Free.
        //====================================================================================
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
    }


    //====================================================================================
    // 7. Update the orbits for next step
    //====================================================================================
    if(status == GSL_SUCCESS)
    {
        //====================================
        // At EML2. The first 4 RCM components are good, as well as the initial time.
        // Hence, we need to update:
        // 1. The last RCM component (unstable part),
        // 2. The final time.
        //====================================
        orbit_EM.si[4] = PROJ_EPSILON;
        orbit_EM.tf    = t_traj_n[traj_grid_size]/SEML.us_em.ns;

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        //====================================
        // At SEML2. The first 5 RCM components are good.
        // Hence, we need to update:
        // 1. The initial time.
        //====================================
        orbit_SEM.t0    = t_traj_n[traj_grid_size];

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);
    }

    cout << "---------------------------------------------"  << endl;
    cout << " crefvtplan. "  << endl;
    cout << "          End of computation.                "  << endl;
    cout << "---------------------------------------------"  << endl;

    //====================================================================================
    // Free
    //====================================================================================
    gnuplot_close(h5);
    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}

/**
 *  \brief Continuation of a single of EML2-to-SEML2 connection, between orbit_EM and orbit_SEM.
 *         The Jacobian of the parameterization of the manifolds are contained in  DCM_EM_TFC and DCMS_SEM_TFC.
 *         The patch and final times are fixed. Hence, the continuation procedure is "forced" to go transversally to the original solution, producing a family of connections.
 **/
void crefftplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        int cont_steps_MAX,
        int isSaved,
        int isUserDefined,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 100;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;        //TOF on EML2 orbit
    double tof_seml_SEM = 30*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction with fixed time    "  << endl;
    cout << "-------------------------------------------"  << endl;
    //Free variables
    int nfv = 4*(man_grid_size-1)+5;
    double *nullvector = dvector(0, nfv-1);


    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");
    gnuplot_cmd(h3, "set grid");

    gnuplot_ctrl *h4;
    h4 = gnuplot_init();
    gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");
    gnuplot_cmd(h4, "set grid");

    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
    gnuplot_cmd(h5, "set grid");

    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = floor((double)cont_steps_MAX/floor(0.1*GP_MAX_TMP_FILES))+1;

    //======================================================================
    //GNUPLOT FOR FINAL TRAJECTORY
    //======================================================================
    gnuplot_ctrl *h6;
    h6 = gnuplot_init();
    gnuplot_cmd(h6,  "set title \"Final trajectory\" ");
    gnuplot_cmd(h6, "set grid");

    double **semP_coord = dmatrix(0, 6, 0, 2);
    switch(coord_type)
    {
    case PSEM:
        semPoints(0.0, semP_coord);
        semPlot(h6, semP_coord);
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h6, semP_coord);
        break;
    case PEM:
        emPoints(0.0, semP_coord);
        emPlot(h6, semP_coord);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h6, semP_coord);
        break;
    }


    //======================================================================
    //Continuation, step one : Differential correction procedure with fixed times
    //======================================================================
    int status = 0;
    double yv[42];
    double **ymc   = dmatrix(0, 5, 0, mPlot);
    double *tmc    = dvector(0, mPlot);

    cout << "Old. y_traj[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][0] << endl;

    cout << "Old. y_traj[end-1] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][man_grid_size-1] << endl;


    cout << "Old. y_traj[end] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj[i][man_grid_size] << endl;


    //Differential correction procedure
    status = msftplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                      traj_grid_size, coord_type, 5e-8, true,
                      DCM_EM_TFC, DCMS_SEM_TFC,
                      CCM_R_RCM_EM,
                      CCM_R_RCM_SEM,
                      orbit_EM, orbit_SEM,
                      h2, isPlotted, isUserDefined, false);

    cout << "Old. y_traj_n[0] = " << endl;
    for(int i = 0; i <6; i++) cout << y_traj_n[i][0] << endl;

    //======================================================================
    // Go on if success
    //======================================================================
    if(status == GSL_SUCCESS)
    {
        //======================================================================
        // Save first entry
        //======================================================================
        string filename      = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF, SEML.li_SEM, orbit_EM.t0/SEML.us_em.T);
        string filename_traj = filenameCUM(OFTS_ORDER, TYPE_CONT_ATF_TRAJ, SEML.li_SEM, orbit_EM.t0/SEML.us_em.T);
        fstream filestream;
        if(isSaved)
        {
            //=================================
            // Main parameters
            //=================================
            filestream.open (filename.c_str(), ios::out);
            //Title
            filestream << "t0_CMU_EM  s1_CMU_EM  s2_CMU_EM  s3_CMU_EM  s4_CMU_EM s5_CMU_EM ";
            filestream << "tf_CMU_EM  s1_CMS_SEM s2_CMS_SEM s3_CMS_SEM s4_CMS_SEM s5_CMS_SEM ";
            filestream << "x0_CMU_NCEM  y0_CMU_NCEM z0_CMU_NCEM px0_CMU_NCEM py0_CMU_NCEM pz0_CMU_NCEM ";
            filestream << "x0_CMS_NCSEM y0_CMS_NCSEM z0_CMS_NCSEM px0_CMS_NCSEM py0_CMS_NCSEM pz0_CMS_NCSEM ";
            filestream << endl;
            //Data
            filestream << setprecision(15) <<  setiosflags(ios::scientific);
            filestream << orbit_EM.t0 << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_EM.si[i]  << "  ";
            filestream << orbit_EM.tf << "  ";
            for(int i = 0; i <5; i++) filestream << orbit_SEM.si[i] << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_EM.z0[i]  << "  ";
            for(int i = 0; i <6; i++) filestream << orbit_SEM.z0[i] << "  ";
            filestream << endl;
            filestream.close();


            //=================================
            // Entire trajectory
            //=================================
            filestream.open (filename_traj.c_str(), ios::out | ios::binary);
            double res;
            //Final trajectory on lines, segment by segment
            for(int k = 0; k < traj_grid_size; k++)
            {
                //Integration segment by segment
                for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        }


        cout << "-------------------------------------------"  << endl;
        cout << "Continuation with fixed time               "  << endl;
        cout << "-------------------------------------------"  << endl;
        //======================================================================
        // Continuation procedure
        //======================================================================
        double ye[42];
        double ds = 5e-1;
        double dkn = 0.0;
        int kn = 0;
        do
        {
            //======================================================================
            //Updating the free variables
            //======================================================================
            //Updating CM_EM_RCM coordinates
            orbit_EM.si[0] += ds*nullvector[0];
            orbit_EM.si[2] += ds*nullvector[1];

            //Updating CM_EM_NCEM coordinates
            orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
            //To CM_EM_NCSEM coordinates
            for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
            qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
            for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


            //The middle (patch points) is classical cartesian coordinates at patch points
            for(int k = 1; k < man_grid_size; k++)
            {
                y_traj_n[0][k] += ds*nullvector[0 + 4*k-2];
                y_traj_n[1][k] += ds*nullvector[1 + 4*k-2];
                y_traj_n[3][k] += ds*nullvector[2 + 4*k-2];
                y_traj_n[4][k] += ds*nullvector[3 + 4*k-2];
            }

            //Last 3 correction variables is orbit.si
            //Updating CM_SEM_RCM coordinates
            orbit_SEM.si[0] += ds*nullvector[4*man_grid_size-2];
            orbit_SEM.si[2] += ds*nullvector[4*man_grid_size-1];
            orbit_SEM.si[4] += ds*nullvector[4*man_grid_size-0];

            //Updating CM_SEM_NCSEM coordinates
            orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

            //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
            for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];

            //======================================================================
            //Diff corr
            //======================================================================
            isPlotted = (kn % plotfreq == 0) ? 1:0;
            status = msftplan(y_traj_n, t_traj_n, y_traj_n, t_traj_n, nullvector, 42,
                              traj_grid_size, coord_type, 5e-8, false,
                              DCM_EM_TFC, DCMS_SEM_TFC,
                              CCM_R_RCM_EM,
                              CCM_R_RCM_SEM,
                              orbit_EM, orbit_SEM,
                              h2, isPlotted, isUserDefined, false);


            //======================================================================
            //Save
            //======================================================================
            if(isSaved && status == GSL_SUCCESS)
            {
                filestream.open (filename.c_str(), ios::out | ios::app);
                //Data
                filestream << setprecision(15) <<  setiosflags(ios::scientific);
                filestream << orbit_EM.t0 << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_EM.si[i]  << "  ";
                filestream << orbit_EM.tf << "  ";
                for(int i = 0; i <5; i++) filestream << orbit_SEM.si[i] << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_EM.z0[i]  << "  ";
                for(int i = 0; i <6; i++) filestream << orbit_SEM.z0[i] << "  ";
                filestream << endl;
                filestream.close();

                //=================================
                // Entire trajectory
                //=================================
                filestream.open (filename_traj.c_str(), ios::out | ios::binary | ios::app);
                double res;
                //Final trajectory on lines, segment by segment
                for(int k = 0; k < traj_grid_size; k++)
                {
                    //Integration segment by segment
                    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
                    ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
            //Display
            //======================================================================
            dkn = (double) kn;
            if(kn % plotfreq == 0)
            {
                gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
                gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
                gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);
            }

            //======================================================================
            //Advance one step
            //======================================================================
            kn++;
            cout << "Step n°" << kn << "/" << cont_steps_MAX << " completed.               "  << endl;
            cout << "---------------------------------------------"  << endl;
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);

        }
        while(kn < cont_steps_MAX && status == GSL_SUCCESS);

        //====================================================================================
        // 5. Compute the initial orbit
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Compute the initial EML2 orbit             "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //---------------------------------------------------------------------
        // Reset the unstable direction
        //---------------------------------------------------------------------
        orbit_EM.si[4] = 0.0;
        orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
        gnuplot_plot_xyz(h6, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"EML2 orbit", "lines", "1", "2", 2);

        //================================================================================
        // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
        //================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Compute the initial SEML2 orbit            "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        vector<Oftsc> Fh;
        Fh.reserve(5);
        for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
        readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

        //---------------------------------------------------------------------
        //For dot(s) = fh(s)
        //---------------------------------------------------------------------
        RVF rvf;
        rvf.ofs_order = SEML.eff_nf_SEM;
        Ofsc AUX(rvf.ofs_order);
        rvf.fh         = &Fh;
        rvf.ofs        = &AUX;
        rvf.order      = OFTS_ORDER;
        rvf.n          = orbit_SEM.n;
        rvf.reduced_nv = 5;

        gsl_odeiv2_system sys_fh;
        sys_fh.function  = qbfbp_fh;
        sys_fh.jacobian  = NULL;
        sys_fh.dimension = 2*rvf.reduced_nv;
        sys_fh.params    = &rvf;

        const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
        gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


        //---------------------------------------------------------------------
        // Temp variables
        //---------------------------------------------------------------------
        double t0_SEM = t_traj_n[traj_grid_size];
        double t1_SEM = t0_SEM+tof_seml_SEM;
        double z[6];
        double t2 = t0_SEM;
        int k  = 0;
        double  s1ccm8[2*rvf.reduced_nv]; //CCM8

        //---------------------------------------------------------------------
        // Initial state in CCM8 form
        //---------------------------------------------------------------------
        RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

        //---------------------------------------------------------------------
        // Loop
        //---------------------------------------------------------------------
        while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
        {
            cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;
            //To NC coordinates
            RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                         orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

            //Save
            for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
            t_man_SEM[k] = t2;

            //Plot
            gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

            //Advance one step
            gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
            CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);
            k++;
        }

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);
        gnuplot_plot_xyz(h6, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"SEML2 orbit", "lines", "1", "2", 3);

        //---------------------------------------------------------------------
        //Old version using classic integrator
        //---------------------------------------------------------------------
        //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
        //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


        //====================================================================================
        // 7. Final trajectory, on a grid
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Final trajectory, on a grid                "  << endl;
        cout << "-------------------------------------------"  << endl;
        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);
            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
            //Plot on h6
            if(k == 0) gnuplot_plot_xyz(h6, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Transfer leg", "lines", "1", "2", 0);
            else gnuplot_plot_xyz(h6, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
        }

        //====================================================================================
        //Plot the final trajectory
        //====================================================================================
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);

        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
    }

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //====================================================================================
    // Free
    //====================================================================================
    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);
    gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}

/**
 *  \brief Same as crefftplan, but in the 3D case.
 **/
void crefft3d(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        int cont_steps_MAX,
        int isUserDefined,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;        //TOF on EML2 orbit
    double tof_seml_SEM = 30*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction & final trajectory "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int nfv = 6*man_grid_size+3;  //free variables
    double *nullvector = dvector(0, nfv-1);

    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");

    gnuplot_ctrl *h4;
    h4 = gnuplot_init();
    gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");

    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s4_EM vs s2_EM\" ");

    //We keep the number of points below 10% of the maximum gnuplot temp file.
    int plotfreq = 1;//floor((double)cont_steps_MAX/floor(0.1*GP_MAX_TMP_FILES))+1;


    //======================================================================
    //Continuation, step one : Differential correction procedure with fixed times
    //======================================================================
    int status = 0;
    double yv[42];
    //Differential correction procedure
    status = msft3d(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                    traj_grid_size, coord_type, 5e-8, true,
                    DCM_EM_TFC, DCMS_SEM_TFC,
                    CCM_R_RCM_EM,
                    CCM_R_RCM_SEM,
                    orbit_EM, orbit_SEM,
                    h2, isPlotted, isUserDefined, false);

    cout << "nullvector = " << endl;
    vector_printf_prec(nullvector, nfv);


    //======================================================================
    // Go on if success
    //======================================================================
    if(status == GSL_SUCCESS)
    {
        cout << "-------------------------------------------"  << endl;
        cout << "Continuation with fixed time               "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //======================================================================
        // Continuation procedure
        //======================================================================
        double ye[42];
        double ds = 1e-2;
        //double dkn = 0.0;
        int kn = 0;
        do
        {
            //======================================================================
            //Updating the free variables
            //======================================================================
            //Updating CM_EM_RCM coordinates
            for(int i = 0; i < 4; i++) orbit_EM.si[i] += ds*nullvector[i];
            //Updating CM_EM_NCEM coordinates
            orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
            //To CM_EM_NCSEM coordinates
            for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
            qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
            for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

            //The middle (patch points) is classical cartesian coordinates at patch points
            for(int k = 1; k < man_grid_size; k++)
            {
                for(int i = 0; i < 6; i++) y_traj_n[i][k] += ds*nullvector[i + 6*k-2];
            }

            //Last 5 correction variables is orbit.si
            //Updating CM_SEM_RCM coordinates
            for(int i = 0; i < 5; i++) orbit_SEM.si[i] += ds*nullvector[i + 6*man_grid_size-2];

            //Updating CM_SEM_NCSEM coordinates
            orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

            //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
            for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];


            cout << "orbit_EM.si (1) = " << endl;
            vector_printf_prec(orbit_EM.si, 5);

            cout << "orbit_SEM.si (1) = " << endl;
            vector_printf_prec(orbit_SEM.si, 5);

            //======================================================================
            //Diff corr
            //======================================================================
            isPlotted = (kn % plotfreq == 0) ? 1:0;
            status = msft3d(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
                            traj_grid_size, coord_type, 5e-8, false,
                            DCM_EM_TFC, DCMS_SEM_TFC,
                            CCM_R_RCM_EM,
                            CCM_R_RCM_SEM,
                            orbit_EM, orbit_SEM,
                            h2, isPlotted, isUserDefined, false);

            cout << "orbit_EM.si (2) = " << endl;
            vector_printf_prec(orbit_EM.si, 5);

            cout << "orbit_SEM.si (2) = " << endl;
            vector_printf_prec(orbit_SEM.si, 5);

            //======================================================================
            //Save
            //======================================================================

            //======================================================================
            //Display
            //======================================================================
            //dkn = (double) kn;
            if(kn % plotfreq == 0)
            {
                gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
                gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
                gnuplot_plot_xy(h5, &orbit_EM.si[1],  &orbit_EM.si[3], 1, (char*)"", "points", "1", "2", 0);
            }


            //======================================================================
            //Advance one step
            //======================================================================
            kn++;
            cout << "Step n°" << kn << " completed.               "  << endl;
            cout << "---------------------------------------------"  << endl;
            //printf("Press ENTER to go on\n");
            //scanf("%c",&ch);

        }
        while(kn < cont_steps_MAX && status == GSL_SUCCESS);

        //====================================================================================
        // 5. Compute the initial orbit
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Compute the initial EML2 orbit             "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //---------------------------------------------------------------------
        // Reset the unstable direction
        //---------------------------------------------------------------------
        orbit_EM.si[4] = 0.0;
        orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

        //---------------------------------------------------------------------
        // Gnuplot (NCEM)
        //---------------------------------------------------------------------
        gnuplot_ctrl *h6;
        h6 = gnuplot_init();
        gnuplot_cmd(h6,  "set title \"EML2 orbit  in NCEM coordinates\" ");

        cout << "orbit_EM.si = " << endl;
        vector_printf_prec(orbit_EM.si, 5);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

        cout << "orbit_EM.z0 = " << endl;
        vector_printf_prec(orbit_EM.z0 , 6);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
        gnuplot_plot_xyz(h6, y_man_NCEM[0], y_man_NCEM[1],  y_man_NCEM[2], mPlot+1, (char*)"", "lines", "1", "1", 2);



        //====================================================================================
        // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Compute the initial SEML2 orbit            "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        vector<Oftsc> Fh;
        Fh.reserve(5);
        for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
        readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

        //---------------------------------------------------------------------
        //For dot(s) = fh(s)
        //---------------------------------------------------------------------
        RVF rvf;
        rvf.ofs_order = SEML.eff_nf_SEM;
        Ofsc AUX(rvf.ofs_order);
        rvf.fh         = &Fh;
        rvf.ofs        = &AUX;
        rvf.order      = OFTS_ORDER;
        rvf.n          = orbit_SEM.n;
        rvf.reduced_nv = 5;

        gsl_odeiv2_system sys_fh;
        sys_fh.function  = qbfbp_fh;
        sys_fh.jacobian  = NULL;
        sys_fh.dimension = 2*rvf.reduced_nv;
        sys_fh.params    = &rvf;

        const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
        gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


        //---------------------------------------------------------------------
        // Temp variables
        //---------------------------------------------------------------------
        double t0_SEM = t_traj_n[traj_grid_size];
        double t1_SEM = t0_SEM+tof_seml_SEM;
        double z[6];
        double t2 = t0_SEM;
        int k  = 0;
        double  s1ccm8[2*rvf.reduced_nv]; //CCM8

        //---------------------------------------------------------------------
        // Initial state in CCM8 form
        //---------------------------------------------------------------------
        RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

        //---------------------------------------------------------------------
        // Loop
        //---------------------------------------------------------------------
        while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
        {
            cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;

            //To NC coordinates
            RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                         orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

            //Save
            for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
            t_man_SEM[k] = t2;

            //Plot
            gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

            //Advance one step
            gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
            CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);

            k++;
        }

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);

        //---------------------------------------------------------------------
        //Old version using classic integrator
        //---------------------------------------------------------------------
        //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
        //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


        //====================================================================================
        // 7. Final trajectory, on a grid
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Final trajectory, on a grid                "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);

        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
        }

        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
        gnuplot_close(h6);


    }
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //====================================================================================
    // Free
    //====================================================================================
    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);


    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}



//====================================
// Single solutions
//====================================
/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D, in the planar case.
 **/
void srefvtplan(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 20*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction & final trajectory "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    double yv[42];
    int nfv = 7*(man_grid_size+1)-4;  //free variables
    double *nullvector = dvector(0, nfv-1);

    //======================================================================
    // NO CONTINUATION:
    // @todo Update msdvt_CMS_RCM to meet msvtplan
    // standards in terms of inputs and plotting. msvtplan is used for now
    // BUT the computation of the null vector is NOT used...
    //======================================================================
    //    msdvt_CMS_RCM(y_traj, t_traj, y_traj_n, t_traj_n, 42,
    //                  traj_grid_size, coord_type,
    //                  DCM_EM_TFC, DCMS_SEM_TFC,
    //                  CCM_R_RCM_EM,
    //                  CCM_R_RCM_SEM,
    //                  orbit_EM, orbit_SEM,
    //                  h2, isPlotted, false);

    //With fixed time
    msftplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
             traj_grid_size, coord_type, 5e-8, true,
             DCM_EM_TFC, DCMS_SEM_TFC,
             CCM_R_RCM_EM,
             CCM_R_RCM_SEM,
             orbit_EM, orbit_SEM,
             h2, isPlotted, true,  false);

    // With unfixed time
    //    msvtplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
    //                                  traj_grid_size, coord_type, 5e-8, true,
    //                                  DCM_EM_TFC, DCMS_SEM_TFC,
    //                                  CCM_R_RCM_EM,
    //                                  CCM_R_RCM_SEM,
    //                                  orbit_EM, orbit_SEM,
    //                                  h2, isPlotted, false);


    //====================================================================================
    // 5. Compute the initial orbit
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial EML2 orbit             "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Reset the unstable direction
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;
    orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

    //---------------------------------------------------------------------
    // Gnuplot (NCEM)
    //---------------------------------------------------------------------
    gnuplot_ctrl *h6;
    h6 = gnuplot_init();
    gnuplot_cmd(h6,  "set title \"EML2 orbit  in NCEM coordinates\" ");

    cout << "orbit_EM.si = " << endl;
    vector_printf_prec(orbit_EM.si, 5);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    cout << "orbit_EM.z0 = " << endl;
    vector_printf_prec(orbit_EM.z0 , 6);

    //---------------------------------------------------------------------
    //Integration on mPlot+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
    gnuplot_plot_xyz(h6, y_man_NCEM[0], y_man_NCEM[1],  y_man_NCEM[2], mPlot+1, (char*)"", "lines", "1", "1", 2);



    //====================================================================================
    // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial SEML2 orbit            "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    vector<Oftsc> Fh;
    Fh.reserve(5);
    for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

    //---------------------------------------------------------------------
    //For dot(s) = fh(s)
    //---------------------------------------------------------------------
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf_SEM;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh         = &Fh;
    rvf.ofs        = &AUX;
    rvf.order      = OFTS_ORDER;
    rvf.n          = orbit_SEM.n;
    rvf.reduced_nv = 5;

    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*rvf.reduced_nv;
    sys_fh.params    = &rvf;

    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


    //---------------------------------------------------------------------
    // Temp variables
    //---------------------------------------------------------------------
    double t0_SEM = t_traj_n[traj_grid_size];
    double t1_SEM = t0_SEM+tof_seml_SEM;
    double z[6];
    double t2 = t0_SEM;
    int k  = 0;
    double  s1ccm8[2*rvf.reduced_nv]; //CCM8

    //---------------------------------------------------------------------
    // Initial state in CCM8 form
    //---------------------------------------------------------------------
    RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

    //---------------------------------------------------------------------
    // Loop
    //---------------------------------------------------------------------
    while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
    {
        cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;

        //To NC coordinates
        RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                     orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

        //Save
        for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
        t_man_SEM[k] = t2;

        //Plot
        gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

        //Advance one step
        gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
        CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);

        k++;
    }

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);

    //---------------------------------------------------------------------
    //Old version using classic integrator
    //---------------------------------------------------------------------
    //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
    //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


    //====================================================================================
    // 7. Final trajectory, on a grid
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Final trajectory, on a grid                "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    double **ymc   = dmatrix(0, 5, 0, mPlot);
    double *tmc    = dvector(0, mPlot);

    //Final trajectory on lines, segment by segment
    for(int k = 0; k < traj_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
    }

    //---------------------------------------------------------------------
    // Final trajectory, whith single shooting integration
    //---------------------------------------------------------------------
    int kstart = 0; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
    ode78(y_traj, t_traj, t_traj_n[kstart], t2, yv, 6, traj_grid_size, dcs, coord_type, coord_type);

    //Plot on h2
    //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}



//========================================================================================
//
//         Refinement of solutions: with additionnal search on the initial parameters (time, or s3_EM...). DEPRECATED
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution. DEPRECATED
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_T(int proj_grid_size,
        int man_grid_size,
        int coord_type,
        double st_EM[],
        double t0_EM,
        double tmax_on_manifold_EM,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        double ynormMax,
        double snormMax,
        int isPlanar)
{
    //====================================================================================
    // 0. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 1. Time on each orbit
    //====================================================================================
    double tof_seml_SEM = 5*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP_coord = dmatrix(0, 6, 0, 2);
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


    //====================================================================================
    // 5. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCMS_SEM_TFC(6, 5, 5, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCMS_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 7. Loop. Only one solution for now!
    //====================================================================================

    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;
    cout << "Begin the search for the min distance of projection..." << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // 6.1 Initialize local variables: EM
    //====================================================================================
    //---------------------------------------------------------------------
    //Initialize the initial conditions (both NC and RCM coordinates)
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    Ofsc orbit_ofs_EM(OFS_ORDER);
    OdeStruct driver_EM;
    SingleOrbit orbit_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
    //Init routine
    init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
               &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
               1, t0_EM, t0_EM+tmax_on_manifold_EM, SEML.us_em.T/5, &driver_EM, &SEML_EM);


    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, st_EM, t0_EM);


    //====================================================================================
    // 6.3 Estimate the distance of projection
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_NCSEM = dmatrix(0, 5, 0, proj_grid_size);
    double *t_man_SEM    = dvector(0, proj_grid_size);


    //---------------------------------------------------------------------
    //Integration on proj_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int status = ode78(y_man_NCSEM, t_man_SEM, t0_EM, t0_EM+tmax_on_manifold_EM, orbit_EM.z0, 6, proj_grid_size, I_NCEM, NCEM, NCSEM);

    //----------------------------------------------------------
    // projection by default is arbitrarily big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //----------------------------------------------------------
    // We need first to check that the integration went well
    //----------------------------------------------------------
    Ofsc ofs(OFS_ORDER);
    double yv[6], tv;
    double proj_dist_SEM, min_proj_dist_SEM = ePdef, y_man_norm_NCSEM = 0.0;
    double yvproj_NCSEM[6], sproj[5], yv_SEM[6], yvproj_SEM[6], yv_VSEM[6], yvproj_VSEM[6];

    //Optimal variables
    double sprojmin[5];
    double tmin = t0_EM+tmax_on_manifold_EM;

    if(status != -1)
    {
        //----------------------------------------------------------
        //Loop on trajectory
        //----------------------------------------------------------
        for(int kman = 0; kman <= proj_grid_size; kman++)
        {
            //Current state
            for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
            tv = t_man_SEM[kman];

            //Current distance from SEMLi in NCSEM units
            y_man_norm_NCSEM = 0.0;
            for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
            y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

            //If the current state is close enough to SEMLi
            if(y_man_norm_NCSEM < ynormMax)
            {
                // Projection on the center manifold
                NCprojCCMtoCUS(yv, tv, SEML_SEM.us.n, sproj, CMS_SEM_TFC, 0.0, MIcoc_SEM, Vcoc_SEM);

                //If the projection is close enough to SEMLi
                if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
                {
                    //yvproj_NCSEM = W(sproj, tv)
                    RCMtoNCbyTFC(sproj, tv, SEML_SEM.us.n, OFTS_ORDER, OFS_ORDER, 5, CMS_SEM_TFC, Mcoc_SEM, Vcoc_SEM, yvproj_NCSEM, 1);

                    //Back to SEM coordinates
                    NCSEMmtoSEMm(tv, yv, yv_SEM, &SEML_SEM);
                    NCSEMmtoSEMm(tv, yvproj_NCSEM, yvproj_SEM, &SEML_SEM);

                    //Back to SEM coordinates with velocity
                    SEMmtoSEMv(tv, yv_SEM, yv_VSEM, &SEML_SEM);
                    SEMmtoSEMv(tv, yvproj_SEM, yvproj_VSEM, &SEML_SEM);

                    //Distance of projection in SEM coordinates
                    proj_dist_SEM = 0.0;
                    for(int i = 0; i < 6; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
                    proj_dist_SEM = sqrt(proj_dist_SEM);
                }
                else proj_dist_SEM = ePdef;
            }
            else proj_dist_SEM = ePdef;

            //Update distance min if necessary
            if(proj_dist_SEM < min_proj_dist_SEM)
            {
                min_proj_dist_SEM = proj_dist_SEM;

                //Update optimal variables
                for(int i = 0; i < 5; i++) sprojmin[i] = sproj[i];
                tmin  = tv;
            }

        }

    }
    else proj_dist_SEM = ePdef;


    cout << "-------------------------------------------"  << endl;
    cout << "End of the search for the min distance of projection..." << endl;
    cout << "The minimum is : " << min_proj_dist_SEM << endl;
    cout << "Obtained for sproj = " << endl;
    vector_printf_prec(sproj, 5);
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    //Update the time in orbit_EM!!
    //---------------------------------------------------------------------
    orbit_EM.tf = tmin/SEML.us_em.ns;

    //====================================================================================
    // 6.2 Initialize local variables: SEM
    //====================================================================================

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    Ofsc orbit_ofs_SEM(OFS_ORDER);
    OdeStruct driver_SEM;
    SingleOrbit orbit_SEM;

    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
    //Init routine
    init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
               &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
               1, tmin, tmin+tof_seml_SEM, SEML.us_sem.T, &driver_SEM, &SEML_SEM);


    //====================================================================================
    // @todo: ADD AN ESTIMATE OF THE DISTANCE OF PROJECTION + by reproducing the PROJECTION METHOD
    // THIS TIME ON THE CMS of SEML2! (set s5 of SEML2 to zero, or maybe a rough search...)
    // In this way, only 3 inputs: s1, s3 and t0.
    //====================================================================================

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_SEM, sprojmin, tmin);


    //====================================================================================
    // 6.2 Multiple shooting procedure
    //====================================================================================
    if(isPlanar)
    {
        srefvtplan(orbit_EM, orbit_SEM,
                                              DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }
    else
    {
        ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(orbit_EM, orbit_SEM,
                                          DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    return 0;
}

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution. DEPRECATED
 **/
int ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_OPT(int proj_grid_size,
        int man_grid_size,
        int coord_type,
        double st_EM[],
        double t0_EM,
        double tmax_on_manifold_EM,
        double s3_MIN_CMU_RCM,
        double s3_MAX_CMU_RCM,
        int s3_grid_size,
        vector<Oftsc> &CM_EM_NC,
        vector<Oftsc> &CM_EM_TFC,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Ofsc>  &Mcoc_EM,
        matrix<Ofsc>  &MIcoc_EM,
        vector<Ofsc>  &Vcoc_EM,
        double ynormMax,
        double snormMax,
        int isPlanar,
        int isPar)
{
    //====================================================================================
    // 0. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 1. Time on each orbit
    //====================================================================================
    double tof_seml_SEM = 5*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 3. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP_coord = dmatrix(0, 6, 0, 2);
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


    //====================================================================================
    // 5. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
    CMS_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CMS_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CMS_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CMS_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCMS_SEM_TFC(6, 5, 5, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCMS_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 7. Loop. Only one solution for now!
    //====================================================================================

    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;
    cout << "Begin the search for the min distance of projection..." << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // 6.1 Initialize variables: EM
    //====================================================================================
    //---------------------------------------------------------------------
    //Initialize the initial conditions (both NC and RCM coordinates)
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    Ofsc orbit_ofs_EM(OFS_ORDER);
    OdeStruct driver_EM;
    SingleOrbit orbit_EM;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
    //Init routine
    init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
               &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
               1, t0_EM, t0_EM+tmax_on_manifold_EM, SEML.us_em.T/5, &driver_EM, &SEML_EM);


    //----------------------------------------------------------
    // projection by default is arbitrarily big
    //----------------------------------------------------------
    double ePdef = 1e5;

    //====================================================================================
    // 6.2 Initialize variables: Optimal variables
    //====================================================================================
    double smin_proj_dist_SEM = ePdef;
    double sprojsmin[5];
    double tsmin = t0_EM+tmax_on_manifold_EM;
    int ks3min = 0;

    //------------------------------------------
    //Building the working grids
    //------------------------------------------
    double *grid_s3_CMU_RCM = dvector(0,  s3_grid_size);
    init_grid(grid_s3_CMU_RCM, s3_MIN_CMU_RCM, s3_MAX_CMU_RCM, s3_grid_size);

    int s5_grid_size = 0;
    double *grid_s5_CMU_RCM = dvector(0,  s5_grid_size);
    init_grid(grid_s5_CMU_RCM, 0.0, +5e-2, s5_grid_size);

    //====================================================================================
    // 6.3 Estimate the distance of projection
    //====================================================================================
    int iter = 0;
    COMPLETION = 0;
    #pragma omp parallel for if(isPar)  shared(iter, sprojsmin, smin_proj_dist_SEM, tsmin)
    for(int ks3 = 0; ks3 <= s3_grid_size; ks3++)
    {
        //---------------------------------------------------------------------
        // Update the initial state in the orbit_EM, with the RCM coordinates
        //---------------------------------------------------------------------
        st_EM[2] = grid_s3_CMU_RCM[ks3];
        st_EM[4] = PROJ_EPSILON;
        orbit_update_ic(orbit_EM, st_EM, t0_EM);

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_NCSEM = dmatrix(0, 5, 0, proj_grid_size);
        double *t_man_SEM    = dvector(0, proj_grid_size);


        //---------------------------------------------------------------------
        //Integration on proj_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        int status = ode78(y_man_NCSEM, t_man_SEM, t0_EM, t0_EM+tmax_on_manifold_EM, orbit_EM.z0, 6, proj_grid_size, I_NCEM, NCEM, NCSEM);

        //----------------------------------------------------------
        // We need first to check that the integration went well
        //----------------------------------------------------------
        Ofsc ofs(OFS_ORDER);
        double yv[6], tv;
        double proj_dist_SEM, proj_dist_SEM_s5, min_proj_dist_SEM = ePdef, y_man_norm_NCSEM = 0.0;
        double yvproj_NCSEM[6], sproj[5], sproj_s5[5], yv_SEM[6], yvproj_SEM[6], yv_VSEM[6], yvproj_VSEM[6];


        double sprojmin[5];
        double tmin = t0_EM+tmax_on_manifold_EM;

        if(status != -1)
        {
            //----------------------------------------------------------
            //Loop on trajectory
            //----------------------------------------------------------
            for(int kman = 0; kman <= proj_grid_size; kman++)
            {
                //Current state
                for(int i = 0; i < 6; i++) yv[i] = y_man_NCSEM[i][kman];
                tv = t_man_SEM[kman];

                //Current distance from SEMLi in NCSEM units
                y_man_norm_NCSEM = 0.0;
                for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
                y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

                //If the current state is close enough to SEMLi
                if(y_man_norm_NCSEM < ynormMax)
                {
                    // Projection on the center manifold
                    NCprojCCMtoCUS(yv, tv, SEML_SEM.us.n, sproj, CMS_SEM_TFC, 0.0, MIcoc_SEM, Vcoc_SEM);

                    //If the projection is close enough to SEMLi
                    if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< snormMax)
                    {
                        proj_dist_SEM = ePdef;
                        for(int ks5 = 0; ks5 <= s5_grid_size ; ks5++)
                        {

                            // Projection on the center manifold
                            NCprojCCMtoCUS(yv, tv, SEML_SEM.us.n, sproj_s5, CMS_SEM_TFC, grid_s5_CMU_RCM[ks5], MIcoc_SEM, Vcoc_SEM);

                            //yvproj_NCSEM = W(sproj, tv)
                            RCMtoNCbyTFC(sproj_s5, tv, SEML_SEM.us.n, OFTS_ORDER, OFS_ORDER, 5, CMS_SEM_TFC, Mcoc_SEM, Vcoc_SEM, yvproj_NCSEM, 1);

                            //Back to SEM coordinates
                            NCSEMmtoSEMm(tv, yv, yv_SEM, &SEML_SEM);
                            NCSEMmtoSEMm(tv, yvproj_NCSEM, yvproj_SEM, &SEML_SEM);

                            //Back to SEM coordinates with velocity
                            SEMmtoSEMv(tv, yv_SEM, yv_VSEM, &SEML_SEM);
                            SEMmtoSEMv(tv, yvproj_SEM, yvproj_VSEM, &SEML_SEM);

                            //Distance of projection in SEM coordinates
                            proj_dist_SEM_s5 = 0.0;
                            for(int i = 0; i < 6; i++) proj_dist_SEM_s5 += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
                            proj_dist_SEM_s5 = sqrt(proj_dist_SEM_s5);

                            //Update distance min if necessary
                            if(proj_dist_SEM_s5 < proj_dist_SEM)
                            {
                                proj_dist_SEM = proj_dist_SEM_s5;
                                //Update optimal variables
                                sproj[0] = sproj_s5[0];
                                sproj[1] = sproj_s5[1];
                                sproj[2] = sproj_s5[2];
                                sproj[3] = sproj_s5[3];
                                sproj[4] = sproj_s5[4];
                            }
                        }
                    }
                    else proj_dist_SEM = ePdef;
                }
                else proj_dist_SEM = ePdef;

                //Update distance min if necessary
                if(proj_dist_SEM < min_proj_dist_SEM)
                {
                    min_proj_dist_SEM = proj_dist_SEM;
                    //Update optimal variables
                    sprojmin[0] = sproj[0];
                    sprojmin[1] = sproj[1];
                    sprojmin[2] = sproj[2];
                    sprojmin[3] = sproj[3];
                    sprojmin[4] = sproj[4];
                    tmin  = tv;
                }

            }

        }
        else proj_dist_SEM = ePdef;

        //Save
        #pragma omp critical
        {

            //Update distance min if necessary
            if(min_proj_dist_SEM < smin_proj_dist_SEM)
            {
                smin_proj_dist_SEM = min_proj_dist_SEM;
                //Update optimal variables
                for(int i = 0; i < 5; i++) sprojsmin[i] = sprojmin[i];
                tsmin  = tmin;
                ks3min = ks3;
            }


            //Display
            displayCompletion("ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE_PROJ_OPT", (double) iter++/max(s3_grid_size+1,0)*100);
            cout << "tmin = " << tmin << endl;
            cout << "min_proj_dist_SEM = " << min_proj_dist_SEM << endl;

        }

        free_dmatrix(y_man_NCSEM, 0, 5, 0, proj_grid_size);
        free_dvector(t_man_SEM, 0, proj_grid_size);

    }


    cout << "-------------------------------------------"  << endl;
    cout << "End of the search for the min distance of projection..." << endl;
    cout << "The minimum is : " << smin_proj_dist_SEM << endl;
    cout << "Obtained for sproj = " << endl;
    vector_printf_prec(sprojsmin, 5);
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    //Update the time in orbit_EM!!
    //---------------------------------------------------------------------
    orbit_EM.tf = tsmin/SEML.us_em.ns;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_EM, with the RCM coordinates
    //---------------------------------------------------------------------
    st_EM[2] = grid_s3_CMU_RCM[ks3min];
    st_EM[4] = PROJ_EPSILON;
    orbit_update_ic(orbit_EM, st_EM, t0_EM);

    //====================================================================================
    // 6.2 Initialize local variables: SEM
    //====================================================================================

    //---------------------------------------------------------------------
    // Initialisation of the orbit structure
    //---------------------------------------------------------------------
    Ofsc orbit_ofs_SEM(OFS_ORDER);
    OdeStruct driver_SEM;
    SingleOrbit orbit_SEM;

    //Init ode structure
    init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
    //Init routine
    init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
               &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
               1, tsmin, tsmin+tof_seml_SEM, SEML.us_sem.T, &driver_SEM, &SEML_SEM);


    //====================================================================================
    // @todo: ADD AN ESTIMATE OF THE DISTANCE OF PROJECTION + by reproducing the PROJECTION METHOD
    // THIS TIME ON THE CMS of SEML2! (set s5 of SEML2 to zero, or maybe a rough search...)
    // In this way, only 3 inputs: s1, s3 and t0.
    //====================================================================================

    //---------------------------------------------------------------------
    // Update the initial state in the orbit_SEM, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_SEM, sprojsmin, tsmin);


    //====================================================================================
    // 6.2 Multiple shooting procedure
    //====================================================================================
    if(isPlanar)
    {
        srefvtplan(orbit_EM, orbit_SEM,
                                              DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }
    else
    {
        ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(orbit_EM, orbit_SEM,
                                          DCM_EM_TFC, DCMS_SEM_TFC, dcs, coord_type, man_grid_size, h2);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    return 0;
}

/**
 *  \brief Computes ONE trajectory from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CMS of SEML2. The time at each point
 *         except the first one is allowed to vary. A continuation procedure can be performed to get more than one refined solution.
 **/
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D(SingleOrbit &orbit_EM,
                                       SingleOrbit &orbit_SEM,
                                       matrix<Oftsc> &DCM_EM_TFC,
                                       matrix<Oftsc> &DCMS_SEM_TFC,
                                       int dcs,
                                       int coord_type,
                                       int man_grid_size,
                                       gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 20*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction & final trajectory "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int nfv = 7*(man_grid_size+1)-4;  //free variables
    double *nullvector = dvector(0, nfv-1);


    //======================================================================
    //NO CONTINUATION
    //======================================================================
    //    msdvt_CMS_RCM(y_traj, t_traj, y_traj_n, t_traj_n, 42,
    //                  traj_grid_size, coord_type,
    //                  DCM_EM_TFC, DCMS_SEM_TFC,
    //                  CCM_R_RCM_EM,
    //                  CCM_R_RCM_SEM,
    //                  orbit_EM, orbit_SEM,
    //                  h2, isPlotted, false);

    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");

    gnuplot_ctrl *h4;
    h4 = gnuplot_init();
    gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");

    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s4_SEM vs s2_SEM\" ");


    //======================================================================
    //Continuation : Update the variables with the null vector
    //======================================================================
    msvt3d(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
           traj_grid_size, coord_type, 1e-8, true,
           DCM_EM_TFC, DCMS_SEM_TFC,
           CCM_R_RCM_EM,
           CCM_R_RCM_SEM,
           orbit_EM, orbit_SEM,
           h2, isPlotted, false);


    cout << "-------------------------------------------"  << endl;
    cout << "Continuation"  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    double ds = 1e-2;
    double yv[42], ye[42];

    for(int kn = 0; kn < 50; kn++)
    {

        //-----------------------------------------------------
        //Updating CM_EM_RCM coordinates
        //-----------------------------------------------------
        for(int i = 0; i < 4; i++) orbit_EM.si[i] += ds*nullvector[i];
        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
        //To CM_EM_NCSEM coordinates
        for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

        //-----------------------------------------------------
        //The middle (patch points) is classical cartesian coordinates at patch points
        //-----------------------------------------------------
        for(int k = 1; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) y_traj_n[i][k] += ds*nullvector[i + 7*k-3];
            t_traj_n[k] += ds*nullvector[7*(k+1)-4];
        }
        //Last time:
        t_traj_n[man_grid_size] += ds*nullvector[ 7*man_grid_size+2];


        //-----------------------------------------------------
        //Last 4 correction variables is orbit.si
        //-----------------------------------------------------
        //Updating CM_SEM_RCM coordinates
        for(int i = 0; i < 5; i++) orbit_SEM.si[i] += ds*nullvector[i + 7*man_grid_size-3];

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];

        msvt3d(y_traj_n, t_traj_n, y_traj_n, t_traj_n, nullvector, 42,
               traj_grid_size, coord_type, 1e-8, false,
               DCM_EM_TFC, DCMS_SEM_TFC,
               CCM_R_RCM_EM,
               CCM_R_RCM_SEM,
               orbit_EM, orbit_SEM,
               h2, isPlotted, false);

        cout << "orbit_SEM.si = " << endl;
        vector_printf_prec(orbit_SEM.si, 5);

        gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h5, &orbit_SEM.si[1], &orbit_SEM.si[3], 1, (char*)"", "points", "1", "2", 0);


        cout << "t_traj_n[end] - t_traj_n[end--1] = " << t_traj_n[man_grid_size] - t_traj_n[man_grid_size-1] << endl;
        //printf("Press ENTER to go on\n");
        //scanf("%c",&ch);

    }

    //====================================================================================
    // 5. Compute the initial orbit
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial EML2 orbit             "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Reset the unstable direction
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;
    orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

    //---------------------------------------------------------------------
    // Gnuplot (NCEM)
    //---------------------------------------------------------------------
    gnuplot_ctrl *h6;
    h6 = gnuplot_init();
    gnuplot_cmd(h6,  "set title \"EML2 orbit  in NCEM coordinates\" ");

    cout << "orbit_EM.si = " << endl;
    vector_printf_prec(orbit_EM.si, 5);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    cout << "orbit_EM.z0 = " << endl;
    vector_printf_prec(orbit_EM.z0 , 6);

    //---------------------------------------------------------------------
    //Integration on mPlot+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
    gnuplot_plot_xyz(h6, y_man_NCEM[0], y_man_NCEM[1],  y_man_NCEM[2], mPlot+1, (char*)"", "lines", "1", "1", 2);



    //====================================================================================
    // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial SEML2 orbit            "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    vector<Oftsc> Fh;
    Fh.reserve(5);
    for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

    //---------------------------------------------------------------------
    //For dot(s) = fh(s)
    //---------------------------------------------------------------------
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf_SEM;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh         = &Fh;
    rvf.ofs        = &AUX;
    rvf.order      = OFTS_ORDER;
    rvf.n          = orbit_SEM.n;
    rvf.reduced_nv = 5;

    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*rvf.reduced_nv;
    sys_fh.params    = &rvf;

    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


    //---------------------------------------------------------------------
    // Temp variables
    //---------------------------------------------------------------------
    double t0_SEM = t_traj_n[traj_grid_size];
    double t1_SEM = t0_SEM+tof_seml_SEM;
    double z[6];
    double t2 = t0_SEM;
    int k  = 0;
    double  s1ccm8[2*rvf.reduced_nv]; //CCM8

    //---------------------------------------------------------------------
    // Initial state in CCM8 form
    //---------------------------------------------------------------------
    RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

    //---------------------------------------------------------------------
    // Loop
    //---------------------------------------------------------------------
    while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
    {
        cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;

        //To NC coordinates
        RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                     orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

        //Save
        for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
        t_man_SEM[k] = t2;

        //Plot
        gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

        //Advance one step
        gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
        CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);

        k++;
    }

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);

    //---------------------------------------------------------------------
    //Old version using classic integrator
    //---------------------------------------------------------------------
    //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
    //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


    //====================================================================================
    // 7. Final trajectory, on a grid
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Final trajectory, on a grid                "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    double **ymc   = dmatrix(0, 5, 0, mPlot);
    double *tmc    = dvector(0, mPlot);

    //Final trajectory on lines, segment by segment
    for(int k = 0; k < traj_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
    }

    //---------------------------------------------------------------------
    // Final trajectory, whith single shooting integration
    //---------------------------------------------------------------------
    int kstart = 0; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
    ode78(y_traj, t_traj, t_traj_n[kstart], t2, yv, 6, traj_grid_size, dcs, coord_type, coord_type);

    //Plot on h2
    //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);

    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);
    gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}


//========================================================================================
//
//         Refinement with Pseudo-Arclenght constraint: NOT WORKING RIGHT NOW (BUT NOT REALLY NECESSARY!)
//
//========================================================================================

/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D, in the planar case.
 **/
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT_PAC_AFT(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        int cont_steps_MAX,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    rotmat_CC_R_RCM_CENTER(CCM_R_RCM_EM);

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM_SEM);

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 20*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction & final trajectory "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int nfv = 4*(man_grid_size-1)+5;       //free variables
    double *nullvector = dvector(0, nfv-1);
    double *conv_free_var = dvector(0, nfv-1);

    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");
    gnuplot_cmd(h3, "set grid");

    gnuplot_ctrl *h4;
    h4 = gnuplot_init();
    gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");
    gnuplot_cmd(h4, "set grid");

    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
    gnuplot_cmd(h5, "set grid");


    //======================================================================
    //GNUPLOT FOR FINAL TRAJECTORY
    //======================================================================
    gnuplot_ctrl *h6;
    h6 = gnuplot_init();
    gnuplot_cmd(h6,  "set title \"Final trajectory\" ");
    gnuplot_cmd(h6, "set grid");

    double **semP_coord = dmatrix(0, 6, 0, 2);
    switch(coord_type)
    {
    case PSEM:
        semPoints(0.0, semP_coord);
        semPlot(h6, semP_coord);
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h6, semP_coord);
        break;
    case PEM:
        emPoints(0.0, semP_coord);
        emPlot(h6, semP_coord);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h6, semP_coord);
        break;
    }


    //======================================================================
    //Continuation : Update the variables with the null vector
    //======================================================================
    msftplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
             traj_grid_size, coord_type, 5e-8, true,
             DCM_EM_TFC, DCMS_SEM_TFC,
             CCM_R_RCM_EM,
             CCM_R_RCM_SEM,
             orbit_EM, orbit_SEM,
             h2, isPlotted, true, false);

    cout << "-------------------------------------------"  << endl;
    cout << "Continuation"  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    double ds = 1e-3;
    double yv[42], ye[42];
    double dkn = 0.0;
    for(int kn = 0; kn < cont_steps_MAX; kn++)
    {
        //====================================================================
        //Updating the (converged) free variable vector
        //====================================================================
        // CM of EML2
        conv_free_var[0] = orbit_EM.si[0];
        conv_free_var[1] = orbit_EM.si[2];

        //Patch point
        for(int p = 1; p < man_grid_size; p++)
        {
            conv_free_var[0 + 4*p-2] = y_traj_n[0][p];
            conv_free_var[1 + 4*p-2] = y_traj_n[1][p];
            conv_free_var[2 + 4*p-2] = y_traj_n[3][p];
            conv_free_var[3 + 4*p-2] = y_traj_n[4][p];
        }

        //CMS of SEML2
        conv_free_var[4*man_grid_size-2] = orbit_SEM.si[0];
        conv_free_var[4*man_grid_size-1] = orbit_SEM.si[2];
        conv_free_var[4*man_grid_size-0] = orbit_SEM.si[4];


        //======================================================================
        //Updating the free variables
        //======================================================================
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += ds*nullvector[0];
        orbit_EM.si[2] += ds*nullvector[1];

        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
        //To CM_EM_NCSEM coordinates
        for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];


        //The middle (patch points) is classical cartesian coordinates at patch points
        for(int k = 1; k < man_grid_size; k++)
        {
            y_traj_n[0][k] += ds*nullvector[0 + 4*k-2];
            y_traj_n[1][k] += ds*nullvector[1 + 4*k-2];
            y_traj_n[3][k] += ds*nullvector[2 + 4*k-2];
            y_traj_n[4][k] += ds*nullvector[3 + 4*k-2];
        }

        //Last 4 correction variables is orbit.si
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += ds*nullvector[4*man_grid_size-2];
        orbit_SEM.si[2] += ds*nullvector[4*man_grid_size-1];
        orbit_SEM.si[4] += ds*nullvector[4*man_grid_size-0];

        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];

        //======================================================================
        //Diff corr
        //======================================================================
        msdvt_CMS_RCM_deps_planar_pac_ATF(y_traj_n, t_traj_n, y_traj_n, t_traj_n, nullvector, conv_free_var, ds, 42,
                                          traj_grid_size, coord_type, 5e-8, false,
                                          DCM_EM_TFC, DCMS_SEM_TFC,
                                          CCM_R_RCM_EM,
                                          CCM_R_RCM_SEM,
                                          orbit_EM, orbit_SEM,
                                          h2, isPlotted, false);

        //======================================================================
        //Display
        //======================================================================
        dkn = (double) kn;
        gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);


        cout << "t_traj_n[end]" << t_traj_n[man_grid_size]<< endl;

        cout << "orbit_EM.si = " << endl;
        vector_printf_prec(orbit_EM.si, 5);

        cout << "orbit_SEM.si = " << endl;
        vector_printf_prec(orbit_SEM.si, 5);
        cout << "-------------------------------------------"  << endl;
        //printf("Press ENTER to go on\n");
        //scanf("%c",&ch);

    }

    //====================================================================================
    // 5. Compute the initial orbit
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial EML2 orbit             "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Reset the unstable direction
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;
    orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

    //---------------------------------------------------------------------
    // Gnuplot (NCEM)
    //---------------------------------------------------------------------
    cout << "orbit_EM.si = " << endl;
    vector_printf_prec(orbit_EM.si, 5);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    cout << "orbit_EM.z0 = " << endl;
    vector_printf_prec(orbit_EM.z0 , 6);

    //---------------------------------------------------------------------
    //Integration on mPlot+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
    gnuplot_plot_xyz(h6, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"EML2 orbit", "lines", "1", "2", 2);


    //====================================================================================
    // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial SEML2 orbit            "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    vector<Oftsc> Fh;
    Fh.reserve(5);
    for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

    //---------------------------------------------------------------------
    //For dot(s) = fh(s)
    //---------------------------------------------------------------------
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf_SEM;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh         = &Fh;
    rvf.ofs        = &AUX;
    rvf.order      = OFTS_ORDER;
    rvf.n          = orbit_SEM.n;
    rvf.reduced_nv = 5;

    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*rvf.reduced_nv;
    sys_fh.params    = &rvf;

    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


    //---------------------------------------------------------------------
    // Temp variables
    //---------------------------------------------------------------------
    double t0_SEM = t_traj_n[traj_grid_size];
    double t1_SEM = t0_SEM+tof_seml_SEM;
    double z[6];
    double t2 = t0_SEM;
    int k  = 0;
    double  s1ccm8[2*rvf.reduced_nv]; //CCM8

    //---------------------------------------------------------------------
    // Initial state in CCM8 form
    //---------------------------------------------------------------------
    RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

    //---------------------------------------------------------------------
    // Loop
    //---------------------------------------------------------------------
    while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
    {
        cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;

        //To NC coordinates
        RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                     orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

        //Save
        for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
        t_man_SEM[k] = t2;

        //Plot
        gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

        //Advance one step
        gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
        CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);

        k++;
    }

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);
    gnuplot_plot_xyz(h6, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"SEML2 orbit", "lines", "1", "2", 3);

    //---------------------------------------------------------------------
    //Old version using classic integrator
    //---------------------------------------------------------------------
    //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
    //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


    //====================================================================================
    // 7. Final trajectory, on a grid
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Final trajectory, on a grid                "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    double **ymc   = dmatrix(0, 5, 0, mPlot);
    double *tmc    = dvector(0, mPlot);

    //Final trajectory on lines, segment by segment
    for(int k = 0; k < traj_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);

        //Plot on h6
        if(k == 0) gnuplot_plot_xyz(h6, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Transfer leg", "lines", "1", "2", 0);
        else gnuplot_plot_xyz(h6, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
    }

    //---------------------------------------------------------------------
    // Final trajectory, whith single shooting integration
    //---------------------------------------------------------------------
    int kstart = 0; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
    ode78(y_traj, t_traj, t_traj_n[kstart], t2, yv, 6, traj_grid_size, dcs, coord_type, coord_type);

    //Plot on h2
    //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //====================================================================================
    //Plot the final trajectory
    //====================================================================================
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);

    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);
    gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}


/**
 *  \brief Same as ref1_CMU_EM_to_CMS_SEM_MSD_RCM_3D, in the planar case.
 **/
void ref1_CMU_EM_to_CMS_SEM_MSD_RCM_PLANAR_CONT_PAC(SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM,
        matrix<Oftsc> &DCM_EM_TFC,
        matrix<Oftsc> &DCMS_SEM_TFC,
        int dcs,
        int coord_type,
        int man_grid_size,
        gnuplot_ctrl *h2)
{
    //====================================================================================
    // 1. Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Local variables for plotting
    //---------------------------------------------------------------------
    int mPlot = 1000;
    double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
    double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
    double *t_man_EM           = dvector(0, mPlot);
    double *t_man_SEM          = dvector(0, mPlot);
    double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
    double *t_man_coord_plot   = dvector(0, mPlot);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    // Initialize the COC matrices
    //---------------------------------------------------------------------
    //CCM_R_RCM_EM: RCM to CCM in R(4,4), about EML2
    gsl_matrix_complex *CCM_R_RCM_EM  = gsl_matrix_complex_calloc(4, 4);
    for(int i = 0; i < 4; i++) gsl_matrix_complex_set(CCM_R_RCM_EM, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_EM, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_EM, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_EM, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_EM, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));

    //CCM_R_RCM_SEM: RCM to CCM in R(5,5), about SEML2
    gsl_matrix_complex *CCM_R_RCM_SEM  = gsl_matrix_complex_calloc(5, 5);
    for(int i = 0; i < 4; i++) gsl_matrix_complex_set(CCM_R_RCM_SEM, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_SEM, 4, 4, gslc_complex(1.0, 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_SEM, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_SEM, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_SEM, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_SEM, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));

    //---------------------------------------------------------------------
    // Time on each orbit
    //---------------------------------------------------------------------
    double tof_eml_EM   = 5*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 20*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2.  Compute the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, orbit_EM.tf, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }


    //====================================================================================
    // 3. Compute the first point of the final SEML2 orbit and add it to the manifold leg
    //====================================================================================
    //---------------------------------------------------------------------
    // Save the last point at position man_grid_size
    //---------------------------------------------------------------------
    double yout[6], tout;
    qbcp_coc(orbit_SEM.t0, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

    for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
    t_traj[man_grid_size] = tout;

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


    //====================================================================================
    // 4.  Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    char ch;

    cout << "-------------------------------------------"  << endl;
    cout << "Differential correction & final trajectory "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    int nfv = 5*man_grid_size+1;    //free variables
    double *nullvector    = dvector(0, nfv-1);
    double *conv_free_var = dvector(0, nfv-1);

    //======================================================================
    //GNUPLOT
    //======================================================================
    gnuplot_ctrl *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3,  "set title \"s3_EM vs s1_EM\" ");
    gnuplot_cmd(h3, "set grid");

    gnuplot_ctrl *h4;
    h4 = gnuplot_init();
    gnuplot_cmd(h4,  "set title \"s3_SEM vs s1_SEM\" ");
    gnuplot_cmd(h4, "set grid");

    gnuplot_ctrl *h5;
    h5 = gnuplot_init();
    gnuplot_cmd(h5,  "set title \"s5_SEM vs steps\" ");
    gnuplot_cmd(h5, "set grid");

    gnuplot_ctrl *h1;
    h1 = gnuplot_init();
    gnuplot_cmd(h1,  "set title \"NCEM focus\" ");
    gnuplot_cmd(h1,  "set grid");

    gnuplot_ctrl *h7;
    h7 = gnuplot_init();
    gnuplot_cmd(h7,  "set title \"Dt vs steps\" ");
    gnuplot_cmd(h7,  "set grid");


    //======================================================================
    //Continuation : Update the variables with the null vector
    //======================================================================
    msvtplan(y_traj, t_traj, y_traj_n, t_traj_n, nullvector, 42,
             traj_grid_size, coord_type, 5e-8, true,
             DCM_EM_TFC, DCMS_SEM_TFC,
             CCM_R_RCM_EM,
             CCM_R_RCM_SEM,
             orbit_EM, orbit_SEM,
             h2, isPlotted, false);

    cout << "-------------------------------------------"  << endl;
    cout << "Continuation"  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    double ds = 2e-1;
    double yv[42], ye[42];
    double dkn = 0.0, dts = 0.0;
    for(int kn = 0; kn < 100; kn++)
    {
        //====================================================================
        //Updating the (converged) free variable vector
        //====================================================================
        // CM of EML2
        conv_free_var[0] = orbit_EM.si[0];
        conv_free_var[1] = orbit_EM.si[2];

        //Patch point
        for(int p = 1; p < man_grid_size; p++)
        {
            conv_free_var[0 + 5*p-3] = y_traj_n[0][p];
            conv_free_var[1 + 5*p-3] = y_traj_n[1][p];
            conv_free_var[2 + 5*p-3] = y_traj_n[3][p];
            conv_free_var[3 + 5*p-3] = y_traj_n[4][p];
            conv_free_var[5*p+1]     = t_traj_n[p];
        }

        //Last time:
        conv_free_var[5*man_grid_size] = t_traj_n[man_grid_size];

        //CMS of SEML2
        conv_free_var[5*man_grid_size-3] = orbit_SEM.si[0];
        conv_free_var[5*man_grid_size-2] = orbit_SEM.si[2];
        conv_free_var[5*man_grid_size-1] = orbit_SEM.si[4];


        //====================================================================
        //Updating the free variables
        //====================================================================
        //Updating CM_EM_RCM coordinates
        orbit_EM.si[0] += ds*nullvector[0];
        orbit_EM.si[2] += ds*nullvector[1];

        //Updating CM_EM_NCEM coordinates
        orbit_update_ic(orbit_EM, orbit_EM.si, t_traj_n[0]/SEML.us_em.ns);
        //To CM_EM_NCSEM coordinates
        for(int i = 0; i < 42; i++) yv[i] = orbit_EM.z0[i];
        qbcp_coc(t_traj_n[0]/SEML.us_em.ns, yv, ye, NCEM, coord_type);
        for(int i = 0; i < 6; i++) y_traj_n[i][0] = ye[i];

        //The middle (patch points) is classical cartesian coordinates at patch points
        for(int k = 1; k < man_grid_size; k++)
        {
            y_traj_n[0][k] += ds*nullvector[0 + 5*k-3];
            y_traj_n[1][k] += ds*nullvector[1 + 5*k-3];
            y_traj_n[3][k] += ds*nullvector[2 + 5*k-3];
            y_traj_n[4][k] += ds*nullvector[3 + 5*k-3];

            t_traj_n[k] += ds*nullvector[5*k+1];
        }

        //Last time:
        t_traj_n[man_grid_size] += ds*nullvector[5*man_grid_size];

        //Last 4 correction variables is orbit.si
        //Updating CM_SEM_RCM coordinates
        orbit_SEM.si[0] += ds*nullvector[5*man_grid_size-3];
        orbit_SEM.si[2] += ds*nullvector[5*man_grid_size-2];
        orbit_SEM.si[4] += ds*nullvector[5*man_grid_size-1];


        //Updating CM_SEM_NCSEM coordinates
        orbit_update_ic(orbit_SEM, orbit_SEM.si, t_traj_n[man_grid_size]);

        //Updating CM_SEM_NCSEM coordinates (identity transformation is performed in qbcp_coc.
        for(int i = 0; i < 6; i++) y_traj_n[i][man_grid_size] = orbit_SEM.z0[i];


        //====================================================================
        //Diff corr
        //====================================================================
        msdvt_CMS_RCM_deps_planar_pac(y_traj_n, t_traj_n, y_traj_n, t_traj_n, nullvector, conv_free_var, ds, 42,
                                      traj_grid_size, coord_type, 5e-8, false,
                                      DCM_EM_TFC, DCMS_SEM_TFC,
                                      CCM_R_RCM_EM,
                                      CCM_R_RCM_SEM,
                                      orbit_EM, orbit_SEM,
                                      h2, isPlotted, false);

        //====================================================================
        //NCEM focus
        //====================================================================
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][0];
        qbcp_coc(t_traj_n[0], yv, ye, coord_type, NCEM);


        //====================================================================
        //Display
        //====================================================================
        dkn = (double) kn;
        dts = (t_traj_n[traj_grid_size] - t_traj_n[0])/SEML.us_sem.T;
        gnuplot_plot_xy(h3, &orbit_EM.si[0],  &orbit_EM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h4, &orbit_SEM.si[0], &orbit_SEM.si[2], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h5, &dkn, &orbit_SEM.si[4], 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h7, &dkn, &dts, 1, (char*)"", "points", "1", "2", 0);
        gnuplot_plot_xy(h1, &ye[0], &ye[1], 1, (char*)"", "points", "1", "2", 0);

        cout << "orbit_EM.si = " << endl;
        vector_printf_prec(orbit_EM.si, 5);

        cout << "orbit_SEM.si = " << endl;
        vector_printf_prec(orbit_SEM.si, 5);

        cout << "t_traj_n[end] - t_traj_n[end-1] = " << t_traj_n[man_grid_size] - t_traj_n[man_grid_size-1] << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

    }

    //====================================================================================
    // 5. Compute the initial orbit
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial EML2 orbit             "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Reset the unstable direction
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;
    orbit_EM.tf    = orbit_EM.t0-tof_eml_EM;

    //---------------------------------------------------------------------
    // Gnuplot (NCEM)
    //---------------------------------------------------------------------
    gnuplot_ctrl *h6;
    h6 = gnuplot_init();
    gnuplot_cmd(h6,  "set title \"EML2 orbit  in NCEM coordinates\" ");

    cout << "orbit_EM.si = " << endl;
    vector_printf_prec(orbit_EM.si, 5);

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    cout << "orbit_EM.z0 = " << endl;
    vector_printf_prec(orbit_EM.z0 , 6);

    //---------------------------------------------------------------------
    //Integration on mPlot+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 2);
    gnuplot_plot_xyz(h6, y_man_NCEM[0], y_man_NCEM[1],  y_man_NCEM[2], mPlot+1, (char*)"", "lines", "1", "1", 2);



    //====================================================================================
    // 6. Compute the final SEML2 orbit, via integration in the reduced coordinates
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Compute the initial SEML2 orbit            "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    vector<Oftsc> Fh;
    Fh.reserve(5);
    for(int i = 0; i < 5; i++) Fh.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
    readVOFTS_bin(Fh, SEML_SEM.cs.F_PMS+"rvf/fh");

    //---------------------------------------------------------------------
    //For dot(s) = fh(s)
    //---------------------------------------------------------------------
    RVF rvf;
    rvf.ofs_order = SEML.eff_nf_SEM;
    Ofsc AUX(rvf.ofs_order);
    rvf.fh         = &Fh;
    rvf.ofs        = &AUX;
    rvf.order      = OFTS_ORDER;
    rvf.n          = orbit_SEM.n;
    rvf.reduced_nv = 5;

    gsl_odeiv2_system sys_fh;
    sys_fh.function  = qbfbp_fh;
    sys_fh.jacobian  = NULL;
    sys_fh.dimension = 2*rvf.reduced_nv;
    sys_fh.params    = &rvf;

    const gsl_odeiv2_step_type *T_fh = gsl_odeiv2_step_rk8pd;
    gsl_odeiv2_driver *d_fh = gsl_odeiv2_driver_alloc_y_new (&sys_fh, T_fh, 1e-12, 1e-10, 1e-10);


    //---------------------------------------------------------------------
    // Temp variables
    //---------------------------------------------------------------------
    double t0_SEM = t_traj_n[traj_grid_size];
    double t1_SEM = t0_SEM+tof_seml_SEM;
    double z[6];
    double t2 = t0_SEM;
    int k  = 0;
    double  s1ccm8[2*rvf.reduced_nv]; //CCM8

    //---------------------------------------------------------------------
    // Initial state in CCM8 form
    //---------------------------------------------------------------------
    RCMtoCCM8(orbit_SEM.si, s1ccm8, 5);

    //---------------------------------------------------------------------
    // Loop
    //---------------------------------------------------------------------
    while(k <= mPlot && orbit_SEM.si[4] > 1e-5)
    {
        cout << "orbit_SEM.si[4] = " << orbit_SEM.si[4] << endl;

        //To NC coordinates
        RCMtoNCbyTFC(orbit_SEM.si, t2, orbit_SEM.n, orbit_SEM.order, orbit_SEM.ofs_order,
                     orbit_SEM.reduced_nv, *orbit_SEM.Wh, *orbit_SEM.PC, *orbit_SEM.V, z, true);

        //Save
        for(int i = 0; i < 6; i++) y_man_NCSEM[i][k] = z[i];
        t_man_SEM[k] = t2;

        //Plot
        gnuplot_plot_xyz(h2, &z[0], &z[1],  &z[2], 1, (char*)"", "points", "1", "1", 6);

        //Advance one step
        gsl_odeiv2_evolve_apply (d_fh->e, d_fh->c, d_fh->s, d_fh->sys, &t2, t1_SEM, &d_fh->h, s1ccm8);
        CCM8toRCM(s1ccm8, orbit_SEM.si, orbit_SEM.reduced_nv);

        k++;
    }

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, max(k-1, 0), NCSEM, coord_type);

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], max(k, 1), (char*)"", "lines", "1", "1", 6);

    //---------------------------------------------------------------------
    //Old version using classic integrator
    //---------------------------------------------------------------------
    //orbit_update_ic(orbit_SEM, orbit_SEM.si, t0_SEM);
    //trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);


    //====================================================================================
    // 7. Final trajectory, on a grid
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "Final trajectory, on a grid                "  << endl;
    cout << "-------------------------------------------"  << endl;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    double **ymc   = dmatrix(0, 5, 0, mPlot);
    double *tmc    = dvector(0, mPlot);

    //Final trajectory on lines, segment by segment
    for(int k = 0; k < traj_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 0);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 0);
    }

    //---------------------------------------------------------------------
    // Final trajectory, whith single shooting integration
    //---------------------------------------------------------------------
    int kstart = 0; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
    for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
    ode78(y_traj, t_traj, t_traj_n[kstart], t2, yv, 6, traj_grid_size, dcs, coord_type, coord_type);

    //Plot on h2
    //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // Free
    //====================================================================================
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);

    gnuplot_close(h3);
    gnuplot_close(h4);
    gnuplot_close(h5);
    gnuplot_close(h6);

    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, mPlot);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, mPlot);
    free_dvector(t_man_EM, 0, mPlot);
    free_dvector(t_man_SEM, 0, mPlot);
    free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
    free_dvector(t_man_coord_plot, 0, mPlot);

    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);

    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);

    gsl_matrix_complex_free(CCM_R_RCM_EM);
    gsl_matrix_complex_free(CCM_R_RCM_SEM);
}


//========================================================================================
//
//         Refinement of solutions: Complete trajectory
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP(int ofts_order,
                                  int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM)
{
    //====================================================================================
    //Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);
    //int coordsys = default_framework(coord_type);

    //====================================================================================
    // 0 . Define the time of flight on each orbit
    //====================================================================================
    double tof_eml_EM   = 5*SEML.us.T;        //TOF on EML2 orbit
    double tof_seml_SEM = 10*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 1. Get the size of the data from the sorted solutions
    //====================================================================================
    int number_of_sol0;
    getLengthIntSortedCU_bin(&number_of_sol0, ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    //---------------------------------------------------------------------
    //In fact, here, we select only ONE solution (!)
    //---------------------------------------------------------------------
    int number_of_sol = 0;

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

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
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        //Comp
        emNCPoints(0.0, semP_comp);
        emPlot(h3, semP_comp);
        break;

    case PEM:
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semNCPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
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
    //---------------------------------------------------------------------
    //To store data from the sorted solutions
    //---------------------------------------------------------------------
    double *label               = dvector(0, number_of_sol);
    double *t0_EM               = dvector(0, number_of_sol);
    double *tf_EM               = dvector(0, number_of_sol);
    double *s1_CMU_EM           = dvector(0, number_of_sol);
    double *s3_CMU_EM           = dvector(0, number_of_sol);
    double *s1_CM_SEM           = dvector(0, number_of_sol);
    double *s3_CM_SEM           = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_1 = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_2 = dvector(0, number_of_sol);
    double **final_state_CMU_SEM     = dmatrix(0, 5, 0, number_of_sol);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, number_of_sol);
    double **init_state_CMU_NCEM     = dmatrix(0, 5, 0, number_of_sol);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size = 3*man_grid_size;
    double **y_traj  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj  = dvector(0, traj_grid_size);

    double **y_traj_comp  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj_comp  = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);
    //TBC
    double **yma        = dmatrix(0, 41, 0, traj_grid_size);

    //====================================================================================
    // 4. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc>  CM_SEM_NC;     ///center manifold in NC coordinates
    vector<Oftsc> CM_SEM_TFC;     ///center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    ///COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    ///COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   ///COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   ///COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      ///COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);


    //====================================================================================
    // 5. Getting back the data
    //====================================================================================
    string filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    readIntProjSortCU_bin(filename, label, t0_EM, tf_EM, s1_CMU_EM, s3_CMU_EM, s1_CM_SEM, s3_CM_SEM,
                          init_state_CMU_NCEM, final_state_CMU_SEM, projected_state_CMU_SEM,
                          min_proj_dist_SEM_1, min_proj_dist_SEM_2, number_of_sol);


    //====================================================================================
    // 6. Loop. Only one solution for now!
    //====================================================================================
    int index = 0;
    COMPLETION = 0;
    for(int kpos = number_of_sol; kpos <= number_of_sol; kpos++)
    {

        cout << "-------------------------------------------"  << endl;
        cout << "  Refinement of EML2-SEML2 arc             "  << endl;
        cout << "-------------------------------------------"  << endl;
        cout << "Estimated error at patch point (km):       "  << endl;
        cout <<  min_proj_dist_SEM_1[kpos]*SEML.cs_sem.cr3bp.L << endl;
        cout << "-------------------------------------------"  << endl;

        //====================================================================================
        // 6.1 Initialize local variables
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        double yv[6], st[5];

        for(int i = 0; i < 6; i ++) yv[i] = init_state_CMU_NCEM[i][kpos];

        st[0] = s1_CMU_EM[kpos];
        st[1] = 0.0;
        st[2] = s3_CMU_EM[kpos];
        st[3] = 0.0;
        st[4] = 0.0;

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_NCEM   = dmatrix(0, 5, 0, man_grid_size);
        double **y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_size);
        double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
        double *t_man_EM      = dvector(0, man_grid_size);
        double *t_man_SEM     = dvector(0, man_grid_size);
        double *t_man_coord   = dvector(0, man_grid_size);

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs(OFS_ORDER);
        OdeStruct driver;
        SingleOrbit orbit;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver, T, T_root, 6, qbcp_vfn, &SEML_EM);
        //Init routine
        init_orbit(orbit, &CM_EM_NC, &CM_EM_TFC,
                   &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs, OFTS_ORDER, OFS_ORDER,
                   1, t0_EM[kpos], tf_EM[kpos], SEML.us.T/5, &driver, &SEML_EM);


        //====================================================================================
        // 6.2 Compute the initial orbit
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, t0_EM[kpos]);

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0_EM[kpos], t0_EM[kpos]-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_size, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, man_grid_size, NCEM, coord_type);

        //---------------------------------------------------------------------
        // Save All BUT the last point
        //---------------------------------------------------------------------
        for(int kman = 0; kman < man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][man_grid_size-kman];
            t_traj[kman] = t_man_coord[man_grid_size-kman];
        }


        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 5);



        //====================================================================================
        // 6.3 Compute the manifold leg
        //====================================================================================
        //---------------------------------------------------------------------
        // Integration
        //---------------------------------------------------------------------
        ode78(y_man_coord, t_man_coord, t0_EM[kpos], tf_EM[kpos], yv, 6, man_grid_size, dcs, NCEM, coord_type);

        //---------------------------------------------------------------------
        // Save All BUT the last point
        //---------------------------------------------------------------------
        for(int kman = 0; kman < man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][man_grid_size + kman] = y_man_coord[i][kman];
            t_traj[man_grid_size + kman] = t_man_coord[kman];
        }

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);

        //Final and projected points: only in NCSEM and VNCSEM for now!
        if(coord_type == VNCSEM || coord_type == NCSEM)
        {
            gnuplot_plot_xyz(h2, &final_state_CMU_SEM[0][kpos], &final_state_CMU_SEM[1][kpos], &final_state_CMU_SEM[2][kpos],  1, (char*)"", "points", "1", "3", 2);
            gnuplot_plot_xyz(h2, &projected_state_CMU_SEM[0][kpos], &projected_state_CMU_SEM[1][kpos], &projected_state_CMU_SEM[2][kpos], 1, (char*)"", "points", "1", "3", 8);
        }


        //====================================================================================
        // 6.4 Compute the final SEML2 orbit
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        st[0] = s1_CM_SEM[kpos];
        st[1] = 0.0;
        st[2] = s3_CM_SEM[kpos];
        st[3] = 0.0;

        //Initial time in SEM units
        double t0_SEM = tf_EM[kpos]*SEML.us_em.ns;

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        //Init ode structure
        init_ode_structure(&driver, T, T_root, 6, qbcp_vfn, &SEML_SEM);
        //Init routine
        init_orbit(orbit, &CM_SEM_NC, &CM_SEM_TFC,
                   &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs, OFTS_ORDER, OFS_ORDER,
                   1, t0_SEM, t0_SEM+SEML.us_sem.T, SEML.us_sem.T/5, &driver, &SEML_SEM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, t0_SEM);

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_size, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, man_grid_size, NCSEM, coord_type);


        //---------------------------------------------------------------------
        // Save ALL points
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][2*man_grid_size + kman] = y_man_coord[i][kman];
            t_traj[2*man_grid_size + kman] = t_man_coord[kman];
        }

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 6);

        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        int isPlotted   = 0;
        multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, traj_grid_size, coord_type, isPlotted, h2);
        //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, traj_grid_size, coord_type, isPlotted, isTimeFixed, h2);


        //---------------------------------------------------------------------
        // Final trajectory, on a grid
        //---------------------------------------------------------------------
        int mPlot = 500;
        int color;
        double **ymc      = dmatrix(0, 5, 0, mPlot);
        double *tmc       = dvector(0, mPlot);
        double **ymc_comp = dmatrix(0, 5, 0, mPlot);
        double *tmc_comp  = dvector(0, mPlot);

        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Color
            //color = (int) 2*t_traj_n[k]/SEML.us_sem.T+1;
            color = 0;

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);


            //Plot on h3
            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
        }

        //Final trajectory on points
        //gnuplot_plot_xyz(h2, ymdn[0], ymdn[1],  ymdn[2], traj_grid_size+1, (char*)"Final trajectory", "points", "1", "3", 7);


        //---------------------------------------------------------------------
        // Final trajectory, whith single shooting integration
        //---------------------------------------------------------------------
        int kstart = man_grid_size; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
        ode78(y_traj, t_traj, t_traj_n[kstart], t_traj_n[traj_grid_size], yv, 6, traj_grid_size, dcs, coord_type, coord_type);

        //Plot on h2
        //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //Plot on h3
        qbcp_coc_vec(y_traj, t_traj, y_traj_comp, t_traj_comp, traj_grid_size, coord_type, comp_type);
        //gnuplot_plot_xyz(h3, y_traj_comp[0], y_traj_comp[1],  y_traj_comp[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //---------------------------------------------------------------------
        // Display
        //---------------------------------------------------------------------
        displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*++index/(number_of_sol+1));


        //---------------------------------------------------------------------
        // Free
        //---------------------------------------------------------------------
        free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
        free_dvector(t_man_EM, 0, man_grid_size);
        free_dvector(t_man_SEM, 0, man_grid_size);
        free_dvector(t_man_coord, 0, man_grid_size);
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
        free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
        free_dvector(tmc_comp, 0, mPlot);

    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dvector(label, 0, number_of_sol);
    free_dvector(t0_EM, 0, number_of_sol);
    free_dvector(tf_EM, 0, number_of_sol);
    free_dvector(s1_CMU_EM, 0, number_of_sol);
    free_dvector(s3_CMU_EM, 0, number_of_sol);
    free_dvector(s1_CM_SEM, 0, number_of_sol);
    free_dvector(s3_CM_SEM, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_1, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_2, 0, number_of_sol);
    free_dmatrix(final_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(projected_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, number_of_sol);
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_comp, 0, traj_grid_size);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);
    free_dmatrix(yma, 0, 41, 0, traj_grid_size);

    return 0;
}

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP(int man_grid_size,
                                  int coord_type,
                                  SingleOrbit &orbit_EM,
                                  SingleOrbit &orbit_SEM)
{
    //====================================================================================
    //Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 1 . Define the time of flight on each orbit
    //====================================================================================
    double tof_eml_EM   = 0.1*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 0.1*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3, *h4;
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();

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

        //h3 displays the system in EM coordinates
        gnuplot_set_xlabel(h3, (char*) "Xem");
        gnuplot_set_ylabel(h3, (char*) "Yem");
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

        //h3 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h3, (char*) "xem");
        gnuplot_set_ylabel(h3, (char*) "yem");
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

        //h2 displays the system in EM coordinates
        gnuplot_set_xlabel(h2, (char*) "Xem");
        gnuplot_set_ylabel(h2, (char*) "Yem");
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

        //h2 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h2, (char*) "xem");
        gnuplot_set_ylabel(h2, (char*) "yem");
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
    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size = 3*man_grid_size;
    double **y_traj  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj  = dvector(0, traj_grid_size);

    double **y_traj_comp  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj_comp  = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);
    //TBC
    double **yma        = dmatrix(0, 41, 0, traj_grid_size);

    //====================================================================================
    // 6. Loop. Only one solution for now!
    //====================================================================================
    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;

    //====================================================================================
    // 6.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_NCEM   = dmatrix(0, 5, 0, man_grid_size);
    double **y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_size);
    double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
    double **y_man_comp   = dmatrix(0, 5, 0, man_grid_size);
    double *t_man_EM      = dvector(0, man_grid_size);
    double *t_man_SEM     = dvector(0, man_grid_size);
    double *t_man_coord   = dvector(0, man_grid_size);
    double *t_man_comp    = dvector(0, man_grid_size);


    //====================================================================================
    // 6.2 Compute the initial orbit
    //====================================================================================
    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_size, 1);

    //---------------------------------------------------------------------
    //To coord_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, man_grid_size, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp, t_man_comp, man_grid_size, NCEM, comp_type);

    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][man_grid_size-kman];
        t_traj[kman] = t_man_coord[man_grid_size-kman];
    }


    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 5);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], man_grid_size+1, (char*)"", "points", "1", "1", 5);

    //====================================================================================
    // 6.3 Compute the manifold leg
    //====================================================================================
    double tf_EM = orbit_SEM.t0/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    ode78(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

    //---------------------------------------------------------------------
    //To comp_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_coord, t_man_coord, y_man_comp, t_man_comp, man_grid_size, coord_type, comp_type);


    //---------------------------------------------------------------------
    // Save All BUT the last point
    //---------------------------------------------------------------------
    for(int kman = 0; kman < man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][man_grid_size + kman] = y_man_coord[i][kman];
        t_traj[man_grid_size + kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], man_grid_size+1, (char*)"", "points", "1", "1", 4);

    //====================================================================================
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_size, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, man_grid_size, NCSEM, coord_type);
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_comp, t_man_comp, man_grid_size, NCSEM, comp_type);


    //---------------------------------------------------------------------
    // Save ALL points
    //---------------------------------------------------------------------
    for(int kman = 0; kman <= man_grid_size; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][2*man_grid_size + kman] = y_man_coord[i][kman];
        t_traj[2*man_grid_size + kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 6);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], man_grid_size+1, (char*)"", "points", "1", "1", 6);


    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    int mPlot = 500;
    int color;
    double yv[6];
    double **ymc      = dmatrix(0, 5, 0, mPlot);
    double *tmc       = dvector(0, mPlot);
    double **ymc_comp = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp  = dvector(0, mPlot);

    //Initial trajectory on lines, segment by segment
    for(int k = 0; k < traj_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

        //Color
        //color = (int) 2*t_traj[k]/SEML.us_sem.T+1;
        color = 2;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);


        //Plot on h3
        qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
        if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, traj_grid_size, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, traj_grid_size, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    //    //Final trajectory on lines, segment by segment
    //    for(int k = 0; k < traj_grid_size; k++)
    //    {
    //        //Integration segment by segment
    //        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
    //        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);
    //
    //        //Color
    //        //color = (int) 2*t_traj_n[k]/SEML.us_sem.T+1;
    //        color = 0;
    //
    //        //Plot on h2
    //        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
    //        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //
    //
    //        //Plot on h3
    //        qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
    //        if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
    //        else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //    }

    char ch;
    printf("Press ENTER to go to JPL on\n");
    scanf("%c",&ch);

    int shift = man_grid_size;
    cout << "Search for best fit in JPL DE430..." << endl;
    //====================================================================================
    // Search for best fit in  JPL ephemerides
    //====================================================================================
    OdeStruct driver_JPL;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf, &SEML);

    //----------------------------------------------------------
    //Get the best fit at t = t_traj_n[0]. et0 is in seconds
    //----------------------------------------------------------
    double et0, et;
    qbcp2jpl(t_traj_n[shift], &et0, coord_type);


    cout << "Change of coordinates syn -> ecliptic..." << endl;
    //====================================================================================
    // Change of coordinates
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_jpl   = dmatrix(0, 41, 0, man_grid_size);
    double *et_traj_jpl    = dvector(0, man_grid_size);
    double **y_traj_jpl_n = dmatrix(0, 41, 0, man_grid_size);
    double *t_traj_jpl_n  = dvector(0, man_grid_size);

    double **y_jpl_temp = dmatrix(0, 41, 0, man_grid_size);
    double *t_jpl_temp  = dvector(0, man_grid_size);

    //---------------------------------------------------------------------
    //Compute the change of coord.: syn -> ecliptic
    //---------------------------------------------------------------------
    double B[3];
    double C[3][3], k, kCprim[3][3], Bprim[3];
    double yvs[6], YV[6];

    for(int p = 0; p <= man_grid_size; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj_n[p+shift]-t_traj_n[shift])/SEML.n_sem;

        //COC at time t
        init_ecl2synstate(et, B, C, &k, Bprim, kCprim, coord_type);

        //Take the COC
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][p+shift];
        qbcp_coc(t_traj_n[p], yv, yvs, coord_type, VSEM); //in VSEM coordinates!
        syn2eclstate(yvs, YV, B, C, k, Bprim, kCprim, SEML.n_sem);

        //Store
        for(int i = 0; i < 6; i++) y_traj_jpl[i][p] = YV[i];
        et_traj_jpl[p] = et;
    }


    cout << "Plotting Moon +  Earth + L2..." << endl;
    //====================================================================================
    // Plotting Moon +  Earth + L2
    //====================================================================================
    //---------------------------------------------------------------------
    //Position of the Moon +  Earth + L2
    //---------------------------------------------------------------------
    double **y_earth_spice = dmatrix(0, 5, 0, man_grid_size);
    double **y_moon_spice  = dmatrix(0, 5, 0, man_grid_size);
    double **y_l2_spice    = dmatrix(0, 5, 0, man_grid_size);
    double lt;

    for(int p = 0; p <= man_grid_size; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj_n[p+shift]-t_traj_n[shift])/SEML.n_sem;

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
    //Position of the Sun +  Earth + Initial guess in ecliptic coordinates
    //---------------------------------------------------------------------
    gnuplot_set_xlabel(h4, (char*) "X [km]");
    gnuplot_set_ylabel(h4, (char*) "Y [km]");
    gnuplot_set_zlabel(h4, (char*) "Z [km]");

    gnuplot_plot_xy(h4, y_earth_spice[0], y_earth_spice[1],  man_grid_size+1, (char*) "EARTH", "lines", "2", "2", 3);
    gnuplot_plot_xy(h4, y_moon_spice[0],  y_moon_spice[1],   man_grid_size+1, (char*) "MOON", "lines", "2", "2", 4);
    gnuplot_plot_xy(h4, y_traj_jpl[0],    y_traj_jpl[1],     man_grid_size+1, (char*) "Initial Guess", "lines", "2", "2", 5);
    gnuplot_plot_xy(h4, y_l2_spice[0],    y_l2_spice[1],     man_grid_size+1, (char*) "L2", "lines", "2", "2", 6);

    //---------------------------------------------------------------------
    //Position of the Sun +  Earth + Initial guess in NCSEM coordinates
    //---------------------------------------------------------------------
    ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, man_grid_size, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  man_grid_size+1, (char*) "EARTH", "points", "3", "2", 8);

    ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, man_grid_size, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  man_grid_size+1, (char*) "MOON", "points", "4", "2", 8);

    ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, man_grid_size, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  man_grid_size+1, (char*) "L2", "points", "5", "2", 8);

    printf("Press ENTER to plot the Initial Guess\n");
    scanf("%c",&ch);
    cout << "Plotting Initial Guess..." << endl;
    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    //---------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < man_grid_size; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl[i][k];
        ode78_jpl(ymc_comp, tmc_comp, et_traj_jpl[k], et_traj_jpl[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

        //Back to NCSEM coordinates
        ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, t_traj[shift], eph_coord(coord_type));

        //Color
        color = 3;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    printf("Press ENTER to refine\n");
    scanf("%c",&ch);

    cout << "Refine trajectory in JPL ephemerides..." << endl;
    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    isPlotted   = 1;
    multiple_shooting_direct(y_traj_jpl, et_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, man_grid_size, VECLI, isPlotted, h4);
    //Plot
    gnuplot_plot_xy(h4, y_traj_jpl_n[0], y_traj_jpl_n[1],  man_grid_size+1, (char*) "Final trajectory", "lines", "2", "2", 7);

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_comp, 0, traj_grid_size);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);
    free_dmatrix(yma, 0, 41, 0, traj_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
    free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
    free_dvector(t_man_EM, 0, man_grid_size);
    free_dvector(t_man_SEM, 0, man_grid_size);
    free_dvector(t_man_coord, 0, man_grid_size);
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    return 0;
}

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID(int coord_type,
        SingleOrbit &orbit_EM,
        SingleOrbit &orbit_SEM)
{
    //====================================================================================
    //Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 1 . Define the time of flight on each orbit
    //====================================================================================
    double tof_eml_EM   =  1*SEML.us.T;       //TOF on EML2 orbit
    double tof_seml_SEM = 10*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3, *h4;
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();

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

        //h3 displays the system in EM coordinates
        gnuplot_set_xlabel(h3, (char*) "Xem");
        gnuplot_set_ylabel(h3, (char*) "Yem");
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

        //h3 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h3, (char*) "xem");
        gnuplot_set_ylabel(h3, (char*) "yem");
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

        //h2 displays the system in EM coordinates
        gnuplot_set_xlabel(h2, (char*) "Xem");
        gnuplot_set_ylabel(h2, (char*) "Yem");
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

        //h2 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h2, (char*) "xem");
        gnuplot_set_ylabel(h2, (char*) "yem");
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
    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int max_grid_size = 1000;
    double **y_traj       = dmatrix(0, 41, 0, max_grid_size);
    double  *t_traj       = dvector(0, max_grid_size);
    double **y_traj_comp  = dmatrix(0, 41, 0, max_grid_size);
    double  *t_traj_comp  = dvector(0, max_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, max_grid_size);
    double *t_traj_n    = dvector(0, max_grid_size);
    //TBC
    double **yma        = dmatrix(0, 41, 0, max_grid_size);


    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of EML2-SEML2 arc             "  << endl;
    cout << "-------------------------------------------"  << endl;

    //====================================================================================
    // 6.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_NCEM   = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_NCSEM  = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_coord  = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_comp   = dmatrix(0, 5, 0, max_grid_size);
    double *t_man_EM      = dvector(0, max_grid_size);
    double *t_man_SEM     = dvector(0, max_grid_size);
    double *t_man_coord   = dvector(0, max_grid_size);
    double *t_man_comp    = dvector(0, max_grid_size);


    //====================================================================================
    // 6.2 Compute the initial orbit
    //====================================================================================
    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on variable grid
    //---------------------------------------------------------------------
    int em_index = trajectory_integration_variable_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, max_grid_size, 1);

    //---------------------------------------------------------------------
    //To coord_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, em_index, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp, t_man_comp, em_index, NCEM, comp_type);


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
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], em_index+1, (char*)"EM first guess", "lines", "1", "1", 5);
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], em_index+1, (char*)"EM first guess", "points", "1", "1", 5);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], em_index+1, (char*)"", "points", "1", "1", 5);

    //====================================================================================
    // 6.3 Compute the manifold leg
    //====================================================================================
    double tf_EM = orbit_SEM.t0/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    int man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, max_grid_size, dcs, NCEM, coord_type, -1);

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
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_index+1, (char*)"", "lines", "1", "1", 4);
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_index+1, (char*)"", "points", "1", "1", 4);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], man_index+1, (char*)"", "points", "1", "1", 4);

    //====================================================================================
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, max_grid_size, 1);
    cout << "sem_index = " << sem_index << endl;

    //---------------------------------------------------------------------
    //To SEM coordinates
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
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"", "lines", "1", "1", 6);
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"", "points", "1", "1", 6);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], sem_index+1, (char*)"", "points", "1", "1", 6);

    //====================================================================================
    // Entire size is: no more than a hundred points or so.
    //====================================================================================
    int final_index = em_index + man_index + sem_index;
    int freq = final_index/60;
    cout << "----------------------------------"<< endl;
    cout << "Number of MSD points:        "     << endl;
    cout << "At EML:           " << em_index    << endl;
    cout << "During coast arc: " << man_index   << endl;
    cout << "At SEML:          " << sem_index   << endl;
    cout << "Total:            " << final_index << endl;
    cout << "----------------------------------"<< endl;
    cout << " Selecting 1/" << freq << " points" << endl;
    cout << "Total:            " << floor(final_index/freq) << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // Subset of points
    //====================================================================================
    for(int kman = 0; kman <= floor(final_index/freq); kman++)
    {
        for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
        t_traj_comp[kman] = t_traj[freq*kman];
    }

    final_index = floor(final_index/freq);

    for(int kman = 0; kman <= final_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
        t_traj[kman] = t_traj_comp[kman];
    }

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    int mPlot = 100;
    int color;
    double yv[6];
    double **ymc      = dmatrix(0, 5, 0, mPlot);
    double *tmc       = dvector(0, mPlot);
    double **ymc_comp = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp  = dvector(0, mPlot);

    //    //Initial trajectory on lines, segment by segment
    //    for(int k = 0; k < final_index; k++)
    //    {
    //        //Integration segment by segment
    //        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
    //        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);
    //
    //        //Color
    //        //color = (int) 2*t_traj[k]/SEML.us_sem.T+1;
    //        color = 2;
    //
    //        //Plot on h2
    //        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess", "lines", "1", "2", color);
    //        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //
    //
    //        //Plot on h3
    //        qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
    //        if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess", "lines", "1", "2", color);
    //        else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //    }

    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 0;
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, traj_grid_size, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    //        //Final trajectory on lines, segment by segment
    //        for(int k = 0; k < final_index; k++)
    //        {
    //            //Integration segment by segment
    //            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
    //            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);
    //
    //            //Color
    //            //color = (int) 2*t_traj_n[k]/SEML.us_sem.T+1;
    //            color = 0;
    //
    //            //Plot on h2
    //            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
    //            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //
    //
    //            //Plot on h3
    //            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
    //            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
    //            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    //        }

    printf("Press ENTER to go to JPL on\n");
    scanf("%c",&ch);


    //====================================================================================
    //Load kernels
    //====================================================================================
    furnsh_c ( "spice/kernels/metakernel.furnsh" );

    int shift = 0;
    cout << "Search for best fit in JPL DE430..." << endl;
    //====================================================================================
    // Search for best fit in  JPL ephemerides
    //====================================================================================
    OdeStruct driver_JPL;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf, &SEML);

    //----------------------------------------------------------
    //Get the best fit at t = t_traj_n[0]. et0 is in seconds
    //----------------------------------------------------------
    double et0, et;
    qbcp2jpl(t_traj_n[shift], &et0, coord_type);


    cout << "Change of coordinates syn -> ecliptic..." << endl;
    //====================================================================================
    // Change of coordinates
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_jpl   = dmatrix(0, 41, 0, final_index);
    double *et_traj_jpl    = dvector(0, final_index);
    double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double *t_traj_jpl_n  = dvector(0, final_index);

    double **y_jpl_temp = dmatrix(0, 41, 0, final_index);
    double *t_jpl_temp  = dvector(0, final_index);

    //---------------------------------------------------------------------
    //Compute the change of coord.: syn -> ecliptic
    //---------------------------------------------------------------------
    double B[3];
    double C[3][3], k, kCprim[3][3], Bprim[3];
    double yvs[6], YV[6];

    for(int p = 0; p <= final_index; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj_n[p+shift]-t_traj_n[shift])/SEML.n_sem;

        //COC at time t
        init_ecl2synstate(et, B, C, &k, Bprim, kCprim, coord_type);

        //Take the COC
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][p+shift];
        qbcp_coc(t_traj_n[p], yv, yvs, coord_type, VSEM); //in VSEM coordinates!
        syn2eclstate(yvs, YV, B, C, k, Bprim, kCprim, SEML.n_sem);

        //Store
        for(int i = 0; i < 6; i++) y_traj_jpl[i][p] = YV[i];
        et_traj_jpl[p] = et;
    }


    cout << "Plotting Moon +  Earth + L2..." << endl;
    //====================================================================================
    // Plotting Moon +  Earth + L2
    //====================================================================================
    //---------------------------------------------------------------------
    //Position of the Moon +  Earth + L2
    //---------------------------------------------------------------------
    double **y_earth_spice = dmatrix(0, 5, 0, final_index);
    double **y_moon_spice  = dmatrix(0, 5, 0, final_index);
    double **y_l2_spice    = dmatrix(0, 5, 0, final_index);
    double lt;

    for(int p = 0; p <= final_index; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj_n[p+shift]-t_traj_n[shift])/SEML.n_sem;

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
    //Position of the Sun +  Earth + Initial guess in ecliptic coordinates
    //---------------------------------------------------------------------
    gnuplot_set_xlabel(h4, (char*) "X [km]");
    gnuplot_set_ylabel(h4, (char*) "Y [km]");
    gnuplot_set_zlabel(h4, (char*) "Z [km]");

    //    gnuplot_plot_xy(h4, y_earth_spice[0], y_earth_spice[1],  final_index+1, (char*) "EARTH", "lines", "2", "2", 3);
    //    gnuplot_plot_xy(h4, y_moon_spice[0],  y_moon_spice[1],   final_index+1, (char*) "MOON", "lines", "2", "2", 4);
    //    gnuplot_plot_xy(h4, y_traj_jpl[0],    y_traj_jpl[1],     final_index+1, (char*) "Initial Guess", "lines", "2", "2", 5);
    //    gnuplot_plot_xy(h4, y_l2_spice[0],    y_l2_spice[1],     final_index+1, (char*) "L2", "lines", "2", "2", 6);

    //---------------------------------------------------------------------
    //Position of the Sun +  Earth + Initial guess in NCSEM coordinates
    //---------------------------------------------------------------------
    ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "EARTH", "points", "3", "2", 8);

    ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "MOON", "points", "4", "2", 8);

    ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "L2", "points", "5", "2", 8);


    printf("Press ENTER to plot the Initial Guess\n");
    scanf("%c",&ch);
    cout << "Plotting Initial Guess..." << endl;
    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    //---------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl[i][k];
        ode78_jpl(ymc_comp, tmc_comp, et_traj_jpl[k], et_traj_jpl[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

        //Back to NCSEM coordinates
        ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, t_traj[shift], eph_coord(coord_type));

        //Color
        color = 3;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    printf("Press ENTER to refine\n");
    scanf("%c",&ch);

    cout << "Refine trajectory in JPL ephemerides..." << endl;
    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    isPlotted   = 1;
    multiple_shooting_direct(y_traj_jpl, et_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, VECLI, isPlotted, h4);
    //Plot
    gnuplot_plot_xy(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


    //---------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
        ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

        //Back to NCSEM coordinates
        ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, t_traj[shift], eph_coord(coord_type));

        //Color
        color = 6;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final guess in JPL", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, max_grid_size);
    free_dvector(t_traj, 0, max_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, max_grid_size);
    free_dvector(t_traj_comp, 0, max_grid_size);
    free_dmatrix(y_traj_n, 0, 41, 0, max_grid_size);
    free_dvector(t_traj_n, 0, max_grid_size);
    free_dmatrix(yma, 0, 41, 0, max_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, max_grid_size);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, max_grid_size);
    free_dmatrix(y_man_coord, 0, 5, 0, max_grid_size);
    free_dvector(t_man_EM, 0, max_grid_size);
    free_dvector(t_man_SEM, 0, max_grid_size);
    free_dvector(t_man_coord, 0, max_grid_size);
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    return 0;
}

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_COMP_VARIABLE_GRID_TEST(int coord_type, SingleOrbit &orbit_SEM)
{
    //====================================================================================
    //Get the default coordinates system from the coord_type
    //====================================================================================
    //int dcs  = default_coordinate_system(coord_type);

    //====================================================================================
    // 1 . Define the time of flight on each orbit
    //====================================================================================
    double tof_seml_SEM = 10*SEML.us_sem.T;   //TOF on SEML2 orbit

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3, *h4;
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();

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

        //h3 displays the system in EM coordinates
        gnuplot_set_xlabel(h3, (char*) "Xem");
        gnuplot_set_ylabel(h3, (char*) "Yem");
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

        //h3 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h3, (char*) "xem");
        gnuplot_set_ylabel(h3, (char*) "yem");
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

        //h2 displays the system in EM coordinates
        gnuplot_set_xlabel(h2, (char*) "Xem");
        gnuplot_set_ylabel(h2, (char*) "Yem");
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

        //h2 displays the system in NCEM coordinates
        gnuplot_set_xlabel(h2, (char*) "xem");
        gnuplot_set_ylabel(h2, (char*) "yem");
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
    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int max_grid_size = 1000;
    double **y_traj       = dmatrix(0, 41, 0, max_grid_size);
    double  *t_traj       = dvector(0, max_grid_size);
    double **y_traj_comp  = dmatrix(0, 41, 0, max_grid_size);
    double  *t_traj_comp  = dvector(0, max_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    //TBC
    double **yma        = dmatrix(0, 41, 0, max_grid_size);


    cout << "-------------------------------------------"  << endl;
    cout << "  Refinement of SEML2 orbit                "  << endl;
    cout << "-------------------------------------------"  << endl;

    //====================================================================================
    // 6.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
    double **y_man_NCEM   = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_NCSEM  = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_coord  = dmatrix(0, 5, 0, max_grid_size);
    double **y_man_comp   = dmatrix(0, 5, 0, max_grid_size);
    double *t_man_EM      = dvector(0, max_grid_size);
    double *t_man_SEM     = dvector(0, max_grid_size);
    double *t_man_coord   = dvector(0, max_grid_size);
    double *t_man_comp    = dvector(0, max_grid_size);


    //====================================================================================
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, max_grid_size, 1);

    //---------------------------------------------------------------------
    //To SEM coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, sem_index, NCSEM, coord_type);
    qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_comp, t_man_comp, sem_index, NCSEM, comp_type);


    //---------------------------------------------------------------------
    // Save ALL points
    //---------------------------------------------------------------------
    for(int kman = 0; kman <= sem_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"", "lines", "1", "1", 6);
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"", "points", "1", "1", 6);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], sem_index+1, (char*)"", "points", "1", "1", 6);

    //====================================================================================
    // Entire size is: no more than a hundred points or so.
    //====================================================================================
    int final_index = sem_index;
    int freq = final_index/50;
    cout << "----------------------------------"<< endl;
    cout << "Number of MSD points:        "     << endl;
    cout << "At SEML:          " << sem_index   << endl;
    cout << "Total:            " << final_index << endl;
    cout << "----------------------------------"<< endl;
    cout << " Selecting 1/" << freq << " points" << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //====================================================================================
    // Subset of points
    //====================================================================================
    for(int kman = 0; kman <= floor(final_index/freq); kman++)
    {
        for(int i = 0; i < 6; i++) y_traj_comp[i][kman] = y_traj[i][freq*kman];
        t_traj_comp[kman] = t_traj[freq*kman];
    }

    final_index = floor(final_index/freq);

    for(int kman = 0; kman <= final_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_traj_comp[i][kman];
        t_traj[kman] = t_traj_comp[kman];
    }

    //====================================================================================
    // Initial trajectory, on a grid
    //====================================================================================
    int mPlot = 100;
    int color;
    double yv[6];
    double **ymc      = dmatrix(0, 5, 0, mPlot);
    double *tmc       = dvector(0, mPlot);
    double **ymc_comp = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp  = dvector(0, mPlot);

    printf("Press ENTER to go to JPL on\n");
    scanf("%c",&ch);

    //====================================================================================
    //Load kernels
    //====================================================================================
    furnsh_c ( "spice/kernels/metakernel.furnsh" );

    int shift = 0;
    cout << "Search for best fit in JPL DE430..." << endl;
    //====================================================================================
    // Search for best fit in  JPL ephemerides
    //====================================================================================
    OdeStruct driver_JPL;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf, &SEML);

    //----------------------------------------------------------
    //Get the best fit at t = t_traj_n[0]. et0 is in seconds
    //----------------------------------------------------------
    double et0, et;
    qbcp2jpl(t_traj[shift], &et0, coord_type);


    cout << "Change of coordinates syn -> ecliptic..." << endl;
    //====================================================================================
    // Change of coordinates
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_jpl   = dmatrix(0, 41, 0, final_index);
    double *et_traj_jpl    = dvector(0, final_index);
    double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
    double *t_traj_jpl_n  = dvector(0, final_index);

    double **y_jpl_temp = dmatrix(0, 41, 0, final_index);
    double *t_jpl_temp  = dvector(0, final_index);

    //---------------------------------------------------------------------
    //Compute the change of coord.: syn -> ecliptic
    //---------------------------------------------------------------------
    double B[3];
    double C[3][3], k, kCprim[3][3], Bprim[3];
    double yvs[6], YV[6];

    for(int p = 0; p <= final_index; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj[p+shift]-t_traj[shift])/SEML.n_sem;

        //COC at time t
        init_ecl2synstate(et, B, C, &k, Bprim, kCprim, coord_type);

        //Take the COC
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][p+shift];
        qbcp_coc(t_traj[p], yv, yvs, coord_type, VSEM); //in VSEM coordinates!
        syn2eclstate(yvs, YV, B, C, k, Bprim, kCprim, SEML.n_sem);

        //Store
        for(int i = 0; i < 6; i++) y_traj_jpl[i][p] = YV[i];
        et_traj_jpl[p] = et;
    }


    cout << "Plotting Moon +  Earth + L2..." << endl;
    //====================================================================================
    // Plotting Moon +  Earth + L2
    //====================================================================================
    //---------------------------------------------------------------------
    //Position of the Moon +  Earth + L2
    //---------------------------------------------------------------------
    double **y_earth_spice = dmatrix(0, 5, 0, final_index);
    double **y_moon_spice  = dmatrix(0, 5, 0, final_index);
    double **y_l2_spice    = dmatrix(0, 5, 0, final_index);
    double lt;

    for(int p = 0; p <= final_index; p++)
    {
        //Current time in seconds
        et = et0 + (t_traj[p+shift]-t_traj[shift])/SEML.n_sem;

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
    //Position of the Sun +  Earth + Initial guess in ecliptic coordinates
    //---------------------------------------------------------------------
    gnuplot_set_xlabel(h4, (char*) "X [km]");
    gnuplot_set_ylabel(h4, (char*) "Y [km]");
    gnuplot_set_zlabel(h4, (char*) "Z [km]");

    //    gnuplot_plot_xy(h4, y_earth_spice[0], y_earth_spice[1],  final_index+1, (char*) "EARTH", "lines", "2", "2", 3);
    //    gnuplot_plot_xy(h4, y_moon_spice[0],  y_moon_spice[1],   final_index+1, (char*) "MOON", "lines", "2", "2", 4);
    //    gnuplot_plot_xy(h4, y_traj_jpl[0],    y_traj_jpl[1],     final_index+1, (char*) "Initial Guess", "lines", "2", "2", 5);
    //    gnuplot_plot_xy(h4, y_l2_spice[0],    y_l2_spice[1],     final_index+1, (char*) "L2", "lines", "2", "2", 6);

    //---------------------------------------------------------------------
    //Position of the Sun +  Earth + Initial guess in NCSEM coordinates
    //---------------------------------------------------------------------
    ecl2coordstate_vec(y_earth_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "EARTH", "points", "3", "2", 8);

    ecl2coordstate_vec(y_moon_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "MOON", "points", "4", "2", 8);

    ecl2coordstate_vec(y_l2_spice, et_traj_jpl, y_jpl_temp, t_jpl_temp, final_index, coord_type, et0, t_traj[shift], eph_coord(coord_type));
    gnuplot_plot_xy(h2, y_jpl_temp[0], y_jpl_temp[1],  final_index+1, (char*) "L2", "points", "5", "2", 8);


    printf("Press ENTER to plot the Initial Guess\n");
    scanf("%c",&ch);
    cout << "Plotting Initial Guess..." << endl;
    //====================================================================================
    // Plotting Initial Guess
    //====================================================================================
    //---------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl[i][k];
        ode78_jpl(ymc_comp, tmc_comp, et_traj_jpl[k], et_traj_jpl[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

        //Back to NCSEM coordinates
        ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, t_traj[shift], eph_coord(coord_type));

        //Color
        color = 3;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    printf("Press ENTER to refine\n");
    scanf("%c",&ch);

    cout << "Refine trajectory in JPL ephemerides..." << endl;
    //====================================================================================
    // 6.5 Differential correction & final trajectory
    //====================================================================================
    int isPlotted   = 1;
    multiple_shooting_direct(y_traj_jpl, et_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, VECLI, isPlotted, h4);
    //Plot
    gnuplot_plot_xy(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


    //---------------------------------------------------------------------
    //Initial trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
        ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

        //Back to NCSEM coordinates
        ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, t_traj[shift], eph_coord(coord_type));

        //Color
        color = 6;

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final guess in JPL", "lines", "1", "2", color);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, max_grid_size);
    free_dvector(t_traj, 0, max_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, max_grid_size);
    free_dvector(t_traj_comp, 0, max_grid_size);
    free_dmatrix(yma, 0, 41, 0, max_grid_size);

    free_dmatrix(y_man_NCEM, 0, 5, 0, max_grid_size);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, max_grid_size);
    free_dmatrix(y_man_coord, 0, 5, 0, max_grid_size);
    free_dvector(t_man_EM, 0, max_grid_size);
    free_dvector(t_man_SEM, 0, max_grid_size);
    free_dvector(t_man_coord, 0, max_grid_size);
    free_dmatrix(ymc, 0, 5, 0, mPlot);
    free_dvector(tmc, 0, mPlot);
    free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
    free_dvector(tmc_comp, 0, mPlot);

    return 0;
}

/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the WHOLE trajectory (EML2 orbit + manifold leg + SEML2 orbit).
 *         It is supposed that orbit_EM and orbit_SEM has been refined with a previous computation, e.g. from ref_CMU_EM_to_CMS_SEM_MSD_RCM_SINGLE.
 **/
int comprefft3d(int man_grid_size_t,
                int coord_type,
                SingleOrbit &orbit_EM,
                SingleOrbit &orbit_SEM,
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
    double tof_eml_EM   = refst.tspan_EM;       //TOF on EML2 orbit
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
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int em_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);
    else //variable grid
        em_index = trajectory_integration_variable_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);

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

    double tf_EM = orbit_SEM.t0/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    int man_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        ode78(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type);
    else //variable grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type, -1);

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
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
    else //variable grid
        sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

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
    int mPlot = 10;

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
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        int fwrk_int  = choice == 0? I_ECLI:I_J2000;

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
            coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));//, SEML.ss);
            break;
        }

        cout << "y_traj_jpl = " << endl;
        cout << y_traj_jpl[0][0] << "  " << y_traj_jpl[0][1] << "  "  << y_traj_jpl[0][2] << "  "  << y_traj_jpl[0][3] << endl;


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
            coord2ecistate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0,  tsys0_comp, eph_coord(comp_type));//, SEML.ss);
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
                eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));//, SEML.ss);
                break;
            }

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);


            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //ECLIPTIC -> comp_type
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));
                break;

                case J2000:
                eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));//, SEML.ss);
                break;
            }


            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

        }


        //================================================================================
        // 6.5 Differential correction & final trajectory
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct_variable_time(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int, 5e-9, isPlotted, h4);
        //multiple_shooting_direct(y_traj_jpl, t_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, coord_int, isPlotted, h4);
        //Plot
        gnuplot_plot_xyz(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], y_traj_jpl_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //---------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, fwrk_int, coord_int, coord_int);

            //Back to coord_type coordinates
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));
                break;

                case J2000:
                eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));//, SEML.ss);
                break;
            }


            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //ECLIPTIC -> comp_type
            switch(coord_int)
            {
                case VECLI:
                ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));
                break;

                case J2000:
                eci2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));//, SEML.ss);
                break;
            }

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

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

/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_seml_dimjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_SEM,
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

    gnuplot_cmd(h2, "set view equal xyz");
    gnuplot_cmd(h3, "set view equal xyz");
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
    int max_grid    = 3000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;

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
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    cout << " comprefft3d. Compute the final SEML2 orbit..."  << endl;
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
    else //variable grid
        sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

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
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"SEML2", "lines", "1", "3", 6);
    //gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], sem_index+1, (char*)"", "points", "1", "3", 6);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], sem_index+1, (char*)"", "lines", "1", "3", 6);

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
    int mPlot = 10;

    double **ymc        = dmatrix(0, 5, 0, mPlot);
    double *tmc         = dvector(0, mPlot);
    double **ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp    = dvector(0, mPlot);
    double **ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_v       = dvector(0, mPlot*final_index);
    double **ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    double yv[6];
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        // Initialize SPICE kernerls & VF
        //================================================================================
        gnuplot_ctrl *h4 = gnuplot_init();
        cout << " comprefft3d. Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        cout << " comprefft3d. Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf, &SEML);

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

        //================================================================================
        // Change of coordinates
        //================================================================================
        cout << " comprefft3d. Change of coordinates syn -> ecliptic..." << endl;
        //---------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //---------------------------------------------------------------------
        double **y_traj_jpl   = dmatrix(0, 41, 0, final_index);
        double *et_traj_jpl    = dvector(0, final_index);
        double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl_n  = dvector(0, final_index);

        double **y_jpl_temp   = dmatrix(0, 41, 0, final_index);
        double *t_jpl_temp    = dvector(0, final_index);

        //---------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //---------------------------------------------------------------------
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, et_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));

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
            ode78_jpl(ymc_comp, tmc_comp, et_traj_jpl[k], et_traj_jpl[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

            //Back to coord_type coordinates
            ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);


            //Back to comp_type coordinates
            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 1);
        }


        //================================================================================
        // 6.5 Differential correction & final trajectory
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_jpl, et_traj_jpl, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, VECLI, isPlotted, h4);
        //Plot
        gnuplot_plot_xy(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //---------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, I_ECLI, VECLI, VECLI);

            //Back to coord_type coordinates
            ecl2coordstate_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_type, et0, tsys0, eph_coord(coord_type));

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

            //Back to comp_type coordinates
            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

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

/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_seml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_SEM,
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
    int max_grid    = 3000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;

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
    // 6.4 Compute the final SEML2 orbit
    //====================================================================================
    cout << " comprefft3d. Compute the final SEML2 orbit..."  << endl;
    //---------------------------------------------------------------------
    // Initialize the initial conditions (both NC and RCM coordinates)
    // We also need to kill the stable part.
    //---------------------------------------------------------------------
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
    else //variable grid
        sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

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
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
        t_traj[kman] = t_man_coord[kman];
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
        final_index = sem_index;
        int freq = final_index/(3*man_grid_size_t);
        cout << "----------------------------------"<< endl;
        cout << "Number of MSD points:        "     << endl;
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
    int mPlot = 10;
    double **ymc        = dmatrix(0, 5, 0, mPlot);
    double *tmc         = dvector(0, mPlot);
    double **ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp    = dvector(0, mPlot);
    double **ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_v       = dvector(0, mPlot*final_index);
    double **ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    double yv[6];
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        // Initialize SPICE kernerls
        //================================================================================
        gnuplot_ctrl *h4 = gnuplot_init();
        cout << " comprefft3d. Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //================================================================================
        // Initialize VF
        //================================================================================
        cout << " comprefft3d. Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &SEML);

        //================================================================================
        // Search for best fit in  JPL ephemerides
        //================================================================================
        cout << " comprefft3d. Search for best fit in JPL DE430..." << endl;

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

        //================================================================================
        // Change of coordinates
        //================================================================================
        cout << " comprefft3d. Change of coordinates syn -> ecliptic..." << endl;
        //---------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //---------------------------------------------------------------------
        double **y_traj_syn   = dmatrix(0, 41, 0, final_index);
        double *t_traj_syn    = dvector(0, final_index);
        double *et_traj_jpl   = dvector(0, final_index);

        double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl_n  = dvector(0, final_index);

        //---------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //---------------------------------------------------------------------
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
        cout << " comprefft3d. Plotting Initial Guess..." << endl;
        //---------------------------------------------------------------------
        //Initial trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 1);
        }

        //================================================================================
        // 6.5 Differential correction & final trajectory
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
        //Plot
        gnuplot_plot_xyz(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], y_traj_jpl_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //---------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 2);
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

/**
 *  \brief Computes only an EML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_eml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_EM,
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
    double tof_eml_EM = refst.tspan_EM;   //TOF on EML2 orbit
    //    double tof_seml_SEM = refst.tspan_SEM;   //TOF on SEML2 orbit

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
    int max_grid    = 3000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;

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
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int em_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);
    else //variable grid
        em_index = trajectory_integration_variable_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);

    //---------------------------------------------------------------------
    //To coord_type coordinates
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord, t_man_coord, em_index, NCEM, coord_type);
    qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_comp,  t_man_comp,  em_index, NCEM, comp_type);

    //---------------------------------------------------------------------
    // Save All points
    //---------------------------------------------------------------------
    for(int kman = 0; kman <= em_index; kman++)
    {
        for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][em_index-kman];
        t_traj[kman] = t_man_coord[em_index-kman];
    }

    //---------------------------------------------------------------------
    //Plot
    //---------------------------------------------------------------------
    gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], em_index+1, (char*)"", "lines", "5", "3", 5);
    //gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], em_index+1, (char*)"EML2", "points", "5", "3", 5);
    gnuplot_plot_xyz(h3, y_man_comp[0], y_man_comp[1],  y_man_comp[2], em_index+1, (char*)"", "lines", "5", "3", 5);


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
    int mPlot = 10;
    double **ymc        = dmatrix(0, 5, 0, mPlot);
    double *tmc         = dvector(0, mPlot);
    double **ymc_comp   = dmatrix(0, 5, 0, mPlot);
    double *tmc_comp    = dvector(0, mPlot);
    double **ymc_v      = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_v       = dvector(0, mPlot*final_index);
    double **ymc_comp_v = dmatrix(0, 5, 0, mPlot*final_index);
    double *tmc_comp_v  = dvector(0, mPlot*final_index);

    //Initial trajectory on lines, segment by segment
    double yv[6];
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n,  42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        // Initialize SPICE kernerls
        //================================================================================
        gnuplot_ctrl *h4 = gnuplot_init();
        cout << " comprefft3d. Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //================================================================================
        // Initialize VF
        //================================================================================
        cout << " comprefft3d. Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &SEML);

        //================================================================================
        // Search for best fit in  JPL ephemerides
        //================================================================================
        cout << " comprefft3d. Search for best fit in JPL DE430..." << endl;

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

        //================================================================================
        // Change of coordinates
        //================================================================================
        cout << " comprefft3d. Change of coordinates syn -> ecliptic..." << endl;
        //---------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //---------------------------------------------------------------------
        double **y_traj_syn   = dmatrix(0, 41, 0, final_index);
        double *t_traj_syn    = dvector(0, final_index);
        double *et_traj_jpl   = dvector(0, final_index);

        double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl_n  = dvector(0, final_index);

        //---------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //---------------------------------------------------------------------
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
        cout << " comprefft3d. Plotting Initial Guess..." << endl;
        //---------------------------------------------------------------------
        //Initial trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

        }

        //================================================================================
        // 6.5 Differential correction & final trajectory
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_jpl_n, t_traj_jpl_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
        //Plot
        gnuplot_plot_xyz(h4, y_traj_jpl_n[0], y_traj_jpl_n[1], y_traj_jpl_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //---------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //----------------------------------------------------
            //Integration segment by segment
            //----------------------------------------------------
            for(int i = 0; i < 6; i++) yv[i] = y_traj_jpl_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_jpl_n[k], t_traj_jpl_n[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //----------------------------------------------------
            //Back to coord_type coordinates
            //----------------------------------------------------
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

            //----------------------------------------------------
            //Back to comp_type coordinates
            //----------------------------------------------------
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            //OR
            //coord_type  -> ECLIPTIC -> comp_type
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 2);
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


/**
 *  \brief Computes only a SEML2 orbit and test a JPL refinement.
 **/
int comprefft3d_test_eml2seml_synjpl(int man_grid_size_t,
                     int coord_type,
                     SingleOrbit &orbit_EM,
                     SingleOrbit &orbit_SEM,
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
    double tof_eml_EM = refst.tspan_EM;      //TOF on EML2 orbit
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
    double** semP_coord = dmatrix(0, 6, 0, 2);
    double** semP_comp  = dmatrix(0, 6, 0, 2);
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
    int max_grid    = 10000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? 3*man_grid_size_t:max_grid;

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj  = dvector(0, traj_grid_2);

    double** y_traj_comp  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj_comp  = dvector(0, traj_grid_2);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, traj_grid_2);
    double* t_traj_n    = dvector(0, traj_grid_2);
    //TBC
    double** yma        = dmatrix(0, 41, 0, traj_grid_2);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================

    //====================================================================================
    // 4.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
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
    cout << " comprefft3d. Compute the initial orbit..."  << endl;
    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int em_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);
    else //variable grid
        em_index = trajectory_integration_variable_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);

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

    double tf_EM = orbit_SEM.t0/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    int man_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        ode78(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type);
    else //variable grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type, -1);

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
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
    else //variable grid
        sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

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
        final_index = em_index + man_index+ sem_index;
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
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
        // Initialize SPICE kernerls
        //================================================================================
        gnuplot_ctrl* h4 = gnuplot_init();
        cout << " comprefft3d. Initialize SPICE kernerls..." << endl;
        furnsh_c("spice/kernels/metakernel.furnsh");

        //================================================================================
        // Initialize VF
        //================================================================================
        cout << " comprefft3d. Initialize associated VF..." << endl;
        int shift = 0;
        OdeStruct driver_JPL;
        //Root-finding
        const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_JPL, T, T_root, 6, jpl_vf_syn, &SEML);

        //================================================================================
        // Search for best fit in  JPL ephemerides
        //================================================================================
        cout << " comprefft3d. Search for best fit in JPL DE430..." << endl;

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

        //================================================================================
        // Change of coordinates
        //================================================================================
        cout << " comprefft3d. Change of coordinates syn -> ecliptic..." << endl;
        //---------------------------------------------------------------------
        //Local variables to store the refined trajectory
        //---------------------------------------------------------------------
        double** y_traj_syn   = dmatrix(0, 41, 0, final_index);
        double* t_traj_syn    = dvector(0, final_index);
        double* et_traj_jpl   = dvector(0, final_index);
        double** y_traj_jpl   = dmatrix(0, 41, 0, final_index);
        double* t_traj_jpl    = dvector(0, final_index);
        double** y_traj_syn_n = dmatrix(0, 41, 0, final_index);
        double* t_traj_syn_n  = dvector(0, final_index);

        //---------------------------------------------------------------------
        //Compute the change of coord: syn -> ecliptic
        //---------------------------------------------------------------------
        //qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_syn, t_traj_syn, final_index, coord_type, eph_coord(coord_type));
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl, t_traj_jpl, final_index, coord_type, et0,  tsys0, eph_coord(coord_type));

        //================================================================================
        //Half of the trajectory is set in the other plane
        //================================================================================
        double **y_traj_jpl_n = dmatrix(0, 41, 0, final_index);
        double *t_traj_jpl_n  = dvector(0, final_index);
        coord2eclstate_vec(y_traj_n, t_traj_n, y_traj_jpl_n, t_traj_jpl_n, final_index, coord_type, et0, tsys0_comp, eph_coord(comp_type));

        //----------------------------------------------------------
        //Store the second half: at the beginning or at the end, depending on the computation coordinates
        //----------------------------------------------------------
        int begind = 0, endind = 0;
        switch(eph_coord(coord_type))
        {
                case VEM:
                begind = floor(final_index/2)-3;
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
        //Whole trajectory back in syn coordinates
        //================================================================================
        ecl2synstate_vec(y_traj_jpl, t_traj_jpl, y_traj_syn, t_traj_syn, final_index, et0, tsys0, eph_coord(coord_type));


        //---------------------------------------------------------------------
        //Time in seconds
        //---------------------------------------------------------------------
        for(int p = 0; p <= final_index; p++)
        {
            et_traj_jpl[p] = et0 + (t_traj_n[p]-tsys0)/mean_motion(eph_coord(coord_type));
        }

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
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_syn[k], t_traj_syn[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));

            //To ecl, then to NCSEM
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

            //Back to comp_type coordinates
            // @todo: put it NCEM  -> ECLIPTIC -> NCSEM, using ecl2coordstate_vec and its counterpart
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess in JPL", "lines", "1", "2", 1);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 1);
        }

        //================================================================================
        // 6.5 Differential correction & final trajectory
        //================================================================================
        printf("Press ENTER to refine\n");
        scanf("%c",&ch);
        cout << " comprefft3d. Refine trajectory in JPL ephemerides..." << endl;

        isPlotted   = 1;
        multiple_shooting_direct(y_traj_syn, t_traj_syn, y_traj_syn_n, t_traj_syn_n, 42, final_index, eph_coord(coord_type), isPlotted, h4);
        //Plot
        gnuplot_plot_xyz(h4, y_traj_syn_n[0], y_traj_syn_n[1], y_traj_syn_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


        //---------------------------------------------------------------------
        //Final trajectory on lines, segment by segment
        //---------------------------------------------------------------------
        for(int k = 0; k < final_index; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_syn_n[i][k];
            ode78_jpl(ymc_comp, tmc_comp, t_traj_syn_n[k], t_traj_syn_n[k+1], yv, 6, mPlot, eph_fwrk(coord_type), eph_coord(coord_type), eph_coord(coord_type));


            //Back to coord_type coordinates
            qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, eph_coord(coord_type), coord_type);

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

            //Back to comp_type coordinates
            //qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            syn2eclstate_vec(ymc_comp, tmc_comp, ymc, tmc,  mPlot, et0, tsys0, eph_coord(coord_type));
            ecl2coordstate_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, comp_type, et0, tsys0_comp, eph_coord(comp_type));

            //Plot on h3
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final traj in JPL", "lines", "1", "2", 2);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 2);
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

/**
 *  \brief Computes the whole trajectory in INSEM coordinates
 **/
int comprefft3d_test_eml2seml_insem(int man_grid_size_t,
                                    int coord_type,
                                    SingleOrbit& orbit_EM,
                                    SingleOrbit& orbit_SEM,
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
    double tof_eml_EM = refst.tspan_EM;      //TOF on EML2 orbit
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
    gnuplot_cmd(h2, "set grid");
    gnuplot_cmd(h3, "set grid");

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double** semP_coord = dmatrix(0, 6, 0, 2);
    double** semP_comp  = dmatrix(0, 6, 0, 2);
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
    int max_grid    = 10000;
    int man_grid_2  = (refst.grid == REF_FIXED_GRID) ? man_grid_size_t:max_grid;
    int traj_grid_2 = (refst.grid == REF_FIXED_GRID) ? 3*man_grid_size_t:max_grid;

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    double** y_traj  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj  = dvector(0, traj_grid_2);

    double** y_traj_comp  = dmatrix(0, 41, 0, traj_grid_2);
    double*  t_traj_comp  = dvector(0, traj_grid_2);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double** y_traj_n   = dmatrix(0, 41, 0, traj_grid_2);
    double* t_traj_n    = dvector(0, traj_grid_2);
    //TBC
    double** yma        = dmatrix(0, 41, 0, traj_grid_2);

    //====================================================================================
    // 4. Build the trajectory
    //====================================================================================

    //====================================================================================
    // 4.1 Initialize local variables
    //====================================================================================
    //---------------------------------------------------------------------
    //Local variables to store the manifold leg
    //---------------------------------------------------------------------
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
    cout << " comprefft3d. Compute the initial orbit..."  << endl;
    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we kill the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = 0.0;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int em_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);
    else //variable grid
        em_index = trajectory_integration_variable_grid(orbit_EM, orbit_EM.t0, orbit_EM.t0-tof_eml_EM, y_man_NCEM, t_man_EM, man_grid_2, 1);

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

    double tf_EM = orbit_SEM.t0/SEML.us_em.ns; //the end time is the starting time of the final orbit, in EM units

    //---------------------------------------------------------------------
    //For the computation of the initial orbit, we update the unstable part.
    //---------------------------------------------------------------------
    orbit_EM.si[4] = PROJ_EPSILON;

    //---------------------------------------------------------------------
    // Update the initial state in the orbit, with the RCM coordinates
    //---------------------------------------------------------------------
    orbit_update_ic(orbit_EM, orbit_EM.si, orbit_EM.t0);

    //---------------------------------------------------------------------
    // Integration
    //---------------------------------------------------------------------
    int man_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        ode78(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type);
    else //variable grid
        man_index = ode78_qbcp_vg(y_man_coord, t_man_coord, orbit_EM.t0, tf_EM, orbit_EM.z0, 6, man_grid_2, dcs, NCEM, coord_type, -1);

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
    orbit_SEM.si[4] = 0.0;
    orbit_update_ic(orbit_SEM, orbit_SEM.si, orbit_SEM.t0);


    //---------------------------------------------------------------------
    //Integration on man_grid_size+1 fixed grid
    //---------------------------------------------------------------------
    int sem_index = man_grid_2;
    if(refst.grid == REF_FIXED_GRID) //fixed grid
        trajectory_integration_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);
    else //variable grid
        sem_index = trajectory_integration_variable_grid(orbit_SEM, orbit_SEM.t0, orbit_SEM.t0+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_2, 1);

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
        final_index = em_index + man_index+ sem_index;
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
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj[i][k];
        ode78(ymc, tmc, t_traj[k], t_traj[k+1], yv, 6, mPlot, NCEM, coord_type, coord_type);

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
    multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, h2);
    //multiple_shooting_direct_variable_time(y_traj, t_traj, y_traj_n, t_traj_n, 42, final_index, coord_type, isPlotted, isTimeFixed, h2);

    //---------------------------------------------------------------------
    // Final trajectory, on a grid
    //---------------------------------------------------------------------
    cout << " comprefft3d. Final trajectory, on a grid..."  << endl;
    //Final trajectory on lines, segment by segment
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
        ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

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
    // QBCP INSEM
    //====================================================================================
    //----------------------------------------------------------
    //Go on
    //----------------------------------------------------------
    printf("Press ENTER to go on with the INSEM/INEM ref");
    scanf("%c",&ch);
    gnuplot_ctrl* h4 = gnuplot_init();

    //================================================================================
    // Select between INEM and INSEM
    //================================================================================
    int choice = 0;
    cout << "Enter 0 for INSEM, 1 for INEM: ";
    cin >> choice;
    int coord_int = choice == 0? INSEM:INEM;
    int fwrk_int  = choice == 0? I_INSEM:I_INEM;

    //================================================================================
    // Change of coordinates
    //================================================================================
    cout << " comprefft3d. Change of coordinates syn -> insem..." << endl;

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double** y_traj_insem   = dmatrix(0, 41, 0, final_index);
    double* t_traj_insem    = dvector(0, final_index);
    double** y_traj_insem_n = dmatrix(0, 41, 0, final_index);
    double* t_traj_insem_n  = dvector(0, final_index);

    //---------------------------------------------------------------------
    //Compute the change of coord: syn -> ecliptic
    //---------------------------------------------------------------------
    qbcp_coc_vec(y_traj_n, t_traj_n, y_traj_insem, t_traj_insem, final_index, coord_type, coord_int);


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
        for(int i = 0; i < 6; i++) yv[i] = y_traj_insem[i][k];
        ode78(ymc_comp, tmc_comp, t_traj_insem[k], t_traj_insem[k+1], yv, 6, mPlot, fwrk_int, coord_int, coord_int);

        //To coord_type
        qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_int, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial guess in INSEM/INEM", "lines", "1", "2", 1);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 1);

        //Back to comp_type coordinates
        qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);

        //Plot on h3
        if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Initial guess in INSEM/INEM", "lines", "1", "2", 1);
        else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 1);
    }

    //================================================================================
    // 6.5 Differential correction & final trajectory
    //================================================================================
    printf("Press ENTER to refine\n");
    scanf("%c",&ch);
    cout << " comprefft3d. Refine trajectory in INSEM/INEM..." << endl;

    isPlotted   = 1;
    multiple_shooting_direct(y_traj_insem, t_traj_insem, y_traj_insem_n, t_traj_insem_n, 42, final_index, coord_int, isPlotted, h4);
    //Plot
    gnuplot_plot_xyz(h4, y_traj_insem_n[0], y_traj_insem_n[1], y_traj_insem_n[2], final_index+1, (char*) "Final trajectory", "lines", "2", "2", 7);


    //---------------------------------------------------------------------
    //Final trajectory on lines, segment by segment
    //---------------------------------------------------------------------
    for(int k = 0; k < final_index; k++)
    {
        //Integration segment by segment
        for(int i = 0; i < 6; i++) yv[i] = y_traj_insem_n[i][k];
        ode78(ymc_comp, tmc_comp, t_traj_insem_n[k], t_traj_insem_n[k+1], yv, 6, mPlot, fwrk_int, coord_int, coord_int);

        //To coord_type
        qbcp_coc_vec(ymc_comp, tmc_comp, ymc, tmc, mPlot, coord_int, coord_type);

        //Plot on h2
        if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final traj in INSEM/INEM", "lines", "1", "2", 2);
        else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 2);

        //Back to comp_type coordinates
        qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);

        //Plot on h3
        if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final traj in INSEM/INEM", "lines", "1", "2", 2);
        else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", 2);
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

//========================================================================================
//
//         Refinement of solutions: CMU to CM
//
//========================================================================================
/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the PARTIAL trajectory (manifold leg + SEML2 orbit).
 *         The initial conditions are EML2 are allowed to move in the CMU of EML2, via its Fourier-Taylor parameterization. Note that the time at each point is NOT allowed to vary.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_PART(int ofts_order,
                                  int man_grid_size,
                                  int coord_type,
                                  vector<Oftsc> &CM_EM_NC,
                                  vector<Oftsc> &CM_EM_TFC,
                                  matrix<Oftsc> &DCM_EM_TFC,
                                  matrix<Ofsc>  &Mcoc_EM,
                                  matrix<Ofsc>  &MIcoc_EM,
                                  vector<Ofsc>  &Vcoc_EM)
{
    //====================================================================================
    // -1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);
    //int coordsys = default_framework(coord_type);

    //====================================================================================
    // 0. Time on each orbit
    //====================================================================================
    double tof_eml_EM   = 5*SEML.us.T;        //TOF on EML2 orbit
    double tof_seml_SEM = 50*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 1. Get the size of the data from the sorted solutions
    //====================================================================================
    int number_of_sol0;
    getLengthIntSortedCU_bin(&number_of_sol0, ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    //---------------------------------------------------------------------
    //In fact, here, we select only ONE solution (!)
    //---------------------------------------------------------------------
    int number_of_sol = 0;

    //====================================================================================
    // 1. Initialise the Jacobian matrices
    //====================================================================================
    //CCM_R_RCM: RCM to CCM
    gsl_matrix_complex *CCM_R_RCM  = gsl_matrix_complex_calloc(4, 4);
    for(int i = 0; i < 4; i++) gsl_matrix_complex_set(CCM_R_RCM, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    gsl_matrix_complex_set(CCM_R_RCM, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));


    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

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
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        //Comp
        emNCPoints(0.0, semP_comp);
        emPlot(h3, semP_comp);
        break;

    case PEM:
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semNCPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
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
    //---------------------------------------------------------------------
    //To store data from the sorted solutions
    //---------------------------------------------------------------------
    double *label               = dvector(0, number_of_sol);
    double *t0_EM               = dvector(0, number_of_sol);
    double *tf_EM               = dvector(0, number_of_sol);
    double *s1_CMU_EM           = dvector(0, number_of_sol);
    double *s3_CMU_EM           = dvector(0, number_of_sol);
    double *s1_CM_SEM           = dvector(0, number_of_sol);
    double *s3_CM_SEM           = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_1 = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_2 = dvector(0, number_of_sol);

    double **final_state_CMU_SEM     = dmatrix(0, 5, 0, number_of_sol);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, number_of_sol);
    double **init_state_CMU_NCEM     = dmatrix(0, 5, 0, number_of_sol);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = 2*man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);
    double **y_traj_comp  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj_comp  = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);
    double **yma        = dmatrix(0, 41, 0, traj_grid_size); //TBC: useful?

    //====================================================================================
    // 4. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CM_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_SEM_TFC;     //center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");

    cout << "PMS = " << SEML_SEM.cs.F_PMS << endl;
    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCM_SEM_TFC(6, 4, 4, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 5. Getting back the data
    //====================================================================================
    string filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    readIntProjSortCU_bin(filename, label, t0_EM, tf_EM, s1_CMU_EM, s3_CMU_EM, s1_CM_SEM, s3_CM_SEM,
                          init_state_CMU_NCEM, final_state_CMU_SEM, projected_state_CMU_SEM,
                          min_proj_dist_SEM_1, min_proj_dist_SEM_2, number_of_sol);


    //====================================================================================
    // 6. Loop. Only one solution for now!
    //====================================================================================
    //    int index = 0;
    for(int kpos = number_of_sol; kpos <= number_of_sol; kpos++)
    {

        cout << "-------------------------------------------"  << endl;
        cout << "  Refinement of EML2-SEML2 arc             "  << endl;
        cout << "-------------------------------------------"  << endl;
        cout << "Estimated error at patch point (km):       "  << endl;
        cout <<  min_proj_dist_SEM_1[kpos]*SEML.cs_sem.cr3bp.L << endl;
        cout << "Estimated error at patch point (SEMSU):    "  << endl;
        cout <<  min_proj_dist_SEM_1[kpos]                     << endl;
        cout << "-------------------------------------------"  << endl;

        char ch;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //====================================================================================
        // 6.1 Initialize local variables
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        double yv[6], st[5], st_SEM[4];

        for(int i = 0; i < 6; i ++) yv[i] = init_state_CMU_NCEM[i][kpos];

        st[0] = s1_CMU_EM[kpos];
        st[1] = 0.0;
        st[2] = s3_CMU_EM[kpos];
        st[3] = 0.0;
        st[4] = PROJ_EPSILON;

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
        double *t_man_coord   = dvector(0, man_grid_size);

        //---------------------------------------------------------------------
        //Local variables for plotting
        //---------------------------------------------------------------------
        int mPlot = 2000;
        double **y_man_NCEM        = dmatrix(0, 5, 0, mPlot);
        double **y_man_NCSEM       = dmatrix(0, 5, 0, mPlot);
        double *t_man_EM           = dvector(0, mPlot);
        double *t_man_SEM          = dvector(0, mPlot);
        double **y_man_coord_plot  = dmatrix(0, 5, 0, mPlot);
        double *t_man_coord_plot   = dvector(0, mPlot);

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs_EM(OFS_ORDER);
        OdeStruct driver_EM;
        SingleOrbit orbit_EM;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
        //Init routine
        init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
                   &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
                   1, t0_EM[kpos], tf_EM[kpos], SEML.us.T/5, &driver_EM, &SEML_EM);


        //====================================================================================
        // 6.2 Compute the manifold leg
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, st, t0_EM[kpos]);

        //---------------------------------------------------------------------
        // Integration
        //---------------------------------------------------------------------
        ode78(y_man_coord, t_man_coord, t0_EM[kpos], tf_EM[kpos], orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

        //---------------------------------------------------------------------
        // Save All BUT the last point
        //---------------------------------------------------------------------
        for(int kman = 0; kman < man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
            t_traj[kman] = t_man_coord[kman];
        }

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);

        //====================================================================================
        // 6.4 Compute the final SEML2 orbit
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        st_SEM[0] = s1_CM_SEM[kpos];
        st_SEM[1] = 0.0;
        st_SEM[2] = s3_CM_SEM[kpos];
        st_SEM[3] = 0.0;

        //Initial time in SEM units
        double t0_SEM = tf_EM[kpos]*SEML.us_em.ns;

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs_SEM(OFS_ORDER);
        OdeStruct driver_SEM;
        SingleOrbit orbit_SEM;

        //Init ode structure
        init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
        //Init routine
        init_orbit(orbit_SEM, &CM_SEM_NC, &CM_SEM_TFC,
                   &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
                   1, t0_SEM, t0_SEM+SEML.us_sem.T, SEML.us_sem.T/5, &driver_SEM, &SEML_SEM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, st_SEM, t0_SEM);

        cout << "s1_EM = " << st[0] << endl;
        cout << "s3_EM = " << st[2] << endl;

        cout << "s1 = " << st_SEM[0] << endl;
        cout << "s3 = " << st_SEM[2] << endl;
        cout << "projected state = " << projected_state_CMU_SEM[0][kpos] << projected_state_CMU_SEM[1][kpos] << endl;

        double ycc[6];
        NCSEMmtoSEMm(t0_SEM, orbit_SEM.z0, ycc, &SEML_SEM);
        cout << "projected state = " << ycc[0] << ycc[1] << endl;

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_SEM, t0_SEM, t0_SEM+tof_seml_SEM, y_man_NCSEM, t_man_SEM, man_grid_size, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord, t_man_coord, man_grid_size, NCSEM, coord_type);


        //---------------------------------------------------------------------
        // Save ALL points
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][man_grid_size + kman] = y_man_coord[i][kman];
            t_traj[man_grid_size + kman] = t_man_coord[kman];
        }

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 6);

        printf("Press ENTER to go on\n");
        scanf("%c",&ch);


        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        int isPlotted   = 1;
        cout << "-------------------------------------------"  << endl;
        cout << "Differential correction & final trajectory "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        msd_CM_RCM(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42,
                   traj_grid_size, coord_type, isPlotted, h2, Mcoc_EM,
                   DCM_EM_TFC, CCM_R_RCM, orbit_EM);



        //====================================================================================
        // 6.2 Compute the initial orbit
        //====================================================================================
        cout << "-------------------------------------------"  << endl;
        cout << "Compute the EML2 initial orbit:            "  << endl;
        cout << "-------------------------------------------"  << endl;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, t0_EM[kpos]);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, t0_EM[kpos], t0_EM[kpos]-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 5);


        //====================================================================================
        cout << "si prior to diffcorr = " << endl;
        for(int i = 0; i <5; i++) cout << st[i] << endl;

        cout << "si after diffcorr = " << endl;
        for(int i = 0; i <5; i++) cout << orbit_EM.si[i] << endl;

        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
        //====================================================================================

        //====================================================================================
        // 6.4 Compute the final SEML2 orbit
        //====================================================================================
        int ind_projection = floor(1.5*man_grid_size);
        //---------------------------------------------------------------------
        // Projection on the center manifold
        //---------------------------------------------------------------------
        for(int i = 0; i <6; i++) yv[i] = y_traj_n[i][ind_projection];
        double tv = t_traj_n[ind_projection];
        NCprojCCMtoCM(yv, tv, SEML_SEM.us.n, orbit_SEM.si, CM_SEM_TFC, MIcoc_SEM, Vcoc_SEM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tv);

        //---------------------------------------------------------------------
        // distance of projection
        //---------------------------------------------------------------------
        double dist_of_projection = 0.0;
        for(int i = 0; i <6; i++) dist_of_projection += (y_traj_n[i][ind_projection] - orbit_SEM.z0[i])*(y_traj_n[i][ind_projection] - orbit_SEM.z0[i]);
        dist_of_projection = sqrt(dist_of_projection);
        dist_of_projection *= SEML.cs_sem.gamma;

        cout << "New dist_of_projection (SEM) = " << dist_of_projection << endl;

        gnuplot_plot_xyz(h2, &y_traj_n[0][ind_projection], &y_traj_n[1][ind_projection],  &y_traj_n[2][ind_projection], 1, (char*)"yf", "points", "1", "5", 7);
        gnuplot_plot_xyz(h2, &orbit_SEM.z0[0],  &orbit_SEM.z0[1],   &orbit_SEM.z0[2], 1, (char*)"yp", "points", "1", "5", 8);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_SEM, tv, tv+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, mPlot, NCSEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"Initial SEML2 orbit", "lines", "1", "2", 7);

        cout << "si prior to diffcorr = " << endl;
        for(int i = 0; i <4; i++) cout << st_SEM[i] << endl;

        cout << "si after diffcorr = " << endl;
        for(int i = 0; i <4; i++) cout << orbit_SEM.si[i] << endl;


        printf("Press ENTER to go on\n");
        scanf("%c",&ch);


        //====================================================================================
        //====================================================================================


        //---------------------------------------------------------------------
        // Final trajectory, on a grid
        //---------------------------------------------------------------------
        int color;
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);

        double **ymc_comp = dmatrix(0, 5, 0, mPlot);
        double *tmc_comp  = dvector(0, mPlot);

        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Color
            //color = (int) 2*t_traj_n[k]/SEML.us_sem.T+1;
            color = 0;

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);


            //Plot on h3
            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
        }

        //Final trajectory on points
        //gnuplot_plot_xyz(h2, ymdn[0], ymdn[1],  ymdn[2], traj_grid_size+1, (char*)"Final trajectory", "points", "1", "3", 7);


        //---------------------------------------------------------------------
        // Final trajectory, whith single shooting integration
        //---------------------------------------------------------------------
        int kstart = man_grid_size; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
        ode78(y_traj, t_traj, t_traj_n[kstart], t_traj_n[traj_grid_size], yv, 6, traj_grid_size, dcs, coord_type, coord_type);

        //Plot on h2
        //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //Plot on h3
        qbcp_coc_vec(y_traj, t_traj, y_traj_comp, t_traj_comp, traj_grid_size, coord_type, comp_type);
        //gnuplot_plot_xyz(h3, y_traj_comp[0], y_traj_comp[1],  y_traj_comp[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //---------------------------------------------------------------------
        // Display
        //---------------------------------------------------------------------
        //displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*++index/(number_of_sol+1));

        //---------------------------------------------------------------------
        // Free
        //---------------------------------------------------------------------
        free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
        free_dvector(t_man_EM, 0, man_grid_size);
        free_dvector(t_man_SEM, 0, man_grid_size);
        free_dvector(t_man_coord, 0, man_grid_size);
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
        free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
        free_dvector(tmc_comp, 0, mPlot);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    char ch;
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dvector(label, 0, number_of_sol);
    free_dvector(t0_EM, 0, number_of_sol);
    free_dvector(tf_EM, 0, number_of_sol);
    free_dvector(s1_CMU_EM, 0, number_of_sol);
    free_dvector(s3_CMU_EM, 0, number_of_sol);
    free_dvector(s1_CM_SEM, 0, number_of_sol);
    free_dvector(s3_CM_SEM, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_1, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_2, 0, number_of_sol);
    free_dmatrix(final_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(projected_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, number_of_sol);

    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_comp, 0, traj_grid_size);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);
    free_dmatrix(yma, 0, 41, 0, traj_grid_size);
    return 0;
}


/**
 *  \brief Refine the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct scheme is applied on the manifold leg ONLY.
 **/
int ref_CMU_EM_to_CM_SEM_MSD(int ofts_order,
                             int man_grid_size,
                             int isPar)
{
    //====================================================================================
    // 1. Get the data from the sorted solutions
    //====================================================================================
    int number_of_sol0;
    getLengthIntSortedCU_bin(&number_of_sol0, ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    //---------------------------------------------------------------------
    //In fact, here, we select only ONE solution (!)
    //---------------------------------------------------------------------
    int number_of_sol = 0;

    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semNCPoints(0.0, semP);
    semPlot(h2, semP);

    //====================================================================================
    // 3. Init the data containers
    //====================================================================================
    double *label               = dvector(0, number_of_sol);
    double *t0_EM               = dvector(0, number_of_sol);
    double *tf_EM               = dvector(0, number_of_sol);
    double *s1_CMU_EM           = dvector(0, number_of_sol);
    double *s3_CMU_EM           = dvector(0, number_of_sol);
    double *s1_CM_SEM           = dvector(0, number_of_sol);
    double *s3_CM_SEM           = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_1 = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_2 = dvector(0, number_of_sol);


    double **final_state_CMU_SEM     = dmatrix(0, 5, 0, number_of_sol);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, number_of_sol);
    double **init_state_CMU_NCEM     = dmatrix(0, 5, 0, number_of_sol);


    //====================================================================================
    // 4. Getting back the data
    //====================================================================================
    string filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    readIntProjSortCU_bin(filename, label, t0_EM, tf_EM, s1_CMU_EM, s3_CMU_EM, s1_CM_SEM, s3_CM_SEM,
                          init_state_CMU_NCEM, final_state_CMU_SEM, projected_state_CMU_SEM,
                          min_proj_dist_SEM_1, min_proj_dist_SEM_2, number_of_sol);

    //====================================================================================
    // 5. Loop. Only one solution for now!
    //====================================================================================
    int index = 0;
    COMPLETION = 0;
    #pragma omp parallel for if(isPar) shared(index)
    for(int kpos = number_of_sol; kpos <= number_of_sol; kpos++)
    {
        //====================================================================================
        // 4.1 Initialize the initial conditions (both NC and RCM coordinates)
        //====================================================================================
        double yvv[6], yv[6];
        for(int i = 0; i < 6; i++) yv[i] = init_state_CMU_NCEM[i][kpos];

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **ymd   = dmatrix(0, 41, 0, man_grid_size);
        double *tmd    = dvector(0, man_grid_size);

        //---------------------------------------------------------------------
        //Local variables to store the refined manifold leg
        //---------------------------------------------------------------------
        double **ymdn   = dmatrix(0, 41, 0, man_grid_size);
        double **yma    = dmatrix(0, 41, 0, man_grid_size);
        double *tmdn    = dvector(0, man_grid_size);

        //====================================================================================
        // 4.2 Compute the manifold leg on a given grid
        //====================================================================================
        //---------------------------------------------------------------------
        //Computation with ode78, on the first 6 variables (the state)
        //---------------------------------------------------------------------
        ode78(ymd, tmd, t0_EM[kpos], tf_EM[kpos], yv, 6, man_grid_size, I_VNCSEM, NCEM, VNCSEM);

        //---------------------------------------------------------------------
        //Again in SEM coordinates?
        //Not necessary, same precision at the end.
        //---------------------------------------------------------------------
        //for(int i = 0; i < 6; i++) yv[i] = ymd[i][0];
        //ode78(ymd, tmd, tmd[0], tmd[man_grid_size], yv, 6, man_grid_size, I_VNCSEM, VNCSEM, VNCSEM);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        //Final position
        for(int i = 0; i <6; i++) yv[i] = final_state_CMU_SEM[i][kpos];
        qbcp_coc(tf_EM[kpos]*SEML.us_em.ns, yv, yvv, VSEM, VNCSEM);
        gnuplot_plot_xyz(h2, &yvv[0], &yvv[1], &yvv[2],  1, (char*)"yf", "points", "2", "5", 2);

        //Projected position
        for(int i = 0; i <6; i++) yv[i] = projected_state_CMU_SEM[i][kpos];
        qbcp_coc(tf_EM[kpos]*SEML.us_em.ns, yv, yvv, VSEM, VNCSEM);
        gnuplot_plot_xyz(h2, &yvv[0], &yvv[1], &yvv[2], 1, (char*)"yp", "points", "2", "5", 1);



        //====================================================================================
        // 4.3 Differential correction
        //====================================================================================
        //---------------------------------------------------------------------
        // Processing ymd. We include the target point at the end:
        // ymd(end) = yproj
        //---------------------------------------------------------------------
        for(int i = 0; i <6; i++) yv[i] = projected_state_CMU_SEM[i][kpos];
        qbcp_coc(tf_EM[kpos]*SEML.us_em.ns, yv, yvv, VSEM, VNCSEM); //to VNCSEM coordinates
        for(int i = 0; i <6; i++) ymd[i][man_grid_size] = yvv[i];

        //---------------------------------------------------------------------
        // After this point, the last point in ymd has been replaced
        // by the 6-dimensionnal target.
        //---------------------------------------------------------------------
        int isPlotted   = 1;
        //differential_correction_level_I(ymd, tmd, ymdn, tmdn, yma, 42, man_grid_size, isPlotted, isTimeFixed, h2);
        tic();
        //multiple_shooting_gomez(ymd, tmd, ymdn, tmdn, 42, man_grid_size, isPlotted, isTimeFixed, h2);
        multiple_shooting_direct(ymd, tmd, ymdn, tmdn, 42, man_grid_size, VNCSEM, isPlotted, h2);
        //multiple_shooting_direct_variable_time(ymd, tmd, ymdn, tmdn, yma, 42, man_grid_size, VNCSEM, isPlotted, isTimeFixed, h2);
        cout << "End of diffcorr in " << toc();

        //====================================================================================
        // 4.4 Costs
        //====================================================================================
        //--------------------------------------------------------------------------
        // Compute DV in m/s
        //--------------------------------------------------------------------------
        double DV[man_grid_size+1];
        for(int k = 0; k <= man_grid_size; k++)
        {
            DV[k] = 0.0;
            for(int i = 3; i < 6; i++)
            {
                DV[k] += (yma[i][k] - ymdn[i][k])*(yma[i][k] - ymdn[i][k]);
            }
            DV[k]  = sqrt(DV[k]);
            DV[k] *= 1e3*SEML.cs_sem.cr3bp.L*2*M_PI/SEML.cs_sem.cr3bp.T*SEML.cs_sem.cr3bp.l2.gamma_i;
        }


        //--------------------------------------------------------------------------
        // Compute DT in hours
        //--------------------------------------------------------------------------
        double DT[man_grid_size+1];
        for(int k = 0; k <= man_grid_size; k++)
        {
            DT[k] = SEML.cs_sem.cr3bp.T/(2*M_PI*3600)*fabs(tmdn[k] - tmd[k]);
        }

        cout << "------------------------------------" << endl;
        cout << "DV (m/s) = " << endl;
        vector_printf_prec(DV, man_grid_size+1);

        cout << "------------------------------------" << endl;
        cout << "DT (hours) = " << endl;
        vector_printf_prec(DT, man_grid_size+1);


        //====================================================================================
        // 4.5 Final trajectory:
        //====================================================================================
        int mPlot = 500;
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);

        for(int k = 0; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) yv[i] = ymdn[i][k];
            ode78(ymc, tmc, tmdn[k], tmdn[k+1], yv, 6, mPlot, I_VNCSEM, VNCSEM, VNCSEM);
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", 8);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", 8);
        }
        gnuplot_plot_xyz(h2, ymdn[0], ymdn[1],  ymdn[2], man_grid_size+1, (char*)"Final trajectory", "points", "1", "2", 8);



        //====================================================================================
        // 4.6 Initial trajectory:
        //====================================================================================
        for(int k = 0; k < man_grid_size; k++)
        {
            for(int i = 0; i < 6; i++) yv[i] = ymd[i][k];
            ode78(ymc, tmc, tmd[k], tmd[k+1], yv, 6, mPlot, I_VNCSEM, VNCSEM, VNCSEM);
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Initial trajectory", "lines", "dashed", "3", 7);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "dashed", "3", 7);
        }


        //---------------------------------------------------------------------
        // Free
        //---------------------------------------------------------------------
        free_dmatrix(ymd, 0, 5, 0, man_grid_size);
        free_dvector(tmd, 0, man_grid_size);
        free_dmatrix(ymdn, 0, 41, 0, man_grid_size);
        free_dmatrix(yma, 0, 41, 0, man_grid_size);
        free_dvector(tmdn, 0, man_grid_size);
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);

        //---------------------------------------------------------------------
        // Display completion
        //---------------------------------------------------------------------
        displayCompletion("readProjMan", 100.0*++index/(number_of_sol+1));
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dvector(label, 0, number_of_sol);
    free_dvector(t0_EM, 0, number_of_sol);
    free_dvector(tf_EM, 0, number_of_sol);
    free_dvector(s1_CMU_EM, 0, number_of_sol);
    free_dvector(s3_CMU_EM, 0, number_of_sol);
    free_dvector(s1_CM_SEM, 0, number_of_sol);
    free_dvector(s3_CM_SEM, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_1, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_2, 0, number_of_sol);
    free_dmatrix(final_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(projected_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, number_of_sol);


    return 0;
}


/**
 *  \brief Computes the best trajectories from int_proj_CMU_EM_on_CM_SEM. A multiple_shooting_direct is applied on the MANIFOLD trajectory (manifold leg).
 *         DEPRECATED. The initial conditions vary in the paramerization of the CMU of EML2. The final conditions vary in the paramerization of the CM of SEML2.
 *         The fact that the CM (and NOT the CMS of SEML2 is used) prevents this routine from converging. Abandoned for now.
 **/
int ref_CMU_EM_to_CM_SEM_MSD_DEP(int ofts_order,
                                 int man_grid_size,
                                 int coord_type,
                                 vector<Oftsc> &CM_EM_NC,
                                 vector<Oftsc> &CM_EM_TFC,
                                 matrix<Oftsc> &DCM_EM_TFC,
                                 matrix<Ofsc>  &Mcoc_EM,
                                 matrix<Ofsc>  &MIcoc_EM,
                                 vector<Ofsc>  &Vcoc_EM,
                                 int isDebug)
{
    //====================================================================================
    // -1. Get the default coordinates system from the coord_type
    //====================================================================================
    int dcs  = default_coordinate_system(coord_type);
    //int coordsys = default_framework(coord_type);

    //====================================================================================
    // 0. Time on each orbit
    //====================================================================================
    double tof_eml_EM   = 5*SEML.us.T;        //TOF on EML2 orbit
    double tof_seml_SEM = 10*SEML.us_sem.T;   //TOF on SEML2 orbit


    //====================================================================================
    // 1. Get the size of the data from the sorted solutions
    //====================================================================================
    int number_of_sol0;
    getLengthIntSortedCU_bin(&number_of_sol0, ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    //---------------------------------------------------------------------
    //In fact, here, we select only ONE solution (!)
    //---------------------------------------------------------------------
    int number_of_sol = 0;

    //====================================================================================
    // 1. Initialise the Jacobian matrices
    //====================================================================================
    //CCM_R_RCM: RCM to CCM
    gsl_matrix_complex *CCM_R_RCM  = gsl_matrix_complex_calloc(4, 4);
    for(int i = 0; i < 4; i++) gsl_matrix_complex_set(CCM_R_RCM, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    gsl_matrix_complex_set(CCM_R_RCM, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));


    //====================================================================================
    // 2. Init the gnuplot
    //====================================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl *h2, *h3;
    h2 = gnuplot_init();
    h3 = gnuplot_init();

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
        break;
    case NCSEM:
    case VNCSEM:
        semNCPoints(0.0, semP_coord);
        semPlot(h2, semP_coord);
        //Comp
        emNCPoints(0.0, semP_comp);
        emPlot(h3, semP_comp);
        break;

    case PEM:
        emPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
        break;
    case NCEM:
    case VNCEM:
        emNCPoints(0.0, semP_coord);
        emPlot(h2, semP_coord);
        //Comp
        semNCPoints(0.0, semP_comp);
        semPlot(h3, semP_comp);
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
    //---------------------------------------------------------------------
    //To store data from the sorted solutions
    //---------------------------------------------------------------------
    double *label               = dvector(0, number_of_sol);
    double *t0_EM               = dvector(0, number_of_sol);
    double *tf_EM               = dvector(0, number_of_sol);
    double *s1_CMU_EM           = dvector(0, number_of_sol);
    double *s3_CMU_EM           = dvector(0, number_of_sol);
    double *s1_CM_SEM           = dvector(0, number_of_sol);
    double *s3_CM_SEM           = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_1 = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_2 = dvector(0, number_of_sol);

    double **final_state_CMU_SEM     = dmatrix(0, 5, 0, number_of_sol);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, number_of_sol);
    double **init_state_CMU_NCEM     = dmatrix(0, 5, 0, number_of_sol);

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    int traj_grid_size    = man_grid_size;
    double **y_traj       = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj       = dvector(0, traj_grid_size);
    double **y_traj_comp  = dmatrix(0, 41, 0, traj_grid_size);
    double  *t_traj_comp  = dvector(0, traj_grid_size);

    //---------------------------------------------------------------------
    //Local variables to store the refined trajectory
    //---------------------------------------------------------------------
    double **y_traj_n   = dmatrix(0, 41, 0, traj_grid_size);
    double *t_traj_n    = dvector(0, traj_grid_size);
    double **yma        = dmatrix(0, 41, 0, traj_grid_size); //TBC: useful?

    //====================================================================================
    // 4. Structures to compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc> CM_SEM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_SEM_TFC;     //center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      //COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);

    //--------------------------------------
    //Jacobian of the center manifold
    //--------------------------------------
    matrix<Oftsc> DCM_SEM_TFC(6, 4, 4, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_SEM_TFC, SEML_SEM.cs.F_PMS+"DWf/DWhc");


    //====================================================================================
    // 5. Getting back the data
    //====================================================================================
    string filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    readIntProjSortCU_bin(filename, label, t0_EM, tf_EM, s1_CMU_EM, s3_CMU_EM, s1_CM_SEM, s3_CM_SEM,
                          init_state_CMU_NCEM, final_state_CMU_SEM, projected_state_CMU_SEM,
                          min_proj_dist_SEM_1, min_proj_dist_SEM_2, number_of_sol);


    //====================================================================================
    // 6. Loop. Only one solution for now!
    //====================================================================================
    //    int index = 0;
    for(int kpos = number_of_sol; kpos <= number_of_sol; kpos++)
    {

        cout << "-------------------------------------------"  << endl;
        cout << "  Refinement of EML2-SEML2 arc             "  << endl;
        cout << "-------------------------------------------"  << endl;
        cout << "Estimated error at patch point (km):       "  << endl;
        cout <<  min_proj_dist_SEM_1[kpos]*SEML.cs_sem.cr3bp.L << endl;
        cout << "Estimated error at patch point (SEMSU):    "  << endl;
        cout <<  min_proj_dist_SEM_1[kpos]                     << endl;
        cout << "-------------------------------------------"  << endl;

        char ch;
        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //====================================================================================
        // 6.1 Initialize local variables
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        double yv[6], st[5], st_SEM[4];

        for(int i = 0; i < 6; i ++) yv[i] = init_state_CMU_NCEM[i][kpos];

        st[0] = s1_CMU_EM[kpos];
        st[1] = 0.0;
        st[2] = s3_CMU_EM[kpos];
        st[3] = 0.0;
        st[4] = PROJ_EPSILON;

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_coord  = dmatrix(0, 5, 0, man_grid_size);
        double *t_man_coord   = dvector(0, man_grid_size);

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
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs_EM(OFS_ORDER);
        OdeStruct driver_EM;
        SingleOrbit orbit_EM;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver_EM, T, T_root, 6, qbcp_vfn, &SEML_EM);
        //Init routine
        init_orbit(orbit_EM, &CM_EM_NC, &CM_EM_TFC,
                   &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs_EM, OFTS_ORDER, OFS_ORDER,
                   1, t0_EM[kpos], tf_EM[kpos], SEML.us.T/5, &driver_EM, &SEML_EM);


        //====================================================================================
        // 6.2 Compute the manifold leg
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, st, t0_EM[kpos]);

        //---------------------------------------------------------------------
        // Integration
        //---------------------------------------------------------------------
        ode78(y_man_coord, t_man_coord, t0_EM[kpos], tf_EM[kpos], orbit_EM.z0, 6, man_grid_size, dcs, NCEM, coord_type);

        //---------------------------------------------------------------------
        // Save All BUT the last point
        //---------------------------------------------------------------------
        for(int kman = 0; kman < man_grid_size; kman++)
        {
            for(int i = 0; i < 6; i++) y_traj[i][kman] = y_man_coord[i][kman];
            t_traj[kman] = t_man_coord[kman];
        }

        //====================================================================================
        // 6.4 Compute the final SEML2 orbit
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        st_SEM[0] = s1_CM_SEM[kpos];
        st_SEM[1] = 0.0;
        st_SEM[2] = s3_CM_SEM[kpos];
        st_SEM[3] = 0.0;

        //Initial time in SEM units
        double t0_SEM = tf_EM[kpos]*SEML.us_em.ns;

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs_SEM(OFS_ORDER);
        OdeStruct driver_SEM;
        SingleOrbit orbit_SEM;

        //Init ode structure
        init_ode_structure(&driver_SEM, T, T_root, 6, qbcp_vfn, &SEML_SEM);
        //Init routine
        init_orbit(orbit_SEM, &CM_SEM_NC, &CM_SEM_TFC,
                   &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
                   1, t0_SEM, t0_SEM+SEML.us_sem.T, SEML.us_sem.T/5, &driver_SEM, &SEML_SEM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, st_SEM, t0_SEM);

        //---------------------------------------------------------------------
        // Save the last point at position man_grid_size
        //---------------------------------------------------------------------
        double yout[6], tout;
        qbcp_coc(t0_SEM, orbit_SEM.z0, yout, &tout, NCSEM, coord_type);

        for(int i = 0; i < 6; i++) y_traj[i][man_grid_size] = yout[i];
        t_traj[man_grid_size] = tout;

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord[0], y_man_coord[1],  y_man_coord[2], man_grid_size+1, (char*)"", "lines", "1", "1", 4);


        //====================================================================================
        // 6.5 Differential correction & final trajectory
        //====================================================================================
        int isPlotted   = 1;
        //multiple_shooting_direct(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42, traj_grid_size, coord_type, isPlotted, isTimeFixed, h2);

//                msd_CM_RCM(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42,
//                                                    traj_grid_size, coord_type, isPlotted, isTimeFixed, h2, Mcoc_EM,
//                                                    DCM_EM_TFC, CCM_R_RCM, orbit_EM);



        msdvt_CM_RCM(y_traj, t_traj, y_traj_n, t_traj_n, yma, 42,
                     traj_grid_size, coord_type, isPlotted, h2,
                     Mcoc_EM, Mcoc_SEM, Vcoc_SEM,
                     DCM_EM_TFC, DCM_SEM_TFC,
                     CM_SEM_TFC,
                     CCM_R_RCM,
                     orbit_EM, orbit_SEM, isDebug);


        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //====================================================================================
        // 6.2 Compute the initial orbit
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit, with the RCM coordinates
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_EM, orbit_EM.si, t0_EM[kpos]);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_EM, t0_EM[kpos], t0_EM[kpos]-tof_eml_EM, y_man_NCEM, t_man_EM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCEM, t_man_EM, y_man_coord_plot, t_man_coord_plot, mPlot, NCEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 5);


        //====================================================================================
        cout << "si prior to diffcorr = " << endl;
        for(int i = 0; i <5; i++) cout << st[i] << endl;

        cout << "si after diffcorr = " << endl;
        for(int i = 0; i <5; i++) cout << orbit_EM.si[i] << endl;

        printf("Press ENTER to go on\n");
        scanf("%c",&ch);
        //====================================================================================

        //====================================================================================
        // 6.4 Compute the final SEML2 orbit
        //====================================================================================
        //---------------------------------------------------------------------
        // Projection on the center manifold
        //---------------------------------------------------------------------
        for(int i = 0; i <6; i++) yv[i] = y_traj_n[i][traj_grid_size];
        double tv = t_traj_n[traj_grid_size];
        NCprojCCMtoCM(yv, tv, SEML_SEM.us.n, orbit_SEM.si, CM_SEM_TFC, MIcoc_SEM, Vcoc_SEM);

        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit_SEM, orbit_SEM.si, tv);

        //---------------------------------------------------------------------
        //Integration on mPlot+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit_SEM, tv, tv+tof_seml_SEM, y_man_NCSEM, t_man_SEM, mPlot, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_coord_plot, t_man_coord_plot, mPlot, NCSEM, coord_type);

        //---------------------------------------------------------------------
        //Plot
        //---------------------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_coord_plot[0], y_man_coord_plot[1],  y_man_coord_plot[2], mPlot+1, (char*)"", "lines", "1", "1", 6);

        cout << "si prior to diffcorr = " << endl;
        for(int i = 0; i <4; i++) cout << st_SEM[i] << endl;

        cout << "si after diffcorr = " << endl;
        for(int i = 0; i <4; i++) cout << orbit_SEM.si[i] << endl;

        //--------------------------------------
        // Center-stable manifold
        //--------------------------------------
        //        double stm[5];
        //        for(int i = 0; i < 4; i++) stm[i] = orbit_SEM.si[i];
        //
        //
        //        vector<Oftsc> CMS_SEM_NC;      //center manifold in NC coordinates
        //        vector<Oftsc> CMS_SEM_TFC;     //center manifold in TFC coordinates
        //        CMS_SEM_NC.reserve(6);
        //        for(int i = 0; i <6; i++) CMS_SEM_NC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
        //        CMS_SEM_TFC.reserve(6);
        //        for(int i = 0; i <6; i++) CMS_SEM_TFC.push_back(Oftsc(5, OFTS_ORDER, OFS_NV, OFS_ORDER));
        //        //Update CM
        //        readVOFTS_bin(CMS_SEM_NC,  "../OOFTDA/data/CS/QBCP/SEM/L2/W/W",  OFS_ORDER);
        //        readVOFTS_bin(CMS_SEM_TFC, "../OOFTDA/data/CS/QBCP/SEM/L2/W/Wh", OFS_ORDER);
        //
        //        //Init routine
        //        init_orbit(orbit_SEM, &CMS_SEM_NC, &CMS_SEM_TFC,
        //                   &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs_SEM, OFTS_ORDER, OFS_ORDER,
        //                   5, 1, t0_SEM, t0_SEM+SEML.us_sem.T, SEML.us_sem.T/5, &driver_SEM, &SEML_SEM);
        //
        //        cout << SEML_SEM.cs.F_PMS << endl;
        //        for(int i = 2; i < 15; i++)
        //        {
        //            for(int j = 0; j < 10; j++)
        //            {
        //                stm[4] = (double) j*pow(10, -i);
        //                orbit_update_ic(orbit_SEM, stm, tv);
        //                gnuplot_plot_xyz(h2, &orbit_SEM.z0[0], &orbit_SEM.z0[1],  &orbit_SEM.z0[2], 1, (char*)"", "points", "1", "3", 1);
        //
        //                stm[4] = (double) -j*pow(10, -i);
        //                orbit_update_ic(orbit_SEM, stm, tv);
        //                gnuplot_plot_xyz(h2, &orbit_SEM.z0[0], &orbit_SEM.z0[1],  &orbit_SEM.z0[2], 1, (char*)"", "points", "1", "3", 1);
        //            }
        //        }
        //
        //        printf("Press ENTER to go on\n");
        //        scanf("%c",&ch);


        //====================================================================================
        //====================================================================================
        //---------------------------------------------------------------------
        // Final trajectory, on a grid
        //---------------------------------------------------------------------
        int color;
        double **ymc   = dmatrix(0, 5, 0, mPlot);
        double *tmc    = dvector(0, mPlot);

        double **ymc_comp = dmatrix(0, 5, 0, mPlot);
        double *tmc_comp  = dvector(0, mPlot);

        //Final trajectory on lines, segment by segment
        for(int k = 0; k < traj_grid_size; k++)
        {
            //Integration segment by segment
            for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][k];
            ode78(ymc, tmc, t_traj_n[k], t_traj_n[k+1], yv, 6, mPlot, dcs, coord_type, coord_type);

            //Color
            //color = (int) 2*t_traj_n[k]/SEML.us_sem.T+1;
            color = 0;

            //Plot on h2
            if(k == 0) gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h2, ymc[0], ymc[1],  ymc[2], mPlot+1, (char*)"", "lines", "1", "2", color);


            //Plot on h3
            qbcp_coc_vec(ymc, tmc, ymc_comp, tmc_comp, mPlot, coord_type, comp_type);
            if(k == 0) gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"Final trajectory", "lines", "1", "2", color);
            else gnuplot_plot_xyz(h3, ymc_comp[0], ymc_comp[1],  ymc_comp[2], mPlot+1, (char*)"", "lines", "1", "2", color);
        }

        //Final trajectory on points
        //gnuplot_plot_xyz(h2, ymdn[0], ymdn[1],  ymdn[2], traj_grid_size+1, (char*)"Final trajectory", "points", "1", "3", 7);


        //---------------------------------------------------------------------
        // Final trajectory, whith single shooting integration
        //---------------------------------------------------------------------
        int kstart = man_grid_size; //starts at the CMU entry, and NOT at the EML2 orbit entry (kstart != 0)
        for(int i = 0; i < 6; i++) yv[i] = y_traj_n[i][kstart];
        ode78(y_traj, t_traj, t_traj_n[kstart], t_traj_n[traj_grid_size], yv, 6, traj_grid_size, dcs, coord_type, coord_type);

        //Plot on h2
        //gnuplot_plot_xyz(h2, y_traj[0], y_traj[1],  y_traj[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //Plot on h3
        qbcp_coc_vec(y_traj, t_traj, y_traj_comp, t_traj_comp, traj_grid_size, coord_type, comp_type);
        //gnuplot_plot_xyz(h3, y_traj_comp[0], y_traj_comp[1],  y_traj_comp[2], traj_grid_size+1, (char*)"Final trajectory (single shooting)", "lines", "1", "3", 8);

        //---------------------------------------------------------------------
        // Display
        //---------------------------------------------------------------------
        //displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*++index/(number_of_sol+1));

        //---------------------------------------------------------------------
        // Free
        //---------------------------------------------------------------------
        free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_coord, 0, 5, 0, man_grid_size);
        free_dvector(t_man_EM, 0, man_grid_size);
        free_dvector(t_man_SEM, 0, man_grid_size);
        free_dvector(t_man_coord, 0, man_grid_size);
        free_dmatrix(ymc, 0, 5, 0, mPlot);
        free_dvector(tmc, 0, mPlot);
        free_dmatrix(ymc_comp, 0, 5, 0, mPlot);
        free_dvector(tmc_comp, 0, mPlot);
        free_dmatrix(y_man_coord_plot, 0, 5, 0, mPlot);
        free_dvector(t_man_coord_plot, 0, mPlot);
    }

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    char ch;
    printf("Press ENTER to go close\n");
    scanf("%c",&ch);
    gnuplot_close(h2);

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dvector(label, 0, number_of_sol);
    free_dvector(t0_EM, 0, number_of_sol);
    free_dvector(tf_EM, 0, number_of_sol);
    free_dvector(s1_CMU_EM, 0, number_of_sol);
    free_dvector(s3_CMU_EM, 0, number_of_sol);
    free_dvector(s1_CM_SEM, 0, number_of_sol);
    free_dvector(s3_CM_SEM, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_1, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_2, 0, number_of_sol);
    free_dmatrix(final_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(projected_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, number_of_sol);

    free_dmatrix(semP_coord, 0, 6, 0, 2);
    free_dmatrix(semP_comp, 0, 6, 0, 2);
    free_dmatrix(y_traj, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj, 0, traj_grid_size);
    free_dmatrix(y_traj_comp, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_comp, 0, traj_grid_size);
    free_dmatrix(y_traj_n, 0, 41, 0, traj_grid_size);
    free_dvector(t_traj_n, 0, traj_grid_size);
    free_dmatrix(yma, 0, 41, 0, traj_grid_size);
    return 0;
}


//========================================================================================
//
//         Integration of sorted solutions
//
//========================================================================================
/**
 *  \brief Integrates the best trajectories from int_proj_CMU_EM_on_CM_SEM.
 **/
int int_sorted_sol_CMU_EM_to_CM_SEM(int ofts_order,
                                    int man_grid_size,
                                    int isPar,
                                    vector<Oftsc> &CM_EM_NC,
                                    vector<Oftsc> &CM_EM_TFC,
                                    matrix<Ofsc>  &Mcoc_EM,
                                    matrix<Ofsc>  &MIcoc_EM,
                                    vector<Ofsc>  &Vcoc_EM)
{
    //====================================================================================
    // 1. Get the data from the sorted solutions
    //====================================================================================
    int number_of_sol;
    getLengthIntSortedCU_bin(&number_of_sol, ofts_order, TYPE_MAN_SORT, SEML.li_SEM);
    cout << "number_of_sol = " << number_of_sol << endl;
    cout << "man_grid_size = " << man_grid_size << endl;

    //---------------------------------------------------------------------
    //To store final data
    //---------------------------------------------------------------------
    double ***y_orb_eml_SEM = d3tensor(0, 5, 0, number_of_sol, 0, man_grid_size+2);
    double **t_orb_eml_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);
    double **h_orb_eml_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);

    double ***y_orb_seml_SEM = d3tensor(0, 5, 0, number_of_sol, 0, man_grid_size+2);
    double **t_orb_seml_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);
    double **h_orb_seml_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);

    double ***y_cmu_SEM = d3tensor(0, 5, 0, number_of_sol, 0, man_grid_size+2);
    double **t_cmu_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);
    double **h_cmu_SEM  = dmatrix(0, number_of_sol, 0, man_grid_size+2);


    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);

    //====================================================================================
    // 2. Init the initial data containers
    //====================================================================================
    double *label               = dvector(0, number_of_sol);
    double *t0_EM               = dvector(0, number_of_sol);
    double *tf_EM               = dvector(0, number_of_sol);
    double *s1_CMU_EM           = dvector(0, number_of_sol);
    double *s3_CMU_EM           = dvector(0, number_of_sol);
    double *s1_CM_SEM           = dvector(0, number_of_sol);
    double *s3_CM_SEM           = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_1 = dvector(0, number_of_sol);
    double *min_proj_dist_SEM_2 = dvector(0, number_of_sol);

    double **final_state_CMU_SEM     = dmatrix(0, 5, 0, number_of_sol);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, number_of_sol);
    double **init_state_CMU_NCEM     = dmatrix(0, 5, 0, number_of_sol);


    //====================================================================================
    // 3. Getting back the data
    //====================================================================================
    string filename = filenameCUM(ofts_order, TYPE_MAN_SORT, SEML.li_SEM);

    readIntProjSortCU_bin(filename, label, t0_EM, tf_EM, s1_CMU_EM, s3_CMU_EM, s1_CM_SEM, s3_CM_SEM,
                          init_state_CMU_NCEM, final_state_CMU_SEM, projected_state_CMU_SEM,
                          min_proj_dist_SEM_1, min_proj_dist_SEM_2, number_of_sol);

    //====================================================================================
    // 4. Loop
    //====================================================================================
    int index = 0;
    COMPLETION = 0;
    #pragma omp parallel for if(isPar) shared(index)
    for(int kpos = 0; kpos <= number_of_sol; kpos++)
    {
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        double yvv[6], yv[6], st[5];

        yv[0] = init_state_CMU_NCEM[0][kpos];
        yv[1] = init_state_CMU_NCEM[1][kpos];
        yv[2] = init_state_CMU_NCEM[2][kpos];
        yv[3] = init_state_CMU_NCEM[3][kpos];
        yv[4] = init_state_CMU_NCEM[4][kpos];
        yv[5] = init_state_CMU_NCEM[5][kpos];

        st[0] = s1_CMU_EM[kpos];
        st[1] = 0.0;
        st[2] = s3_CMU_EM[kpos];
        st[3] = 0.0;
        st[4] = 0.0;

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_NCEM  = dmatrix(0, 5, 0, man_grid_size);
        double **y_man_SEM   = dmatrix(0, 5, 0, man_grid_size);
        double *t_man_EM     = dvector(0, man_grid_size);
        double *t_man_SEM    = dvector(0, man_grid_size);
        double *HManSEM    = dvector(0, man_grid_size);

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs(OFS_ORDER);
        OdeStruct driver;
        SingleOrbit orbit;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver, T, T_root, 6, qbcp_vfn, &SEML_EM);
        //Init routine
        init_orbit(orbit, &CM_EM_NC, &CM_EM_TFC,
                   &Mcoc_EM, &MIcoc_EM, &Vcoc_EM, &orbit_ofs, OFTS_ORDER, OFS_ORDER,
                   1, t0_EM[kpos], tf_EM[kpos], SEML.us.T/5, &driver, &SEML_EM);


        //====================================================================================
        // 4.1 Compute the initial orbit
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, t0_EM[kpos]);

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0_EM[kpos], t0_EM[kpos]-SEML.us.T, y_man_NCEM, t_man_EM, man_grid_size, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        NCEMmtoSEMm_vec(y_man_NCEM, t_man_EM, y_man_SEM, t_man_SEM, man_grid_size, &SEML_EM);

        //---------------------------------------------------------------------
        //Hamiltonian
        //---------------------------------------------------------------------
        HSEM_vec(t_man_SEM, y_man_SEM, HManSEM, man_grid_size, &SEML_SEM);

        //---------------------------------------------------------------------
        //Store for future save
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= man_grid_size; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) y_orb_eml_SEM[i][kpos][kman] = y_man_SEM[i][kman];
            t_orb_eml_SEM[kpos][kman] = t_man_SEM[kman];
            h_orb_eml_SEM[kpos][kman] = HManSEM[kman];
        }

        //====================================================================================
        // 4.2 Compute the manifold leg
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, yv);

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0_EM[kpos], tf_EM[kpos], y_man_NCEM, t_man_EM, man_grid_size, 0);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        NCEMmtoSEMm_vec(y_man_NCEM, t_man_EM, y_man_SEM, t_man_SEM, man_grid_size, &SEML_EM);

        //---------------------------------------------------------------------
        //Hamiltonian
        //---------------------------------------------------------------------
        HSEM_vec(t_man_SEM, y_man_SEM, HManSEM, man_grid_size, &SEML_SEM);

        //---------------------------------------------------------------------
        //Store for future save
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= man_grid_size; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) y_cmu_SEM[i][kpos][kman] = y_man_SEM[i][kman];
            t_cmu_SEM[kpos][kman] = t_man_SEM[kman];
            h_cmu_SEM[kpos][kman] = HManSEM[kman];
        }

        //----------------------------------------------------------
        // Update final data
        //----------------------------------------------------------
        // Careful: we save the final state at the end of y_cmu_SEM
        // Careful: we save the projection state at the end of y_cmu_SEM
        //----------------------------------------------------------
        for(int i = 0; i < 6; i++) y_cmu_SEM[i][kpos][man_grid_size+1] = final_state_CMU_SEM[i][kpos];
        for(int i = 0; i < 6; i++) y_cmu_SEM[i][kpos][man_grid_size+2] = projected_state_CMU_SEM[i][kpos];

        //----------------------------------------------------------
        // Careful: we save the projection time at the end of t_cmu_SEM
        //----------------------------------------------------------
        t_cmu_SEM[kpos][man_grid_size+1] = t_man_SEM[man_grid_size];
        t_cmu_SEM[kpos][man_grid_size+2] = t_man_SEM[man_grid_size];
        //----------------------------------------------------------
        // Careful: we save the Hamiltonian of the projected state at the end of h_cmu_SEM
        //----------------------------------------------------------
        h_cmu_SEM[kpos][man_grid_size+1] = HManSEM[man_grid_size];
        for(int i = 0; i < 6; i++) yvv[i] = projected_state_CMU_SEM[i][kpos];
        SEMvtoSEMm(t_man_SEM[man_grid_size], yvv, yv, &SEML_SEM);
        h_cmu_SEM[kpos][man_grid_size+2] = qbfbp_H_SEM(t_man_SEM[man_grid_size], yv, &SEML_SEM);

        #pragma omp critical
        {
            //---------------------------------------------------------------------
            // Display
            //---------------------------------------------------------------------
            displayCompletion("int_proj_CMU_EM_on_CM_SEM", 100.0*++index/(number_of_sol+1));
        }

        //---------------------------------------------------------------------
        // Free
        //---------------------------------------------------------------------
        free_dmatrix(y_man_NCEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_SEM, 0, 5, 0, man_grid_size);
        free_dvector(t_man_EM, 0, man_grid_size);
        free_dvector(t_man_SEM, 0, man_grid_size);
        free_dvector(HManSEM, 0, man_grid_size);
    }


    //====================================================================================
    // 5. Compute the final orbit about SEML2
    //====================================================================================
    //--------------------------------------
    // Center-manifold
    //--------------------------------------
    vector<Oftsc>  CM_SEM_NC;     ///center manifold in NC coordinates
    vector<Oftsc> CM_SEM_TFC;     ///center manifold in TFC coordinates
    CM_SEM_NC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_NC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CM_SEM_TFC.reserve(6);
    for(int i = 0; i <6; i++) CM_SEM_TFC.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    //Update CM
    readVOFTS_bin(CM_SEM_NC,  SEML_SEM.cs.F_PMS+"W/W");
    readVOFTS_bin(CM_SEM_TFC, SEML_SEM.cs.F_PMS+"W/Wh");


    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc_SEM(6,6);    ///COC matrix
    matrix<Ofsc>  Pcoc_SEM(6,6);    ///COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc_SEM(6,6);   ///COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc_SEM(6,6);   ///COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc_SEM(6);      ///COC vector
    //Update COC
    initCOC(Pcoc_SEM, Mcoc_SEM, PIcoc_SEM, MIcoc_SEM, Vcoc_SEM, SEML_SEM);


    //--------------------------------------
    // Loop
    //--------------------------------------
    #pragma omp parallel for if(isPar) shared(index)
    for(int kpos = 0; kpos <= number_of_sol; kpos++)
    {
        //====================================================================================
        // 4.1 Initialization
        //====================================================================================
        //---------------------------------------------------------------------
        //Initialize the initial conditions (both NC and RCM coordinates)
        //---------------------------------------------------------------------
        double st[4];
        st[0] = s1_CM_SEM[kpos];
        st[1] = 0.0;
        st[2] = s3_CM_SEM[kpos];
        st[3] = 0.0;

        //Initial time in SEM units
        double t0_SEM = tf_EM[kpos]*SEML.us_em.ns;

        //---------------------------------------------------------------------
        //Local variables to store the manifold leg
        //---------------------------------------------------------------------
        double **y_man_NCSEM  = dmatrix(0, 5, 0, man_grid_size);
        double **y_man_SEM    = dmatrix(0, 5, 0, man_grid_size);
        double *t_man_SEM     = dvector(0, man_grid_size);
        double *HManSEM       = dvector(0, man_grid_size);

        //---------------------------------------------------------------------
        // Initialisation of the orbit structure
        //---------------------------------------------------------------------
        Ofsc orbit_ofs(OFS_ORDER);
        OdeStruct driver;
        SingleOrbit orbit;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver, T, T_root, 6, qbcp_vfn, &SEML_SEM);
        //Init routine
        init_orbit(orbit, &CM_SEM_NC, &CM_SEM_TFC,
                   &Mcoc_SEM, &MIcoc_SEM, &Vcoc_SEM, &orbit_ofs, OFTS_ORDER, OFS_ORDER,
                   1, t0_SEM, t0_SEM+SEML.us_sem.T, SEML.us_sem.T/5, &driver, &SEML_SEM);


        //====================================================================================
        // 4.2 Compute the orbit
        //====================================================================================
        //---------------------------------------------------------------------
        // Update the initial state in the orbit
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, t0_SEM);

        //---------------------------------------------------------------------
        //Integration on man_grid_size+1 fixed grid
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0_SEM, t0_SEM+3*SEML.us_sem.T, y_man_NCSEM, t_man_SEM, man_grid_size, 1);

        //---------------------------------------------------------------------
        //To SEM coordinates
        //---------------------------------------------------------------------
        //NCtoSYS_vec(y_man_NCSEM, t_man_SEM, y_man_SEM, man_grid_size, &SEML_SEM); //old version
        qbcp_coc_vec(y_man_NCSEM, t_man_SEM, y_man_SEM, man_grid_size, NCSEM, PSEM);


        //---------------------------------------------------------------------
        //Hamiltonian
        //---------------------------------------------------------------------
        HSEM_vec(t_man_SEM, y_man_SEM, HManSEM, man_grid_size, &SEML_SEM);

        //---------------------------------------------------------------------
        //Store for future save
        //---------------------------------------------------------------------
        for(int kman = 0; kman <= man_grid_size; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) y_orb_seml_SEM[i][kpos][kman] = y_man_SEM[i][kman];
            t_orb_seml_SEM[kpos][kman] = t_man_SEM[kman];
            h_orb_seml_SEM[kpos][kman] = HManSEM[kman];
        }

        //---------------------------------------------------------------------
        //Free
        //---------------------------------------------------------------------
        free_dmatrix(y_man_NCSEM, 0, 5, 0, man_grid_size);
        free_dmatrix(y_man_SEM, 0, 5, 0, man_grid_size);
        free_dvector(t_man_SEM, 0, man_grid_size);
        free_dvector(HManSEM, 0, man_grid_size);
    }


    //====================================================================================
    // 6. Save
    //====================================================================================
    filename = filenameCUM(ofts_order, TYPE_MAN_SORT_IN, SEML.li_SEM);

    //---------------------
    //Open datafile
    //---------------------
    fstream filestream(filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        //---------------------
        //Loop
        //---------------------
        for(int kpos = 0; kpos <= number_of_sol; kpos++)
        {
            for(int kman = 0; kman <= man_grid_size+2; kman++)
            {
                //1. label
                res = kpos;
                filestream.write((char*) &res, sizeof(double));

                //2. s1 (EM)
                res = s1_CMU_EM[kpos];
                filestream.write((char*) &res, sizeof(double));

                //3. s3 (EM)
                res = s3_CMU_EM[kpos];
                filestream.write((char*) &res, sizeof(double));

                //4. s1 (SEM)
                res = s1_CM_SEM[kpos];
                filestream.write((char*) &res, sizeof(double));

                //5. s3 (SEM)
                res = s3_CM_SEM[kpos];
                filestream.write((char*) &res, sizeof(double));

                //6. min_proj_dist_SEM (1)
                res = min_proj_dist_SEM_1[kpos];
                filestream.write((char*) &res, sizeof(double));

                //7. The current time
                res = t_cmu_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));

                //8-13. The SE state
                for (int k = 0; k < 6; k++)
                {
                    res = y_cmu_SEM[k][kpos][kman];
                    filestream.write((char*) &res, sizeof(double));
                }

                //14. Hamiltonian
                res = h_cmu_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));

                //15. The current time on the EMLi orbit
                res = t_orb_eml_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));

                //16-21. The SE on the EMLi orbit
                for (int k = 0; k < 6; k++)
                {
                    res = y_orb_eml_SEM[k][kpos][kman];
                    filestream.write((char*) &res, sizeof(double));
                }

                //22. Hamiltonian on the EMLi orbit
                res = h_orb_eml_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));

                //23. The current time on the SEMLi orbit
                res = t_orb_seml_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));

                //24-29. The SE on the SEMLi orbit
                for (int k = 0; k < 6; k++)
                {
                    res = y_orb_seml_SEM[k][kpos][kman];
                    filestream.write((char*) &res, sizeof(double));
                }

                //30. Hamiltonian on the SEMLi orbit
                res = h_orb_seml_SEM[kpos][kman];
                filestream.write((char*) &res, sizeof(double));
            }
        }
        filestream.close();
    }

    //====================================================================================
    // 7. Free the data containers
    //====================================================================================
    free_dvector(label, 0, number_of_sol);
    free_dvector(t0_EM, 0, number_of_sol);
    free_dvector(tf_EM, 0, number_of_sol);
    free_dvector(s1_CMU_EM, 0, number_of_sol);
    free_dvector(s3_CMU_EM, 0, number_of_sol);
    free_dvector(s1_CM_SEM, 0, number_of_sol);
    free_dvector(s3_CM_SEM, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_1, 0, number_of_sol);
    free_dvector(min_proj_dist_SEM_2, 0, number_of_sol);
    free_dmatrix(final_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(projected_state_CMU_SEM, 0, 5, 0, number_of_sol);
    free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, number_of_sol);
    free_d3tensor(y_cmu_SEM, 0, 5, 0, number_of_sol, 0, man_grid_size+2);
    free_dmatrix(t_cmu_SEM, 0, number_of_sol, 0, man_grid_size+2);
    free_dmatrix(h_cmu_SEM, 0, number_of_sol, 0, man_grid_size+2);

    return 1;
}
