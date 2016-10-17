//=======================================================================================================================================
//
//          Main routines (3): old versions of connections
//
//=======================================================================================================================================
/**
 * Exponential moving average with alpha = 1/N.
 * See http://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average
 **/
double approxRollingAverage (double avg, double input, int N)
{
    avg -= avg/N;
    avg += input/N;
    return avg;
}

/**
 * Get Exponential moving average with alpha = 1/N.
 * See http://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average
 **/
double getDeltaMovingAverage(double delta, list<double> &listDeltaMA, unsigned int N)
{
    listDeltaMA.push_back(delta);
    if (listDeltaMA.size() > N) listDeltaMA.pop_front();
    double sum = 0;
    for (std::list<double>::iterator p = listDeltaMA.begin(); p != listDeltaMA.end(); ++p)
        sum += (double)*p;
    return sum / listDeltaMA.size();
}

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tmax_on_manifold_EM], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli2(double st0[],
                 double t0,
                 double tf,
                 double tmax_on_manifold_EM,
                 double tplot,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CM_TFC,
                 matrix<Ofsc>  &Mcoc,
                 matrix<Ofsc>  &Pcoc,
                 matrix<Ofsc>  &MIcoc,
                 matrix<Ofsc>  &PIcoc,
                 vector<Ofsc>  &Vcoc)
{
    //List
    std::list<double> listDeltaMA;
    //---------------------------------------------------------------------
    // Initialisation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,       //precision: root
                       PREC_DIFF,       //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,     //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;
    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Notable points at t = 0
    //------------------
    //in SEM system
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);
    //in EM system
    double **emP = dmatrix(0, 6, 0, 2);
    emPoints(0.0, emP);

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_ctrl  *h1, *h2, *h3, *h4;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    h3 = gnuplot_init();
    h4 = gnuplot_init();

    //------------------------------------------
    // Energy variations
    //------------------------------------------
    gnuplot_cmd(h4, "set title \"Variations of the energy\" ");
    gnuplot_cmd(h4, "set grid");
    gnuplot_set_xlabel(h4, (char*)"t [SEM units]");
    gnuplot_set_ylabel(h4, (char*)"H [SEM units]");

    //---------------------------------------------------------------------
    // Computation of the orbit leg
    //---------------------------------------------------------------------
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int N = floor(fabs(orbit.tf - orbit.t0)/tplot);
    double **yNCE = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------
    //Plotting
    //------------------------------------------
    double **yEM = dmatrix(0, 5, 0, N);
    NCtoSYS_vec(yNCE, tNCE, yEM, N, &SEML);
    gnuplot_plot_xyz(h1, yEM[0], yEM[1], yEM[2], N+1, (char*)"EM coordinates", "lines", "1", "1", 1);
    //Notable points
    emPlot(h1, emP);

    //------------------------------------------
    //Plotting in SEM coordinates
    //------------------------------------------
    double **ySEM = dmatrix(0, 5, 0, N);
    NCEMmtoSEMm_vec(yNCE, tNCE, ySEM, N, &SEML);
    gnuplot_plot_xyz(h2, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "1", 1);
    //Notable points
    semPlot(h2, semP);

    char ch;
    printf("Press ENTER to go on\n");
    //scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------
    //------------------------------------------
    //Loop parameters
    //------------------------------------------
    int kmin = 0;//1280;
    int kmax = N;//1300;
    int kstep = 1;
    int kn = floor((kmax-kmin)/kstep)+1;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int Nman = 200*floor(fabs(tmax_on_manifold_EM)/SEML.us.T);
    double **y_man_NCEM  = dmatrix(0, 5, 0, Nman);
    double **y_man_SEM   = dmatrix(0, 5, 0, Nman);
    double *t_man_EM     = dvector(0, Nman);
    double *t_man_SEM    = dvector(0, Nman);

    double *HMan     = dvector(0, Nman);
    double *HSEMLi   = dvector(0, Nman);
    double **y_man_NCSEM  = dmatrix(0, 5, 0, Nman);
    double DH = 0, DHmax = 1e20;

    //------------------------------------------
    // Find the minimum distance to Li
    //------------------------------------------
    int kopt = 0;
    double proj_dist_SEM, ePmOpt, ePav, min_proj_dist_SEM[kn], kvec[kn], yvmin[6];
    ePmOpt = 1e6;
    int Naverage = floor(SEML.us.T/tplot);
    cout << "Naverage = " << Naverage << endl;


    //------------------------------------------
    //Loop on the positions on the manifold
    //------------------------------------------
    double yv[6], yvc[6], tv;
    double sti[5];
    int indix  = 0;
    int indOpt = 0;
    int ntm;
    for(int k = kmin; k <= kmax; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yv[i] = yNCE[i][k];
        //The current time is the time on the orbit
        tv = tNCE[k];
        NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, +PROJ_EPSILON, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        ntm = trajectory_integration_variable_grid(orbit, tv, tv+tmax_on_manifold_EM, y_man_NCEM, t_man_EM, Nman, 0);

        if(ntm > 0)
        {
            //------------------------------------------
            //NC EM to SEM
            //------------------------------------------
            NCEMmtoSEMm_vec(y_man_NCEM, t_man_EM, y_man_SEM, t_man_SEM, ntm, &SEML);

            //------------------------------------------
            //Integrated distance to Li
            //------------------------------------------
            kvec[indix] = indix;
            min_proj_dist_SEM[indix] = 1e6;
            ePav = 0.0;
            for(int p = 0; p <= ntm; p++)
            {
                //Store in yv
                for(int i = 0; i < 6; i++) yv[i] = y_man_SEM[i][p];
                //From momenta to velocities
                SEMmtoSEMv(t_man_SEM[p], yv, yvc, &SEML);
                //distance between yvc and the Li point
                proj_dist_SEM = 0;
                for(int i = 0; i < 3; i++) proj_dist_SEM += (semP[SEML.li_SEM+4][i] - yvc[i])*(semP[SEML.li_SEM+4][i] - yvc[i]);
                //Adding velocities
                for(int i = 4; i < 6; i++) proj_dist_SEM += (yvc[i])*(yvc[i]);
                proj_dist_SEM = sqrt(proj_dist_SEM);
                //Add to the current moving average
                ePav = getDeltaMovingAverage(proj_dist_SEM, listDeltaMA, Naverage);
                if(ePav < min_proj_dist_SEM[indix]) min_proj_dist_SEM[indix] = ePav;
            }

            //------------------------------------------
            //Energy
            //------------------------------------------
            //Along the manifold
            HSEM_vec(t_man_SEM, y_man_SEM, HMan, ntm, &SEML);
            //Along SEMLi on the same time span
            HSEMLi_vec(t_man_SEM, HSEMLi, ntm, &SEML);
            //Delta
            DH = 0.0;
            for(int i = 0; i <= ntm; i++) DH += fabs(HMan[i] - HSEMLi[i]);
            DH = sqrt(DH);

            //Update the min distance
            if(DH < DHmax)
            {
                DHmax = DH;
                //If we want the energy to be the selection parameter, uncomment:
                //-------------------------------------------------------------------
                //for(int i = 0; i < 6; i++) yvmin[i] = y_man_SEM[i][ntm];     //yvmin = yv
                //kopt = k;
                //indOpt = indix;
                //Orbit
                //gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], ntm+1, (char*)"", "lines", "1", "1", 2);
                //gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 7);
                //Energy
                //gnuplot_plot_xy(h4, t_man_SEM, HMan, ntm+1, (char*)"", "lines", "1", "2", 1);
                //gnuplot_plot_xy(h4, t_man_SEM, HSEMLi, ntm+1, (char*)"", "lines", "1", "2", 2);
                //-------------------------------------------------------------------
            }

            //Update the min distance
            if(ePmOpt > min_proj_dist_SEM[indix])
            {
                ePmOpt = min_proj_dist_SEM[indix];
                for(int i = 0; i < 6; i++) yvmin[i] = y_man_SEM[i][ntm];     //yvmin = yv
                kopt = k;
                indOpt = indix;
                gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], ntm+1, (char*)"", "lines", "1", "1", 2);
                gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 7);
                //Energy
                gnuplot_plot_xy(h4, t_man_SEM, HMan, ntm+1, (char*)"", "lines", "1", "2", 1);
                gnuplot_plot_xy(h4, t_man_SEM, HSEMLi, ntm+1, (char*)"", "lines", "1", "2", 2);
            }
            //else gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], ntm+1, (char*)"", "lines", "1", "1", 7);
        }
        indix++;
    }

    //----------------------
    // Projection on the center-stable manifold
    //----------------------
    //The current state is the point on the orbit.
    for(int i = 0; i <6; i++) yv[i] = yNCE[i][kopt];
    //The current time is the time on the orbit
    tv = tNCE[kopt];
    NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, +PROJ_EPSILON, MIcoc, Vcoc);

    //----------------------
    // Update the state
    //----------------------
    orbit_update_ic(orbit, sti, tv);

    //------------------------------------------
    //Integration
    //------------------------------------------
    trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, y_man_NCEM, t_man_EM, Nman, 0);

    //------------------------------------------
    //NC EM to SEM
    //------------------------------------------
    NCEMmtoSEMm_vec(y_man_NCEM, t_man_EM, y_man_SEM, t_man_SEM, Nman, &SEML);

    //------------------------------------------
    //Energy
    //------------------------------------------
    //Along the manifold
    HSEM_vec(t_man_SEM, y_man_SEM, HMan, Nman, &SEML);
    //Along SEMLi on the same time span
    HSEMLi_vec(t_man_SEM, HSEMLi, Nman, &SEML);

    //Misc
    gnuplot_plot_xy(h4, t_man_SEM, HMan, Nman+1, (char*)"orbit", "lines", "1", "3", 4);
    gnuplot_plot_xy(h4, t_man_SEM, HSEMLi, Nman+1, (char*)"SEMLi", "lines", "1", "2", 2);


    //Plot the min distance
    gnuplot_cmd(h3, "set grid");
    gnuplot_set_xlabel(h3, (char*)"kvec [-]");
    gnuplot_set_ylabel(h3, (char*)"proj_dist_SEM [-]");
    gnuplot_plot_xy(h3, kvec, min_proj_dist_SEM, kn, (char*)"", "points", "1", "2", 1);
    gnuplot_plot_xy(h3, &kvec[indOpt], &min_proj_dist_SEM[indOpt], 1, (char*)"", "points", "2", "2", 4);
    //Plot the corresponding optimum
    gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 5);
    //Plot the corresponding manifold leg
    gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);


    cout << "Best solution has the indix " << kopt << endl;
    cout << "Which corresponds to the initial time " << tNCE[kopt] << "in EM units" << endl;
    printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);

    //-----------------------------------------------------------------------------------------------------------
    // Change to center manifold
    //-----------------------------------------------------------------------------------------------------------
    cout << "Changing the default parameters..." << endl;
    OFTS_ORDER=15;
    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc>  CMc;     ///center manifold in NC coordinates
    vector<Oftsc> CENTER_MANIFOLD_TFCc;     ///center manifold in TFC coordinates
    CMc.reserve(6);
    for(int i = 0; i <6; i++) CMc.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CENTER_MANIFOLD_TFCc.reserve(6);
    for(int i = 0; i <6; i++) CENTER_MANIFOLD_TFCc.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));

    //-----------------------------------------------------------------------------------------------------------
    //Changing the scope to SEM...
    //-----------------------------------------------------------------------------------------------------------
    cout << "Changing the scope to SEM..." << endl;
    //------------------------------------------
    //Change FOCUS
    //------------------------------------------
    changeDCS(SEML, F_SEM);

    //------------------------------------------
    // Update of the central manifold
    //------------------------------------------
    //Update CM
    readVOFTS_bin(CMc,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CENTER_MANIFOLD_TFCc, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //--------------------------------------
    //Init routine for the orbit
    //--------------------------------------
    init_orbit(orbit, &CMc, &CENTER_MANIFOLD_TFCc, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 4, 1, t0, tf, tproj, &driver, &SEML);

    //------------------------------------------
    //Get the minimum distance to SEML manifold
    //------------------------------------------
    //To NC SEM coordinates
    SYStoNC_vec(y_man_SEM, t_man_SEM, y_man_NCSEM, Nman, &SEML);
    double zproj[6], zprojmin[6], stiproj[5];
    double tvproj = 0.0;
    ePmOpt = 2;
    double ePvec[6];
    double ePvecOpt[6], ePmat[Nman+1], pmat[Nman+1];
    int pmin = 0;
    for(int p = 0; p <= Nman; p++)
    {
        //----------------------
        // Projection on the center manifold
        //----------------------
        //The current state is the point on the orbit, in NC SEM coordinates
        for(int i = 0; i <6; i++) yv[i] = y_man_NCSEM[i][p];
        //The current time is the time on the orbit, in NC SEM coordinates
        tv = t_man_SEM[p];
        //Projection on the center manifold
        NCprojCCMtoCUS(yv, tv, SEML.us_sem.n, sti, CENTER_MANIFOLD_TFCc, 0.0, MIcoc, Vcoc);
        //----------------------
        //Minimum distance?
        //----------------------
        //zproj = W(sti, tv)
        RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, 4, CENTER_MANIFOLD_TFCc, *orbit.ofs, Mcoc, Vcoc, zproj, 1);
        //distance between zproj and yv
        proj_dist_SEM = 0;
        for(int i = 0; i < 6; i++)
        {
            proj_dist_SEM += (zproj[i] - yv[i])*(zproj[i] - yv[i]);
            ePvec[i] = fabs(zproj[i] - yv[i]);
        }
        proj_dist_SEM = sqrt(proj_dist_SEM);
        ePmat[p] = proj_dist_SEM;
        pmat[p] = p;

        //Update the min distance
        if(ePmOpt > proj_dist_SEM)
        {
            pmin = p;
            ePmOpt = proj_dist_SEM;
            for(int i = 0; i < 6; i++) ePvecOpt[i] = ePvec[i];
            //Update yvmin
            for(int i = 0; i < 6; i++) yvmin[i] = yv[i];     //yvmin = yv
            NCtoSEM(tv, yvmin, yv, &SEML);                   //yv = yvmin in SEM
            for(int i = 0; i < 6; i++) yvmin[i] = yv[i];     //yvmin = yv
            //Update zprojmin
            for(int i = 0; i < 6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
            NCtoSEM(tv, zprojmin, zproj, &SEML);               //zproj = zprojmin in SEM
            for(int i = 0; i < 6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
            tvproj = tv;
            for(int i = 0; i <5; i++) stiproj[i] = sti[i];
        }

    }

    string error = "Error min = "+numTostring(ePmOpt);

    //Orbit
    gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)error.c_str(), "points", "1", "5", 5);
    gnuplot_plot_xyz(h2, zprojmin, zprojmin+1, zprojmin+2, 1, (char*)"", "points", "1", "5", 7);

    //eP
    gnuplot_plot_xy(h3, pmat, ePmat, Nman+1, (char*)"eP", "lines", "1", "2", 7);
    gnuplot_plot_xy(h3, pmat+pmin, ePmat+pmin, 1, (char*)"min_proj_dist_SEM", "points", "1", "2", 4);

    cout << "Error vector = " << endl;
    vector_printf(ePvecOpt, 6);
    cout << "Time at projection = " << tvproj << endl;
    cout << "RCM state at projection = " << endl;
    vector_printf(stiproj, 5);

    printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);

    //-----------------------------------------------------------------------------------------------------------
    // Plot the orbit starting from the projected state
    //-----------------------------------------------------------------------------------------------------------
    //Relax the condition on the projection
    orbit.ePmax = 1e-2;
    cout << "orbit.ePmax is relaxed." << endl;

    //The default interval of projection is set to Tproj = T/5
    tproj = SEML.us.T/5.0;

    //------------------------------------------
    //Update the IC
    //------------------------------------------
    orbit_update_ic(orbit, stiproj, tvproj);

    printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);


    //Integration
    int Nsem = Nman;
    double **yncsem = dmatrix(0, 5, 0, Nsem);
    double **ysem   = dmatrix(0, 5, 0, Nsem);
    double *tsem    = dvector(0, Nsem);
    double *Hsem    = dvector(0, Nsem);
    trajectory_integration_grid(orbit, tvproj, tvproj-1*SEML.us.T, yncsem, tsem, Nsem, 1);
    NCtoSYS_vec(yncsem, tsem, ysem, Nsem, &SEML);
    gnuplot_plot_xyz(h2, ysem[0], ysem[1], ysem[2], Nsem+1, (char*)"", "lines", "1", "1", 7);

    trajectory_integration_grid(orbit, tvproj, tvproj+1*SEML.us.T, yncsem, tsem, Nsem, 1);
    NCtoSYS_vec(yncsem, tsem, ysem, Nsem, &SEML);
    gnuplot_plot_xyz(h2, ysem[0], ysem[1], ysem[2], Nsem+1, (char*)"", "lines", "1", "1", 8);

    //------------------------------------------
    //Energy
    //------------------------------------------
    //Along the manifold
    HSEM_vec(tsem, ysem, Hsem, Nsem, &SEML);
    gnuplot_plot_xy(h4, tsem, Hsem, Nsem+1, (char*)"SEM orbit", "lines", "1", "2", 3);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(yEM, 0, 5, 0, N);
    free_dmatrix(ySEM, 0, 5, 0, N);
    free_dmatrix(y_man_NCEM, 0, 5, 0, Nman);
    free_dmatrix(y_man_NCSEM, 0, 5, 0, Nman);
    free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
    free_dvector(tNCE, 0, N);
    free_dvector(t_man_EM, 0, Nman);
    free_dvector(t_man_SEM, 0, Nman);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);
    gnuplot_close(h4);

    return status;
}


/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tmax_on_manifold_EM], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli_torb_vs_tman(double st0[],
                             double t0,
                             double tf,
                             double tmax_on_manifold_EM,
                             vector<Oftsc> &CM,
                             vector<Oftsc> &CM_TFC,
                             matrix<Ofsc>  &Mcoc,
                             matrix<Ofsc>  &Pcoc,
                             matrix<Ofsc>  &MIcoc,
                             matrix<Ofsc>  &PIcoc,
                             vector<Ofsc>  &Vcoc)
{
    //---------------------------------------------------------------------
    // Initialisation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,       //precision: root
                       PREC_DIFF,       //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,     //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;
    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Notable points at t = 0
    //------------------
    //in SEM system
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);
    //in EM system
    double **emP = dmatrix(0, 6, 0, 2);
    emPoints(0.0, emP);

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_ctrl  *h1, *h2, *h3;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    //------------------------------------------
    // Projection distance
    //------------------------------------------
    gnuplot_cmd(h3, "set title \"Projection distance\" ");
    gnuplot_cmd(h3, "set grid");
    gnuplot_set_xlabel(h3, (char*)"torb [ x Tsun]");
    gnuplot_set_ylabel(h3, (char*)"tmax_on_manifold_EM [ x Tsun]");
    gnuplot_set_zlabel(h3, (char*)"proj_dist_SEM  [SEM units]");

    //---------------------------------------------------------------------
    // Computation of the orbit leg
    //---------------------------------------------------------------------
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    double tplot = 1e-2; //we plot every tplot = 0.05
    int N = floor(fabs(orbit.tf - orbit.t0)/tplot);
    double **yNCE = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------
    //Plotting in SEM coordinates
    //------------------------------------------
    double **ySEM = dmatrix(0, 5, 0, N);
    NCEMmtoSEMm_vec(yNCE, tNCE, ySEM, N, &SEML);
    gnuplot_plot_xyz(h2, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "1", 1);
    //Notable points
    semPlot(h2, semP);

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------
    //------------------------------------------
    //Loop parameters
    //------------------------------------------
    int kmin  = 0;
    int kmax  = 500;//N;
    int kstep = 1;
    int kn    = floor((kmax-kmin)/kstep)+1;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int Nman = 200*floor(fabs(tmax_on_manifold_EM)/SEML.us.T);
    double **yMan = dmatrix(0, 5, 0, Nman);
    double **y_man_SEM = dmatrix(0, 5, 0, Nman);
    double *tMan  = dvector(0, Nman);

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, Nman, 0, kn);
    double **tTens  = dmatrix(0, Nman, 0, kn);

    //------------------------------------------
    //Loop on the positions on the manifold
    //------------------------------------------
    double yv[6], tv;
    double sti[5];
    int indix = 0;
    for(int k = kmin; k <= kmax; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yv[i] = yNCE[i][k];
        //The current time is the time on the orbit
        tv = tNCE[k];
        NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, +PROJ_EPSILON, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, yMan, tMan, Nman, 0);

        //------------------------------------------
        //Store data
        //------------------------------------------
        //NC EM to NC SEM
        NCEMmtoNCSEMm_vec(yMan, tMan, y_man_SEM, Nman, &SEML);
        for(int i = 0; i <= Nman; i++)
        {
            for(int k = 0; k <6; k++) yTens[k][i][indix] = y_man_SEM[k][i];
            tTens[i][indix] = tMan[i]*SEML.us_em.ns; //time is also given in SEM units!
        }
        indix++;

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        //NC EM to SEM (no NC for the output here)
        //NCEMmtoSEMm_vec(yMan, tMan, y_man_SEM, Nman, &SEML);
        //gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "1", 2);
        //printf("Press ENTER to go on\n");
        //scanf("%c",&ch);
    }

    //---------------------------------------------------------------------
    // Find the minimum distance to the center manifold of SEMLi
    //---------------------------------------------------------------------
    cout << "Changing the scope to SEM..." << endl;
    //------------------------------------------
    //Change FOCUS
    //------------------------------------------
    changeDCS(SEML, F_SEM);
    cout << "SEML.cs.F_PMS = " << SEML.cs.F_PMS << endl;
    cout << "SEML.cs.F_COC = " << SEML.cs.F_COC << endl;

    //------------------------------------------
    // Update of the central manifold
    //------------------------------------------
    //Update CM
    readVOFTS_bin(CM,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CM_TFC, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //------------------------------------------
    // Projection loop
    //------------------------------------------
    double *distMan  = dvector(0, Nman);
    double *disttorb = dvector(0, Nman);
    double *disttman = dvector(0, Nman);

    double zproj[6], yvSEM[6], zprojSEM[6], proj_dist_SEM;
    for(int k = 0; k < kn; k++)
    {
        proj_dist_SEM = 2;
        for(int p = 0; p <= Nman; p++)
        {
            //----------------------
            // Projection on the center manifold
            //----------------------
            //The current state is the point on the orbit, in NC SEM coordinates
            for(int i = 0; i <6; i++) yv[i] = yTens[i][p][k];
            //The current time is the time on the orbit, in NC SEM coordinates
            tv = tTens[p][k];
            //Projection on the center manifold
            NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, 0.0, MIcoc, Vcoc);

            //----------------------
            //Minimum distance?
            //----------------------
            //zproj = W(sti, tv)
            RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, 4, CM_TFC, *orbit.ofs, Mcoc, Vcoc, zproj, 1);

            //Back to SEM coordinates
            NCtoSYS(tv, yv, yvSEM, &SEML);
            NCtoSYS(tv, zproj, zprojSEM, &SEML);

            //distance between zproj and yv
            proj_dist_SEM = 0;
            for(int i = 0; i < 6; i++) proj_dist_SEM += (zprojSEM[i] - yvSEM[i])*(zprojSEM[i] - yvSEM[i]);
            proj_dist_SEM = sqrt(proj_dist_SEM);

            //----------------------
            //Save value
            //----------------------
            distMan[p]  = proj_dist_SEM;
            disttorb[p] = tTens[0][k]/SEML.us_sem.T;//-SEML.cs_sem.gamma*(yv[0] - SEML.cs_sem.c1);// //in SEML.us.T units
            disttman[p] = tTens[p][k]/SEML.us_sem.T;//-SEML.cs_sem.gamma*yv[1];//; //in SEML.us.T units
        }
        gnuplot_plot_xyz(h3, disttorb, disttman, distMan, Nman+1, (char*)"", "points", "7", "1", 2);
        gnuplot_cmd(h3, "set zrange [0:1e-2]");
        gnuplot_plot_xy(h1, disttorb, distMan, Nman+1, (char*)"", "points", "7", "1", 2);
        gnuplot_cmd(h1, "set yrange [0:1e-2]");
        //gnuplot_cmd(h3, "set xrange [-0.02:0.02]");
    }



    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);



    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(ySEM, 0, 5, 0, N);
    free_dmatrix(yMan, 0, 5, 0, Nman);
    free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
    free_dmatrix(tTens, 0, Nman, 0, kn);
    free_d3tensor(yTens, 0, 5, 0, Nman, 0, kn);
    free_dvector(tNCE, 0, N);
    free_dvector(tMan, 0, Nman);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);

    return status;
}


/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tmax_on_manifold_EM], with initial conditions st0, on the center-(un)stable manifold CM.
 **/
int orbit_cus(double st0[],
              double t0,
              double tf,
              double tmax_on_manifold_EM,
              vector<Oftsc> &CM,
              vector<Oftsc> &CM_TFC,
              matrix<Ofsc>  &Mcoc,
              matrix<Ofsc>  &Pcoc,
              matrix<Ofsc>  &MIcoc,
              matrix<Ofsc>  &PIcoc,
              vector<Ofsc>  &Vcoc)
{
    //---------------------------------------------------------------------
    // Initialisation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,           //precision: root
                       PREC_DIFF,           //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,            //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;

    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Notable points at t = 0
    //------------------
    //in SEM system
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);
    //in EM system
    double **emP = dmatrix(0, 6, 0, 2);
    emPoints(0.0, emP);

    //------------------
    //Plotting handlers
    //------------------
    gnuplot_ctrl  *h1, *h2, *h3;
    h1 = gnuplot_init();
    h2 = gnuplot_init();
    h3 = gnuplot_init();

    gnuplot_cmd(h3, "set title \"Variations of the energy\" ");
    gnuplot_cmd(h3, "set grid");
    gnuplot_set_xlabel(h3, (char*)"t [SEM units]");
    gnuplot_set_ylabel(h3, (char*)"H [SEM units]");

    //---------------------------------------------------------------------
    // Computation of the orbit leg
    //---------------------------------------------------------------------
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    double tplot = 1e-2; //we plot every tplot = 0.05
    int N = floor(fabs(orbit.tf - orbit.t0)/tplot);
    double **yNCE = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);
    //int status = trajectory_integration_variable_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------
    //Plotting in EM coordinates
    //------------------------------------------
    double **yEM = dmatrix(0, 5, 0, N);
    NCtoSYS_vec(yNCE, tNCE, yEM, N, &SEML);
    gnuplot_plot_xyz(h1, yEM[0], yEM[1], yEM[2], N+1, (char*)"EM coordinates", "lines", "1", "1", 1);

    //Notable points
    emPlot(h1, emP);

    //------------------------------------------
    //Plotting in SEM coordinates
    //------------------------------------------
    double **ySEM = dmatrix(0, 5, 0, N);
    NCEMmtoSEMm_vec(yNCE, tNCE, ySEM, N, &SEML);
    gnuplot_plot_xyz(h2, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "1", 1);

    //Notable points
    semPlot(h2, semP);

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------
    // Loop parameters
    int kmin  = 0;
    int kmax  = N;
    int kstep = 5;

    //Plotting parameters.
    // We need ~ 200 points per Sun period.
    int Nman = 200*floor(fabs(tmax_on_manifold_EM)/SEML.us.T);
    int color = (SEML.cs_em.manType == MAN_CENTER_S)?  2:7;
    cout << "color = " << color << endl;

    //Manifold
    double **yManNC  = dmatrix(0, 5, 0, Nman);
    double **y_man_SEM = dmatrix(0, 5, 0, Nman);
    double **yManEM  = dmatrix(0, 5, 0, Nman);
    double *tMan     = dvector(0, Nman);
    double *t_man_SEM  = dvector(0, Nman);
    double *HMan     = dvector(0, Nman);
    double *HSEMLi   = dvector(0, Nman);
    double DH = 0, DHmax = 1e20;

    //Current state
    double yv[6], tv;
    double sti[5];
    int nt = 0;
    for(int k = kmin; k <= kmax; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yv[i] = yNCE[i][k];
        //The current time is the time on the orbit
        tv = tNCE[k];
        NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, +PROJ_EPSILON, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        //trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, yManNC, tMan, Nman, 0);
        nt = trajectory_integration_variable_grid(orbit, tv, tv+tmax_on_manifold_EM, yManNC, tMan, Nman, 0);

        //------------------------------------------
        //Plotting
        //------------------------------------------
        //NCtoSYS_vec(yManNC, tMan, yManEM, nt, &SEML);
        //gnuplot_plot_xyz(h1, yManEM[0], yManEM[1], yManEM[2], nt+1, (char*)"", "lines", "1", "1", color);

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        NCEMmtoSEMm_vec(yManNC, tMan, y_man_SEM, t_man_SEM, nt, &SEML);

        //Compute the maximum distance from the Earth
        double maxDistToEarth = 0.0;
        double res = 0.0;
        for(int i = 0; i <= nt; i++)
        {
            //Earth is situated at Pe
            res = 0.0;
            for(int p = 0; p < 3; p++) res += (y_man_SEM[p][i] - semP[0][p])*(y_man_SEM[p][i] - semP[0][p]);
            res = sqrt(res);
            //See if maximum distance is topped
            if(res > maxDistToEarth)
            {
                maxDistToEarth = res;
            }
        }

        //------------------------------------------
        //Energy
        //------------------------------------------
        //Along the manifold
        HSEM_vec(t_man_SEM, y_man_SEM, HMan, nt, &SEML);
        //Along SEMLi on the same time span
        HSEMLi_vec(t_man_SEM, HSEMLi, nt, &SEML);
        //Delta
        DH = 0;
        for(int i = 0; i <= nt; i++) DH += fabs(HMan[i] - HSEMLi[i]);
        DH = sqrt(DH);

        if(maxDistToEarth < 0.02 && DH < DHmax)
        {
            DHmax = DH;
            //Orbit
            gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], nt+1, (char*)"", "lines", "1", "1", color);
            //Energy
            gnuplot_plot_xy(h3, t_man_SEM, HMan, nt+1, (char*)"", "lines", "1", "2", 1);
        }
    }

    //------------------------------------------
    //Level of energy of SEML2 on the last solution
    //------------------------------------------
    gnuplot_plot_xy(h3, t_man_SEM, HSEMLi, nt+1, (char*)"SEML2", "lines", "1", "2", 2);

    //------------------------------------------
    //Level of energy of EML2
    // TODO: - do the same in EM units (Energy)
    //       - plot the energy along the initial orbit
    //------------------------------------------
    HEMLi_in_SEM_vec(tMan, HMan, nt, &SEML);
    gnuplot_plot_xy(h3, t_man_SEM, HMan, nt+1, (char*)"EML2", "lines", "1", "2", 3);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    free_dmatrix(yManNC, 0, 5, 0, Nman);
    free_dmatrix(yManEM, 0, 5, 0, Nman);
    free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
    free_dvector(tMan, 0, Nman);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);

    return status;
}


/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tmax_on_manifold_EM], with initian conditions st0, on the center-unstable manifold CM (SEM framework).
 **/
int orbit_cu_sem(double st0[],
                 double t0,
                 double tf,
                 double tmax_on_manifold_EM,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CM_TFC,
                 matrix<Ofsc>  &Mcoc,
                 matrix<Ofsc>  &Pcoc,
                 matrix<Ofsc>  &MIcoc,
                 matrix<Ofsc>  &PIcoc,
                 vector<Ofsc>  &Vcoc)
{
    //---------------------------------------------------------------------
    // Initialisation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,           //precision: root
                       PREC_DIFF,           //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,            //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;

    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Primaries and geometrical libration points
    //------------------
    double Lsem = SEML.cs_sem.cr3bp.L;
    //Earth position in SEM coordinates
    double Pe[3], PeL[3], Pm[3];
    evaluateCoef(Pe, 0.0, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
    for(int i = 0; i <3; i++) PeL[i]= Pe[i];

    double SEML1[3], SEML2[3];
    for(int i = 0; i <3; i++) SEML1[i] = SEML.cs_sem.cr3bp.l1.position[i];
    for(int i = 0; i <3; i++) SEML2[i] = SEML.cs_sem.cr3bp.l2.position[i];
    //Get the right sign
    SEML1[0] *= -1.0;
    SEML2[0] *= -1.0;

    //---------------------------------------------------------------------
    // Computation of the orbit leg
    //---------------------------------------------------------------------
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    double tplot = 1e-2; //we plot every tplot = 0.01
    int N = floor(fabs(orbit.tf - orbit.t0)/tplot);
    double **yNCE = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------
    // SEM coordinates
    //------------------------------------------
    double **ySEM = dmatrix(0, 5, 0, N);
    NCtoSYS_vec(yNCE, tNCE, ySEM, N, &SEML);


    //------------------------------------------
    // size (Az)
    //------------------------------------------
    double Az = 0.0;
    for(int i = 0; i <= N; i+=10)
    {
        if(Az < fabs(ySEM[2][i])) Az = fabs(ySEM[2][i]);
    }
    Az *= Lsem;
    cout << "Az = " << Az << " km" << endl;


    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_ctrl  *hsem;
    hsem = gnuplot_init();
    gnuplot_cmd(hsem, "set grid");
    gnuplot_set_xlabel(hsem, (char*)"X [-]");
    gnuplot_set_ylabel(hsem, (char*)"Y [-]");
    gnuplot_set_zlabel(hsem, (char*)"Z [-]");
    gnuplot_plot_xyz(hsem, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "1", 1);
    gnuplot_plot_xyz(hsem, SEML1, SEML1+1,  SEML1+2, 1, (char*)"SEML1", "points", "1", "2", 8);
    gnuplot_plot_xyz(hsem, SEML2, SEML2+1,  SEML2+2, 1, (char*)"SEML2", "points", "2", "2", 8);
    gnuplot_plot_xyz(hsem, PeL, PeL+1,  PeL+2, 1, (char*)"Earth", "points", "3", "2", 3);

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);


    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------

    //------------------------------------------
    // Loop parameters
    //------------------------------------------
    int kmin  = 0;
    int kmax  = N;
    int kstep = 1;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    tplot = fabs(tmax_on_manifold_EM)/5000;
    int Nman = floor(fabs(tmax_on_manifold_EM)/tplot);

    double **yMan = dmatrix(0, 5, 0, Nman);
    double **y_man_SEM = dmatrix(0, 5, 0, Nman);
    double *tMan     = dvector(0, Nman);
    //Storing optimum values
    double **yManSEMOpt = dmatrix(0, 5, 0, Nman);

    //------------------------------------------
    //Loop on the position on the orbit
    //------------------------------------------
    //Keeping the position of the primaries
    double yMoon[3], yEarth[3], yMoonOpt[3];
    //Moon is situated at Pm
    evaluateCoef(yMoon, 0.0, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
    //Earth is situated at Pe
    evaluateCoef(yEarth, 0.0, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
    //Current state
    double yv[6], tv;
    double sti[5];

    //Initial min distance to Earth set arbitrarily large
    double minDistToEarthOpt = 1e6;
    double minDistToEarthMoonOpt = 1e6;
    for(int k = kmin; k <= kmax; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yv[i] = yNCE[i][k];
        //The current time is the time on the orbit
        tv = tNCE[k];
        NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, +PROJ_EPSILON, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, yMan, tMan, Nman, 0);

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        NCtoSYS_vec(yMan, tMan, y_man_SEM, Nman, &SEML);

        //------------------------------------------
        //Compute the minimum distance from the Earth
        //------------------------------------------
        double minDistToEarth = 1e6;
        double res = 0.0;
        for(int i = 0; i <= Nman; i++)
        {
            //Distance
            res = 0.0;
            for(int k = 0; k < 3; k++) res += (y_man_SEM[k][i] - Pe[k])*(y_man_SEM[k][i] - Pe[k]);
            res = sqrt(res);
            //See if maximum distance is topped
            if(res < minDistToEarth)
            {
                minDistToEarth = res;

            }
        }

        //------------------------------------------
        //Compute the minimum distance from the Moon
        //------------------------------------------
        double minDistToMoon = 1e6;
        res = 0.0;
        for(int i = 0; i <= Nman; i++)
        {

            //Moon is situated at Pm
            evaluateCoef(Pm, tMan[i], SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
            //Distance
            res = 0.0;
            for(int k = 0; k < 3; k++) res += (y_man_SEM[k][i] - Pm[k])*(y_man_SEM[k][i] - Pm[k]);
            res = sqrt(res);
            //See if maximum distance is topped
            if(res < minDistToMoon)
            {
                minDistToMoon = res;
                for(int k = 0; k < 3; k++) yMoon[k]  = Pm[k];

            }
        }

        if(minDistToEarthMoonOpt > minDistToEarth+minDistToMoon)
        {
            minDistToEarthMoonOpt = minDistToEarth+minDistToMoon;
            minDistToEarthOpt = minDistToEarth;
            for(int k = 0; k < 3; k++) yMoonOpt[k]  = yMoon[k];
            for(int i = 0; i <=Nman; i++) for(int k = 0; k < 6; k++) yManSEMOpt[k][i] = y_man_SEM[k][i];
            //gnuplot_plot_xyz(hsem, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "1", 1);
            //gnuplot_plot_xyz(hsem, yMoon, yMoon+1, yMoon+2, 1, (char*)"", "points", "6", "2", 1);
        }
    }

    //Scaling in km
    minDistToEarthOpt *= Lsem;
    string title = "Minimum distance to Earth = " + numTostring(minDistToEarthOpt) + " km";

    gnuplot_plot_xyz(hsem, yManSEMOpt[0], yManSEMOpt[1], yManSEMOpt[2], Nman+1, (char*)title.c_str(), "lines", "1", "2", 6);
    gnuplot_plot_xyz(hsem, yMoonOpt, yMoonOpt+1, yMoonOpt+2, 1, (char*)"", "points", "6", "2", 6);

    //------------------------------------------
    //Additional plots
    //------------------------------------------
    gnuplot_ctrl  *hxz, *hyz, *hxy;
    hxz = gnuplot_init();
    hyz = gnuplot_init();
    hxy = gnuplot_init();

    gnuplot_cmd(hxz, "set grid");
    gnuplot_set_xlabel(hxz, (char*)"X [-]");
    gnuplot_set_ylabel(hxz, (char*)"Z [-]");

    gnuplot_cmd(hyz, "set grid");
    gnuplot_set_xlabel(hyz, (char*)"Y [-]");
    gnuplot_set_ylabel(hyz, (char*)"Z [-]");

    gnuplot_cmd(hxy, "set grid");
    gnuplot_set_xlabel(hxy, (char*)"X [-]");
    gnuplot_set_ylabel(hxy, (char*)"Y [-]");

    gnuplot_plot_xy(hxy, yManSEMOpt[0], yManSEMOpt[1], Nman+1, (char*)title.c_str(), "lines", "1", "2", 6);
    gnuplot_plot_xy(hxy, ySEM[0], ySEM[1], N+1, (char*)"", "lines", "1", "1", 1);
    gnuplot_plot_xy(hxy, yMoonOpt, yMoonOpt+1, 1, (char*)"", "points", "6", "2", 6);
    gnuplot_plot_xy(hxy, PeL, PeL+1, 1, (char*)"Earth", "points", "3", "2", 3);

    gnuplot_plot_xy(hxz, yManSEMOpt[0], yManSEMOpt[2], Nman+1, (char*)"", "lines", "1", "2", 6);
    gnuplot_plot_xy(hxz, ySEM[0], ySEM[2], N+1, (char*)"", "lines", "1", "1", 1);
    gnuplot_plot_xy(hxz, yMoonOpt, yMoonOpt+2, 1, (char*)"", "points", "6", "2", 6);
    gnuplot_plot_xy(hxz, PeL, PeL+2, 1, (char*)"Earth", "points", "3", "2", 3);

    gnuplot_plot_xy(hyz, yManSEMOpt[1], yManSEMOpt[2], Nman+1, (char*)"", "lines", "1", "2", 6);
    gnuplot_plot_xy(hyz, ySEM[1], ySEM[2], N+1, (char*)"", "lines", "1", "1", 1);
    gnuplot_plot_xy(hyz, yMoonOpt+1, yMoonOpt+2, 1, (char*)"", "points", "6", "2", 6);
    gnuplot_plot_xy(hyz, PeL+1,  PeL+2, 1, (char*)"Earth", "points", "3", "2", 3);


    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);



    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    free_dmatrix(yMan, 0, 5, 0, Nman);
    free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
    free_dvector(tMan, 0, Nman);

    return status;
}

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tmax_on_manifold_EM], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int semliToseml2(double st0[],
                 double t0,
                 double tf,
                 double tmax_on_manifold_EM,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CM_TFC,
                 matrix<Ofsc>  &Mcoc,
                 matrix<Ofsc>  &Pcoc,
                 matrix<Ofsc>  &MIcoc,
                 matrix<Ofsc>  &PIcoc,
                 vector<Ofsc>  &Vcoc)
{
    //---------------------------------------------------------------------
    // Initialisation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,       //precision: root
                       PREC_DIFF,       //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,     //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5.0;

    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs,
               OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Primaries and geometrical libration points
    // €€TODO: put into a seperate function
    //------------------
    //Earth position in SEM coordinates
    double Pe[3];
    evaluateCoef(Pe, 0.0, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
    //Moon position in EM coordinates
    double Pm[3];
    evaluateCoef(Pm, 0.0, SEML.us_em.n, SEML.nf, SEML.cs_em.Pm, 3);
    //SEML1,2
    double SEML1[3], SEML2[3];
    for(int i = 0; i <3; i++) SEML1[i] = SEML.cs_sem.cr3bp.l1.position[i];
    for(int i = 0; i <3; i++) SEML2[i] = SEML.cs_sem.cr3bp.l2.position[i];
    //Get the right sign
    SEML1[0] *= -1.0;
    SEML2[0] *= -1.0;
    //EML2
    double EML2[3];
    for(int i = 0; i <3; i++) EML2[i] = SEML.cs_em.cr3bp.l2.position[i];
    //Get the right sign
    EML2[0] *= -1.0;


    //---------------------------------------------------------------------
    // Computation of the orbit leg
    //---------------------------------------------------------------------
    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);
    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    double tplot = 1e-2; //we plot every tplot = 0.01
    int N = floor(fabs(orbit.tf - orbit.t0)/tplot);
    //State in SEM coordinates
    double **yNCSEM = dmatrix(0, 5, 0, N);
    double **ySEM = dmatrix(0, 5, 0, N);
    double *tNCSEM  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCSEM, tNCSEM, N, 1);

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_ctrl  *hem, *hse, *h3;
    hem = gnuplot_init();
    hse = gnuplot_init();
    h3  = gnuplot_init();

    //------------------------------------------
    //Plotting in SEM coordinates
    //------------------------------------------
    NCtoSYS_vec(yNCSEM, tNCSEM, ySEM, N, &SEML);
    gnuplot_setstyle(hse, (char*)"lines");
    gnuplot_cmd(hse, "set grid");
    gnuplot_set_xlabel(hse, (char*)"X [-]");
    gnuplot_set_ylabel(hse, (char*)"Y [-]");
    gnuplot_set_zlabel(hse, (char*)"Z [-]");
    gnuplot_plot_xyz(hse, ySEM[0], ySEM[1], ySEM[2], N+1, (char*)"SEM coordinates", "lines", "1", "1", 1);
    gnuplot_plot_xyz(hse, SEML1, SEML1+1, SEML1+2, 1, (char*)"SEML1", "points", "1", "2", 4);
    gnuplot_plot_xyz(hse, SEML2, SEML2+1, SEML2+2, 1, (char*)"SEML2", "points", "2", "2", 4);
    gnuplot_plot_xyz(hse, Pe, Pe+1,  Pe+2, 1, (char*)"Earth", "points", "6", "2", 8);

    //------------------------------------------
    //Plotting in EM coordinates
    //------------------------------------------
    double **yEM = dmatrix(0, 5, 0, N);
    NCSEMmtoEMm_vec(yNCSEM, tNCSEM, yEM, N, &SEML);
    gnuplot_setstyle(hem, (char*)"lines");
    gnuplot_cmd(hem, "set grid");
    gnuplot_set_xlabel(hem, (char*)"x [-]");
    gnuplot_set_ylabel(hem, (char*)"y [-]");
    gnuplot_set_zlabel(hem, (char*)"z [-]");
    gnuplot_plot_xyz(hem, yEM[0], yEM[1], yEM[2], N+1, (char*)"EM coordinates", "lines", "1", "1", 1);
    gnuplot_plot_xyz(hem, EML2, EML2+1,  EML2+2, 1, (char*)"EML2", "points", "1", "2", 4);
    gnuplot_plot_xyz(hem, Pm, Pm+1,  Pm+2, 1, (char*)"Moon", "points", "6", "2", 4);

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

    //---------------------------------------------------------------------
    // Computation of the manifold leg
    //---------------------------------------------------------------------
    //------------------------------------------
    //Loop parameters
    //------------------------------------------
    int kmin  = 0;
    int kmax  = N;
    int kstep = 20;
    int kn = floor((kmax-kmin)/kstep)+1;

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    tplot = tmax_on_manifold_EM/5000;
    int Nman = floor(fabs(tmax_on_manifold_EM)/tplot);
    double **yMan  = dmatrix(0, 5, 0, Nman);
    double **y_man_SEM  = dmatrix(0, 5, 0, Nman);
    double **yManEM   = dmatrix(0, 5, 0, Nman);
    double *tMan      = dvector(0, Nman);

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, Nman, 0, kn);
    double **tTens  = dmatrix(0, Nman, 0, kn);

    //------------------------------------------
    //Loop on the positions on the manifold
    //------------------------------------------
    double yv[6], tv;
    double sti[5];
    int indix = 0;
    for(int k = kmin; k <= kmax/*N*/; k+=kstep)
    {
        //----------------------
        // Projection on the center-stable manifold
        //----------------------
        //The current state is the point on the orbit.
        for(int i = 0; i <6; i++) yv[i] = yNCSEM[i][k];
        //The current time is the time on the orbit
        tv = tNCSEM[k];

        //Direction for the (un)stable manifold
        double epsilon = PROJ_EPSILON;
        if(SEML.li_SEM == 1) epsilon = +PROJ_EPSILON;
        else if(SEML.li_SEM == 2) epsilon = -PROJ_EPSILON;

        //Projection
        NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, epsilon, MIcoc, Vcoc);
        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, yMan, tMan, Nman, 0);

        //------------------------------------------
        //Store data
        //------------------------------------------
        //NC SEM to NC EM
        NCSEMmtoNCEMm_vec(yMan, tMan, yManEM, Nman, &SEML);
        for(int i = 0; i <= Nman; i++)
        {
            for(int k = 0; k <6; k++) yTens[k][i][indix] = yManEM[k][i];
            tTens[i][indix] = tMan[i]/SEML.us_em.ns; //time is also given in EM units!
        }
        indix++;

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        //NC SEM to SEM (no NC for the output here)
        NCtoSYS_vec(yMan, tMan, y_man_SEM, Nman, &SEML);
        gnuplot_plot_xyz(hse, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "1", 2);

        //------------------------------------------
        //Plotting in EM coordinates
        //------------------------------------------
        //NC SEM to EM (no NC for the output here)
        //NCSEMmtoEMm_vec(yMan, tMan, yManEM, Nman, &SEML);
        //gnuplot_plot_xyz(hem, yManEM[0], yManEM[1], yManEM[2], Nman+1, (char*)"", "lines", "1", "1", 2);

    }

    //---------------------------------------------------------------------
    // Find the minimum distance to the center manifold of SEMLi
    //---------------------------------------------------------------------
    cout << "Changing the scope to EM..." << endl;
    //------------------------------------------
    //Change FOCUS
    //------------------------------------------
    changeDCS(SEML, F_EM);

    //------------------------------------------
    // Update of the central manifold
    //------------------------------------------
    //Update CM
    readVOFTS_bin(CM,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CM_TFC, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //------------------------------------------
    // Projection loop
    //------------------------------------------
    int kopt = 0;
    double zproj[6], proj_dist_SEM, min_proj_dist_SEM[kn], kvec[kn];
    double yvmin[6], zprojmin[6];
    double yvopt[6], zprojopt[6];
    double tvproj = 0.0;
    double tvopt = 0.0;
    double ePopt = 0.0;
    for(int k = 0; k < kn; k++)
    {
        kvec[k] = k;
        min_proj_dist_SEM[k] = 0.0;
        for(int p = 0; p <= Nman; p+=50)
        {
            //----------------------
            // Projection on the center manifold
            //----------------------
            //The current state is the point on the orbit, in NC EM coordinates
            for(int i = 0; i <6; i++) yv[i] = yTens[i][p][k];
            //The current time is the time on the orbit, in NC EM coordinates
            tv = tTens[p][k];
            //Projection on the center manifold
            NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, 0.0, MIcoc, Vcoc);
            //----------------------
            //Minimum distance?
            //----------------------
            //zproj = W(sti, tv)
            RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, 4, CM_TFC, *orbit.ofs, Mcoc, Vcoc, zproj, 1);
            //distance between zproj and yv
            proj_dist_SEM = 0;
            for(int i = 0; i < 6; i++) proj_dist_SEM += (zproj[i] - yv[i])*(zproj[i] - yv[i]);
            proj_dist_SEM = sqrt(proj_dist_SEM);

            //cout << "proj_dist_SEM = " << proj_dist_SEM << endl;
            //Update the min distance
            if(p == 0)
            {
                min_proj_dist_SEM[k] = proj_dist_SEM;
                //Update yvmin
                for(int i = 0; i <6; i++) yvmin[i] = yv[i];     //yvmin = yv
                NCtoEM(tv, yvmin, yv, &SEML);                   //yv = yvmin in EM
                for(int i = 0; i <6; i++) yvmin[i] = yv[i];     //yvmin = yv
                //Update zprojmin
                for(int i = 0; i <6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
                NCtoEM(tv, zprojmin, zproj, &SEML);               //zproj = zprojmin in EM
                for(int i = 0; i <6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
                tvproj = tv;
            }
            else
            {
                if(min_proj_dist_SEM[k] > proj_dist_SEM)
                {
                    min_proj_dist_SEM[k] = proj_dist_SEM;
                    //Update yvmin
                    for(int i = 0; i <6; i++) yvmin[i] = yv[i];     //yvmin = yv
                    NCtoEM(tv, yvmin, yv, &SEML);                   //yv = yvmin in EM
                    for(int i = 0; i <6; i++) yvmin[i] = yv[i];     //yvmin = yv
                    //Update zprojmin
                    for(int i = 0; i <6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
                    NCtoEM(tv, zprojmin, zproj, &SEML);           //zproj = zprojmin in EM
                    for(int i = 0; i <6; i++) zprojmin[i] = zproj[i]; //zprojmin = zproj
                    tvproj = tv;
                }
            }
        }


        //----------------------------
        //Super optimum
        //----------------------------
        if(k==0)
        {
            kopt = k;
            tvopt = tvproj;
            for(int i = 0; i <6; i++)
            {
                yvopt[i] = yvmin[i];
                zprojopt[i] = zprojmin[i];
            }
            ePopt = min_proj_dist_SEM[k];
        }
        else
        {
            if(min_proj_dist_SEM[k] < ePopt)
            {
                kopt = k;
                tvopt = tvproj;
                for(int i = 0; i <6; i++)
                {
                    yvopt[i] = yvmin[i];
                    zprojopt[i] = zprojmin[i];
                }
                ePopt = min_proj_dist_SEM[k];
            }
        }
    }

    //Plot the corresponding points
    gnuplot_plot_xyz(hem, yvopt, yvopt+1, yvopt+2, 1, (char*)"yv", "points", "1", "3", 5);
    gnuplot_plot_xyz(hem, zprojopt, zprojopt+1, zprojopt+2, 1, (char*)"zproj", "points", "1", "3", 6);

    //Plot the min distance
    gnuplot_setstyle(h3, (char*)"points");
    gnuplot_cmd(h3, "set grid");
    gnuplot_set_xlabel(h3, (char*)"kvec [-]");
    gnuplot_set_ylabel(h3, (char*)"proj_dist_SEM [-]");
    gnuplot_plot_xy(h3, kvec, min_proj_dist_SEM, kn, (char*)"", "points", "1", "1", 1);
    gnuplot_plot_xy(h3, &kvec[kopt], &min_proj_dist_SEM[kopt], 1, (char*)"", "points", "3", "1", 4);

    //Plot the corresponding manifold leg in EM coordinates
    for(int i = 0; i <= Nman; i++)
    {
        for(int k = 0; k < 6; k++) yMan[k][i] = yTens[k][i][kopt];
        tMan[i] = tTens[i][kopt];
    }
    NCtoSYS_vec(yMan, tMan, yManEM, Nman, &SEML);
    gnuplot_plot_xyz(hem, yManEM[0], yManEM[1], yManEM[2], Nman+1, (char*)"", "lines", "1", "1", 4);

    //Plot the manifold leg in SEM coordinates
    NCEMmtoSEMm_vec(yMan, tMan, y_man_SEM, Nman, &SEML);
    gnuplot_plot_xyz(hse, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "1", 4);
    //Moon position in SEM coordinates
    evaluateCoef(Pm, tvopt*SEML.us_em.ns, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
    gnuplot_plot_xyz(hse, Pm, Pm+1,  Pm+2, 1, (char*)"Moon", "points", "6", "2", 4);

    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    free_dmatrix(yNCSEM, 0, 5, 0, N);
    free_dmatrix(yEM, 0, 5, 0, N);
    free_dmatrix(ySEM, 0, 5, 0, N);
    free_dmatrix(yMan, 0, 5, 0, Nman);
    free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
    free_dmatrix(yManEM, 0, 5, 0, Nman);
    free_dmatrix(tTens, 0, Nman, 0, kn);
    free_d3tensor(yTens, 0, 5, 0, Nman, 0, kn);
    free_dvector(tNCSEM, 0, N);
    free_dvector(tMan, 0, Nman);
    gnuplot_close(hse);
    gnuplot_close(hem);
    gnuplot_close(h3);

    return status;
}


//=======================================================================================================================================
//
//          I/O (2). Used in manOrbit/manOrbitProj/manOrbitPostProcess/manOrbitSEM
//
//=======================================================================================================================================
/**
 * \brief Store the manifold (tTens, yTens) in the  data file filenameOrbit(ofts_order, sizeOrbit, type). The manifodl is composed of N+1 branches, each of them
 * described by Nman+1 points.
 **/
int writeManifold_bin(double **tTens, double ***yTens, int ofts_order, int sizeOrbit, int N, int Nman, int type)
{
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);
    fstream filestream;


    filestream.open (filename.c_str(), ios::binary | ios::out);
    if (filestream.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of data
        //---------------------
        resi = N;
        filestream.write((char*) &resi, sizeof(int));
        resi = Nman;
        filestream.write((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int n1 = 0; n1 < N; n1++)
        {
            //Current homogeneous polynomial
            for (int n2 = 0; n2 < Nman; n2++)
            {
                res = tTens[n1][n2];
                filestream.write((char*) &res, sizeof(double));
                for (int k = 0; k < 6; k++)
                {
                    res = yTens[k][n1][n2];
                    filestream.write((char*) &res, sizeof(double));
                }
            }
        }
        filestream.close();
    }
    else return 0;

    return 1;
}

/**
* \brief Read the manifold (tTens, yTens) in the  data file filenameOrbit(ofts_order, sizeOrbit, type). The manifodl is composed of N+1 branches, each of them
* described by Nman+1 points.
**/
int readManifold_bin(double **tTens, double ***yTens, int ofts_order, int sizeOrbit, int N, int Nman, int type)
{
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);
    fstream filestream;


    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        int N0, Nman0;
        //---------------------
        //Number of data
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        N0 = resi;
        filestream.read((char*) &resi, sizeof(int));
        Nman0 = resi;

        if(N0 < N || Nman0 < Nman)
        {
            cout << "readManifold_bin: wrong inputs" << endl;
            cout << "N = " << N << ", but N0 = " << N0 << endl;
            cout << "Nman = " << Nman << ", but Nman0 = " << Nman0 << endl;
            return 0;
        }

        double res;
        //Loop
        for(int n1 = 0; n1 < N; n1++)
        {
            //Current homogeneous polynomial
            for (int n2 = 0; n2 < Nman; n2++)
            {
                filestream.read((char*) &res, sizeof(double));
                tTens[n1][n2] = res;

                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                    yTens[k][n1][n2] = res;
                }
            }
        }
        filestream.close();
    }
    else return 0;

    return 1;
}

/**
 * \brief Get the length of the data file filenameOrbit(ofts_order, sizeOrbit, type) containing the manifold (tTens, yTens). The manifodl is composed of N+1 branches, each of them
 * described by Nman+1 points.
 **/
int getLenghtManifold_bin(int ofts_order, int sizeOrbit, int *N, int *Nman, int type)
{
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);
    fstream filestream;


    filestream.open (filename.c_str(), ios::binary | ios::in);
    if (filestream.is_open())
    {
        int resi;
        //---------------------
        //Number of data
        //---------------------
        filestream.read((char*) &resi, sizeof(int));
        *N = resi;
        filestream.read((char*) &resi, sizeof(int));
        *Nman = resi;
        filestream.close();
    }
    else return 0;

    return 1;
}


//=======================================================================================================================================
//
//          Main routines (2): connections with a fixed orbit
//
//=======================================================================================================================================
/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initial conditions st0, on the center manifold CM.
 **/
int gridOrbit(double st0[],
              double t0,
              double tf,
              double dt,
              vector<Oftsc> &CM,
              vector<Oftsc> &CM_TFC,
              matrix<Ofsc>  &Mcoc,
              matrix<Ofsc>  &Pcoc,
              matrix<Ofsc>  &MIcoc,
              matrix<Ofsc>  &PIcoc,
              vector<Ofsc>  &Vcoc)
{
    //---------------------------------------------------------------------
    // Computation
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver,
                       T,               //stepper: int
                       T_root,          //stepper: root
                       PREC_ABS,        //precision: int_abs
                       PREC_REL,        //precision: int_rel
                       PREC_ROOT,       //precision: root
                       PREC_DIFF,       //precision: diffcorr
                       6,               //dimension
                       PREC_HSTART,     //initial int step
                       qbfbp_vfn_novar, //vector field
                       NULL,            //jacobian
                       &SEML);          //current four-body system

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us.T/5;

    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 4, 1, t0, tf, tproj, &driver, &SEML);

    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int N         = floor(fabs(orbit.tf - orbit.t0)/dt);
    double **yNCE = dmatrix(0, 5, 0, N);
    double **ySYS = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------
    //Integration
    //------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------
    //To SYS coordinates
    //------------------------------------------
    NCtoSYS_vec(yNCE, tNCE, ySYS,N, &SEML);

    //------------------------------------------
    //Plotting
    //------------------------------------------
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, ySYS[0], ySYS[1], ySYS[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);


    //------------------------------------------
    //Save in file
    //------------------------------------------
    writeOrbit(tNCE, yNCE, OFTS_ORDER, floor(fabs(st0[0])), N+1, TYPE_ORBIT);

    //------------------------------------------
    //User check
    //------------------------------------------
    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(h1);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_orbit(&orbit);
    //Free
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(ySYS, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    return status;
}


/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file (obtained from gridOrbit).
 **/
int cusOrbit(int ofts_order,
             int sizeOrbit,
             double epsilon,
             vector<Oftsc> &CM_TFC,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &MIcoc,
             vector<Ofsc>  &Vcoc,
             bool isPar)
{
    //----------------------------------------------------------
    // Get the number of data lines from file
    //----------------------------------------------------------
    int N = getLineNumber(ofts_order, sizeOrbit, TYPE_ORBIT)-1;

    //----------------------------------------------------------
    //If the file exists and contains data, go on
    //----------------------------------------------------------
    if(N > 0)
    {
        //----------------------
        // Read data from file
        //----------------------
        double **yNCE = dmatrix(0, 5, 0, N);
        double *tNCE  = dvector(0, N);
        int status    = readOrbit(tNCE, yNCE, ofts_order, sizeOrbit, N+1, TYPE_ORBIT);

        //----------------------
        // Loop on all elements
        //----------------------
        if(status)
        {
            double **init_state_CMU_NCEM = dmatrix(0, 5, 0, N);
            double **init_state_CMU_RCM = dmatrix(0, 4, 0, N);

            #pragma omp parallel for if(isPar)
            for(int k = 0; k <= N; k++)
            {
                Ofsc ofs(OFS_ORDER);
                double *yv  = dvector(0,5);
                double *yvu = dvector(0,5);
                double *sti = dvector(0,4);
                double tv;
                //----------------------
                // Projection on the center-stable manifold
                //----------------------
                //The current state is the point on the orbit.
                for(int i = 0; i <6; i++) yv[i] = yNCE[i][k];
                //The current time is the time on the orbit
                tv = tNCE[k];
                //Projection along the center-stable direction
                NCprojCCMtoCUS(yv, tv, SEML.us.n, sti, CM_TFC, epsilon, MIcoc, Vcoc);
                //Equivalent state
                RCMtoNCbyTFC(sti, tv, SEML.us.n, OFTS_ORDER, OFS_ORDER, 5, CM_TFC, ofs, Mcoc, Vcoc, yvu, 1);

                #pragma omp critical
                {
                    //Save
                    init_state_CMU_NCEM[0][k] = yvu[0];
                    init_state_CMU_NCEM[1][k] = yvu[1];
                    init_state_CMU_NCEM[2][k] = yvu[2];
                    init_state_CMU_NCEM[3][k] = yvu[3];
                    init_state_CMU_NCEM[4][k] = yvu[4];
                    init_state_CMU_NCEM[5][k] = yvu[5];

                    //Save
                    init_state_CMU_RCM[0][k] = sti[0];
                    init_state_CMU_RCM[1][k] = sti[1];
                    init_state_CMU_RCM[2][k] = sti[2];
                    init_state_CMU_RCM[3][k] = sti[3];
                    init_state_CMU_RCM[4][k] = sti[4];
                }

                free_dvector(yv, 0, 5);
                free_dvector(yvu, 0, 5);
                free_dvector(sti, 0, 4);
            }
            //Store values
            writeOrbit(tNCE, init_state_CMU_NCEM, init_state_CMU_RCM, ofts_order, sizeOrbit, N+1, TYPE_CU);
            //Free
            free_dmatrix(init_state_CMU_NCEM, 0, 5, 0, N);
            free_dmatrix(init_state_CMU_RCM, 0, 4, 0, N);
        }
        else
        {
            cout << "cusOrbit: error during data reading." << endl;
            return 0;
        }

        //Free
        free_dmatrix(yNCE, 0, 5, 0, N);
        free_dvector(tNCE, 0, N);
    }
    else
    {
        cout << "cusOrbit: no datafile." << endl;
        return 0;

    }

    return 1;
}

/**
 *  \brief Computes the manifold branches from a discrete set of unstable directions along a given LPO (obtained from cusOrbit)
 **/
int manOrbit(double tmax_on_manifold_EM,
             int ofts_order,
             int sizeOrbit,
             int number_of_sol,
             int Nman,
             int isPar,
             vector<Oftsc> &CM,
             vector<Oftsc> &CM_TFC,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc)
{
    //----------------------------------------------------------
    // Get the number of data lines from file
    //----------------------------------------------------------
    int N;
    int N0 = getLineNumber(ofts_order, sizeOrbit, TYPE_CU)-1;

    //----------------------------------------------------------
    //If the file exists and contains data, go on
    //----------------------------------------------------------
    if(N0 > 0)
    {
        if(number_of_sol <= 0) N = N0;
        else N = number_of_sol;

        if(N > N0)
        {
            cout << "manOrbit: wrong input, number_of_sol = " << number_of_sol << ", but N0 = " << N0 << endl;
            return 0;
        }

        //----------------------
        // Read data from file
        //----------------------
        double **yNCU = dmatrix(0, 5, 0, N);
        double **sNCU = dmatrix(0, 4, 0, N);
        double *tNCU  = dvector(0, N);
        int status    = readOrbit(tNCU, yNCU, sNCU, ofts_order, sizeOrbit, N+1, TYPE_CU);

        //----------------------
        // Loop on all elements
        //----------------------
        if(status)
        {
            //To store all data
            double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
            double **tTens  = dmatrix(0, N, 0, Nman);


            //The default interval of projection is set to Tproj = T/5
            double tproj = SEML.us.T/5.0;

            int index = 0;
            #pragma omp parallel for if(isPar)
            for(int k = 0; k <= N; k++)
            {
                //---------------------------------------------------------------------
                //Temp variables
                //---------------------------------------------------------------------
                double yv[6], st[5], tv;
                yv[0] = yNCU[0][k];
                yv[1] = yNCU[1][k];
                yv[2] = yNCU[2][k];
                yv[3] = yNCU[3][k];
                yv[4] = yNCU[4][k];
                yv[5] = yNCU[5][k];
                st[0] = sNCU[0][k];
                st[1] = sNCU[1][k];
                st[2] = sNCU[2][k];
                st[3] = sNCU[3][k];
                st[4] = sNCU[4][k];
                tv    = tNCU[k];

                //---------------------------------------------------------------------
                //Local variables
                //---------------------------------------------------------------------
                double **y_man_NCEM  = dmatrix(0, 5, 0, Nman);
                double **y_man_NCSEM = dmatrix(0, 5, 0, Nman);
                double *t_man_EM     = dvector(0, Nman);

                //---------------------------------------------------------------------
                // Initialisation
                //---------------------------------------------------------------------
                //------------------
                // ODE
                //------------------
                OdeStruct driver;
                //Root-finding
                const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
                //Stepper
                const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
                //Init ode structure
                init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF, 6, PREC_HSTART, qbfbp_vfn_novar, NULL, &SEML);

                //------------------
                //Orbit structure
                //------------------
                Ofsc orbit_ofs(OFS_ORDER);
                SingleOrbit orbit;
                //Init routine
                init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, tv, tv+tmax_on_manifold_EM, tproj, &driver, &SEML);

                //---------------------------------------------------------------------
                // Update the state
                //---------------------------------------------------------------------
                orbit_update_ic(orbit, st, yv, tv);

                //---------------------------------------------------------------------
                //Integration
                //---------------------------------------------------------------------
                trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, y_man_NCEM, t_man_EM, Nman, 0);

                //                //---------------------------------------------------------------------
                //                //Integration check
                //                //---------------------------------------------------------------------
                //                reset_ode_structure(orbit.driver);
                //                double t0 = t_man_EM[0];
                //                double yv0[6];
                //                for(int i = 0; i < 6; i++) yv0[i] = yv[i];
                //                int kk = Nman;
                //
                //                cout << "Time 0" << endl;
                //                cout << t_man_EM[0] << "  " << t0 << "   " << endl;
                //
                //                cout << "Init" << endl;
                //                for(int i = 0; i < 6; i++) cout << y_man_NCEM[i][0] << endl;
                //
                //                //Integration
                //                gsl_odeiv2_driver_apply (driver.d, &t0, t_man_EM[kk], yv);
                //                cout << "Time" << endl;
                //                cout << t_man_EM[kk] << "  " << t0 << "   " << endl;
                //                cout << "Result" << endl;
                //                for(int i = 0; i <6; i++) cout << (y_man_NCEM[i][kk] - yv[i])/max(1.0, fabs(yv[i])) << endl;
                //
                //
                //                t0 = t_man_EM[0];
                //                reset_ode_structure(&driver);
                //                gsl_odeiv2_driver_apply (driver.d, &t0, t_man_EM[kk], yv0);
                //                cout << "Result 2: " << endl;
                //                for(int i = 0; i < 6; i++) cout << (yv[i] - yv0[i])/max(1.0, fabs(yv[i])) << endl;


                //---------------------------------------------------------------------
                //Storage
                //---------------------------------------------------------------------
                #pragma omp critical
                {
                    //To NCSEM coordinates
                    NCEMmtoNCSEMm_vec(y_man_NCEM, t_man_EM, y_man_NCSEM, Nman, &SEML);

                    //                    //----------------------------------------------------------
                    //                    //Gnuplot window
                    //                    //----------------------------------------------------------
                    //                    gnuplot_ctrl  *h1;
                    //                    h1 = gnuplot_init();
                    //                    gnuplot_plot_xyz(h1, y_man_NCSEM[0], y_man_NCSEM[1], y_man_NCSEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
                    //                    char ch;
                    //                    printf("Press ENTER to go on\n");
                    //                    scanf("%c",&ch);
                    //                    gnuplot_close(h1);

                    for(int p = 0; p <= Nman; p++)
                    {
                        for(int i = 0; i < 6; i++) yTens[i][k][p] = y_man_NCSEM[i][p];
                        tTens[k][p] = t_man_EM[p]*SEML.us_em.ns; //time is also given in SEM units!
                    }
                }

                //---------------------------------------------------------------------
                //Free
                //---------------------------------------------------------------------
                free_dmatrix(y_man_NCEM, 0, 5, 0, Nman);
                free_dmatrix(y_man_NCSEM, 0, 5, 0, Nman);
                free_dvector(t_man_EM, 0, Nman);

                //----------------------------------------------------------
                //Display completion
                //----------------------------------------------------------
                displayCompletion("manOrbit", 100.0*index++/N);
            }


            //---------------------------------------------------------------------
            //Save
            //---------------------------------------------------------------------
            writeManifold_bin(tTens, yTens, ofts_order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

            //---------------------------------------------------------------------
            //Free
            //---------------------------------------------------------------------
            free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
            free_dmatrix(tTens, 0, N, 0, Nman);
        }
        else
        {
            cout << "cusOrbit: error during data reading." << endl;
            return 0;
        }
        //Free
        free_dmatrix(yNCU, 0, 5, 0, N);
        free_dvector(tNCU, 0, N);
    }
    return 1;
}


/**
 *  \brief Postprocess of the manifold branches obtained with manOrbit. Sort with them with respect to the integrated distance between the current state along
 * the trajectories and the SEMLi point selected in the SEML structure.
 **/
void manOrbitPostProcess(int ofts_order,
                         int sizeOrbit,
                         int Nex,
                         int isPar)
{
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int N, Nman;
    getLenghtManifold_bin(ofts_order, sizeOrbit, &N, &Nman, TYPE_MAN);
    N--;   //need to shift to take into account that vectors starts at 0
    Nman--; //need to shift to take into account that vectors starts at 0

    int number_of_sol = Nex < 0 ? N:Nex;

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //----------------------------------------------------------
    //Read data
    //----------------------------------------------------------
    int status = readManifold_bin(tTens, yTens, ofts_order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

    //----------------------------------------------------------
    //Selection
    //----------------------------------------------------------
    if(status)
    {
        //----------------------------------------------------------
        //Notable points in SEM system
        //----------------------------------------------------------
        double **semP = dmatrix(0, 6, 0, 2);
        semPoints(0.0, semP);

        //----------------------------------------------------------
        // To store
        //----------------------------------------------------------
        vector<double> DH(N+1);
        vector<double> DE(N+1);

        //----------------------------------------------------------
        // Loop on all positions
        //----------------------------------------------------------
        #pragma omp parallel for if(isPar)
        for(int k = 0; k <= N; k++)
        {
            //----------------------------------------------------------
            //Temp variable
            //----------------------------------------------------------
            double **y_man_NCSEM = dmatrix(0, 5, 0, Nman);
            double **y_man_SEM   = dmatrix(0, 5, 0, Nman);
            double *t_man_SEM    = dvector(0, Nman);
            double *HMan       = dvector(0, Nman);
            double *HSEMLi     = dvector(0, Nman);


            //-------------------------------
            // Update y_man_NCSEM, tManNCSEM
            //-------------------------------
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) y_man_NCSEM[i][p] = yTens[i][k][p];
                t_man_SEM[p] = tTens[k][p];
            }

            //-------------------------------
            //To SEM coordinates
            //-------------------------------
            NCtoSEM_vec(y_man_NCSEM, t_man_SEM, y_man_SEM, Nman, &SEML);

            //-------------------------------
            //Compute the energy
            //on the last part of the trajectory
            //-------------------------------
            //Along the manifold
            HSEM_vec(t_man_SEM, y_man_SEM, HMan, Nman, &SEML);
            //Along SEMLi on the same time span
            HSEMLi_vec(t_man_SEM, HSEMLi, Nman, &SEML);
            //Delta
            DH[k] = 0.0;
            for(int p = 0; p <= Nman; p++)
            {
                if(fabs(t_man_SEM[Nman] - t_man_SEM[p]) < SEML.us_sem.T) DH[k] += fabs(HMan[p] - HSEMLi[p]);
            }

            //-------------------------------
            //Compute the distance the SEMLi
            //on the last part of the trajectory
            //-------------------------------
            double proj_dist_SEM;
            double yv[6], yvc[6];
            for(int p = 0; p <= Nman; p++)
            {
                if(true /*fabs(t_man_SEM[Nman] - t_man_SEM[p]) < SEML.us_sem.T*/)
                {
                    //Store in yv
                    for(int i = 0; i < 6; i++) yv[i] = y_man_SEM[i][p];
                    //From momenta to velocities
                    SEMmtoSEMv(t_man_SEM[p], yv, yvc, &SEML);
                    //distance between yvc and the Li point
                    proj_dist_SEM = 0;
                    for(int i = 0; i < 3; i++) proj_dist_SEM += (semP[SEML.li_SEM+4][i] - yvc[i])*(semP[SEML.li_SEM+4][i] - yvc[i]);
                    //Adding velocities
                    //for(int i = 4; i < 6; i++) proj_dist_SEM += (yvc[i])*(yvc[i]);
                    //Store
                    DE[k] += sqrt(proj_dist_SEM);
                }
            }

            //----------------------------------------------------------
            //Free
            //----------------------------------------------------------
            free_dmatrix(y_man_NCSEM, 0, 5, 0, Nman);
            free_dmatrix(y_man_SEM, 0, 5, 0, Nman);
            free_dvector(t_man_SEM, 0, Nman);
            free_dvector(HMan, 0, Nman);
            free_dvector(HSEMLi, 0, Nman);
        }

        //----------------------------------------------------------
        //Sorting DH by indexes
        //----------------------------------------------------------
        vector<size_t> sortId = sort_indexes(DE);

        //----------------------------------------------------------
        //Select the first number_of_sol values
        //----------------------------------------------------------
        double **y_man_NCSEM = dmatrix(0, 5, 0, Nman);
        double *t_man_SEM    = dvector(0, Nman);
        double ***y_cmu_SEM = d3tensor(0, 5, 0, number_of_sol, 0, Nman);
        double **t_cmu_SEM  = dmatrix(0, number_of_sol, 0, Nman);
        int sortk;
        for(int k = 0; k <= number_of_sol; k++)
        {
            sortk = sortId[k];
            //-------------------------------
            // Update y_man_NCSEM, tManNCSEM
            //-------------------------------
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) y_man_NCSEM[i][p] = yTens[i][sortk][p];
                t_man_SEM[p] = tTens[sortk][p];

                for(int i = 0; i < 6; i++) y_cmu_SEM[i][k][p] = yTens[i][sortk][p];
                t_cmu_SEM[k][p] = tTens[sortk][p];
            }
        }

        //---------------------------------------------------------------------
        //Save
        //---------------------------------------------------------------------
        writeManifold_bin(t_cmu_SEM, y_cmu_SEM, ofts_order, sizeOrbit, number_of_sol+1, Nman+1, TYPE_MAN_SORT_DR);

        //---------------------------------------------------------------------
        //Free
        //---------------------------------------------------------------------
        free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
        free_d3tensor(y_cmu_SEM, 0, 5, 0, number_of_sol, 0, Nman);
        free_dmatrix(tTens, 0, N, 0, Nman);
        free_dmatrix(t_cmu_SEM, 0, number_of_sol, 0, Nman);
        free_dmatrix(y_man_NCSEM, 0, 5, 0, Nman);
        free_dvector(t_man_SEM, 0, Nman);
    }
    else
    {
        cout << "manOrbitPostProcess: error while loading the file" << endl;
    }
}



/**
 *  \brief Projection of the manifold branches obtained with manOrbit+manOrbitPostProcess. Each state on the manifold branches are projected on the Center Manifold
 *   of the SEMLi point selected in the SEML structure. Then, the corresponding solutions are stored in a data file.
 **/
void manOrbitProj(int order_em,
                  int size_em,
                  int type_em,
                  int order_sem,
                  int Nex,
                  int Nmant,
                  matrix<Ofsc>  &Mcoc,
                  matrix<Ofsc>  &Pcoc,
                  matrix<Ofsc>  &MIcoc,
                  matrix<Ofsc>  &PIcoc,
                  vector<Ofsc>  &Vcoc,
                  int isPar)
{
    //-------------------------------------------------------------------------------
    //Getting back the data from EM files
    //-------------------------------------------------------------------------------
    //Read data sizeOrbit
    int N, Nman;
    getLenghtManifold_bin(order_em, size_em, &N, &Nman, type_em);
    N--;     //need to shift to take into account that vectors starts at 0
    Nman--;  //need to shift to take into account that vectors starts at 0

    int number_of_sol = Nex < 0 ? N:Nex;

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //To store final data
    double ***y_cmu_SEM = d3tensor(0, 5, 0, number_of_sol, 0, Nman+1);
    double **t_cmu_SEM  = dmatrix(0, number_of_sol, 0, Nman+1);

    //Read data
    readManifold_bin(tTens, yTens, order_em, size_em, N+1, Nman+1, type_em);

    //-------------------------------------------------------------------------------
    //Changing the scope to CM SEM...
    //-------------------------------------------------------------------------------
    changeScope(order_sem, F_SEM);

    //--------------------------------------
    // Center- manifold
    //--------------------------------------
    vector<Oftsc>  CMc;     ///center manifold in NC coordinates
    vector<Oftsc> CENTER_MANIFOLD_TFCc;     ///center manifold in TFC coordinates
    CMc.reserve(6);
    for(int i = 0; i <6; i++) CMc.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));
    CENTER_MANIFOLD_TFCc.reserve(6);
    for(int i = 0; i <6; i++) CENTER_MANIFOLD_TFCc.push_back(Oftsc(4, OFTS_ORDER, OFS_NV, OFS_ORDER));

    //------------------------------------------
    // Update of the central manifold
    //------------------------------------------
    //Update CM
    readVOFTS_bin(CMc,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CENTER_MANIFOLD_TFCc, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    vector<int>    indexMin(N+1);
    vector<double> distMin(N+1);
    double **final_state_CMU_SEM = dmatrix(0, 5, 0, N);
    double **projected_state_CMU_SEM = dmatrix(0, 5, 0, N);

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl  *h1, *h2;
    h1 = gnuplot_init();
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);
    semPlot(h2, semP);

    //-------------------------------------------------------------------------------
    // Superloop
    //-------------------------------------------------------------------------------
    int index = 0;
    #pragma omp parallel for if(isPar) shared(index)
    for(int kpos = 0; kpos <= N; kpos++)
    {
        //----------------------------------------------------------
        //Temp variable
        //----------------------------------------------------------
        Ofsc ofs(OFS_ORDER);
        double yv[6], yvproj[6], sproj[4], tv;
        double yv_SEM[6], yvproj_SEM[6];
        double proj_dist_SEM, min_proj_dist_SEM = 1e6;
        double kmin = 0;
        double y_man_norm_NCSEM = 0.0;
        //----------------------------------------------------------
        //Loop on trajectory
        //----------------------------------------------------------
        for(int kman = 0; kman <= Nman; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) yv[i] = yTens[i][kpos][kman];
            tv = tTens[kpos][kman];

            y_man_norm_NCSEM = 0.0;
            for(int i = 0; i < 2; i++) y_man_norm_NCSEM += yv[i]*yv[i];
            y_man_norm_NCSEM = sqrt(y_man_norm_NCSEM);

            if(y_man_norm_NCSEM < 0.5)
            {
                // Projection on the center manifold
                NCprojCCMtoCM(yv, tv, SEML.us.n, sproj, CENTER_MANIFOLD_TFCc, MIcoc, Vcoc);
                //yvproj = W(sproj, tv)
                RCMtoNCbyTFC(sproj, tv, SEML.us.n, OFTS_ORDER, OFS_ORDER, 4, CENTER_MANIFOLD_TFCc, ofs, Mcoc, Vcoc, yvproj, 1);
                //Distance of projection in SEM coordinates
                NCtoSEM(tv, yv, yv_SEM, &SEML);
                NCtoSEM(tv, yvproj, yvproj_SEM, &SEML);
                proj_dist_SEM = 0.0;
                for(int i = 0; i < 3; i++) proj_dist_SEM += (yvproj_SEM[i] - yv_SEM[i])*(yvproj_SEM[i] - yv_SEM[i]);
                proj_dist_SEM = sqrt(proj_dist_SEM);
            }
            else proj_dist_SEM = 1e5;

            //Update distance min if necessary
            if(proj_dist_SEM < min_proj_dist_SEM)
            {
                min_proj_dist_SEM = proj_dist_SEM;
                kmin  = kman;
                for(int i = 0; i < 6; i++) final_state_CMU_SEM[i][kpos] = yvproj_SEM[i];
                for(int i = 0; i < 6; i++) projected_state_CMU_SEM[i][kpos] = yvproj[i];
            }

        }

        //----------------------------------------------------------
        //Store results
        //----------------------------------------------------------
        indexMin[kpos] = (int) kmin;
        distMin[kpos]  = min_proj_dist_SEM;

        //----------------------------------------------------------
        //Display completion
        //----------------------------------------------------------
        displayCompletion("manOrbitProj", 100.0*(index++)/N );
    }


    //----------------------------------------------------------
    //Sorting DH by indexes
    //----------------------------------------------------------
    vector<size_t> sortId = sort_indexes(distMin);


    //----------------------------------------------------------
    //Display results
    //----------------------------------------------------------
    int ksortpos;
    double temp;
    for(int kpos = 0; kpos <= number_of_sol; kpos++)
    {
        ksortpos = sortId[kpos];
        //----------------------------------------------------------
        //Plot
        //----------------------------------------------------------
        temp  = (double) indexMin[ksortpos];
        gnuplot_plot_xy(h1, (double*) &temp, &distMin[ksortpos], 1, (char*)"", "points", "1", "2", 4);
        cout << indexMin[ksortpos] << "   "  << distMin[ksortpos] << endl;

        //----------------------------------------------------------
        //Plot the solution +  the good point
        //----------------------------------------------------------
        double **y_man_NCSEM  = dmatrix(0, 5, 0, Nman);
        double **y_man_SEM   = dmatrix(0, 5, 0, Nman);
        double *t_man_SEM    = dvector(0, Nman);
        for(int kman = 0; kman <= Nman; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) y_man_NCSEM[i][kman] = yTens[i][ksortpos][kman];
            t_man_SEM[kman] = tTens[ksortpos][kman];

            // Current state
            for(int i = 0; i < 6; i++) y_cmu_SEM[i][kpos][kman] = yTens[i][ksortpos][kman];
            t_cmu_SEM[kpos][kman] = tTens[ksortpos][kman];
        }

        //----------------------------------------------------------
        //Update final data
        // Careful: we save the projection state at the end of y_cmu_SEM
        // Careful: we save the indix of the minimum projection distance the end of t_cmu_SEM
        //----------------------------------------------------------
        for(int i = 0; i < 6; i++) y_cmu_SEM[i][kpos][Nman+1] = projected_state_CMU_SEM[i][ksortpos];
        t_cmu_SEM[kpos][Nman+1] = indexMin[ksortpos];

        //To SEM coordinates
        NCtoSEM_vec(y_man_NCSEM, t_man_SEM, y_man_SEM, Nman, &SEML);
        //plot
        gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
        //plot the right solution
        gnuplot_plot_xyz(h2, &y_man_SEM[0][indexMin[ksortpos]], &y_man_SEM[1][indexMin[ksortpos]], &y_man_SEM[2][indexMin[ksortpos]], 1, (char*)"", "points", "1", "2", 7);
        gnuplot_plot_xyz(h2, &final_state_CMU_SEM[0][ksortpos], &final_state_CMU_SEM[1][ksortpos], &final_state_CMU_SEM[2][ksortpos], 1, (char*)"", "points", "1", "2", 8);
    }


    //---------------------------------------------------------------------
    // Save
    // Careful: we save the projection state at the end of y_cmu_SEM
    // Careful: we save the indix of the minimum projection distance at the end of t_cmu_SEM
    //---------------------------------------------------------------------
    cout << "here" << endl;
    writeManifold_bin(t_cmu_SEM, y_cmu_SEM, order_em, size_em, number_of_sol+1, Nman+2, TYPE_MAN_PROJ);

    //---------------------------------------------------------------------
    //Free
    //---------------------------------------------------------------------
    free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
    free_d3tensor(y_cmu_SEM, 0, 5, 0, number_of_sol, 0, Nman+1);
    free_dmatrix(tTens, 0, N, 0, Nman);
    free_dmatrix(t_cmu_SEM, 0, number_of_sol, 0, Nman+1);


    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

}

/**
 *  \brief Differential correction process on the solutions obtained with manOrbit+manOrbitPostProcess+manOrbitProj. WORK IN PROGRESS
 **/
void manOrbitLambert(int order_em,
                     int size_em,
                     int type_em,
                     int order_sem,
                     vector<Oftsc> &CM,
                     vector<Oftsc> &CM_TFC,
                     matrix<Ofsc>  &Mcoc,
                     matrix<Ofsc>  &Pcoc,
                     matrix<Ofsc>  &MIcoc,
                     matrix<Ofsc>  &PIcoc,
                     vector<Ofsc>  &Vcoc)
{
    //Change scope
    changeDCS(SEML, F_SEM);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);


    //-------------------------------------------------------------------------------
    //Getting back the data from SEM files
    //-------------------------------------------------------------------------------
    //Read data sizeOrbit
    int N, Nman;
    getLenghtManifold_bin(order_em, size_em, &N, &Nman, type_em);
    N--;     //need to shift to take into account that vectors starts at 0
    Nman--;  //need to shift to take into account that vectors starts at 0

    cout << "N = " << N << endl;
    cout << "Nman = " << Nman << endl;

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //Read data
    readManifold_bin(tTens, yTens, order_em, size_em, N+1, Nman+1, type_em);

    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl  *h2;
    h2 = gnuplot_init();

    //----------------------------------------------------------
    //Notable points in SEM system
    //----------------------------------------------------------
    double **semP = dmatrix(0, 6, 0, 2);
    semPoints(0.0, semP);
    semPlot(h2, semP);

    //---------------------------------------------------------------------
    // Initialisation of the orbit
    //---------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF, 6, PREC_HSTART, qbfbp_vfn_novar, NULL, &SEML);

    char ch;
    //----------------------------------------------------------
    //Display results
    //----------------------------------------------------------
    double yv0[6], yvf[6], yvfproj[6];
    double yv_SEM[6];
    double t0, tf;
    int ksol;
    for(int kpos = 0; kpos <= 0; kpos++)
    {
        //---------------------------------------------------------------------
        // Plot the valuable points
        //---------------------------------------------------------------------
        //Initial point: ksol = 0
        ksol = 0;
        for(int i = 0; i < 6; i++) yv0[i] = yTens[i][kpos][ksol];
        t0 = tTens[kpos][ksol];
        //To SE coordinates
        NCtoSEM(t0, yv0, yv_SEM, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yv_SEM[0], &yv_SEM[1], &yv_SEM[2], 1, (char*)"", "points", "1", "2", 6);

        //Proj point: ksol = Nman
        ksol = Nman;
        for(int i = 0; i < 6; i++) yvfproj[i] = yTens[i][kpos][ksol];
        tf = tTens[kpos][ksol];
        //To SE coordinates
        NCtoSEM(tf, yvfproj, yv_SEM, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yv_SEM[0], &yv_SEM[1], &yv_SEM[2], 1, (char*)"", "points", "1", "2", 8);

        //Final point: ksol = last indix stored in tTens
        ksol = (int) tTens[kpos][Nman];
        tf = tTens[kpos][ksol];
        for(int i = 0; i < 6; i++) yvf[i] = yTens[i][kpos][ksol];
        //To SE coordinates
        NCtoSEM(tf, yvf, yv_SEM, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yv_SEM[0], &yv_SEM[1], &yv_SEM[2], 1, (char*)"", "points", "1", "2", 7);


        //---------------------------------------------------------------------
        //Equivalent from data file
        //---------------------------------------------------------------------
        //From data
        double **y_man_NCSEM  = dmatrix(0, 5, 0, Nman);
        double *t_man_SEM     = dvector(0, Nman);
        double **y_man_SEM    = dmatrix(0, 5, 0, Nman);
        for(int kman = 0; kman <= Nman; kman++)
        {
            for(int i = 0; i < 6; i++) y_man_NCSEM[i][kman] = yTens[i][kpos][kman];
            t_man_SEM[kman] = tTens[kpos][kman];
        }

        //To SE coordinates
        NCtoSEM_vec(y_man_NCSEM, t_man_SEM, y_man_SEM, Nman, &SEML);

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        gnuplot_plot_xyz(h2, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], ksol+1, (char*)"", "lines", "1", "2", 5);

        printf("Press ENTER to go on\n");
        scanf("%c",&ch);


        //---------------------------------------------------------------------
        //Local variables
        //---------------------------------------------------------------------
        //New vectors
        int Nman0 = 5000;
        double **yManNCSEM0 = dmatrix(0, 5, 0, Nman0);
        double **yManSEM0   = dmatrix(0, 5, 0, Nman0);
        double *tManSEM0    = dvector(0, Nman0);


        //---------------------------------------------------------------------
        // Initialisation
        //---------------------------------------------------------------------
        double st[5];
        //------------------
        // ODE
        //------------------
        OdeStruct driver;
        //Root-finding
        const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
        //Stepper
        const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
        //Init ode structure
        init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF, 6, PREC_HSTART, qbfbp_vfn_novar, NULL, &SEML);
        //The default interval of projection is set to Tproj = T/5
        double tproj = SEML.us.T/5.0;
        //------------------
        //Orbit structure
        //------------------
        Ofsc orbit_ofs(OFS_ORDER);
        SingleOrbit orbit;
        //Init routine
        init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, t0, tf, tproj, &driver, &SEML);

        //---------------------------------------------------------------------
        // Update the state
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, yv0, t0);

        //---------------------------------------------------------------------
        //Integration
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0, tf, yManNCSEM0, tManSEM0, Nman0, 0);

        //----------------------------------------------------------
        //To SE coordinates
        //----------------------------------------------------------
        NCtoSEM_vec(yManNCSEM0, tManSEM0, yManSEM0, Nman0, &SEML);

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        gnuplot_plot_xyz(h2, yManSEM0[0], yManSEM0[1], yManSEM0[2], Nman0+1, (char*)"", "lines", "1", "2", 4);

        cout << "Final state: " << endl;
        for(int i = 0; i <6; i++) cout << yManNCSEM0[i][Nman0] << endl;

        printf("Press ENTER to go on\n");
        scanf("%c",&ch);

        //---------------------------------------------------------------------
        // Lambert arc
        //---------------------------------------------------------------------
        OdeStruct driver2;
        init_ode_structure(&driver2, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF, 42, PREC_HSTART, qbfbp_vfn_varnonlin, NULL, &SEML);

        double ystart[42];
        //Storing initial position and momenta into ystart
        for(int i = 0; i < 6; i++) ystart[i] = yv0[i];
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_calloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(ystart, Id, 6, 6, 6);

        //---------------------------------------------------------------------
        //Differential correction process
        //---------------------------------------------------------------------
        double tof = tf;
        differential_correction_mns(ystart, yvfproj, t0, &tof, 1e-5, driver2.d, 42, 1);
        //differential_correction_ft(ystart, yvfproj, t0, &tof, 1e-4, driver2.d, 42, 1);

        cout << "Final state: " << endl;
        for(int i = 0; i <6; i++) cout << ystart[i] << endl;
        //To SE coordinates
        NCtoSEM(t0, ystart, yv_SEM, &SEML);
        gnuplot_plot_xyz(h2, &yv_SEM[0], &yv_SEM[1], &yv_SEM[2], 1, (char*)"", "points", "1", "5", 7);


        //---------------------------------------------------------------------
        // Update the state
        //---------------------------------------------------------------------
        orbit_update_ic(orbit, st, ystart, t0);

        //---------------------------------------------------------------------
        //Integration
        //---------------------------------------------------------------------
        trajectory_integration_grid(orbit, t0, tof, yManNCSEM0, tManSEM0, Nman0, 0);

        //----------------------------------------------------------
        //To SE coordinates
        //----------------------------------------------------------
        NCtoSEM_vec(yManNCSEM0, tManSEM0, yManSEM0, Nman0, &SEML);

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        gnuplot_plot_xyz(h2, yManSEM0[0], yManSEM0[1], yManSEM0[2], Nman0+1, (char*)"", "lines", "1", "2", 2);


        //----------------------------------------------------------
        // Final Velocity gap
        //----------------------------------------------------------
        double yprojVel[6], yfinalVel[6];
        //Final position on the new trajectory
        for(int i = 0; i <6; i++) yvf[i] = yManSEM0[i][Nman0];
        SEMmtoSEMv(tof, yvf, yfinalVel, &SEML);
        //Projection position
        NCtoSEM(tf, yvfproj, yv_SEM, &SEML);
        SEMmtoSEMv(tf, yv_SEM, yprojVel, &SEML);

        double dV = 0.0;
        for(int i = 3; i < 6; i++) dV += (yprojVel[i] -  yfinalVel[i])*(yprojVel[i] -  yfinalVel[i]);
        dV = sqrt(dV);
        cout << "Final Velocity gap: (m/s)" << endl;
        cout << 1e3*dV*SEML.cs_sem.cr3bp.L*2*M_PI/(SEML.cs_sem.cr3bp.T) << endl;

        //----------------------------------------------------------
        // Starting Velocity gap
        //----------------------------------------------------------
        //Initial position on the old trajectory
        NCtoSEM(t0, yv0, yv_SEM, &SEML);
        SEMmtoSEMv(t0, yv_SEM, yprojVel, &SEML);
        //Initial position on the new trajectory
        for(int i = 0; i <6; i++) yv0[i] = ystart[i];
        NCtoSEM(t0, yv0, yv_SEM, &SEML);
        SEMmtoSEMv(t0, yv_SEM, yfinalVel, &SEML);


        dV = 0.0;
        for(int i = 3; i < 6; i++) dV += (yprojVel[i] -  yfinalVel[i])*(yprojVel[i] -  yfinalVel[i]);
        dV = sqrt(dV);
        cout << "Initial Velocity gap: (m/s)" << endl;
        cout << 1e3*dV*SEML.cs_sem.cr3bp.L*2*M_PI/(SEML.cs_sem.cr3bp.T) << endl;

    }
    //---------------------------------------------------------------------
    //Free
    //---------------------------------------------------------------------
    free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
    free_dmatrix(tTens, 0, N, 0, Nman);

    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

}



/**
 *  \brief Plot of the manifold branches obtained with manOrbit+manOrbitPostProcess.
 **/
void manOrbitPlot(int ofts_order,
                  int sizeOrbit,
                  int number_of_sol,
                  int Nmant,
                  int type)
{
    //Read data sizeOrbit
    int N0, Nman0, N, Nman;
    getLenghtManifold_bin(ofts_order, sizeOrbit, &N0, &Nman0, type);
    N0--;   //need to shift to take into account that vectors starts at 0
    Nman0--; //need to shift to take into account that vectors starts at 0

    //Use value from file if the user wants it
    if(number_of_sol < 0) N = N0;
    else N = number_of_sol;
    //Use value from file if the user wants it
    if(Nmant < 0) Nman = Nman0;
    else Nman = Nmant;

    //Check the values are consistent
    if(N0 < N || Nman0 < Nman)
    {
        cout << "manOrbitPlot: wrong inputs" << endl;
        cout << "N = " << N << ", but N0 = " << N0 << endl;
        cout << "Nman = " << Nman << ", but Nman0 = " << Nman0 << endl;
        return ;
    }

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //Read data
    int status = readManifold_bin(tTens, yTens, ofts_order, sizeOrbit, N+1, Nman+1, type);

    //Plot loop
    if(status)
    {
        //Gnuplot window
        gnuplot_ctrl  *h1;
        h1 = gnuplot_init();


        //Temp variable
        double **y_man_NCSEM = dmatrix(0, 5, 0, Nman);
        double **y_man_SEM   = dmatrix(0, 5, 0, Nman);
        double *tManNCSEM  = dvector(0, Nman);

        for(int k = 0; k <= N; k++)
        {
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) y_man_NCSEM[i][p] = yTens[i][k][p];
                tManNCSEM[p] = tTens[k][p];
            }

            //To SEM coordinates
            NCtoSEM_vec(y_man_NCSEM, tManNCSEM, y_man_SEM, Nman, &SEML);
            gnuplot_plot_xyz(h1, y_man_SEM[0], y_man_SEM[1], y_man_SEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
        }

    }
    else
    {
        cout << "manOrbitPlot: error while loading the file" << endl;
    }

    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);

}



/**
 *  \brief Same as manOrbit, but with integration in the SEM units and coordinates.
 **/
int manOrbitSEM(double tmax_on_manifold_EM,
                int ofts_order,
                int sizeOrbit,
                int number_of_sol,
                int Nman,
                int isPar,
                vector<Oftsc> &CM,
                vector<Oftsc> &CM_TFC,
                matrix<Ofsc>  &Mcoc,
                matrix<Ofsc>  &Pcoc,
                matrix<Ofsc>  &MIcoc,
                matrix<Ofsc>  &PIcoc,
                vector<Ofsc>  &Vcoc)
{
    //----------------------------------------------------------
    // Get the number of data lines from file
    //----------------------------------------------------------
    int N;
    int N0 = getLineNumber(ofts_order, sizeOrbit, TYPE_CU)-1;

    //----------------------------------------------------------
    //If the file exists and contains data, go on
    //----------------------------------------------------------
    if(N0 > 0)
    {
        if(number_of_sol <= 0) N = N0;
        else N = number_of_sol;

        if(N > N0)
        {
            cout << "manOrbit: wrong input, number_of_sol = " << number_of_sol << ", but N0 = " << N0 << endl;
            return 0;
        }

        //----------------------
        // Read data from file
        //----------------------
        double **yNCU = dmatrix(0, 5, 0, N);
        double **sNCU = dmatrix(0, 4, 0, N);
        double *tNCU  = dvector(0, N);
        int status    = readOrbit(tNCU, yNCU, sNCU, ofts_order, sizeOrbit, N+1, TYPE_CU);

        //----------------------
        // Loop on all elements
        //----------------------
        if(status)
        {
            //CHANGE SCOPE
            changeDCS(SEML, F_SEM);

            //To store all data
            double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
            double **tTens  = dmatrix(0, N, 0, Nman);


            //The default interval of projection is set to Tproj = T/5
            double tproj = SEML.us.T/5.0;

            #pragma omp parallel for if(isPar)
            for(int k = 0; k <= 0; k++)
            {
                //---------------------------------------------------------------------
                //Temp variables
                //---------------------------------------------------------------------
                double yvEM[6], tvEM;
                double yv[6], st[5], tv;
                yvEM[0] = yNCU[0][k];
                yvEM[1] = yNCU[1][k];
                yvEM[2] = yNCU[2][k];
                yvEM[3] = yNCU[3][k];
                yvEM[4] = yNCU[4][k];
                yvEM[5] = yNCU[5][k];

                st[0] = sNCU[0][k];
                st[1] = sNCU[1][k];
                st[2] = sNCU[2][k];
                st[3] = sNCU[3][k];
                st[4] = sNCU[4][k];
                tvEM  = tNCU[k];


                //---------------------------------------------------------------------
                //To NC SEM
                //---------------------------------------------------------------------
                NCEMmtoNCSEMm(tvEM, yvEM, yv, &SEML);
                tv = tvEM * SEML.us_em.ns;


                //---------------------------------------------------------------------
                //Local variables
                //---------------------------------------------------------------------
                double **y_man_NCEM  = dmatrix(0, 5, 0, Nman);
                double **y_man_NCSEM = dmatrix(0, 5, 0, Nman);
                double *t_man_EM     = dvector(0, Nman);

                //---------------------------------------------------------------------
                // Initialisation
                //---------------------------------------------------------------------
                //------------------
                // ODE
                //------------------
                OdeStruct driver;
                //Root-finding
                const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
                //Stepper
                const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
                //Init ode structure
                init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF, 6, PREC_HSTART, qbfbp_vfn_novar, NULL, &SEML);

                //------------------
                //Orbit structure
                //------------------
                Ofsc orbit_ofs(OFS_ORDER);
                SingleOrbit orbit;
                //Init routine
                init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 5, 1, tv, tv+tmax_on_manifold_EM, tproj, &driver, &SEML);

                //---------------------------------------------------------------------
                // Update the state
                //---------------------------------------------------------------------
                orbit_update_ic(orbit, st, yv, tv);

                //---------------------------------------------------------------------
                //Integration
                //---------------------------------------------------------------------
                trajectory_integration_grid(orbit, tv, tv+tmax_on_manifold_EM, y_man_NCEM, t_man_EM, Nman, 0);

                //---------------------------------------------------------------------
                //Integration check
                //---------------------------------------------------------------------
                reset_ode_structure(orbit.driver);
                double yv0[6];
                for(int i = 0; i <6; i++) yv0[i] = yv[i];

                double t0 = t_man_EM[0];
                int kk = Nman;

                cout << "Time" << endl;
                cout << t_man_EM[0] << "  " << t0 << "   " << endl;


                cout << "Init" << endl;
                for(int i = 0; i <6; i++) cout << y_man_NCEM[i][0] - yv[i] << endl;


                //Integration
                gsl_odeiv2_driver_apply (orbit.driver->d, &t0, t_man_EM[kk], yv);

                cout << "Time" << endl;
                cout << t_man_EM[kk] << "  " << t0 << "   " << endl;

                cout << "Result" << endl;
                for(int i = 0; i <6; i++) cout << y_man_NCEM[i][kk] - yv[i] << endl;

                //Integration
                reset_ode_structure(orbit.driver);
                t0 = t_man_EM[0];
                gsl_odeiv2_driver_apply(orbit.driver->d, &t0, t_man_EM[kk], yv0);

                cout << "Result" << endl;
                for(int i = 0; i <6; i++) cout << yv[i] - yv0[i] << endl;

                //---------------------------------------------------------------------
                //Storage
                //---------------------------------------------------------------------
                #pragma omp critical
                {
                    //To NCSEM coordinates
                    NCEMmtoNCSEMm_vec(y_man_NCEM, t_man_EM, y_man_NCSEM, Nman, &SEML);
                    for(int p = 0; p <= Nman; p++)
                    {
                        for(int i = 0; i < 6; i++) yTens[i][k][p] = y_man_NCSEM[i][p];
                        tTens[k][p] = t_man_EM[p]*SEML.us_em.ns; //time is also given in SEM units!
                    }
                }

                //---------------------------------------------------------------------
                //Free
                //---------------------------------------------------------------------
                free_dmatrix(y_man_NCEM, 0, 5, 0, Nman);
                free_dmatrix(y_man_NCSEM, 0, 5, 0, Nman);
                free_dvector(t_man_EM, 0, Nman);
            }


            //---------------------------------------------------------------------
            //Save
            //---------------------------------------------------------------------
            writeManifold_bin(tTens, yTens, ofts_order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

            //---------------------------------------------------------------------
            //Free
            //---------------------------------------------------------------------
            free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
            free_dmatrix(tTens, 0, N, 0, Nman);
        }
        else
        {
            cout << "cusOrbit: error during data reading." << endl;
            return 0;
        }
        //Free
        free_dmatrix(yNCU, 0, 5, 0, N);
        free_dvector(tNCU, 0, N);
    }
    return 1;
}



//---------------------------------------------------------------------------------------------------
//
// OLD VERSION FOR CONNECTIONS: EML2 to SEMLi
//
//---------------------------------------------------------------------------------------------------
//----------------------
                // Finding connections with SEMLi (NEW)
                //----------------------
                //            int Size = floor(fabs(st0[0]));
                //            int Nt = 20;
                //            int MSIZE = 1000;

                //            //1. Building an EML2 orbit on a grid
                //            gridOrbit(st0, 0, 2*SEML.us.T, +1e-2, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
                //            //2. Computing the unstable directions along this orbit
                //            cusOrbit(OFTS_ORDER, Size, +1e-5, CM_TFC, Mcoc, MIcoc, Vcoc, 1);
                //            //3. Integrating forward these unstable directions on a grid
                //            manOrbit(5*SEML.us.T, OFTS_ORDER, Size, -1, MSIZE, 1, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
                //            //4. Sort the manifold branches wrt to a certain criterion
                //            manOrbitPostProcess(OFTS_ORDER, Size, -1, 1);
                //            //5. Project the manifold branches on the center manifold of SEML1,2
                //            manOrbitProj(OFTS_ORDER, Size, TYPE_MAN_SORT_DR, OFTS_ORDER, -1, MSIZE, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc, 1);
                //            //6. Diplay the solution, or...
                //            //manOrbitPlot(OFTS_ORDER, Size, Nt, MSIZE, TYPE_MAN_SORT_DR);
                //            //6. Refine the solutions to actually target the projected state on the center manifold
                //            manOrbitLambert(OFTS_ORDER, Size, TYPE_MAN_PROJ, OFTS_ORDER, CM_NC, CM_TFC, Mcoc, Pcoc, MIcoc, PIcoc, Vcoc);
