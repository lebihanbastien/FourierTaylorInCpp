#include "single_orbit.h"
#include "timec.h"
#include <list>
#include <algorithm>    // std::sort
#include <vector>       // std::vector
#include "diffcorr.h"

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          SingleOrbit structure
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Initialize one SingleOrbit structure
 **/
void init_orbit(SingleOrbit &orbit,
                vector<Oftsc>*  W,
                vector<Oftsc>*  Wh,
                matrix<Ofsc>*  PC,
                matrix<Ofsc>*  CQ,
                vector<Ofsc>*  V,
                Ofsc* orbit_ofs,
                int order,
                int ofs_order,
                int reduced_nv,
                int isGS,
                double t0,
                double tf,
                double tproj,
                OdeStruct *driver,
                QBCP_L *qbcp_l)
{
    //-----------
    //Parameterization (common to all orbits)
    //-----------
    orbit.W          =  W;            //z(t) = W(s(t), t)
    orbit.Wh         =  Wh;           //zh(t) = Wh(s(t), t)

    //-----------
    //COC (common to all orbits)
    //-----------
    orbit.PC  = PC;            //COC matrix
    orbit.CQ  = CQ;            //inv COC matrix
    orbit.V   = V;             //COC vector

    //-----------
    //Orders
    //-----------
    orbit.ofs        =  orbit_ofs;    //Auxiliary Ofs object
    orbit.order      =  order;        //order of the expansions
    orbit.ofs_order  =  ofs_order;    //order of the Fourier coefficients
    orbit.reduced_nv =  reduced_nv;   //reduced number of variables
    orbit.isGS       =  isGS;         //Was the pm obtained through graph style?

    //-----------
    //Characteristics
    //-----------
    orbit.z0  = dvector(0, 5);                 //Initial position in NC coordinates dim = 6
    orbit.si  = dvector(0, REDUCED_NV-1);      //Initial RCM configuration dim = 4
    orbit.s0d = dvector(0, 2*REDUCED_NV-1);    //Initial position in CCM8 coordinates (real+imag part) dim = 8
    orbit.xf  = dvector(0, 5);                 //Final position NC dim = 6
    orbit.s0  = dcvector(0,REDUCED_NV-1);      //Initial position in CCM4 coordinates (real+imag part) dim = 4
    orbit.t0  = t0;                            //Initial time
    orbit.tf  = tf;                            //Final time after computation
    orbit.tproj  = tproj;                      //default time between each projection
    orbit.tprojmin  = 1e-2;                    //minimum time between each projection
    orbit.ePmax = (qbcp_l->fwrk == F_SEM)? 1e-3:1e-5; //maximum projection distance allowed

    //-----------
    //ODE integration
    //-----------
    orbit.driver  = driver;              //NC ode struct

    //-----------
    //Parent
    //-----------
    orbit.qbcp_l = qbcp_l;  //QBCP around a given Li point (parent)

    //-----------
    //Pulsation
    //-----------
    orbit.n = qbcp_l->us.n;
}


/**
 *   \brief Free one orbit
 **/
void free_orbit(SingleOrbit *orbit)
{
    //-----------
    //Characteristics
    //-----------
    free_dvector(orbit->z0,  0, 5);
    free_dvector(orbit->si,  0, orbit->reduced_nv-1);
    free_dvector(orbit->s0d, 0, 2*orbit->reduced_nv-1);
    free_dvector(orbit->xf,  0, 5);
    free_dcvector(orbit->s0, 0, orbit->reduced_nv-1);
    //-----------
    //Ode
    //-----------
    free_ode_structure(orbit->driver);
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Initial conditions
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0)
{
    //------------------------------------------
    // 1. Update si
    //------------------------------------------
    for(int p = 0; p < orbit.reduced_nv; p++) orbit.si[p] = si[p];

    //------------------------------------------
    // 2. Update s0
    //------------------------------------------
    RCMtoCCM(si, orbit.s0, orbit.reduced_nv);

    //------------------------------------------
    // 2. Update s0d
    //------------------------------------------
    RCMtoCCM8(si, orbit.s0d, orbit.reduced_nv);

    //------------------------------------------
    // 4. Update z0
    //------------------------------------------
    //z0 = W(si, 0.0)
    RCMtoNCbyTFC(si,
                 t0,
                 orbit.n,
                 orbit.order,
                 orbit.ofs_order,
                 orbit.reduced_nv,
                 *orbit.Wh,
                 *orbit.ofs,
                 *orbit.PC,
                 *orbit.V,
                 orbit.z0,
                 orbit.isGS);
}

/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 *          and an array of initial NC conditions z0.
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], const double z0[], double t0)
{
    //------------------------------------------
    // 1. Update si
    //------------------------------------------
    for(int p = 0; p < orbit.reduced_nv; p++) orbit.si[p] = si[p];

    //------------------------------------------
    // 2. Update s0
    //------------------------------------------
    RCMtoCCM(si, orbit.s0, orbit.reduced_nv);

    //------------------------------------------
    // 2. Update s0d
    //------------------------------------------
    RCMtoCCM8(si, orbit.s0d, orbit.reduced_nv);

    //------------------------------------------
    // 4. Update z0
    //------------------------------------------
    //z0 = W(si, 0.0)
    for(int p = 0; p < NV; p++) orbit.z0[p] = z0[p];
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Integration
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(SingleOrbit &orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *ePm,
                   int *nreset,
                   int isResetOn)
{
    int status;
    double yvp[6], yvi[6];
    cdouble scp[REDUCED_NV];

    //----------------------
    //Evolve one step of z(t)
    //----------------------
    orbit.driver->h = (t1 >= *t)? fabs(orbit.driver->h):-fabs(orbit.driver->h);
    status = gsl_odeiv2_evolve_apply (orbit.driver->e, orbit.driver->c, orbit.driver->s, &orbit.driver->sys, t, t1, &orbit.driver->h, yv);
    if (status != 0)
    {
        cout << "error in gslc_proj_step: integration of z(t) has gone wrong. break." << endl;
        return -1;
    }

    //----------------------
    //Projection if necessary
    //----------------------
    if(isResetOn && fabs(*t-t0) > fabs(*nreset*orbit.tproj))
    {
        //----------------------
        //Projection tools
        //----------------------
        double omega1 = cimag(orbit.Wh->at(0).getCoef(1,0)->getCoef(0));
        double omega3 = cimag(orbit.Wh->at(2).getCoef(1,1)->getCoef(0));

        //----------------------
        // Projection on the center manifold
        //----------------------
        //Get the closest point on the center manifold
        NCprojCCM(yv, *t, SEML.us_em.n, OFS_ORDER, *orbit.CQ, *orbit.V, omega1, omega3, scp, REDUCED_NV);
        //Update the state
        CCMtoNCbyTFC(scp, *t, orbit.qbcp_l->us.n, orbit.order,  orbit.ofs_order,  *orbit.Wh,  *orbit.ofs, *orbit.PC, *orbit.V,  yvp,  orbit.isGS);
        //For comparison
        for(int i = 0; i <6; i++) yvi[i] = yv[i];
        // Copy of yvp in current state
        for(int i=0; i<6; i++) yv[i]  = yvp[i];

        //-----------------
        // Get the current projection error
        //-----------------
        //Get the current error
        *ePm = fabs(yvi[0] - yv[0]);
        for(int i = 1; i <6 ; i++)
        {
            if(fabs(yvi[i] - yv[i]) > *ePm) *ePm = fabs(yvi[i] - yv[i]);
        }

        if(*ePm > orbit.ePmax)
        {
            cout << "Warning: Reset n° " << *nreset << ". ePm = " << *ePm << endl;
            return -2;
        }

        //-----------------
        //Reset ode structure for next step
        //-----------------
        reset_ode_structure(orbit.driver);

        //-----------------
        //One additional reset
        //-----------------
        *nreset = *nreset +1;
    }


    return 0;
}



/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(SingleOrbit &orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *ePm,
                     int *nreset,
                     int isResetOn)
{
    //reset_ode_structure(orbit.driver);
    int status;
    do
    {
        status = gslc_proj_step(orbit, yv, t, t0, t1, ePm, nreset, isResetOn);
    }
    while(status == 0 && fabs(*t)<fabs(t1));

    return status;
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Projection on (un)stable manifold
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CMh, MIcoc, Vcoc). Then the CCM state is extended by adding a non-null direction
 *        along the hyperbolic direction (sti[4]).
 **/
void NCprojCCMtoCUS(double *yv, double tv, double sti[5], vector<Oftsc> &CMh, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc)
{
    //Projection tools
    double omega1 = cimag(CMh[0].getCoef(1,0)->getCoef(0));
    double omega3 = cimag(CMh[2].getCoef(1,1)->getCoef(0));
    cdouble scp[REDUCED_NV];
    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, SEML.us.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, REDUCED_NV);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, REDUCED_NV);
    //Add a given quantity on the hyperbolic direction
    sti[4] = epsilon;
}

/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CMh, MIcoc, Vcoc).
 **/
void NCprojCCMtoCM(double *yv, double tv, double sti[5], vector<Oftsc> &CMh, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc)
{
    //Projection tools
    double omega1 = cimag(CMh[0].getCoef(1,0)->getCoef(0));
    double omega3 = cimag(CMh[2].getCoef(1,1)->getCoef(0));
    cdouble scp[REDUCED_NV];
    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, SEML.us.n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, REDUCED_NV);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, REDUCED_NV);
}

//---------------------------------------------------------------------------------------------------------------------------------------
//                  Energy on vectors
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Hamiltonian along one orbit in SEM units and SEM coordinates
 **/
void HSEM_vec(double *tSEM, double **ySEM, double *Hvec, int N, QBCP_L *qbcp_l)
{
    double yv[6];
    for(int i = 0; i <= N; i++)
    {
        for(int p = 0; p < 6; p++) yv[p] = ySEM[p][i];
        Hvec[i] = qbfbp_H_SEM(tSEM[i], yv, qbcp_l); //careful: time in SEM units!
    }
}

/**
 *  \brief Hamiltonian along the dyneq of SEMLi for a given time vector in SEM units and SEM coordinates
 **/
void HSEMLi_vec(double *tSEM, double *Hvec, int N, QBCP_L *qbcp_l)
{
    double yv[6], yvSEM[6];
    for(int p = 0; p < 6; p++) yv[p] = 0.0;
    for(int i = 0; i <= N; i++)
    {
        NCtoSEM(tSEM[i], yv, yvSEM, qbcp_l);
        Hvec[i] = qbfbp_H_SEM(tSEM[i], yvSEM, qbcp_l);
    }
}

/**
 *  \brief Hamiltonian along the dyneq of EMLi for a given time vector in SEM units and SEM coordinates
 **/
void HEMLi_in_SEM_vec(double *tEM, double *Hvec, int N, QBCP_L *qbcp_l)
{
    double yv[6], yvSEM[6];
    for(int p = 0; p < 6; p++) yv[p] = 0.0;
    for(int i = 0; i <= N; i++)
    {
        NCEMmtoSEMm(tEM[i], yv, yvSEM, qbcp_l);
        Hvec[i] = qbfbp_H_SEM(tEM[i]*qbcp_l->us_em.ns, yvSEM, qbcp_l); //careful, time must be in SEM units at this step!
    }
}



//---------------------------------------------------------------------------------------------------------------------------------------
//
//         Update points
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Update some key positions of notable points in the EM system
 **/
void emPoints(double t, double **emP)
{
    //----------------------------------
    //All in EM coordinates, at time t in EM units
    //----------------------------------
    //Earth position
    evaluateCoef(emP[0], t, SEML.us_em.n, SEML.nf, SEML.cs_em.Pe, 3);
    //Moon position
    evaluateCoef(emP[1], t, SEML.us_em.n, SEML.nf, SEML.cs_em.Pm, 3);
    //Sun position
    evaluateCoef(emP[2], t, SEML.us_em.n, SEML.nf, SEML.cs_em.Ps, 3);

    //EML1 position
    for(int i = 0; i < 3; i++) emP[3][i] = SEML.cs_em.cr3bp.l1.position[i];
    emP[3][0] *= -1.0;
    //EML2 position
    for(int i = 0; i < 3; i++) emP[4][i] = SEML.cs_em.cr3bp.l2.position[i];
    emP[4][0] *= -1.0;

    //SEML1 position in EM coordinates
    double temp1[6], temp2[6];
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l1.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in EM coordinates, starting with a time in SEM units !
    SEMmtoEMm(t*SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) emP[5][i] = temp2[i];

    //SEML2 position in EM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l2.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in EM coordinates, starting with a time in SEM units !
    SEMmtoEMm(t*SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) emP[6][i] = temp2[i];
}


/**
 *  \brief Update some key positions of notable points in the SEM system
 **/
void semPoints(double t, double **semP)
{
    //----------------------------------
    //All in SEM coordinates, at time t in SEM units
    //----------------------------------
    //Earth position
    evaluateCoef(semP[0], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
    //Moon position
    evaluateCoef(semP[1], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
    //Sun position
    evaluateCoef(semP[2], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Ps, 3);

    //SEML1 position
    for(int i = 0; i < 3; i++) semP[5][i] = SEML.cs_sem.cr3bp.l1.position[i];
    semP[5][0] *= -1.0;
    //SEML2 position
    for(int i = 0; i < 3; i++) semP[6][i] = SEML.cs_sem.cr3bp.l2.position[i];
    semP[6][0] *= -1.0;


    //EML1 position in SEM coordinates
    double temp1[6], temp2[6];
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l1.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in SEM coordinates, starting with a time in EM units !
    EMmtoSEMm(t/SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) semP[3][i] = temp2[i];

    //EML2 position in SEM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l2.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in SEM coordinates, starting with a time in EM units !
    EMmtoSEMm(t/SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) semP[4][i] = temp2[i];
}

//---------------------------------------------------------------------------------------------------------------------------------------
//
//        Plots
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Sets notable points in EM system on the gnuplot window ctrl h1
 **/
void emPlot(gnuplot_ctrl *h1, double **emP)
{
    //Misc
    gnuplot_cmd(h1, "set grid");
    gnuplot_set_xlabel(h1, (char*)"Xem [-]");
    gnuplot_set_ylabel(h1, (char*)"Yem [-]");
    gnuplot_set_zlabel(h1, (char*)"Zem [-]");
    //EML point
    switch(SEML.li_EM)
    {
    case 1:
        gnuplot_plot_xyz(h1, emP[3], emP[3]+1,  emP[3]+2, 1, (char*)"EML1", "points", "1", "2", 8);
        break;
    case 2:
        gnuplot_plot_xyz(h1, emP[4], emP[4]+1,  emP[4]+2, 1, (char*)"EML2", "points", "1", "2", 8);
        break;
    }
    //Primaries
    gnuplot_plot_xyz(h1, emP[1], emP[1]+1,  emP[1]+2, 1, (char*)"Moon", "points", "7", "2", 9);


}

/**
 *  \brief Sets notable points in SEM system on the gnuplot window ctrl h2
 **/
void semPlot(gnuplot_ctrl *h2, double **semP)
{
    //Misc
    gnuplot_cmd(h2, "set grid");
    gnuplot_set_xlabel(h2, (char*)"Xsem [-]");
    gnuplot_set_ylabel(h2, (char*)"Ysem [-]");
    gnuplot_set_zlabel(h2, (char*)"Zsem [-]");

    // SEML points
    switch(SEML.li_SEM)
    {
    case 1:
        gnuplot_plot_xyz(h2, semP[5], semP[5]+1,  semP[5]+2, 1, (char*)"SEML1", "points", "1", "2", 9);
        break;
    case 2:
        gnuplot_plot_xyz(h2, semP[6], semP[6]+1,  semP[6]+2, 1, (char*)"SEML2", "points", "2", "2", 9);
        break;
    }
    //EML points
    gnuplot_plot_xyz(h2, semP[3], semP[3]+1,  semP[3]+2, 1, (char*)"EML1", "points", "1", "2", 8);
    gnuplot_plot_xyz(h2, semP[4], semP[4]+1,  semP[4]+2, 1, (char*)"EML2", "points", "2", "2", 8);
    //Primaries
    gnuplot_plot_xyz(h2, semP[0], semP[0]+1,  semP[0]+2, 1, (char*)"Earth", "points", "7", "2", 6);
}

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          I/O orbit
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computes a data filename as a string, depending on the order, size, and type of the data.
 **/
string filenameOrbit(int order, int sizeOrbit, int type)
{
    switch(type)
    {
    case TYPE_ORBIT:
        return SEML.cs.F_PLOT+"orbit_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_CU:
        return SEML.cs.F_PLOT+"cu_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_CS:
        return SEML.cs.F_PLOT+"cs_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_MAN:
        return SEML.cs.F_PLOT+"man_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".bin";
    case TYPE_MAN_SORT_DR:
        return SEML.cs.F_PLOT+"man_sort_dr_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".bin";
    case TYPE_MAN_SORT_DH:
        return SEML.cs.F_PLOT+"man_sort_dh_order_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".bin";
    case TYPE_MAN_PROJ:
        return SEML.cs.F_PLOT+"man_proj_"+numTostring(order)+"_size_"+numTostring(sizeOrbit)+".bin";
    default:
        cout << "filenameOrbit: unknown type." << endl;
        return "";
    }
}

/**
 * \brief Store the orbit (tNCE, yNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, int order, int sizeOrbit, int N, int type)
{
    string filename = filenameOrbit(order, sizeOrbit, type);
    gnuplot_fplot_txp(tNCE, yNCE[0], yNCE[1], yNCE[2], yNCE[3], yNCE[4], yNCE[5], N, filename.c_str());
    return 1;
}

/**
 * \brief Store the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, double **sNCE, int order, int sizeOrbit, int N, int type)
{
    string filename = filenameOrbit(order, sizeOrbit, type);
    gnuplot_fplot_txps(tNCE, yNCE[0], yNCE[1], yNCE[2], yNCE[3], yNCE[4], yNCE[5],
                       sNCE[0], sNCE[1], sNCE[2], sNCE[3], sNCE[4],
                       N, filename.c_str());
    return 1;
}

/**
 * \brief Get the length of the data file filenameOrbit(order, sizeOrbit, type).
 **/
int getLineNumber(int order, int sizeOrbit, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    string line;
    ifstream myfile (filename.c_str());
    int Nfile = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) ) Nfile++;
        myfile.close();
    }
    else return 0;

    return Nfile;
}

/**
 * \brief Read the orbit (tNCE, yNCE) in the  data file filenameOrbit(order, size, type)
 **/
int readOrbit(double *tNCE, double **yNCE, int order, int sizeOrbit, int N, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    int Nfile = getLineNumber(order, sizeOrbit, type);
    //If N is greater than Nfile, the end of the vectors will not be updated
    if(Nfile < N)
    {
        cout << "readOrbit: wrong input, N = " << N << ", but Nfile = " << Nfile << endl;
        return 0;
    }
    else Nfile = N;

    //------------------------
    //Updating the vectors
    //------------------------
    ifstream myfile (filename.c_str());
    if (myfile.is_open())
    {
        for(int i = 0; i < Nfile; i++)
        {
            myfile >> tNCE[i];
            for(int k = 0; k < 6; k++) myfile >> yNCE[k][i];
        }
    }
    else return 0;

    return 1;
}

/**
 * \brief Read the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(order, sizeOrbit, type)
 **/
int readOrbit(double *tNCE, double **yNCE, double **sNCE, int order, int sizeOrbit, int N, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    int Nfile = getLineNumber(order, sizeOrbit, type);
    //If N is greater than Nfile, the end of the vectors will not be updated
    if(Nfile < N)
    {
        cout << "readOrbit: wrong input, N = " << N << ", but Nfile = " << Nfile << endl;
        return 0;
    }
    else Nfile = N;

    //------------------------
    //Updating the vectors
    //------------------------
    ifstream myfile (filename.c_str());
    if (myfile.is_open())
    {
        for(int i = 0; i < Nfile; i++)
        {
            myfile >> tNCE[i];
            for(int k = 0; k < 6; k++) myfile >> yNCE[k][i];
            for(int k = 0; k < 5; k++) myfile >> sNCE[k][i];
        }
    }
    else return 0;

    return 1;
}



//---------------------------------------------------------------------------------------------------------------------------------------
//
//          I/O manifold
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 * \brief Store the manifold (tTens, yTens) in the  data file filenameOrbit(order, sizeOrbit, type). The manifodl is composed of N+1 branches, each of them
 * described by Nman+1 points.
 **/
int writeManifold_bin(double **tTens, double ***yTens, int order, int sizeOrbit, int N, int Nman, int type)
{
    string filename = filenameOrbit(order, sizeOrbit, type);
    fstream myfile;


    myfile.open (filename.c_str(), ios::binary | ios::out);
    if (myfile.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of data
        //---------------------
        resi = N;
        myfile.write((char*) &resi, sizeof(int));
        resi = Nman;
        myfile.write((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int n1 = 0; n1 < N; n1++)
        {
            //Current homogeneous polynomial
            for (int n2 = 0; n2 < Nman; n2++)
            {
                res = tTens[n1][n2];
                myfile.write((char*) &res, sizeof(double));
                for (int k = 0; k < 6; k++)
                {
                    res = yTens[k][n1][n2];
                    myfile.write((char*) &res, sizeof(double));
                }
            }
        }
        myfile.close();
    }
    else return 0;

    return 1;
}

/**
* \brief Read the manifold (tTens, yTens) in the  data file filenameOrbit(order, sizeOrbit, type). The manifodl is composed of N+1 branches, each of them
* described by Nman+1 points.
**/
int readManifold_bin(double **tTens, double ***yTens, int order, int sizeOrbit, int N, int Nman, int type)
{
    string filename = filenameOrbit(order, sizeOrbit, type);
    fstream myfile;


    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        int N0, Nman0;
        //---------------------
        //Number of data
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        N0 = resi;
        myfile.read((char*) &resi, sizeof(int));
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
                myfile.read((char*) &res, sizeof(double));
                tTens[n1][n2] = res;

                for (int k = 0; k < 6; k++)
                {
                    myfile.read((char*) &res, sizeof(double));
                    yTens[k][n1][n2] = res;
                }
            }
        }
        myfile.close();
    }
    else return 0;

    return 1;
}

/**
 * \brief Get the length of the data file filenameOrbit(order, sizeOrbit, type) containing the manifold (tTens, yTens). The manifodl is composed of N+1 branches, each of them
 * described by Nman+1 points.
 **/
int getLenghtManifold_bin(int order, int sizeOrbit, int *N, int *Nman, int type)
{
    string filename = filenameOrbit(order, sizeOrbit, type);
    fstream myfile;


    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        //---------------------
        //Number of data
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        *N = resi;
        myfile.read((char*) &resi, sizeof(int));
        *Nman = resi;
        myfile.close();
    }
    else return 0;

    return 1;
}

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          I/O CU
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computes a data filename as a string, depending on the order, sizeOrbit, and type of the data.
 **/
string filenameCUM(int order, int type)
{
    switch(type)
    {
    case TYPE_CU:
        return SEML.cs_em.F_PLOT+"cu_order_"+numTostring(order)+".bin";
    case TYPE_MAN:
        return SEML.cs_em.F_PLOT+"intcu_order_"+numTostring(order)+".bin";
    case TYPE_MAN_PROJ:
        return SEML.cs_em.F_PLOT+"projcu_order_"+numTostring(order)+".bin";
    case TYPE_MAN_SORT:
        return SEML.cs_em.F_PLOT+"sortprojcu_order_"+numTostring(order)+".bin";

    default:
        cout << "filenameOrbit: unknown type." << endl;
        return "";
    }
}

int writeCU_bin(double ****yNCE, double ****sNCE, double *tGrid, int gSize, int tSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::out);
    if (myfile.is_open())
    {
        double res;
        int resi;

        //---------------------
        //Number of data on the time grid
        //---------------------
        resi = tSize;
        myfile.write((char*) &resi, sizeof(int));

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        resi = gSize;
        myfile.write((char*) &resi, sizeof(int));

        //---------------------
        //Loop
        //---------------------
        for(int nt = 0; nt <= tSize; nt++)
        {
            //Store the current time
            res = tGrid[nt];
            myfile.write((char*) &res, sizeof(double));

            //Store the data at current time
            for(int n1 = 0; n1 <= gSize; n1++)
            {
                for (int n2 = 0; n2 <= gSize; n2++)
                {
                    //NC state
                    for (int k = 0; k < 6; k++)
                    {
                        res = yNCE[k][nt][n1][n2];
                        myfile.write((char*) &res, sizeof(double));
                    }

                    //RCM state
                    for (int k = 0; k < 5; k++)
                    {
                        res = sNCE[k][nt][n1][n2];
                        myfile.write((char*) &res, sizeof(double));
                    }
                }
            }
        }
        myfile.close();
    }
    else return 0;

    return 1;
}

int readCU_bin(double ****yNCE, double ****sNCE, double *tGrid, int gSize, int tSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        int tSize0, gSize0;

        //---------------------
        //Number of data on the time grid
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        tSize0 = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        gSize0 = resi;


        if(tSize0 < tSize || gSize0 < gSize)
        {
            cout << "readManifold_bin: wrong inputs" << endl;
            cout << "tSize = " << tSize << ", but tSize0 = " << tSize0 << endl;
            cout << "gSize = " << gSize << ", but gSize0 = " << gSize0 << endl;
            return 0;
        }

        //---------------------
        //Loop
        //---------------------
        double res;
        for(int nt = 0; nt <= tSize; nt++)
        {
            //Read the current time
            myfile.read((char*) &res, sizeof(double));
            tGrid[nt]= res;

            //Read the data at current time
            for(int n1 = 0; n1 <= gSize; n1++)
            {
                for (int n2 = 0; n2 <= gSize; n2++)
                {
                    //NC state
                    for (int k = 0; k < 6; k++)
                    {
                        myfile.read((char*) &res, sizeof(double));
                        yNCE[k][nt][n1][n2] = res;
                    }

                    //RCM state
                    for (int k = 0; k < 5; k++)
                    {
                        myfile.read((char*) &res, sizeof(double));
                        sNCE[k][nt][n1][n2] = res;
                    }
                }
            }
        }
        myfile.close();
    }
    else return 0;

    return 1;
}

int getLenghtCU_bin(int *gSize, int *tSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        //---------------------
        //Number of data on the time grid
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        *tSize = resi;

        //---------------------
        //Number of data on the manifold grid
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        *gSize = resi;

        myfile.close();

    }
    else return 0;

    return 1;
}

int writeIntCU_bin(double *****y5Tens, double ****t4Tens, int gSize, int tSize, int mSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::out | ios::app); //append at the end of file
    if (myfile.is_open())
    {

        double res;
        //---------------------
        // Store all data
        //---------------------
        for(int kt = 0; kt <= tSize; kt++)
        {
            for(int ks1 = 0; ks1 <= gSize; ks1++)
            {
                for(int ks3 = 0; ks3 <= gSize; ks3++)
                {
                    for(int km = 0; km <= mSize; km++)
                    {
                        //NC time
                        res = t4Tens[kt][ks1][ks3][km];
                        myfile.write((char*) &res, sizeof(double));
                        //NC state
                        for (int k = 0; k < 6; k++)
                        {
                            res = y5Tens[k][kt][ks1][ks3][km];
                            myfile.write((char*) &res, sizeof(double));
                        }
                    }
                }
            }
        }
        myfile.close();

    }
    else return 0;

    return 1;
}

int initIntCU_bin(int gSize, int tSize, int mSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::out);
    if (myfile.is_open())
    {
        int resi;
        //---------------------
        //Number of data
        //---------------------
        resi = tSize;
        myfile.write((char*) &resi, sizeof(int));
        resi = gSize;
        myfile.write((char*) &resi, sizeof(int));
        resi = mSize;
        myfile.write((char*) &resi, sizeof(int));
    }
    else return 0;

    return 1;
}

int getLengthIntCU_bin(int *gSize, int *tSize, int *mSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        //---------------------
        //Number of data
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        *tSize = resi;
        myfile.read((char*) &resi, sizeof(int));
        *gSize = resi;
        myfile.read((char*) &resi, sizeof(int));
        *mSize = resi;
    }
    else return 0;

    return 1;
}

int readIntCU_bin(double *****y5Tens, double ****t4Tens, int gSize, int tSize, int mSize, int order, int type)
{
    //---------------------
    //Filename
    //---------------------
    string filename = filenameCUM(order, type);

    //---------------------
    //Open datafile
    //---------------------
    fstream myfile;
    myfile.open (filename.c_str(), ios::binary | ios::in);
    if (myfile.is_open())
    {
        int resi;
        int tSize0, gSize0, mSize0;

        //---------------------
        //Number of data
        //---------------------
        myfile.read((char*) &resi, sizeof(int));
        tSize0 = resi;
        myfile.read((char*) &resi, sizeof(int));
        gSize0 = resi;
        myfile.read((char*) &resi, sizeof(int));
        mSize0 = resi;


        if(tSize0 < tSize || gSize0 < gSize || mSize0 < mSize)
        {
            cout << "readManifold_bin: wrong inputs" << endl;
            cout << "tSize = " << tSize << ", but tSize0 = " << tSize0 << endl;
            cout << "gSize = " << gSize << ", but gSize0 = " << gSize0 << endl;
            cout << "mSize = " << mSize << ", but mSize0 = " << mSize0 << endl;
            return 0;
        }

        //---------------------
        //Loop
        //---------------------
        double res;
        for(int nt = 0; nt <= tSize; nt++)
        {
            //Read the data at current time
            for(int n1 = 0; n1 <= gSize; n1++)
            {
                for (int n2 = 0; n2 <= gSize; n2++)
                {
                    for (int nm = 0; nm <= mSize; nm++)
                    {
                        //Read the current time
                        myfile.read((char*) &res, sizeof(double));
                        t4Tens[nt][n1][n2][nm] = res;

                        //NC state
                        for (int k = 0; k < 6; k++)
                        {
                            myfile.read((char*) &res, sizeof(double));
                            y5Tens[k][nt][n1][n2][nm] = res;
                        }
                    }
                }
            }
        }
        myfile.close();
    }
    else return 0;

    return 1;
}


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Main routines (1): connections in the plane
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Initialize the grid on which the unstable manifold will be evaluated
 **/
void init_grid(double *grid, double gmin, double gmax, int gsize)
{
    double di, ds;
    for(int i = 0; i <= gsize; i++)
    {
        di = (double) i;
        ds = (double) gsize;
        if(gsize > 0) grid[i] = gmin +  (gmax - gmin)*di/ds;
        else grid[i] = gmin;
    }
}


void displayCompletion(double percent)
{
    //percent = (double) label/numberOfOrbits*100;
    std::cout << "\r" << percent << "% completed: ";
    std::cout << std::string(floor(0.1*percent), '|') << endl;
    std::cout.flush();
}


/**
 *  \brief Computes the unstable directions of a discrete set of points on an orbit stored in a data file (obtained from gridOrbit).
 **/
int cusMan(int order,
           double epsilon,
           double tMin,
           double tMax,
           int tSize,
           double gMin,
           double gMax,
           int gSize,
           vector<Oftsc> &CMh,
           matrix<Ofsc>  &Mcoc,
           matrix<Ofsc>  &MIcoc,
           vector<Ofsc>  &Vcoc,
           bool isPar)

{
    //------------------------------------------
    //Building the working grid
    //------------------------------------------
    double *grid = dvector(0,  gSize);
    init_grid(grid, gMin, gMax, gSize);

    //------------------------------------------
    //Building the time grid
    //------------------------------------------
    double *tgrid = dvector(0,  tSize);
    init_grid(tgrid, tMin, tMax, tSize);

    //------------------------------------------
    //Number of elements
    //------------------------------------------
    int noe = pow(1+gSize, 2.0)*(1+tSize);
    int iter = 1;

    //------------------------------------------
    // Loop on all elements
    //------------------------------------------
    double ****yNCEU = d4tensor(0, 5, 0, tSize, 0, gSize, 0, gSize);
    double ****sNCEU = d4tensor(0, 4, 0, tSize, 0, gSize, 0, gSize);

    //#pragma omp parallel for if(isPar) shared(iter)
    for(int kt = 0; kt <= tSize; kt++)
    {
        //#pragma omp parallel for if(isPar) shared(iter)
        for(int ks1 = 0; ks1 <= gSize; ks1++)
        {
            #pragma omp parallel for if(isPar)  shared(iter)
            for(int ks3 = 0; ks3 <= gSize; ks3++)
            {
                Ofsc ofs(OFS_ORDER);
                double *yvu = dvector(0,5);
                double *sti = dvector(0,4);

                //----------------------
                // Initialization on the center-unstable manifold
                //----------------------
                //Init sti
                sti[0] = grid[ks1];
                sti[1] = 0.0;
                sti[2] = grid[ks3];
                sti[3] = 0.0;
                sti[4] = epsilon;

                //Equivalent state
                RCMtoNCbyTFC(sti, tgrid[kt], SEML.us.n, OFTS_ORDER, OFS_ORDER, CMh, ofs, Mcoc, Vcoc, yvu, 1);

                #pragma omp critical
                {
                    //Save
                    yNCEU[0][kt][ks1][ks3] = yvu[0];
                    yNCEU[1][kt][ks1][ks3] = yvu[1];
                    yNCEU[2][kt][ks1][ks3] = yvu[2];
                    yNCEU[3][kt][ks1][ks3] = yvu[3];
                    yNCEU[4][kt][ks1][ks3] = yvu[4];
                    yNCEU[5][kt][ks1][ks3] = yvu[5];

                    //Save
                    sNCEU[0][kt][ks1][ks3] = sti[0];
                    sNCEU[1][kt][ks1][ks3] = sti[1];
                    sNCEU[2][kt][ks1][ks3] = sti[2];
                    sNCEU[3][kt][ks1][ks3] = sti[3];
                    sNCEU[4][kt][ks1][ks3] = sti[4];

                    //Display
                    displayCompletion((double) iter++/noe*100);
                }

                free_dvector(yvu, 0, 5);
                free_dvector(sti, 0, 4);
            }
        }
    }

    //------------------------------------------
    //Store values
    //------------------------------------------
    writeCU_bin(yNCEU, sNCEU, tgrid, gSize, tSize, order, TYPE_CU);

    //------------------------------------------
    //Free
    //------------------------------------------
    free_d4tensor(yNCEU, 0, 5, 0, tSize, 0, gSize, 0, gSize);
    free_d4tensor(sNCEU, 0, 4, 0, tSize, 0, gSize, 0, gSize);

    return 1;
}

/**
 *  \brief Computes the manifold branches from a discrete set of unstable directions obtained from cusMan
 **/
int intMan(double tman,
           int order,
           int tSizex,
           int gSizex,
           int mSize,
           int isPar,
           vector<Oftsc> &CM,
           vector<Oftsc> &CMh,
           matrix<Ofsc>  &Mcoc,
           matrix<Ofsc>  &Pcoc,
           matrix<Ofsc>  &MIcoc,
           matrix<Ofsc>  &PIcoc,
           vector<Ofsc>  &Vcoc)
{
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int tSize, gSize;
    getLenghtCU_bin(&gSize, &tSize, order, TYPE_CU);

    cout << "intMan. tSize = " << tSize << endl;
    cout << "intMan. gSize = " << gSize << endl;

    //----------------------------------------------------------
    //Initialize the final size comparing tSize/tSizex
    //----------------------------------------------------------
    int tSizet = tSizex < 0 ? tSize:tSizex;
    int gSizet = gSizex < 0 ? gSize:gSizex;

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double ****yNCEU = d4tensor(0, 5, 0, tSize, 0, gSize, 0, gSize);
    double ****sNCEU = d4tensor(0, 4, 0, tSize, 0, gSize, 0, gSize);
    double *tGrid    = dvector(0, tSize);

    //----------------------------------------------------------
    //Read data
    //----------------------------------------------------------
    int status = readCU_bin(yNCEU, sNCEU, tGrid, gSize, tSize, order, TYPE_CU);

    //----------------------------------------------------------
    // Loop on all elements
    //----------------------------------------------------------
    if(status)
    {
        //----------------------------------------------------------
        // Initialize the final data file
        //----------------------------------------------------------
        initIntCU_bin(gSizet, tSizet, mSize, order, TYPE_MAN);

        //----------------------------------------------------------
        //To store all data
        //----------------------------------------------------------
        double *****y5Tens = d5tensor(0, 5, 0, tSize, 0, gSize, 0, gSize, 0, mSize);
        double ****t4Tens  = d4tensor(0, tSize, 0, gSize, 0, gSize, 0, mSize);

        //----------------------------------------------------------
        //The default interval of projection is set to Tproj = T/5,
        // although it is not used here
        //----------------------------------------------------------
        double tproj = SEML.us.T/5.0;

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        //gnuplot_ctrl  *h1;
        //h1 = gnuplot_init();

        int index = 0;
        //#pragma omp parallel for if(isPar) shared(index)
        for(int kt = 0; kt <= tSizet; kt++)
        {
            //#pragma omp parallel for if(isPar) shared(index)
            for(int ks1 = 0; ks1 <= gSizet; ks1++)
            {
                #pragma omp parallel for if(isPar) shared(index)
                for(int ks3 = 0; ks3 <= gSizet; ks3++)
                {
                    //---------------------------------------------------------------------
                    //Temp variables
                    //---------------------------------------------------------------------
                    double yv[6], st[5], tv;
                    tv    = tGrid[kt];

                    yv[0] = yNCEU[0][kt][ks1][ks3];
                    yv[1] = yNCEU[1][kt][ks1][ks3];
                    yv[2] = yNCEU[2][kt][ks1][ks3];
                    yv[3] = yNCEU[3][kt][ks1][ks3];
                    yv[4] = yNCEU[4][kt][ks1][ks3];
                    yv[5] = yNCEU[5][kt][ks1][ks3];

                    st[0] = sNCEU[0][kt][ks1][ks3];
                    st[1] = sNCEU[1][kt][ks1][ks3];
                    st[2] = sNCEU[2][kt][ks1][ks3];
                    st[3] = sNCEU[3][kt][ks1][ks3];
                    st[4] = sNCEU[4][kt][ks1][ks3];

                    //---------------------------------------------------------------------
                    //Local variables
                    //---------------------------------------------------------------------
                    double **yManNCEM  = dmatrix(0, 5, 0, mSize);
                    double **yManNCSEM = dmatrix(0, 5, 0, mSize);
                    double *tManEM     = dvector(0, mSize);

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
                    init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 6, 1e-6, qbfbp_vfn_novar, NULL, &SEML);

                    //------------------
                    //Orbit structure
                    //------------------
                    Ofsc orbit_ofs(OFS_ORDER);
                    SingleOrbit orbit;
                    //Init routine
                    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, tv, tv+tman, tproj, &driver, &SEML);

                    //---------------------------------------------------------------------
                    // Update the state
                    //---------------------------------------------------------------------
                    orbit_update_ic(orbit, st, yv, tv);

                    //---------------------------------------------------------------------
                    //Integration
                    //---------------------------------------------------------------------
                    trajectory_integration_grid(orbit, tv, tv+tman, yManNCEM, tManEM, mSize, 0);

                    //---------------------------------------------------------------------
                    //Storage
                    //---------------------------------------------------------------------
                    #pragma omp critical
                    {
                        //To NCSEM coordinates
                        NCEMmtoNCSEMm_vec(yManNCEM, tManEM, yManNCSEM, mSize, &SEML);
                        for(int km = 0; km <= mSize; km++)
                        {
                            for(int i = 0; i < 6; i++) y5Tens[i][kt][ks1][ks3][km] = yManNCSEM[i][km];
                            t4Tens[kt][ks1][ks3][km] = tManEM[km]*SEML.us_em.ns; //time is also given in SEM units!
                        }

                        //----------------------------------------------------------
                        //Gnuplot window
                        //----------------------------------------------------------
                        //if(yManNCSEM[2][0] == 0) gnuplot_plot_xy(h1, yManNCSEM[0], yManNCSEM[1], mSize+1, (char*)"", "lines", "1", "2", 4);
                        //else gnuplot_plot_xyz(h1, yManNCSEM[0], yManNCSEM[1], yManNCSEM[2], mSize+1, (char*)"", "lines", "1", "2", 4);
                        cout << "intMan : " << 100.0*index++/((1+gSizet)*(1+gSizet)*(1+tSizet)) << endl;
                    }

                    //---------------------------------------------------------------------
                    //Free
                    //---------------------------------------------------------------------
                    free_dmatrix(yManNCEM, 0, 5, 0, mSize);
                    free_dmatrix(yManNCSEM, 0, 5, 0, mSize);
                    free_dvector(tManEM, 0, mSize);


                }
            }
        }

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        //char ch;
        //printf("Press ENTER to go on\n");
        //scanf("%c",&ch);
        //gnuplot_close(h1);


        //---------------------------------------------------------------------
        //Storage
        //---------------------------------------------------------------------
        writeIntCU_bin(y5Tens, t4Tens, gSizet, tSizet, mSize, order, TYPE_MAN);

        //----------------------------------------------------------
        //Free
        //----------------------------------------------------------
        free_d5tensor(y5Tens, 0, 5, 0, tSize, 0, gSize, 0, gSize, 0, mSize);
        free_d4tensor(t4Tens, 0, tSize, 0, gSize, 0, gSize, 0, mSize);
    }
    else
    {
        cout << "intMan: error during data reading." << endl;
        return 0;
    }

    //---------------------------------------------------------------------
    //Free
    //---------------------------------------------------------------------
    free_d4tensor(yNCEU, 0, 5, 0, tSize, 0, gSize, 0, gSize);
    free_d4tensor(sNCEU, 0, 4, 0, tSize, 0, gSize, 0, gSize);
    free_dvector(tGrid, 0, tSize);


    return 1;
}

/**
 *  \brief Changing the scope of the computation:
 *       1. Set OFTS_ORDER=order, REDUCED_NV = rnv.
 *       2. Focus on the coordinate system defined by focus (F_EM, F_SEM...).
 **/
void changeScope(int order, int rnv, int focus)
{
    cout << "changeScope: changing the default parameters..." << endl;
    REDUCED_NV=rnv;
    OFTS_ORDER=order;
    cout << "changeScope: changing focus..." << endl;
    changeDCS(SEML, focus);
}


void projMan(int order_em,
             int order_sem,
             int tSizex,
             int gSizex,
             int mSizex,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &Pcoc,
             matrix<Ofsc>  &MIcoc,
             matrix<Ofsc>  &PIcoc,
             vector<Ofsc>  &Vcoc,
             int isPar)
{
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int gSize, tSize, mSize;
    getLengthIntCU_bin(&gSize, &tSize, &mSize, order_em, TYPE_MAN);

    //----------------------------------------------------------
    //Initialize the final size comparing tSize/tSizex
    //----------------------------------------------------------
    int tSizet = tSizex < 0 ? tSize:tSizex;
    int gSizet = gSizex < 0 ? gSize:gSizex;
    int mSizet = mSizex < 0 ? mSize:mSizex;

    //----------------------------------------------------------
    //To store all data from the IC
    //----------------------------------------------------------
    double ****yNCEU = d4tensor(0, 5, 0, tSize, 0, gSize, 0, gSize);
    double ****sNCEU = d4tensor(0, 4, 0, tSize, 0, gSize, 0, gSize);
    double *tGrid    = dvector(0, tSize);

    //----------------------------------------------------------
    //Read data from the IC
    //----------------------------------------------------------
    int status = readCU_bin(yNCEU, sNCEU, tGrid, gSize, tSize, order_em, TYPE_CU);

    //----------------------------------------------------------
    //To store all data from the Manifold
    //----------------------------------------------------------
    double *****y5Tens = d5tensor(0, 5, 0, tSize, 0, gSize, 0, gSize, 0, mSize);
    double ****t4Tens  = d4tensor(0, tSize, 0, gSize, 0, gSize, 0, mSize);

    //----------------------------------------------------------
    //To sort data
    //----------------------------------------------------------
    int Nsort = (gSize+1)*(gSize+1)*(tSize+1);
    vector<int>    indexMin(Nsort);
    vector<int>    ks1Min(Nsort);
    vector<int>    ks3Min(Nsort);
    vector<int>    ktMin(Nsort);
    vector<double> distMin(Nsort);

    //----------------------------------------------------------
    //Read data from the Manifold
    //----------------------------------------------------------
    status = readIntCU_bin(y5Tens, t4Tens, gSize, tSize, mSize, order_em, TYPE_MAN);


    //----------------------------------------------------------
    // Loop on all elements
    //----------------------------------------------------------
    if(status)
    {
        //-------------------------------------------------------------------------------
        //Changing the scope to CM SEM...
        //-------------------------------------------------------------------------------
        //changeScope(order_sem, 5, F_SEM); //if CSM is used
        changeScope(order_sem, 4, F_SEM);  //if CM is used

        //--------------------------------------
        // Center- manifold
        //--------------------------------------
        vector<Oftsc>  CMc(6);     ///center manifold in NC coordinates
        vector<Oftsc> CMhc(6);     ///center manifold in TFC coordinates

        //------------------------------------------
        // Update of the central manifold
        //------------------------------------------
        //Update CM
        readVOFTS_bin(CMc,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
        readVOFTS_bin(CMhc, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
        //Update COC
        initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

        //----------------------------------------------------------
        //To store
        //----------------------------------------------------------
        double ****yfSEMin    = d4tensor(0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        double ****yprojSEMin = d4tensor(0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        double ****sprojSEMin = d4tensor(0, 3, 0, tSizet, 0, gSizet, 0, gSizet);
        double ****y0SE       = d4tensor(0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        double ****y0EM       = d4tensor(0, 5, 0, tSizet, 0, gSizet, 0, gSizet);

        //----------------------------------------------------------
        //To save
        //----------------------------------------------------------
        double ***ePmin4Tensor  = d3tensor(0, tSize, 0, gSize, 0, gSize);
        double ***DV4Tensor     = d3tensor(0, tSize, 0, gSize, 0, gSize);
        double ***tvMinTensor   = d3tensor(0, tSize, 0, gSize, 0, gSize);


        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        //gnuplot_ctrl  *h1;
        //h1 = gnuplot_init();
        //gnuplot_ctrl *h2;
        //h2 = gnuplot_init();

        //----------------------------------------------------------
        //Notable points in SEM system
        //----------------------------------------------------------
        double **semP = dmatrix(0, 6, 0, 2);
        semPoints(0.0, semP);
        //semPlot(h2, semP);

        //---------------------------------------------------------------------
        //SuperLoop on all positions
        //---------------------------------------------------------------------
        int index = 0;
        //#pragma omp parallel for if(isPar) shared(index)
        for(int kt = 0; kt <= tSizet; kt++)
        {
            //#pragma omp parallel for if(isPar) shared(index)
            for(int ks1 = 0; ks1 <= gSizet; ks1++)
            {
                #pragma omp parallel for if(isPar) shared(index)
                for(int ks3 = 0; ks3 <= gSizet; ks3++)
                {

                    //----------------------------------------------------------
                    //Temp variable
                    //----------------------------------------------------------
                    Ofsc ofs(OFS_ORDER);
                    double yv[6], yvproj[6], sproj[4], tv;
                    double yvSE[6], yvprojSE[6], yvEM[6];
                    double yvSEVel[6], yvprojSEVel[6];
                    double ePm, ePmin = 1e6, DV = 0.0;
                    double yvnorm = 0.0;
                    int ksort = (gSizet+1)*(gSizet+1)*kt + (gSizet+1)*ks1 + ks3;
                    int kmin = 0;

                    //----------------------------------------------------------
                    //Loop on trajectory
                    //----------------------------------------------------------
                    for(int kman = 0; kman <= mSizet; kman++)
                    {
                        // Current state
                        for(int i = 0; i < 6; i++) yv[i] = y5Tens[i][kt][ks1][ks3][kman];
                        tv = t4Tens[kt][ks1][ks3][kman];

                        yvnorm = 0.0;
                        for(int i = 0; i < 2; i++) yvnorm += yv[i]*yv[i];
                        yvnorm = sqrt(yvnorm);

                        //If close enough to SEMLi
                        if(yvnorm < 0.5)
                        {
                            // Projection on the center manifold
                            NCprojCCMtoCM(yv, tv, sproj, CMhc, MIcoc, Vcoc);

                            //If close enough the SEMLi
                            if(sqrt(sproj[0]*sproj[0]+sproj[2]*sproj[2])< 0.5)
                            {
                                //NCprojCCMtoCUS(yv, tv, sproj, CMhc, 1e-6, MIcoc, Vcoc);
                                //yvproj = W(sproj, tv)
                                RCMtoNCbyTFC(sproj, tv, SEML.us.n, OFTS_ORDER, OFS_ORDER, CMhc, ofs, Mcoc, Vcoc, yvproj, 1);
                                //Back to SEM coordinates
                                NCtoSEM(tv, yv, yvSE, &SEML);
                                NCtoSEM(tv, yvproj, yvprojSE, &SEML);
                                //Back to SEM coordinates with velocity
                                SEMmtoSEMv(tv, yvSE, yvSEVel, &SEML);
                                SEMmtoSEMv(tv, yvprojSE, yvprojSEVel, &SEML);
                                //Distance of projection in SEM coordinates
                                ePm = 0.0;
                                for(int i = 0; i < 3; i++) ePm += (yvprojSE[i] - yvSE[i])*(yvprojSE[i] - yvSE[i]);
                                ePm = sqrt(ePm);
                            }
                            else ePm = 1e5;
                        }
                        else ePm = 1e5;

                        //Update distance min if necessary
                        if(ePm < ePmin)
                        {
                            ePmin = ePm;
                            kmin = kman;
                            for(int i = 0; i < 6; i++) yfSEMin[i][kt][ks1][ks3]    = yvSEVel[i];
                            for(int i = 0; i < 6; i++) yprojSEMin[i][kt][ks1][ks3] = yvprojSEVel[i];
                            for(int i = 0; i < 4; i++) sprojSEMin[i][kt][ks1][ks3] = sproj[i];

                            //Associated DV
                            DV = 0.0;
                            for(int i = 3; i < 6; i++) DV += (yvprojSEVel[i] - yvSEVel[i])*(yvprojSEVel[i] - yvSEVel[i]);
                            DV = sqrt(DV);
                        }

                    }

                    #pragma omp critical
                    {
                        //----------------------------------------------------------
                        //Save
                        //----------------------------------------------------------
                        ePmin4Tensor[kt][ks1][ks3] = ePmin;
                        DV4Tensor[kt][ks1][ks3] = DV;
                        tvMinTensor[kt][ks1][ks3] = t4Tens[kt][ks1][ks3][kmin];
                        //Initial position in SE coordinates
                        for(int i = 0; i < 6; i++) yv[i] = y5Tens[i][kt][ks1][ks3][0];
                        tv = t4Tens[kt][ks1][ks3][0];
                        NCtoSEM(tv, yv, yvSE, &SEML);
                        for(int i = 0; i < 6; i++) y0SE[i][kt][ks1][ks3] = yvSE[i];
                        //Same in NCEM coordinates
                        NCSEMmtoNCEMm(tv, yv, yvEM, &SEML);
                        for(int i = 0; i < 6; i++) y0EM[i][kt][ks1][ks3] = yvEM[i];
                        //For sorting
                        indexMin[ksort] = (int) kmin;
                        ks1Min[ksort]   = (int) ks1;
                        ks3Min[ksort]   = (int) ks3;
                        ktMin[ksort]    = (int) kt;
                        distMin[ksort]  = ePmin;

                        //----------------------------------------------------------
                        //Gnuplot window
                        //----------------------------------------------------------
                        //if(ePmin < 1.0) gnuplot_plot_xyz(h1, &sNCEU[0][kt][ks1][ks3], &sNCEU[2][kt][ks1][ks3], &ePmin, 1, (char*)"", "points", "1", "2", 8);
                        cout << "projMan : " << 100.0*index++/((1+gSizet)*(1+gSizet)*(1+tSizet)) << endl;
                    }
                }
            }
        }


        //------------------------------------------------------------------------------------------------------
        //Sorting the best solutions
        //------------------------------------------------------------------------------------------------------
        vector<size_t> sortId = sort_indexes(distMin);

        //----------------------------------------------------------
        //Saving the 50 best results or all results if less than 50 have been computed
        //----------------------------------------------------------
        int Nt  = min(50, Nsort-2);
        //To store final data
        double ***yTR = d3tensor(0, 5, 0, Nt, 0, mSizet+2);
        double **tTR  = dmatrix(0, Nt, 0, mSizet+2);

        //----------------------------------------------------------
        //Display results
        //----------------------------------------------------------
        double **yManNCSEM  = dmatrix(0, 5, 0, mSizet);
        double **yManSEM    = dmatrix(0, 5, 0, mSizet);
        double *tManSEM     = dvector(0, mSizet);

        int ksortpos, ks1pos, ks3pos, kmpos, ktpos;
        //double temp;
        for(int kpos = 0; kpos <= Nt; kpos++)
        {
            ksortpos = sortId[kpos];
            ks1pos   = ks1Min[ksortpos];
            ks3pos   = ks3Min[ksortpos];
            ktpos    = ktMin[ksortpos];
            kmpos    = indexMin[ksortpos];

            //----------------------------------------------------------
            //Plot the solution +  the good point
            //----------------------------------------------------------
            for(int kman = 0; kman <= mSizet; kman++)
            {
                // Current state
                for(int i = 0; i < 6; i++) yManNCSEM[i][kman] = y5Tens[i][ktpos][ks1pos][ks3pos][kman];
                tManSEM[kman] = t4Tens[ktpos][ks1pos][ks3pos][kman];


            }

            //To SEM coordinates
            NCtoSEM_vec(yManNCSEM, tManSEM, yManSEM, mSizet, &SEML);

            //Store for future save
            for(int kman = 0; kman <= mSizet; kman++)
            {
                // Current state
                for(int i = 0; i < 6; i++) yTR[i][kpos][kman] = yManSEM[i][kman];
                tTR[kpos][kman] = tManSEM[kman];
            }

            //----------------------------------------------------------
            // Update final data
            // Careful: we save the final state at the end of yTR
            // Careful: we save the projection state at the end of yTR
            // Careful: we save the projection time at the end of tTR
            //----------------------------------------------------------
            for(int i = 0; i < 6; i++) yTR[i][kpos][mSizet+1] = yfSEMin[i][ktpos][ks1pos][ks3pos];
            for(int i = 0; i < 6; i++) yTR[i][kpos][mSizet+2] = yprojSEMin[i][ktpos][ks1pos][ks3pos];
            tTR[kpos][mSizet+1] = tManSEM[kmpos];
            tTR[kpos][mSizet+2] = tManSEM[kmpos];


            //plot
            //gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1],  yManSEM[2], mSizet+1, (char*)"", "lines", "1", "3", 4);
            //plot the right solution
            //gnuplot_plot_xyz(h2, &yfSEMin[0][ktpos][ks1pos][ks3pos], &yfSEMin[1][ktpos][ks1pos][ks3pos], &yfSEMin[2][ktpos][ks1pos][ks3pos],  1, (char*)"", "points", "1", "3", 2);
            //gnuplot_plot_xyz(h2, &yprojSEMin[0][ktpos][ks1pos][ks3pos], &yprojSEMin[1][ktpos][ks1pos][ks3pos], &yprojSEMin[2][ktpos][ks1pos][ks3pos], 1, (char*)"", "points", "1", "3", 8);

        }

        //---------------------------------------------------------------------
        // Free sorted version
        //---------------------------------------------------------------------
        free_dmatrix(yManNCSEM, 0, 5, 0, mSizet);
        free_dmatrix(yManSEM, 0, 5, 0, mSizet);
        free_dvector(tManSEM, 0, mSizet);


        //---------------------------------------------------------------------
        // Save sorted version
        //---------------------------------------------------------------------
        //---------------------
        //Filename
        //---------------------
        string filename = filenameCUM(order_em, TYPE_MAN_SORT);

        //---------------------
        //Open datafile
        //---------------------
        fstream myfile;
        myfile.open (filename.c_str(), ios::binary | ios::out);
        if (myfile.is_open())
        {
            double res;
            //---------------------
            //Loop
            //---------------------
            for(int kpos = 0; kpos <= Nt; kpos++)
            {
                //Get sorted indices
                ksortpos = sortId[kpos];
                ks1pos   = ks1Min[ksortpos];
                ks3pos   = ks3Min[ksortpos];
                ktpos    = ktMin[ksortpos];


                //Current homogeneous polynomial
                for(int kman = 0; kman <= mSizet+2; kman++)
                {
                    //1. label
                    res = kpos;
                    myfile.write((char*) &res, sizeof(double));

                    //2. The current time
                    res = tTR[kpos][kman];
                    myfile.write((char*) &res, sizeof(double));

                    //3-8. The SE state
                    for (int k = 0; k < 6; k++)
                    {
                        res = yTR[k][kpos][kman];
                        myfile.write((char*) &res, sizeof(double));
                    }

                    //9. s1 (EM)
                    res = sNCEU[0][ktpos][ks1pos][ks3pos];
                    myfile.write((char*) &res, sizeof(double));

                    //10. s3 (EM)
                    res = sNCEU[2][ktpos][ks1pos][ks3pos];
                    myfile.write((char*) &res, sizeof(double));


                    //11. s1 (SEM)
                    res = sprojSEMin[0][ktpos][ks1pos][ks3pos];
                    myfile.write((char*) &res, sizeof(double));

                    //12. s3 (SEM)
                    res = sprojSEMin[2][ktpos][ks1pos][ks3pos];
                    myfile.write((char*) &res, sizeof(double));

                    //13. ePmin (1)
                    res = ePmin4Tensor[ktpos][ks1pos][ks3pos];
                    myfile.write((char*) &res, sizeof(double));

                    //14. ePmin (2)
                    res = distMin[ksortpos];
                    myfile.write((char*) &res, sizeof(double));
                }
            }
            myfile.close();
        }


        //------------------------------------------------------------------------------------------------------
        //Storing all solutions
        //------------------------------------------------------------------------------------------------------
        //---------------------
        //Filename
        //---------------------
        filename = filenameCUM(order_em, TYPE_MAN_PROJ);

        //---------------------
        //Open datafile
        //---------------------
        myfile.open (filename.c_str(), ios::binary | ios::out);
        if (myfile.is_open())
        {
            double res;
            //---------------------
            // Store all data
            //---------------------
            for(int kt = 0; kt <= tSize; kt++)
            {
                for(int ks1 = 0; ks1 <= gSize; ks1++)
                {
                    for(int ks3 = 0; ks3 <= gSize; ks3++)
                    {
                        // 1. NC time
                        res = tGrid[kt];
                        myfile.write((char*) &res, sizeof(double));

                        // 2-7. yNCEU state in NCEM coordinates
                        for (int k = 0; k < 6; k++)
                        {
                            res = yNCEU[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 8-13. yNCEU state in NCEM coordinates again
                        for (int k = 0; k < 6; k++)
                        {
                            res = y0EM[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 14-19. yNCEU state in SE coordinates
                        for (int k = 0; k < 6; k++)
                        {
                            res = y0SE[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 20-24. sNCEU state
                        for (int k = 0; k < 5; k++)
                        {
                            res = sNCEU[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 25. ePmin4Tensor
                        res = ePmin4Tensor[kt][ks1][ks3];
                        myfile.write((char*) &res, sizeof(double));

                        // 26. DV4Tensor
                        res = DV4Tensor[kt][ks1][ks3];
                        myfile.write((char*) &res, sizeof(double));

                        // 27. tvMinTensor
                        res = tvMinTensor[kt][ks1][ks3];
                        myfile.write((char*) &res, sizeof(double));

                        // 28-33. yfSEMin state in SE coordinates
                        for (int k = 0; k < 6; k++)
                        {
                            res = yfSEMin[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 34-39. yprojSEMin state in SE coordinates
                        for (int k = 0; k < 6; k++)
                        {
                            res = yprojSEMin[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }

                        // 40-43. sprojSEMin state in SE coordinates
                        for (int k = 0; k < 4; k++)
                        {
                            res = sprojSEMin[k][kt][ks1][ks3];
                            myfile.write((char*) &res, sizeof(double));
                        }
                    }
                }
            }
            myfile.close();
        }


        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        //char ch;
        //printf("Press ENTER to go on\n");
        //scanf("%c",&ch);
        //gnuplot_close(h1);
        //gnuplot_close(h2);

        //----------------------------------------------------------
        //Free
        //----------------------------------------------------------
        free_d4tensor(yfSEMin, 0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        free_d4tensor(yprojSEMin, 0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        free_d4tensor(sprojSEMin, 0, 3, 0, tSizet, 0, gSizet, 0, gSizet);
        free_d4tensor(y0SE, 0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        free_d4tensor(y0EM, 0, 5, 0, tSizet, 0, gSizet, 0, gSizet);
        free_d3tensor(ePmin4Tensor, 0, tSize, 0, gSize, 0, gSize);
        free_d3tensor(DV4Tensor, 0, tSize, 0, gSize, 0, gSize);
        free_d3tensor(tvMinTensor, 0, tSize, 0, gSize, 0, gSize);
    }
    else
    {
        cout << "projMan: error during data reading." << endl;
    }



    //----------------------------------------------------------
    //Free
    //----------------------------------------------------------
    free_d5tensor(y5Tens, 0, 5, 0, tSize, 0, gSize, 0, gSize, 0, mSize);
    free_d4tensor(t4Tens, 0, tSize, 0, gSize, 0, gSize, 0, mSize);

    cout << "projMan: done." << endl;

}

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Main routines (2): connections with a fixed orbit
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initial conditions st0, on the center manifold CM.
 **/
int gridOrbit(double st0[],
              double t0,
              double tf,
              double dt,
              vector<Oftsc> &CM,
              vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
int cusOrbit(int order,
             int sizeOrbit,
             double epsilon,
             vector<Oftsc> &CMh,
             matrix<Ofsc>  &Mcoc,
             matrix<Ofsc>  &MIcoc,
             vector<Ofsc>  &Vcoc,
             bool isPar)
{
    //----------------------------------------------------------
    // Get the number of data lines from file
    //----------------------------------------------------------
    int N = getLineNumber(order, sizeOrbit, TYPE_ORBIT)-1;

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
        int status    = readOrbit(tNCE, yNCE, order, sizeOrbit, N+1, TYPE_ORBIT);

        //----------------------
        // Loop on all elements
        //----------------------
        if(status)
        {
            double **yNCEU = dmatrix(0, 5, 0, N);
            double **sNCEU = dmatrix(0, 4, 0, N);

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
                NCprojCCMtoCUS(yv, tv, sti, CMh, epsilon, MIcoc, Vcoc);
                //Equivalent state
                RCMtoNCbyTFC(sti, tv, SEML.us.n, OFTS_ORDER, OFS_ORDER, CMh, ofs, Mcoc, Vcoc, yvu, 1);

                #pragma omp critical
                {
                    //Save
                    yNCEU[0][k] = yvu[0];
                    yNCEU[1][k] = yvu[1];
                    yNCEU[2][k] = yvu[2];
                    yNCEU[3][k] = yvu[3];
                    yNCEU[4][k] = yvu[4];
                    yNCEU[5][k] = yvu[5];

                    //Save
                    sNCEU[0][k] = sti[0];
                    sNCEU[1][k] = sti[1];
                    sNCEU[2][k] = sti[2];
                    sNCEU[3][k] = sti[3];
                    sNCEU[4][k] = sti[4];
                }

                free_dvector(yv, 0, 5);
                free_dvector(yvu, 0, 5);
                free_dvector(sti, 0, 4);
            }
            //Store values
            writeOrbit(tNCE, yNCEU, sNCEU, order, sizeOrbit, N+1, TYPE_CU);
            //Free
            free_dmatrix(yNCEU, 0, 5, 0, N);
            free_dmatrix(sNCEU, 0, 4, 0, N);
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
int manOrbit(double tman,
             int order,
             int sizeOrbit,
             int Nt,
             int Nman,
             int isPar,
             vector<Oftsc> &CM,
             vector<Oftsc> &CMh,
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
    int N0 = getLineNumber(order, sizeOrbit, TYPE_CU)-1;

    //----------------------------------------------------------
    //If the file exists and contains data, go on
    //----------------------------------------------------------
    if(N0 > 0)
    {
        if(Nt <= 0) N = N0;
        else N = Nt;

        if(N > N0)
        {
            cout << "manOrbit: wrong input, Nt = " << Nt << ", but N0 = " << N0 << endl;
            return 0;
        }

        //----------------------
        // Read data from file
        //----------------------
        double **yNCU = dmatrix(0, 5, 0, N);
        double **sNCU = dmatrix(0, 4, 0, N);
        double *tNCU  = dvector(0, N);
        int status    = readOrbit(tNCU, yNCU, sNCU, order, sizeOrbit, N+1, TYPE_CU);

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
                double **yManNCEM  = dmatrix(0, 5, 0, Nman);
                double **yManNCSEM = dmatrix(0, 5, 0, Nman);
                double *tManEM     = dvector(0, Nman);

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
                init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 6, 1e-6, qbfbp_vfn_novar, NULL, &SEML);

                //------------------
                //Orbit structure
                //------------------
                Ofsc orbit_ofs(OFS_ORDER);
                SingleOrbit orbit;
                //Init routine
                init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, tv, tv+tman, tproj, &driver, &SEML);

                //---------------------------------------------------------------------
                // Update the state
                //---------------------------------------------------------------------
                orbit_update_ic(orbit, st, yv, tv);

                //---------------------------------------------------------------------
                //Integration
                //---------------------------------------------------------------------
                trajectory_integration_grid(orbit, tv, tv+tman, yManNCEM, tManEM, Nman, 0);

                //                //---------------------------------------------------------------------
                //                //Integration check
                //                //---------------------------------------------------------------------
                //                reset_ode_structure(orbit.driver);
                //                double t0 = tManEM[0];
                //                double yv0[6];
                //                for(int i = 0; i < 6; i++) yv0[i] = yv[i];
                //                int kk = Nman;
                //
                //                cout << "Time 0" << endl;
                //                cout << tManEM[0] << "  " << t0 << "   " << endl;
                //
                //                cout << "Init" << endl;
                //                for(int i = 0; i < 6; i++) cout << yManNCEM[i][0] << endl;
                //
                //                //Integration
                //                gsl_odeiv2_driver_apply (driver.d, &t0, tManEM[kk], yv);
                //                cout << "Time" << endl;
                //                cout << tManEM[kk] << "  " << t0 << "   " << endl;
                //                cout << "Result" << endl;
                //                for(int i = 0; i <6; i++) cout << (yManNCEM[i][kk] - yv[i])/max(1.0, fabs(yv[i])) << endl;
                //
                //
                //                t0 = tManEM[0];
                //                reset_ode_structure(&driver);
                //                gsl_odeiv2_driver_apply (driver.d, &t0, tManEM[kk], yv0);
                //                cout << "Result 2: " << endl;
                //                for(int i = 0; i < 6; i++) cout << (yv[i] - yv0[i])/max(1.0, fabs(yv[i])) << endl;


                //---------------------------------------------------------------------
                //Storage
                //---------------------------------------------------------------------
                #pragma omp critical
                {
                    //To NCSEM coordinates
                    NCEMmtoNCSEMm_vec(yManNCEM, tManEM, yManNCSEM, Nman, &SEML);

                    //                    //----------------------------------------------------------
                    //                    //Gnuplot window
                    //                    //----------------------------------------------------------
                    //                    gnuplot_ctrl  *h1;
                    //                    h1 = gnuplot_init();
                    //                    gnuplot_plot_xyz(h1, yManNCSEM[0], yManNCSEM[1], yManNCSEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
                    //                    char ch;
                    //                    printf("Press ENTER to go on\n");
                    //                    scanf("%c",&ch);
                    //                    gnuplot_close(h1);

                    for(int p = 0; p <= Nman; p++)
                    {
                        for(int i = 0; i < 6; i++) yTens[i][k][p] = yManNCSEM[i][p];
                        tTens[k][p] = tManEM[p]*SEML.us_em.ns; //time is also given in SEM units!
                    }
                }

                //---------------------------------------------------------------------
                //Free
                //---------------------------------------------------------------------
                free_dmatrix(yManNCEM, 0, 5, 0, Nman);
                free_dmatrix(yManNCSEM, 0, 5, 0, Nman);
                free_dvector(tManEM, 0, Nman);

                cout << "manOrbit : " << 100.0*index++/N << endl;
            }


            //---------------------------------------------------------------------
            //Save
            //---------------------------------------------------------------------
            writeManifold_bin(tTens, yTens, order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

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
void manOrbitPostProcess(int order,
                         int sizeOrbit,
                         int Nex,
                         int isPar)
{
    //----------------------------------------------------------
    //Read data size
    //----------------------------------------------------------
    int N, Nman;
    getLenghtManifold_bin(order, sizeOrbit, &N, &Nman, TYPE_MAN);
    N--;   //need to shift to take into account that vectors starts at 0
    Nman--; //need to shift to take into account that vectors starts at 0

    int Nt = Nex < 0 ? N:Nex;

    //----------------------------------------------------------
    //To store all data
    //----------------------------------------------------------
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //----------------------------------------------------------
    //Read data
    //----------------------------------------------------------
    int status = readManifold_bin(tTens, yTens, order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

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
            double **yManNCSEM = dmatrix(0, 5, 0, Nman);
            double **yManSEM   = dmatrix(0, 5, 0, Nman);
            double *tManSEM    = dvector(0, Nman);
            double *HMan       = dvector(0, Nman);
            double *HSEMLi     = dvector(0, Nman);


            //-------------------------------
            // Update yManNCSEM, tManNCSEM
            //-------------------------------
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) yManNCSEM[i][p] = yTens[i][k][p];
                tManSEM[p] = tTens[k][p];
            }

            //-------------------------------
            //To SEM coordinates
            //-------------------------------
            NCtoSEM_vec(yManNCSEM, tManSEM, yManSEM, Nman, &SEML);

            //-------------------------------
            //Compute the energy
            //on the last part of the trajectory
            //-------------------------------
            //Along the manifold
            HSEM_vec(tManSEM, yManSEM, HMan, Nman, &SEML);
            //Along SEMLi on the same time span
            HSEMLi_vec(tManSEM, HSEMLi, Nman, &SEML);
            //Delta
            DH[k] = 0.0;
            for(int p = 0; p <= Nman; p++)
            {
                if(fabs(tManSEM[Nman] - tManSEM[p]) < SEML.us_sem.T) DH[k] += fabs(HMan[p] - HSEMLi[p]);
            }

            //-------------------------------
            //Compute the distance the SEMLi
            //on the last part of the trajectory
            //-------------------------------
            double ePm;
            double yv[6], yvc[6];
            for(int p = 0; p <= Nman; p++)
            {
                if(true /*fabs(tManSEM[Nman] - tManSEM[p]) < SEML.us_sem.T*/)
                {
                    //Store in yv
                    for(int i = 0; i < 6; i++) yv[i] = yManSEM[i][p];
                    //From momenta to velocities
                    SEMmtoSEMv(tManSEM[p], yv, yvc, &SEML);
                    //distance between yvc and the Li point
                    ePm = 0;
                    for(int i = 0; i < 3; i++) ePm += (semP[SEML.li_SEM+4][i] - yvc[i])*(semP[SEML.li_SEM+4][i] - yvc[i]);
                    //Adding velocities
                    //for(int i = 4; i < 6; i++) ePm += (yvc[i])*(yvc[i]);
                    //Store
                    DE[k] += sqrt(ePm);
                }
            }

            //----------------------------------------------------------
            //Free
            //----------------------------------------------------------
            free_dmatrix(yManNCSEM, 0, 5, 0, Nman);
            free_dmatrix(yManSEM, 0, 5, 0, Nman);
            free_dvector(tManSEM, 0, Nman);
            free_dvector(HMan, 0, Nman);
            free_dvector(HSEMLi, 0, Nman);
        }

        //----------------------------------------------------------
        //Sorting DH by indexes
        //----------------------------------------------------------
        vector<size_t> sortId = sort_indexes(DE);

        //----------------------------------------------------------
        //Select the first Nt values
        //----------------------------------------------------------
        double **yManNCSEM = dmatrix(0, 5, 0, Nman);
        double *tManSEM    = dvector(0, Nman);
        double ***yTR = d3tensor(0, 5, 0, Nt, 0, Nman);
        double **tTR  = dmatrix(0, Nt, 0, Nman);
        int sortk;
        for(int k = 0; k <= Nt; k++)
        {
            sortk = sortId[k];
            //-------------------------------
            // Update yManNCSEM, tManNCSEM
            //-------------------------------
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) yManNCSEM[i][p] = yTens[i][sortk][p];
                tManSEM[p] = tTens[sortk][p];

                for(int i = 0; i < 6; i++) yTR[i][k][p] = yTens[i][sortk][p];
                tTR[k][p] = tTens[sortk][p];
            }
        }

        //---------------------------------------------------------------------
        //Save
        //---------------------------------------------------------------------
        writeManifold_bin(tTR, yTR, order, sizeOrbit, Nt+1, Nman+1, TYPE_MAN_SORT_DR);

        //---------------------------------------------------------------------
        //Free
        //---------------------------------------------------------------------
        free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
        free_d3tensor(yTR, 0, 5, 0, Nt, 0, Nman);
        free_dmatrix(tTens, 0, N, 0, Nman);
        free_dmatrix(tTR, 0, Nt, 0, Nman);
        free_dmatrix(yManNCSEM, 0, 5, 0, Nman);
        free_dvector(tManSEM, 0, Nman);
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

    int Nt = Nex < 0 ? N:Nex;

    //To store all data
    double ***yTens = d3tensor(0, 5, 0, N, 0, Nman);
    double **tTens  = dmatrix(0, N, 0, Nman);

    //To store final data
    double ***yTR = d3tensor(0, 5, 0, Nt, 0, Nman+1);
    double **tTR  = dmatrix(0, Nt, 0, Nman+1);

    //Read data
    readManifold_bin(tTens, yTens, order_em, size_em, N+1, Nman+1, type_em);

    //-------------------------------------------------------------------------------
    //Changing the scope to CM SEM...
    //-------------------------------------------------------------------------------
    changeScope(order_sem, 4, F_SEM);

    //--------------------------------------
    // Center- manifold
    //--------------------------------------
    vector<Oftsc>  CMc(6);     ///center manifold in NC coordinates
    vector<Oftsc> CMhc(6);     ///center manifold in TFC coordinates

    //------------------------------------------
    // Update of the central manifold
    //------------------------------------------
    //Update CM
    readVOFTS_bin(CMc,  SEML.cs.F_PMS+"W/W",  OFS_ORDER);
    readVOFTS_bin(CMhc, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //----------------------------------------------------------
    //To store
    //----------------------------------------------------------
    vector<int>    indexMin(N+1);
    vector<double> distMin(N+1);
    double **yfSEMin = dmatrix(0, 5, 0, N);
    double **yprojSEMin = dmatrix(0, 5, 0, N);

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
        double yvSE[6], yvprojSE[6];
        double ePm, ePmin = 1e6;
        double kmin = 0;
        double yvnorm = 0.0;
        //----------------------------------------------------------
        //Loop on trajectory
        //----------------------------------------------------------
        for(int kman = 0; kman <= Nman; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) yv[i] = yTens[i][kpos][kman];
            tv = tTens[kpos][kman];

            yvnorm = 0.0;
            for(int i = 0; i < 2; i++) yvnorm += yv[i]*yv[i];
            yvnorm = sqrt(yvnorm);

            if(yvnorm < 0.5)
            {
                // Projection on the center manifold
                NCprojCCMtoCM(yv, tv, sproj, CMhc, MIcoc, Vcoc);
                //yvproj = W(sproj, tv)
                RCMtoNCbyTFC(sproj, tv, SEML.us.n, OFTS_ORDER, OFS_ORDER, CMhc, ofs, Mcoc, Vcoc, yvproj, 1);
                //Distance of projection in SEM coordinates
                NCtoSEM(tv, yv, yvSE, &SEML);
                NCtoSEM(tv, yvproj, yvprojSE, &SEML);
                ePm = 0.0;
                for(int i = 0; i < 3; i++) ePm += (yvprojSE[i] - yvSE[i])*(yvprojSE[i] - yvSE[i]);
                ePm = sqrt(ePm);
            }
            else ePm = 1e5;

            //Update distance min if necessary
            if(ePm < ePmin)
            {
                ePmin = ePm;
                kmin  = kman;
                for(int i = 0; i < 6; i++) yfSEMin[i][kpos] = yvprojSE[i];
                for(int i = 0; i < 6; i++) yprojSEMin[i][kpos] = yvproj[i];
            }

        }

        //----------------------------------------------------------
        //Store results
        //----------------------------------------------------------
        indexMin[kpos] = (int) kmin;
        distMin[kpos]  = ePmin;

        cout << "progression : " << 100.0*(index++)/N << endl;
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
    for(int kpos = 0; kpos <= Nt; kpos++)
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
        double **yManNCSEM  = dmatrix(0, 5, 0, Nman);
        double **yManSEM   = dmatrix(0, 5, 0, Nman);
        double *tManSEM    = dvector(0, Nman);
        for(int kman = 0; kman <= Nman; kman++)
        {
            // Current state
            for(int i = 0; i < 6; i++) yManNCSEM[i][kman] = yTens[i][ksortpos][kman];
            tManSEM[kman] = tTens[ksortpos][kman];

            // Current state
            for(int i = 0; i < 6; i++) yTR[i][kpos][kman] = yTens[i][ksortpos][kman];
            tTR[kpos][kman] = tTens[ksortpos][kman];
        }

        //----------------------------------------------------------
        //Update final data
        // Careful: we save the projection state at the end of yTR
        // Careful: we save the indix of the minimum projection distance the end of tTR
        //----------------------------------------------------------
        for(int i = 0; i < 6; i++) yTR[i][kpos][Nman+1] = yprojSEMin[i][ksortpos];
        tTR[kpos][Nman+1] = indexMin[ksortpos];

        //To SEM coordinates
        NCtoSEM_vec(yManNCSEM, tManSEM, yManSEM, Nman, &SEML);
        //plot
        gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
        //plot the right solution
        gnuplot_plot_xyz(h2, &yManSEM[0][indexMin[ksortpos]], &yManSEM[1][indexMin[ksortpos]], &yManSEM[2][indexMin[ksortpos]], 1, (char*)"", "points", "1", "2", 7);
        gnuplot_plot_xyz(h2, &yfSEMin[0][ksortpos], &yfSEMin[1][ksortpos], &yfSEMin[2][ksortpos], 1, (char*)"", "points", "1", "2", 8);
    }


    //---------------------------------------------------------------------
    // Save
    // Careful: we save the projection state at the end of yTR
    // Careful: we save the indix of the minimum projection distance at the end of tTR
    //---------------------------------------------------------------------
    cout << "here" << endl;
    writeManifold_bin(tTR, yTR, order_em, size_em, Nt+1, Nman+2, TYPE_MAN_PROJ);

    //---------------------------------------------------------------------
    //Free
    //---------------------------------------------------------------------
    free_d3tensor(yTens, 0, 5, 0, N, 0, Nman);
    free_d3tensor(yTR, 0, 5, 0, Nt, 0, Nman+1);
    free_dmatrix(tTens, 0, N, 0, Nman);
    free_dmatrix(tTR, 0, Nt, 0, Nman+1);


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
                     vector<Oftsc> &CMh,
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
    init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 6, 1e-6, qbfbp_vfn_novar, NULL, &SEML);

    char ch;
    //----------------------------------------------------------
    //Display results
    //----------------------------------------------------------
    double yv0[6], yvf[6], yvfproj[6];
    double yvSE[6];
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
        NCtoSEM(t0, yv0, yvSE, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yvSE[0], &yvSE[1], &yvSE[2], 1, (char*)"", "points", "1", "2", 6);

        //Proj point: ksol = Nman
        ksol = Nman;
        for(int i = 0; i < 6; i++) yvfproj[i] = yTens[i][kpos][ksol];
        tf = tTens[kpos][ksol];
        //To SE coordinates
        NCtoSEM(tf, yvfproj, yvSE, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yvSE[0], &yvSE[1], &yvSE[2], 1, (char*)"", "points", "1", "2", 8);

        //Final point: ksol = last indix stored in tTens
        ksol = (int) tTens[kpos][Nman];
        tf = tTens[kpos][ksol];
        for(int i = 0; i < 6; i++) yvf[i] = yTens[i][kpos][ksol];
        //To SE coordinates
        NCtoSEM(tf, yvf, yvSE, &SEML);
        //Plot
        gnuplot_plot_xyz(h2, &yvSE[0], &yvSE[1], &yvSE[2], 1, (char*)"", "points", "1", "2", 7);


        //---------------------------------------------------------------------
        //Equivalent from data file
        //---------------------------------------------------------------------
        //From data
        double **yManNCSEM  = dmatrix(0, 5, 0, Nman);
        double *tManSEM     = dvector(0, Nman);
        double **yManSEM    = dmatrix(0, 5, 0, Nman);
        for(int kman = 0; kman <= Nman; kman++)
        {
            for(int i = 0; i < 6; i++) yManNCSEM[i][kman] = yTens[i][kpos][kman];
            tManSEM[kman] = tTens[kpos][kman];
        }

        //To SE coordinates
        NCtoSEM_vec(yManNCSEM, tManSEM, yManSEM, Nman, &SEML);

        //----------------------------------------------------------
        //Gnuplot window
        //----------------------------------------------------------
        gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], ksol+1, (char*)"", "lines", "1", "2", 5);

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
        init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 6, 1e-6, qbfbp_vfn_novar, NULL, &SEML);
        //The default interval of projection is set to Tproj = T/5
        double tproj = SEML.us.T/5.0;
        //------------------
        //Orbit structure
        //------------------
        Ofsc orbit_ofs(OFS_ORDER);
        SingleOrbit orbit;
        //Init routine
        init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
        init_ode_structure(&driver2, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 42, 1e-6, qbfbp_vfn_varnonlin, NULL, &SEML);

        double ystart[42];
        //Storing initial position and momenta into ystart
        for(int i = 0; i < 6; i++) ystart[i] = yv0[i];
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
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
        NCtoSEM(t0, ystart, yvSE, &SEML);
        gnuplot_plot_xyz(h2, &yvSE[0], &yvSE[1], &yvSE[2], 1, (char*)"", "points", "1", "5", 7);


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
        NCtoSEM(tf, yvfproj, yvSE, &SEML);
        SEMmtoSEMv(tf, yvSE, yprojVel, &SEML);

        double dV = 0.0;
        for(int i = 3; i < 6; i++) dV += (yprojVel[i] -  yfinalVel[i])*(yprojVel[i] -  yfinalVel[i]);
        dV = sqrt(dV);
        cout << "Final Velocity gap: (m/s)" << endl;
        cout << 1e3*dV*SEML.cs_sem.cr3bp.L*2*M_PI/(SEML.cs_sem.cr3bp.T) << endl;

        //----------------------------------------------------------
        // Starting Velocity gap
        //----------------------------------------------------------
        //Initial position on the old trajectory
        NCtoSEM(t0, yv0, yvSE, &SEML);
        SEMmtoSEMv(t0, yvSE, yprojVel, &SEML);
        //Initial position on the new trajectory
        for(int i = 0; i <6; i++) yv0[i] = ystart[i];
        NCtoSEM(t0, yv0, yvSE, &SEML);
        SEMmtoSEMv(t0, yvSE, yfinalVel, &SEML);


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
void manOrbitPlot(int order,
                  int sizeOrbit,
                  int Nt,
                  int Nmant,
                  int type)
{
    //Read data sizeOrbit
    int N0, Nman0, N, Nman;
    getLenghtManifold_bin(order, sizeOrbit, &N0, &Nman0, type);
    N0--;   //need to shift to take into account that vectors starts at 0
    Nman0--; //need to shift to take into account that vectors starts at 0

    //Use value from file if the user wants it
    if(Nt < 0) N = N0;
    else N = Nt;
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
    int status = readManifold_bin(tTens, yTens, order, sizeOrbit, N+1, Nman+1, type);

    //Plot loop
    if(status)
    {
        //Gnuplot window
        gnuplot_ctrl  *h1;
        h1 = gnuplot_init();


        //Temp variable
        double **yManNCSEM = dmatrix(0, 5, 0, Nman);
        double **yManSEM   = dmatrix(0, 5, 0, Nman);
        double *tManNCSEM  = dvector(0, Nman);

        for(int k = 0; k <= N; k++)
        {
            for(int p = 0; p <= Nman; p++)
            {
                for(int i = 0; i < 6; i++) yManNCSEM[i][p] = yTens[i][k][p];
                tManNCSEM[p] = tTens[k][p];
            }

            //To SEM coordinates
            NCtoSEM_vec(yManNCSEM, tManNCSEM, yManSEM, Nman, &SEML);
            gnuplot_plot_xyz(h1, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);
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
int manOrbitSEM(double tman,
                int order,
                int sizeOrbit,
                int Nt,
                int Nman,
                int isPar,
                vector<Oftsc> &CM,
                vector<Oftsc> &CMh,
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
    int N0 = getLineNumber(order, sizeOrbit, TYPE_CU)-1;

    //----------------------------------------------------------
    //If the file exists and contains data, go on
    //----------------------------------------------------------
    if(N0 > 0)
    {
        if(Nt <= 0) N = N0;
        else N = Nt;

        if(N > N0)
        {
            cout << "manOrbit: wrong input, Nt = " << Nt << ", but N0 = " << N0 << endl;
            return 0;
        }

        //----------------------
        // Read data from file
        //----------------------
        double **yNCU = dmatrix(0, 5, 0, N);
        double **sNCU = dmatrix(0, 4, 0, N);
        double *tNCU  = dvector(0, N);
        int status    = readOrbit(tNCU, yNCU, sNCU, order, sizeOrbit, N+1, TYPE_CU);

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
                double **yManNCEM  = dmatrix(0, 5, 0, Nman);
                double **yManNCSEM = dmatrix(0, 5, 0, Nman);
                double *tManEM     = dvector(0, Nman);

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
                init_ode_structure(&driver, T, T_root, PREC_ABS, PREC_REL, 1e-13, 1e-12, 6, 1e-6, qbfbp_vfn_novar, NULL, &SEML);

                //------------------
                //Orbit structure
                //------------------
                Ofsc orbit_ofs(OFS_ORDER);
                SingleOrbit orbit;
                //Init routine
                init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, tv, tv+tman, tproj, &driver, &SEML);

                //---------------------------------------------------------------------
                // Update the state
                //---------------------------------------------------------------------
                orbit_update_ic(orbit, st, yv, tv);

                //---------------------------------------------------------------------
                //Integration
                //---------------------------------------------------------------------
                trajectory_integration_grid(orbit, tv, tv+tman, yManNCEM, tManEM, Nman, 0);

                //---------------------------------------------------------------------
                //Integration check
                //---------------------------------------------------------------------
                reset_ode_structure(orbit.driver);
                double yv0[6];
                for(int i = 0; i <6; i++) yv0[i] = yv[i];

                double t0 = tManEM[0];
                int kk = Nman;

                cout << "Time" << endl;
                cout << tManEM[0] << "  " << t0 << "   " << endl;


                cout << "Init" << endl;
                for(int i = 0; i <6; i++) cout << yManNCEM[i][0] - yv[i] << endl;


                //Integration
                gsl_odeiv2_driver_apply (orbit.driver->d, &t0, tManEM[kk], yv);

                cout << "Time" << endl;
                cout << tManEM[kk] << "  " << t0 << "   " << endl;

                cout << "Result" << endl;
                for(int i = 0; i <6; i++) cout << yManNCEM[i][kk] - yv[i] << endl;

                //Integration
                reset_ode_structure(orbit.driver);
                t0 = tManEM[0];
                gsl_odeiv2_driver_apply(orbit.driver->d, &t0, tManEM[kk], yv0);

                cout << "Result" << endl;
                for(int i = 0; i <6; i++) cout << yv[i] - yv0[i] << endl;

                //---------------------------------------------------------------------
                //Storage
                //---------------------------------------------------------------------
                #pragma omp critical
                {
                    //To NCSEM coordinates
                    NCEMmtoNCSEMm_vec(yManNCEM, tManEM, yManNCSEM, Nman, &SEML);
                    for(int p = 0; p <= Nman; p++)
                    {
                        for(int i = 0; i < 6; i++) yTens[i][k][p] = yManNCSEM[i][p];
                        tTens[k][p] = tManEM[p]*SEML.us_em.ns; //time is also given in SEM units!
                    }
                }

                //---------------------------------------------------------------------
                //Free
                //---------------------------------------------------------------------
                free_dmatrix(yManNCEM, 0, 5, 0, Nman);
                free_dmatrix(yManNCSEM, 0, 5, 0, Nman);
                free_dvector(tManEM, 0, Nman);
            }


            //---------------------------------------------------------------------
            //Save
            //---------------------------------------------------------------------
            writeManifold_bin(tTens, yTens, order, sizeOrbit, N+1, Nman+1, TYPE_MAN);

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


//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Main routines (3): old versions of connections
//
//---------------------------------------------------------------------------------------------------------------------------------------
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
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli2(double st0[],
                 double t0,
                 double tf,
                 double tman,
                 double tplot,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
    int Nman = 200*floor(fabs(tman)/SEML.us.T);
    double **yManNCEM  = dmatrix(0, 5, 0, Nman);
    double **yManSEM   = dmatrix(0, 5, 0, Nman);
    double *tManEM     = dvector(0, Nman);
    double *tManSEM    = dvector(0, Nman);

    double *HMan     = dvector(0, Nman);
    double *HSEMLi   = dvector(0, Nman);
    double **yManNCSEM  = dmatrix(0, 5, 0, Nman);
    double DH = 0, DHmax = 1e20;

    //------------------------------------------
    // Find the minimum distance to Li
    //------------------------------------------
    int kopt = 0;
    double ePm, ePmOpt, ePav, ePmin[kn], kvec[kn], yvmin[6];
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
        NCprojCCMtoCUS(yv, tv, sti, CMh, +1e-6, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        ntm = trajectory_integration_variable_grid(orbit, tv, tv+tman, yManNCEM, tManEM, Nman, 0);

        if(ntm > 0)
        {
            //------------------------------------------
            //NC EM to SEM
            //------------------------------------------
            NCEMmtoSEMm_vec(yManNCEM, tManEM, yManSEM, tManSEM, ntm, &SEML);

            //------------------------------------------
            //Integrated distance to Li
            //------------------------------------------
            kvec[indix] = indix;
            ePmin[indix] = 1e6;
            ePav = 0.0;
            for(int p = 0; p <= ntm; p++)
            {
                //Store in yv
                for(int i = 0; i < 6; i++) yv[i] = yManSEM[i][p];
                //From momenta to velocities
                SEMmtoSEMv(tManSEM[p], yv, yvc, &SEML);
                //distance between yvc and the Li point
                ePm = 0;
                for(int i = 0; i < 3; i++) ePm += (semP[SEML.li_SEM+4][i] - yvc[i])*(semP[SEML.li_SEM+4][i] - yvc[i]);
                //Adding velocities
                for(int i = 4; i < 6; i++) ePm += (yvc[i])*(yvc[i]);
                ePm = sqrt(ePm);
                //Add to the current moving average
                ePav = getDeltaMovingAverage(ePm, listDeltaMA, Naverage);
                if(ePav < ePmin[indix]) ePmin[indix] = ePav;
            }

            //------------------------------------------
            //Energy
            //------------------------------------------
            //Along the manifold
            HSEM_vec(tManSEM, yManSEM, HMan, ntm, &SEML);
            //Along SEMLi on the same time span
            HSEMLi_vec(tManSEM, HSEMLi, ntm, &SEML);
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
                //for(int i = 0; i < 6; i++) yvmin[i] = yManSEM[i][ntm];     //yvmin = yv
                //kopt = k;
                //indOpt = indix;
                //Orbit
                //gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], ntm+1, (char*)"", "lines", "1", "1", 2);
                //gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 7);
                //Energy
                //gnuplot_plot_xy(h4, tManSEM, HMan, ntm+1, (char*)"", "lines", "1", "2", 1);
                //gnuplot_plot_xy(h4, tManSEM, HSEMLi, ntm+1, (char*)"", "lines", "1", "2", 2);
                //-------------------------------------------------------------------
            }

            //Update the min distance
            if(ePmOpt > ePmin[indix])
            {
                ePmOpt = ePmin[indix];
                for(int i = 0; i < 6; i++) yvmin[i] = yManSEM[i][ntm];     //yvmin = yv
                kopt = k;
                indOpt = indix;
                gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], ntm+1, (char*)"", "lines", "1", "1", 2);
                gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 7);
                //Energy
                gnuplot_plot_xy(h4, tManSEM, HMan, ntm+1, (char*)"", "lines", "1", "2", 1);
                gnuplot_plot_xy(h4, tManSEM, HSEMLi, ntm+1, (char*)"", "lines", "1", "2", 2);
            }
            //else gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], ntm+1, (char*)"", "lines", "1", "1", 7);
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
    NCprojCCMtoCUS(yv, tv, sti, CMh, +1e-6, MIcoc, Vcoc);

    //----------------------
    // Update the state
    //----------------------
    orbit_update_ic(orbit, sti, tv);

    //------------------------------------------
    //Integration
    //------------------------------------------
    trajectory_integration_grid(orbit, tv, tv+tman, yManNCEM, tManEM, Nman, 0);

    //------------------------------------------
    //NC EM to SEM
    //------------------------------------------
    NCEMmtoSEMm_vec(yManNCEM, tManEM, yManSEM, tManSEM, Nman, &SEML);

    //------------------------------------------
    //Energy
    //------------------------------------------
    //Along the manifold
    HSEM_vec(tManSEM, yManSEM, HMan, Nman, &SEML);
    //Along SEMLi on the same time span
    HSEMLi_vec(tManSEM, HSEMLi, Nman, &SEML);

    //Misc
    gnuplot_plot_xy(h4, tManSEM, HMan, Nman+1, (char*)"orbit", "lines", "1", "3", 4);
    gnuplot_plot_xy(h4, tManSEM, HSEMLi, Nman+1, (char*)"SEMLi", "lines", "1", "2", 2);


    //Plot the min distance
    gnuplot_cmd(h3, "set grid");
    gnuplot_set_xlabel(h3, (char*)"kvec [-]");
    gnuplot_set_ylabel(h3, (char*)"ePm [-]");
    gnuplot_plot_xy(h3, kvec, ePmin, kn, (char*)"", "points", "1", "2", 1);
    gnuplot_plot_xy(h3, &kvec[indOpt], &ePmin[indOpt], 1, (char*)"", "points", "2", "2", 4);
    //Plot the corresponding optimum
    gnuplot_plot_xyz(h2, yvmin, yvmin+1, yvmin+2, 1, (char*)"", "points", "1", "1", 5);
    //Plot the corresponding manifold leg
    gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "2", 4);


    cout << "Best solution has the indix " << kopt << endl;
    cout << "Which corresponds to the initial time " << tNCE[kopt] << "in EM units" << endl;
    printf("Press ENTER to close the gnuplot window(s)\n");
    //scanf("%c",&ch);

    //-----------------------------------------------------------------------------------------------------------
    // Change to center manifold
    //-----------------------------------------------------------------------------------------------------------
    cout << "Changing the default parameters..." << endl;
    REDUCED_NV=4;
    OFTS_ORDER=15;
    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc>  CMc(6);     ///center manifold in NC coordinates
    vector<Oftsc> CMhc(6);     ///center manifold in TFC coordinates

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
    readVOFTS_bin(CMhc, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //--------------------------------------
    //Init routine for the orbit
    //--------------------------------------
    init_orbit(orbit, &CMc, &CMhc, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

    //------------------------------------------
    //Get the minimum distance to SEML manifold
    //------------------------------------------
    //To NC SEM coordinates
    SYStoNC_vec(yManSEM, tManSEM, yManNCSEM, Nman, &SEML);
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
        for(int i = 0; i <6; i++) yv[i] = yManNCSEM[i][p];
        //The current time is the time on the orbit, in NC SEM coordinates
        tv = tManSEM[p];
        //Projection on the center manifold
        NCprojCCMtoCUS(yv, tv, sti, CMhc, 0.0, MIcoc, Vcoc);
        //----------------------
        //Minimum distance?
        //----------------------
        //zproj = W(sti, tv)
        RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, CMhc, *orbit.ofs, Mcoc, Vcoc, zproj, 1);
        //distance between zproj and yv
        ePm = 0;
        for(int i = 0; i < 6; i++)
        {
            ePm += (zproj[i] - yv[i])*(zproj[i] - yv[i]);
            ePvec[i] = fabs(zproj[i] - yv[i]);
        }
        ePm = sqrt(ePm);
        ePmat[p] = ePm;
        pmat[p] = p;

        //Update the min distance
        if(ePmOpt > ePm)
        {
            pmin = p;
            ePmOpt = ePm;
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
    gnuplot_plot_xy(h3, pmat+pmin, ePmat+pmin, 1, (char*)"ePmin", "points", "1", "2", 4);

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
    free_dmatrix(yManNCEM, 0, 5, 0, Nman);
    free_dmatrix(yManNCSEM, 0, 5, 0, Nman);
    free_dmatrix(yManSEM, 0, 5, 0, Nman);
    free_dvector(tNCE, 0, N);
    free_dvector(tManEM, 0, Nman);
    free_dvector(tManSEM, 0, Nman);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);
    gnuplot_close(h4);

    return status;
}


/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int eml2Tosemli_torb_vs_tman(double st0[],
                             double t0,
                             double tf,
                             double tman,
                             vector<Oftsc> &CM,
                             vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
    gnuplot_set_ylabel(h3, (char*)"tman [ x Tsun]");
    gnuplot_set_zlabel(h3, (char*)"ePm  [SEM units]");

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
    int Nman = 200*floor(fabs(tman)/SEML.us.T);
    double **yMan = dmatrix(0, 5, 0, Nman);
    double **yManSEM = dmatrix(0, 5, 0, Nman);
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
        NCprojCCMtoCUS(yv, tv, sti, CMh, +1e-6, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tman, yMan, tMan, Nman, 0);

        //------------------------------------------
        //Store data
        //------------------------------------------
        //NC EM to NC SEM
        NCEMmtoNCSEMm_vec(yMan, tMan, yManSEM, Nman, &SEML);
        for(int i = 0; i <= Nman; i++)
        {
            for(int k = 0; k <6; k++) yTens[k][i][indix] = yManSEM[k][i];
            tTens[i][indix] = tMan[i]*SEML.us_em.ns; //time is also given in SEM units!
        }
        indix++;

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        //NC EM to SEM (no NC for the output here)
        //NCEMmtoSEMm_vec(yMan, tMan, yManSEM, Nman, &SEML);
        //gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "1", 2);
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
    readVOFTS_bin(CMh, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //------------------------------------------
    // Projection loop
    //------------------------------------------
    double *distMan  = dvector(0, Nman);
    double *disttorb = dvector(0, Nman);
    double *disttman = dvector(0, Nman);

    double zproj[6], yvSEM[6], zprojSEM[6], ePm;
    for(int k = 0; k < kn; k++)
    {
        ePm = 2;
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
            NCprojCCMtoCUS(yv, tv, sti, CMh, 0.0, MIcoc, Vcoc);

            //----------------------
            //Minimum distance?
            //----------------------
            //zproj = W(sti, tv)
            RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, CMh, *orbit.ofs, Mcoc, Vcoc, zproj, 1);

            //Back to SEM coordinates
            NCtoSYS(tv, yv, yvSEM, &SEML);
            NCtoSYS(tv, zproj, zprojSEM, &SEML);

            //distance between zproj and yv
            ePm = 0;
            for(int i = 0; i < 6; i++) ePm += (zprojSEM[i] - yvSEM[i])*(zprojSEM[i] - yvSEM[i]);
            ePm = sqrt(ePm);

            //----------------------
            //Save value
            //----------------------
            distMan[p]  = ePm;
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
    free_dmatrix(yManSEM, 0, 5, 0, Nman);
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
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-(un)stable manifold CM.
 **/
int orbit_cus(double st0[],
              double t0,
              double tf,
              double tman,
              vector<Oftsc> &CM,
              vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
    int Nman = 200*floor(fabs(tman)/SEML.us.T);
    int color = (SEML.cs_em.manType == MAN_CENTER_S)?  2:7;
    cout << "color = " << color << endl;

    //Manifold
    double **yManNC  = dmatrix(0, 5, 0, Nman);
    double **yManSEM = dmatrix(0, 5, 0, Nman);
    double **yManEM  = dmatrix(0, 5, 0, Nman);
    double *tMan     = dvector(0, Nman);
    double *tManSEM  = dvector(0, Nman);
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
        NCprojCCMtoCUS(yv, tv, sti, CMh, +1e-6, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        //trajectory_integration_grid(orbit, tv, tv+tman, yManNC, tMan, Nman, 0);
        nt = trajectory_integration_variable_grid(orbit, tv, tv+tman, yManNC, tMan, Nman, 0);

        //------------------------------------------
        //Plotting
        //------------------------------------------
        //NCtoSYS_vec(yManNC, tMan, yManEM, nt, &SEML);
        //gnuplot_plot_xyz(h1, yManEM[0], yManEM[1], yManEM[2], nt+1, (char*)"", "lines", "1", "1", color);

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        NCEMmtoSEMm_vec(yManNC, tMan, yManSEM, tManSEM, nt, &SEML);

        //Compute the maximum distance from the Earth
        double maxDistToEarth = 0.0;
        double res = 0.0;
        for(int i = 0; i <= nt; i++)
        {
            //Earth is situated at Pe
            res = 0.0;
            for(int p = 0; p < 3; p++) res += (yManSEM[p][i] - semP[0][p])*(yManSEM[p][i] - semP[0][p]);
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
        HSEM_vec(tManSEM, yManSEM, HMan, nt, &SEML);
        //Along SEMLi on the same time span
        HSEMLi_vec(tManSEM, HSEMLi, nt, &SEML);
        //Delta
        DH = 0;
        for(int i = 0; i <= nt; i++) DH += fabs(HMan[i] - HSEMLi[i]);
        DH = sqrt(DH);

        if(maxDistToEarth < 0.02 && DH < DHmax)
        {
            DHmax = DH;
            //Orbit
            gnuplot_plot_xyz(h2, yManSEM[0], yManSEM[1], yManSEM[2], nt+1, (char*)"", "lines", "1", "1", color);
            //Energy
            gnuplot_plot_xy(h3, tManSEM, HMan, nt+1, (char*)"", "lines", "1", "2", 1);
        }
    }

    //------------------------------------------
    //Level of energy of SEML2 on the last solution
    //------------------------------------------
    gnuplot_plot_xy(h3, tManSEM, HSEMLi, nt+1, (char*)"SEML2", "lines", "1", "2", 2);

    //------------------------------------------
    //Level of energy of EML2
    // TODO: - do the same in EM units (Energy)
    //       - plot the energy along the initial orbit
    //------------------------------------------
    HEMLi_in_SEM_vec(tMan, HMan, nt, &SEML);
    gnuplot_plot_xy(h3, tManSEM, HMan, nt+1, (char*)"EML2", "lines", "1", "2", 3);


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
    free_dmatrix(yManSEM, 0, 5, 0, Nman);
    free_dvector(tMan, 0, Nman);
    gnuplot_close(h1);
    gnuplot_close(h2);
    gnuplot_close(h3);

    return status;
}


/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initian conditions st0, on the center-unstable manifold CM (SEM framework).
 **/
int orbit_cu_sem(double st0[],
                 double t0,
                 double tf,
                 double tman,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
    tplot = fabs(tman)/5000;
    int Nman = floor(fabs(tman)/tplot);

    double **yMan = dmatrix(0, 5, 0, Nman);
    double **yManSEM = dmatrix(0, 5, 0, Nman);
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
        NCprojCCMtoCUS(yv, tv, sti, CMh, +1e-6, MIcoc, Vcoc);

        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tman, yMan, tMan, Nman, 0);

        //------------------------------------------
        //Plotting in SEM coordinates
        //------------------------------------------
        NCtoSYS_vec(yMan, tMan, yManSEM, Nman, &SEML);

        //------------------------------------------
        //Compute the minimum distance from the Earth
        //------------------------------------------
        double minDistToEarth = 1e6;
        double res = 0.0;
        for(int i = 0; i <= Nman; i++)
        {
            //Distance
            res = 0.0;
            for(int k = 0; k < 3; k++) res += (yManSEM[k][i] - Pe[k])*(yManSEM[k][i] - Pe[k]);
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
            for(int k = 0; k < 3; k++) res += (yManSEM[k][i] - Pm[k])*(yManSEM[k][i] - Pm[k]);
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
            for(int i = 0; i <=Nman; i++) for(int k = 0; k < 6; k++) yManSEMOpt[k][i] = yManSEM[k][i];
            //gnuplot_plot_xyz(hsem, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "1", 1);
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
    free_dmatrix(yManSEM, 0, 5, 0, Nman);
    free_dvector(tMan, 0, Nman);

    return status;
}

/**
 *  \brief Computes an orbit on a given time grid [t0 tf tf+tman], with initial conditions st0, on the center-unstable manifold CM.
 *         Then, computes the minimum distance to the center-manifold of a SEM libration point.
 **/
int semliToseml2(double st0[],
                 double t0,
                 double tf,
                 double tman,
                 vector<Oftsc> &CM,
                 vector<Oftsc> &CMh,
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
                       1e-13,           //precision: root
                       1e-12,           //precision: diffcorr
                       6,               //dimension
                       1e-6,            //initial int step
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
    init_orbit(orbit, &CM, &CMh, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs,
               OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &driver, &SEML);

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
    tplot = tman/5000;
    int Nman = floor(fabs(tman)/tplot);
    double **yMan  = dmatrix(0, 5, 0, Nman);
    double **yManSEM  = dmatrix(0, 5, 0, Nman);
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
        double epsilon = 1e-6;
        if(SEML.li_SEM == 1) epsilon = +1e-6;
        else if(SEML.li_SEM == 2) epsilon = -1e-6;

        //Projection
        NCprojCCMtoCUS(yv, tv, sti, CMh, epsilon, MIcoc, Vcoc);
        //----------------------
        // Update the state
        //----------------------
        orbit_update_ic(orbit, sti, tv);

        //------------------------------------------
        //Integration
        //------------------------------------------
        trajectory_integration_grid(orbit, tv, tv+tman, yMan, tMan, Nman, 0);

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
        NCtoSYS_vec(yMan, tMan, yManSEM, Nman, &SEML);
        gnuplot_plot_xyz(hse, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "1", 2);

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
    readVOFTS_bin(CMh, SEML.cs.F_PMS+"W/Wh", OFS_ORDER);
    //Update COC
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //------------------------------------------
    // Projection loop
    //------------------------------------------
    int kopt = 0;
    double zproj[6], ePm, ePmin[kn], kvec[kn];
    double yvmin[6], zprojmin[6];
    double yvopt[6], zprojopt[6];
    double tvproj = 0.0;
    double tvopt = 0.0;
    double ePopt = 0.0;
    for(int k = 0; k < kn; k++)
    {
        kvec[k] = k;
        ePmin[k] = 0.0;
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
            NCprojCCMtoCUS(yv, tv, sti, CMh, 0.0, MIcoc, Vcoc);
            //----------------------
            //Minimum distance?
            //----------------------
            //zproj = W(sti, tv)
            RCMtoNCbyTFC(sti, tv, orbit.n, orbit.order, orbit.ofs_order, CMh, *orbit.ofs, Mcoc, Vcoc, zproj, 1);
            //distance between zproj and yv
            ePm = 0;
            for(int i = 0; i < 6; i++) ePm += (zproj[i] - yv[i])*(zproj[i] - yv[i]);
            ePm = sqrt(ePm);

            //cout << "ePm = " << ePm << endl;
            //Update the min distance
            if(p == 0)
            {
                ePmin[k] = ePm;
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
                if(ePmin[k] > ePm)
                {
                    ePmin[k] = ePm;
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
            ePopt = ePmin[k];
        }
        else
        {
            if(ePmin[k] < ePopt)
            {
                kopt = k;
                tvopt = tvproj;
                for(int i = 0; i <6; i++)
                {
                    yvopt[i] = yvmin[i];
                    zprojopt[i] = zprojmin[i];
                }
                ePopt = ePmin[k];
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
    gnuplot_set_ylabel(h3, (char*)"ePm [-]");
    gnuplot_plot_xy(h3, kvec, ePmin, kn, (char*)"", "points", "1", "1", 1);
    gnuplot_plot_xy(h3, &kvec[kopt], &ePmin[kopt], 1, (char*)"", "points", "3", "1", 4);

    //Plot the corresponding manifold leg in EM coordinates
    for(int i = 0; i <= Nman; i++)
    {
        for(int k = 0; k < 6; k++) yMan[k][i] = yTens[k][i][kopt];
        tMan[i] = tTens[i][kopt];
    }
    NCtoSYS_vec(yMan, tMan, yManEM, Nman, &SEML);
    gnuplot_plot_xyz(hem, yManEM[0], yManEM[1], yManEM[2], Nman+1, (char*)"", "lines", "1", "1", 4);

    //Plot the manifold leg in SEM coordinates
    NCEMmtoSEMm_vec(yMan, tMan, yManSEM, Nman, &SEML);
    gnuplot_plot_xyz(hse, yManSEM[0], yManSEM[1], yManSEM[2], Nman+1, (char*)"", "lines", "1", "1", 4);
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
    free_dmatrix(yManSEM, 0, 5, 0, Nman);
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

//---------------------------------------------------------------------------------------------------------------------------------------
//
//          Integration
//
//---------------------------------------------------------------------------------------------------------------------------------------
/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t;       //current state and time

    //Projection tools
    double ePm;
    int nreset, nt;

    //Plot
    double ti;

    //Change sign of step if necessary
    if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(orbit.driver);

        //Init the state & time
        for(int i = 0; i < 6; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 1;
        nreset = 1;

        //First position
        for(int k = 0; k < 6; k++) yNCE[k][0] = yv[k];
        tNCE[0] = t0;

        //Loop
        do
        {
            ti = t0 + (double) nt *(tf-t0)/N;
            if(isResetOn) status = gslc_proj_evolve(orbit, yv, &t, t0, ti, &ePm, &nreset, isResetOn);
            else  status = gsl_odeiv2_driver_apply (orbit.driver->d, &t, ti, yv);

            for(int k = 0; k < 6; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = ti;

            //Advance one step
            nt++;
        }
        while((nt<=N) && (status == 0) && (orbit.tproj > orbit.tprojmin));

        //If a new reset is necessary
        if (status == -2 && isResetOn)
        {
            cout << "Warning in trajectory_integration_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == -1)
        {
            cout << "Warning in trajectory_integration_variable_grid: the stepper went wrong. No output is produced.";
            return -1;
        }

    }
    while(status!= 0 && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return -1;
    }

    return 0;
}


/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N. Return the last position that is filled on the grid.
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------
    //Initialization
    //------------------------------------------
    int status;            //current status
    double yv[6], t;       //current state and time

    //Projection tools
    double ePm;
    int nreset, nt;

    //Change sign of step if necessary
    if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;

    //------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(orbit.driver);

        //Init the state & time
        for(int i = 0; i < 6; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 1;
        nreset = 1;

        //First step
        for(int k = 0; k < 6; k++) yNCE[k][0] = orbit.z0[k];
        tNCE[0] = t0;

        //Loop
        do
        {
            status = gslc_proj_step(orbit, yv, &t, t0, tf, &ePm, &nreset, isResetOn);
            for(int k = 0; k < 6; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = t;
            //Advance one step
            nt++;
        }
        while((t < tf) && (nt <= N) && (status == 0) && (orbit.tproj > orbit.tprojmin));

        if (status == -2 && isResetOn) //something went wrong during projection
        {
            cout << "Warning in trajectory_integration_variable_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

        if(status == -1) //something went wrong inside the stepper
        {
            cout << "Warning in trajectory_integration_variable_grid: the stepper went wrong. No output is produced.";
            return 0;
        }

        if(nt == N) //the maximum number of points is reached
        {
            cout << "Warning in trajectory_integration_variable_grid: the final time was not reached because the maximum number of points is reached." << endl;

        }

    }
    while(status!= 0 && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return -1;
    }

    return nt-1;
}
