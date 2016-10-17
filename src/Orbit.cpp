#include "Orbit.h"


//========================================================================================
// Constructors
//========================================================================================
/**
 *  \brief Constructor. Some values are hard coded, such as the default time interval of
 *         projection and the minimum time interval of projection
 **/
Orbit::Orbit(Invman const *invman_, QBCP_L const *qbcp_l_, OdeStruct *driver_, int ofts_order_,
                int ofs_order_, double t0_, double tf_):
//----------------------------------------------------------------------------------------
// Initialization list
//----------------------------------------------------------------------------------------
    order(ofts_order_),
    ofs_order(ofs_order_),
    reduced_nv(invman_->getRnv()),
    fwrk(invman_->getFwrk()),
    ofs(ofs_order_),
    tfx(tf_),
    t0x(t0_),
    tprojx(0.20*invman_->getCS()->us.T),   //The default tproj is set to (period of the model)/5
    tprojminx(1e-4*invman_->getCS()->us.T) //equal to 1e-4*period of the model
//----------------------------------------------------------------------------------------
// Body of the constructor
//----------------------------------------------------------------------------------------
{
    //------------------------------------------------------------------------------------
    // ePmax: maximum projection error allowed during the computation €€TODO
    //------------------------------------------------------------------------------------
    ePmaxx = (fwrk == F_SEM)? 5e-3:5e-5;

    //------------------------------------------------------------------------------------
    // Initial, final and current state
    //------------------------------------------------------------------------------------
    z0x  = dvector(0, 5);              //Initial position in NC coordinates dim = 6
    six  = dvector(0, reduced_nv-1);   //Initial RCM configuration dim = 4
    s0dx = dvector(0, 2*reduced_nv-1); //Initial position in CCM8 coordinates (real+imag part) dim = 8
    s0x  = dcvector(0,reduced_nv-1);   //Initial position in CCM4 coordinates (real+imag part) dim = 4
    xfx  = dvector(0, 5);              //Final position NC dim = 6

    //------------------------------------------------------------------------------------
    // Macro structures and objects
    //------------------------------------------------------------------------------------
    invman = invman_;
    qbcp_l = qbcp_l_;
    driver = driver_;
}


//========================================================================================
// Destructor
//========================================================================================
/**
 *  \brief Destructor.
 **/
Orbit::~Orbit()
{
    free_dvector(z0x,  0, 5);
    free_dvector(six,  0, reduced_nv-1);
    free_dvector(s0dx, 0, 2*reduced_nv-1);
    free_dvector(xfx,  0, 5);
    free_dcvector(s0x, 0, reduced_nv-1);
    free_ode_structure(driver);
}


//========================================================================================
// Update
//========================================================================================
/**
 *  \brief Update the initial conditions: new initial time is t0c, new initial RCM state
 *         is si.
 **/
void Orbit::update_ic(const double si[], double t0c)
{
    //------------------------------------------
    // 1. Update si
    //------------------------------------------
    for(int p = 0; p < reduced_nv; p++) this->six[p] = si[p];

    //------------------------------------------
    // 2. Update s0
    //------------------------------------------
    RCMtoCCM(si, s0x, reduced_nv);

    //------------------------------------------
    // 2. Update s0d
    //------------------------------------------
    RCMtoCCM8(si, s0dx, reduced_nv);

    //------------------------------------------
    // 4. Update z0
    //------------------------------------------
    //z0 = W(si, t0)
    invman->evalRCMtoNC(si, t0c, z0x, order, ofs_order);
}

void Orbit::update_ic(const double si[])
{
    //Same as before, using the time inside the object
    this->update_ic(si, t0x);
}

void Orbit::evalRCMtoNC(double const t, double z1[]) const
{
    this->invman->evalRCMtoNC(this->six, t, z1, this->order, this->ofs_order);
}

void Orbit::ccm8torcm(const double s0d[])
{
    double sit[reduced_nv];
    CCM8toRCM(s0d, sit, reduced_nv);

    for(int i = 0; i < reduced_nv; i++) six[i] = sit[i];
}


//========================================================================================
// Getters
//========================================================================================
const double Orbit::getN() const
{
    return this->invman->getN();
}

const double Orbit::getT0() const
{
    return t0x;
}

const double Orbit::getTf() const
{
    return tfx;
}

const double* Orbit::getZ0() const
{
    return z0x;
}

const double* Orbit::getSi() const
{
    return six;
}

const Invman* Orbit::getInvman() const
{
    return invman;
}

//========================================================================================
// Setters
//========================================================================================
/**
 *  \brief Set tf.
 **/
void Orbit::setTf(double tf)
{
    tfx = tf;
}

/**
 *  \brief Set t0.
 **/
void Orbit::setT0(double t0)
{
    t0x = t0;
}

/**
 *  \brief Set value at the ith dimension of si.
 **/
void Orbit::setSi(double value, int i)
{
    //------------------------------------------------------------------------------------
    // Check sizes
    //------------------------------------------------------------------------------------
    if(i >= reduced_nv)
    {
        cout << "Orbit::addSi. The desired dimension is out of scope" << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Set value at the ith dimensions
    //------------------------------------------------------------------------------------
    six[i] = value;
}


/**
 *  \brief Add value to the ith dimension of si.
 **/
void Orbit::addSi(double value, int i)
{
    //------------------------------------------------------------------------------------
    // Check sizes
    //------------------------------------------------------------------------------------
    if(i >= reduced_nv)
    {
        cout << "Orbit::addSi. The desired dimension is out of scope" << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Add value to the ith dimensions
    //------------------------------------------------------------------------------------
    six[i] += value;
}

//========================================================================================
// Integrate
//========================================================================================
/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary.
 *          This routine returns:
 *                  - 0  if everything went well
 *                  - -1 if the stepper went wrong (via the GSL routine
 *                     gsl_odeiv2_evolve_apply)
 *                  - -2 if the projection procedure went wrong (we are probably out of
 *                     the domain of convergence of the parameterization of the
 *                     center manifold).
 **/
int Orbit::gslc_proj_step(double yv[], double *t, double t0, double t1,
                          double *projdist, int *nreset, int isResetOn)
{
    //------------------------------------------------------------------------------------
    //Evolve one step of z(t)
    //------------------------------------------------------------------------------------
    //Check that the direction of integration is consistent with the time interval
    driver->h = (t1 >= *t)? fabs(driver->h):-fabs(driver->h);
    //Advance one step
    int status = gsl_odeiv2_evolve_apply(driver->e, driver->c, driver->s, &driver->sys, t, t1, &driver->h, yv);
    //Check the status
    if (status != 0)
    {
        cout << "error in Orbit::gslc_proj_step: integration of z(t) has gone wrong. break." << endl;
        return -1;
    }

    //------------------------------------------------------------------------------------
    //Projection if necessary, every "orbit.tprojx" intervals
    //------------------------------------------------------------------------------------
    if(isResetOn && fabs(*t-t0) > fabs(*nreset*tprojx))
    {
        double yvp[6], yvi[6];
        cdouble scp[reduced_nv];

        //--------------------------------------------------------------------------------
        //Projection tools
        //--------------------------------------------------------------------------------
        double omega1 = invman->getOmega1();
        double omega3 = invman->getOmega3();

        //--------------------------------------------------------------------------------
        // Projection on the center manifold
        //--------------------------------------------------------------------------------
        //Get the closest point on the center manifold
        NCprojCCM(yv, *t, invman->getN(), OFS_ORDER, invman->getMIcoc(), invman->getVcoc(), omega1, omega3, scp, reduced_nv);
        //Update the state
        invman->evalCCMtoNC(scp, *t, yvp, order, ofs_order);

        //For comparison
        for(int i = 0; i <6; i++) yvi[i] = yv[i];
        // Copy of yvp in current state
        for(int i=0; i<6; i++) yv[i]  = yvp[i];

        //--------------------------------------------------------------------------------
        // Get the current projection error
        //--------------------------------------------------------------------------------
        //Get the current error
        *projdist = DENorm(yvi, yvp, 6);

        if(*projdist > ePmaxx)
        {
            cout << "-----------------------------------------" << endl;
            cout << "Orbit::gslc_proj_step. Warning: Reset n° " << *nreset << endl;
            cout << "Error (NC)  = " << *projdist << endl;
            cout << "Error (SYS) = " << *projdist*invman->getCS()->gamma << endl;
            cout << "Error (km)  = " << *projdist*invman->getCS()->gamma*invman->getCS()->cr3bp.L << endl;
            cout << "-----------------------------------------" << endl;
            return -2;
        }

        //--------------------------------------------------------------------------------
        //Reset ode structure for next step
        //--------------------------------------------------------------------------------
        reset_ode_structure(driver);

        //--------------------------------------------------------------------------------
        //One additional reset
        //--------------------------------------------------------------------------------
        *nreset = *nreset +1;
    }


    return 0;
}

/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int Orbit::gslc_proj_evolve(double yv[], double *t, double t0, double t1,
                            double *projdist, int *nreset, int isResetOn)
{
    int status;
    do
    {
        status = this->gslc_proj_step(yv, t, t0, t1, projdist, nreset, isResetOn);
    }
    while(status == 0 && fabs(*t - t1) > 1e-16);

    return status;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid.
 *          return 0 if the integration went well, -1 otherwise.
 **/
int Orbit::traj_int_grid(double tfc, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------------
    int nvar = driver->dim;   //number of variables
    int status;               //current status
    double yv[nvar], t;       //current state and time

    //Projection tools
    double projdist;
    int nreset, nt;

    //Plot
    double ti;

    //------------------------------------------------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(driver);

        //Change sign of step if necessary
        if((tfc < t0x && driver->h>0)    || (tfc > t0x && driver->h<0)) driver->h *= -1;
        if((tfc < t0x && driver->d->h>0) || (tfc > t0x && driver->d->h<0)) driver->d->h *= -1;

        //Init the state & time
        for(int i = 0; i < nvar; i++) yv[i] = z0x[i];
        t = t0x;

        //Init the indexes that evolve along with the state
        nt = 1;
        nreset = 1;

        //First position
        for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];
        tNCE[0] = t0x;

        //Loop
        do
        {
            ti = t0x + (double) nt *(tfc-t0x)/(1.0*N);
            if(isResetOn) status = this->gslc_proj_evolve(yv, &t, t0x, ti, &projdist, &nreset, isResetOn);
            else  status = gsl_odeiv2_driver_apply (driver->d, &t, ti, yv);

            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = ti;

            //Advance one step
            nt++;
        }
        while((nt<=N) && (status == 0) && (tprojx > tprojminx));

        //If a new reset is necessary
        if (status == -2 && isResetOn)
        {
            cout << "Orbit::traj_int_grid. Warning: the interval of projection" << endl;
            cout << "has to be reduced: ";
            cout << setprecision(3) << "tproj : " << tprojx << " -> ";
            tprojx *= 0.5;
            cout << tprojx << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == -1)
        {
            cout << "Orbit::traj_int_grid. Warning: the stepper went wrong." << endl;
            cout << "No output is produced.";
            return -1;
        }

    }
    while(status!= 0 && tprojx > tprojminx);

    if(tprojx < tprojminx)
    {
        cout << "Orbit::traj_int_grid. Error: the interval of projection is too small." << endl;
        return -1;
    }

    return 0;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid.
 *          return 0 if the integration went well, -1 otherwise.
 **/
int Orbit::traj_int_var_grid(double tfc, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------------
    int nvar = driver->dim;   //number of variables
    int status;               //current status
    double yv[nvar], t;       //current state and time

    //Projection tools
    double projdist;
    int nreset, nt;

    //------------------------------------------------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(driver);

        //Change sign of step if necessary
        if((tfc < t0x && driver->h>0)    || (tfc > t0x && driver->h<0)) driver->h *= -1;
        if((tfc < t0x && driver->d->h>0) || (tfc > t0x && driver->d->h<0)) driver->d->h *= -1;

        //Init the state & time
        for(int i = 0; i < nvar; i++) yv[i] = z0x[i];
        t = t0x;

        //Init the indexes that evolve along with the state
        nt = 1;
        nreset = 1;

        //First position
        for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];
        tNCE[0] = t0x;

        //Loop
        do
        {
            status = this->gslc_proj_step(yv, &t, t0x, tfc, &projdist, &nreset, isResetOn);
            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = t;
            //Advance one step
            nt++;

        }
        while((t != tfc) && (nt<=N) && (status == 0) && (tprojx > tprojminx));

        //If a new reset is necessary
        if (status == -2 && isResetOn)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the interval of projection" << endl;
            cout << "has to be reduced: ";
            cout << setprecision(3) << "tproj : " << tprojx << " -> ";
            tprojx *= 0.5;
            cout << tprojx << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == -1)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the stepper went wrong." << endl;
            cout << "No output is produced.";
            return -1;
        }

        //the maximum number of points is reached
        if(nt == N)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the final time was not reached because the maximum number of points is reached." << endl;
        }

    }
    while(status!= 0 && tprojx > tprojminx);

    if(tprojx < tprojminx)
    {
        cout << "Orbit::traj_int_var_grid. Error: the interval of projection is too small." << endl;
        return -1;
    }

    return 0;
}

//========================================================================================
//          Projection on (un)stable manifold
//========================================================================================
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given
 *        by this->invman, seen as a center manifold in graph style. If the graph style
 *        is not detected, a warning message is displayed and nothing is done.
 **/
void Orbit::NCprojCCMtoCM(double *yv, double tv, double sti[5])
{
    //------------------------------------------------------------------------------------
    //Check that the graph style is used in this->invman
    //------------------------------------------------------------------------------------
    if(invman->getPmType() != PMS_GRAPH)
    {
        cout << "Orbit::NCprojCCMtoCM. Impossible to use if the invariant manifold ";
        cout << "of the current orbit is not based on the graph style. return." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Projection tools
    //------------------------------------------------------------------------------------
    double omega1 = invman->getOmega1();
    double omega3 = invman->getOmega3();

    //------------------------------------------------------------------------------------
    //Get the closest point on the center manifold, in scp[4]
    //------------------------------------------------------------------------------------
    cdouble scp[4];
    NCprojCCM(yv, tv, invman->getN(), OFS_ORDER, invman->getMIcoc(),
              invman->getVcoc(), omega1, omega3, scp, 4);

    //------------------------------------------------------------------------------------
    //Get the correspondance in RCM coordinates
    //------------------------------------------------------------------------------------
    CCMtoRCM(scp, sti, 4);
}


//========================================================================================
//          Orbit on a grid
//========================================================================================
/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initial conditions st0, on the center manifold CM.
 **/
int oo_gridOrbit(double st0[], double t0, double tf, double dt, Invman &invman)
{
    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    //------------------
    // ODE
    //------------------
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Init ode structure
    init_ode_structure(&driver, T, T_root, 6, qbfbp_vfn_novar, &SEML);

    //------------------
    //Orbit
    //------------------
    Orbit orbit(&invman, &SEML, &driver, OFTS_ORDER, OFS_ORDER, t0, tf);

    //------------------
    //Orbit IC
    //------------------
    orbit.update_ic(st0, t0);

    //------------------------------------------
    //Plotting parameters
    //------------------------------------------
    int N         = floor(fabs(tf - t0)/dt);
    double **yNCE = dmatrix(0, 5, 0, N);
    double **ySYS = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    int status = orbit.traj_int_grid(tf, yNCE, tNCE, N, true);

    //------------------------------------------
    //To SYS coordinates
    //------------------------------------------
    NCtoSYS_vec(yNCE, tNCE, ySYS,N, &SEML);

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    gnuplot_setstyle(h1, (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, ySYS[0], ySYS[1], ySYS[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);


    //------------------------------------------------------------------------------------
    //Save in file
    //------------------------------------------------------------------------------------
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
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(ySYS, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    return status;

    return 0;
}
