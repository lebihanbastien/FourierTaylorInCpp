#include "Orbit.h"
#include "lenconref.h"

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
    tprojminx(1e-6*invman_->getCS()->us.T) //equal to 1e-4*period of the model
//----------------------------------------------------------------------------------------
// Body of the constructor
//----------------------------------------------------------------------------------------
{
    //------------------------------------------------------------------------------------
    // ePmax: maximum projection error allowed during the computation €€TODO
    //------------------------------------------------------------------------------------
    ePmaxx = (fwrk == F_SEM)? 5e-3:1e-4;

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

    //------------------------------------------
    // 5. Update t0
    //------------------------------------------
    t0x = t0c;
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

const double Orbit::getSi(int dim) const
{
    if(dim < reduced_nv) return six[dim];
    else return six[0];
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
 *                  - FTC_SUCCESS = GSL_SUCCESS = 0  if everything went well
 *                  - ORBIT_EINT = -1 if the stepper went wrong (via the GSL routine
 *                     gsl_odeiv2_evolve_apply)
 *                  - ORBIT_EPROJ = -2 if the projection procedure went wrong
 *                     (we are probably out of the domain of convergence of the
 *                     parameterization of the center manifold).
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
    if (status != GSL_SUCCESS)
    {
        cout << "error in Orbit::gslc_proj_step: integration of z(t) has gone wrong. break." << endl;
        return ORBIT_EINT;
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
            return ORBIT_EPROJ;
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

    return FTC_SUCCESS;
}

/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary.
 *          This routine returns:
 *                  - FTC_SUCCESS = GSL_SUCCESS = 0   if everything went well
 *                  - ORBIT_EINT = -1 if the stepper gslc_proj_step went wrong
 *                     (via the GSL routine gsl_odeiv2_evolve_apply)
 *                  - ORBIT_EPROJ = -2 if the projection procedure went wrong
 *                     (we are probably out of the domain of convergence of the
 *                     parameterization of the center manifold).
 **/
int Orbit::gslc_proj_evolve(double yv[], double *t, double t0, double t1,
                            double *projdist, int *nreset, int isResetOn)
{
    int status;
    do
    {
        status = this->gslc_proj_step(yv, t, t0, t1, projdist, nreset, isResetOn);
    }
    while(status == FTC_SUCCESS && fabs(*t - t1) > 1e-16);

    return status;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid.
 *          This routine returns:
 *                  - FTC_SUCCESS = GSL_SUCCESS = 0   if everything went well
 *                  - ORBIT_EINT = -1 if the stepper gslc_proj_step went wrong
 *                     (via the GSL routine gsl_odeiv2_evolve_apply)
 *                  - The index corresponding to the last good state
 *                    if the projection procedure went wrong several
 *                    times inside gslc_proj_step so that the interval of projection
 *                    is reduced below its minimum value.
 *                    (we are probably out of the domain of convergence of the
 *                    parameterization of the center manifold). In the (improbable but
 *                    still possible case that this index is exactly zero, ORBIT_EPROJ
 *                    is returned.
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
        //--------------------------------------------------------------------------------
        // Initialize the integration
        //--------------------------------------------------------------------------------
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

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
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

        //--------------------------------------------------------------------------------
        //Checks
        //--------------------------------------------------------------------------------
        //If a new reset is necessary
        if (status == ORBIT_EPROJ && isResetOn)
        {
            cout << "Orbit::traj_int_grid. Warning: the interval of projection" << endl;
            cout << "has to be reduced: ";
            cout << setprecision(3) << "tproj : " << tprojx << " -> ";
            tprojx *= 0.5;
            cout << tprojx << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == ORBIT_EINT)
        {
            cout << "Orbit::traj_int_grid. Warning: the stepper went wrong." << endl;
            cout << "No output is produced.";
            return ORBIT_EINT;
        }

    }
    while(status!= FTC_SUCCESS && tprojx > tprojminx);

    //------------------------------------------------------------------------------------
    //If the projection procedure failed enough time to reach the minimum projection time,
    //the last "good" indix is returned. If this good indix is exactly zero (which might
    //happen, if the projection procedure is really bad), we have a problem, because then
    //it is equal to FTC_SUCCESS = 0. Since FTC_SUCCESS is set equal to zero
    //to match GSL standard (GSL_SUCCESS := 0), we cannot change its value now.
    //So, if the last indix is exactly zero, we output ORBIT_EPROJ
    //------------------------------------------------------------------------------------
    if(tprojx < tprojminx)
    {
        cout << "Orbit::traj_int_grid. Error: the interval of projection is too small." << endl;
        cout << "The index corresponding to the last good state is returned if it is not zero." << endl;
        if(nt-1 > 0) return nt-1;
        else return ORBIT_EPROJ;
    }

    return FTC_SUCCESS;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid given by tgridNCE.
 *          This routine returns:
 *                  - FTC_SUCCESS = GSL_SUCCESS = 0   if everything went well
 *                  - ORBIT_EINT = -1 if the stepper gslc_proj_step went wrong
 *                     (via the GSL routine gsl_odeiv2_evolve_apply)
 *                  - The index corresponding to the last good state
 *                    if the projection procedure went wrong several
 *                    times inside gslc_proj_step so that the interval of projection
 *                    is reduced below its minimum value.
 *                    (we are probably out of the domain of convergence of the
 *                    parameterization of the center manifold). In the (improbable but
 *                    still possible case that this index is exactly zero, ORBIT_EPROJ
 *                    is returned.
 **/
int Orbit::traj_int_grid(double **yNCE, double *tgridNCE, int N, int isResetOn)
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

    //Initial & final time
    double t0c = tgridNCE[0];
    double tfc = tgridNCE[N];

    //Plot
    double ti;

    //------------------------------------------------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------------------------------------------------
    do
    {
        //--------------------------------------------------------------------------------
        // Initialize the integration
        //--------------------------------------------------------------------------------
        //Reset ode structure.
        reset_ode_structure(driver);

        //Change sign of step if necessary
        if((tfc < t0c && driver->h>0)    || (tfc > t0c && driver->h<0)) driver->h *= -1;
        if((tfc < t0c && driver->d->h>0) || (tfc > t0c && driver->d->h<0)) driver->d->h *= -1;

        //Init the state & time
        for(int i = 0; i < nvar; i++) yv[i] = z0x[i];
        t = t0c;

        //Init the indexes that evolve along with the state
        nt = 1;
        nreset = 1;

        //First position
        for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        do
        {
            ti = tgridNCE[nt];
            if(isResetOn) status = this->gslc_proj_evolve(yv, &t, t0c, ti, &projdist, &nreset, isResetOn);
            else  status = gsl_odeiv2_driver_apply (driver->d, &t, ti, yv);

            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            //Advance one step
            nt++;
        }
        while((nt<=N) && (status == 0) && (tprojx > tprojminx));

        //--------------------------------------------------------------------------------
        //Checks
        //--------------------------------------------------------------------------------
        //If a new reset is necessary
        if (status == ORBIT_EPROJ && isResetOn)
        {
            cout << "Orbit::traj_int_grid. Warning: the interval of projection" << endl;
            cout << "has to be reduced: ";
            cout << setprecision(3) << "tproj : " << tprojx << " -> ";
            tprojx *= 0.5;
            cout << tprojx << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == ORBIT_EINT)
        {
            cout << "Orbit::traj_int_grid. Warning: the stepper went wrong." << endl;
            cout << "No output is produced.";
            return ORBIT_EINT;
        }

    }
    while(status!= FTC_SUCCESS && tprojx > tprojminx);

    //------------------------------------------------------------------------------------
    //If the projection procedure failed enough time to reach the minimum projection time,
    //the last "good" indix is returned. If this good indix is exactly zero (which might
    //happen, if the projection procedure is really bad), we have a problem, because then
    //it is equal to FTC_SUCCESS = 0. Since FTC_SUCCESS is set equal to zero
    //to match GSL standard (GSL_SUCCESS := 0), we cannot change its value now.
    //So, if the last indix is exactly zero, we output ORBIT_EPROJ
    //------------------------------------------------------------------------------------
    if(tprojx < tprojminx)
    {
        cout << "Orbit::traj_int_grid. Error: the interval of projection is too small." << endl;
        cout << "The index corresponding to the last good state is returned if it is not zero." << endl;
        if(nt-1 > 0) return nt-1;
        else return ORBIT_EPROJ;
    }

    return FTC_SUCCESS;
}


/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid.
 *          This routine returns:
 *                  - The last indix that contains data in the data vectors.
 *                  - ORBIT_EINT = -1 if the stepper gslc_proj_step went wrong
 *                     (via the GSL routine gsl_odeiv2_evolve_apply)
 *                  - The index corresponding to the last good state
 *                    if the projection procedure went wrong several
 *                    times inside gslc_proj_step so that the interval of projection
 *                    is reduced below its minimum value.
 *                    (we are probably out of the domain of convergence of the
 *                    parameterization of the center manifold). In the (improbable but
 *                    still possible case that this index is exactly zero, ORBIT_EPROJ
 *                    is returned.
 *
 *          Note that if the maximum number of points N is reached before the end of the
 *          integration procedure, the last indix is still returned as an output, but
 *          a warning is displayed to the user.
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
        //--------------------------------------------------------------------------------
        // Initialize the integration
        //--------------------------------------------------------------------------------
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

        //--------------------------------------------------------------------------------
        //Loop
        //--------------------------------------------------------------------------------
        do
        {
            status = this->gslc_proj_step(yv, &t, t0x, tfc, &projdist, &nreset, isResetOn);
            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = t;
            //Advance one step
            nt++;

        }
        while((t != tfc) && (nt<=N) && (status == 0) && (tprojx > tprojminx));

        //--------------------------------------------------------------------------------
        //Checks
        //--------------------------------------------------------------------------------
        //If a new reset is necessary
        if (status == ORBIT_EPROJ && isResetOn)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the interval of projection" << endl;
            cout << "has to be reduced: ";
            cout << setprecision(3) << "tproj : " << tprojx << " -> ";
            tprojx *= 0.5;
            cout << tprojx << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == ORBIT_EINT)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the stepper went wrong." << endl;
            cout << "No output is produced.";
            return ORBIT_EINT;
        }

        //The maximum number of points is reached. The loop will end and the last indix
        //N-1 will be output.
        if(nt == N)
        {
            cout << "Orbit::traj_int_var_grid. Warning: the final time was not reached because the maximum number of points is reached." << endl;
        }

    }
    while(status!= FTC_SUCCESS && tprojx > tprojminx);

    //------------------------------------------------------------------------------------
    //If the projection procedure failed enough time to reach the minimum projection time,
    //the last "good" indix is returned. If this good indix is exactly zero (which might
    //happen, if the projection procedure is really bad), we have a problem, because then
    //it is equal to FTC_SUCCESS = 0. Since FTC_SUCCESS is set equal to zero
    //to match GSL standard (GSL_SUCCESS := 0), we cannot change its value now.
    //So, if the last indix is exactly zero, we output ORBIT_EPROJ
    //------------------------------------------------------------------------------------
    if(tprojx < tprojminx)
    {
        cout << "Orbit::traj_int_grid. Error: the interval of projection is too small." << endl;
        cout << "The index corresponding to the last good state is returned if it is not zero." << endl;
        if(nt-1 > 0) return nt-1;
        else return ORBIT_EPROJ;
    }

    //------------------------------------------------------------------------------------
    //If the whole integration went well, the last indix is returned.
    //------------------------------------------------------------------------------------
    return nt-1;
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
int oo_gridOrbit(double st0[], double t0, double tf, double dt)
{
    //====================================================================================
    // Initialization
    //====================================================================================
    //------------------------------------------------------------------------------------
    // Invariant manifold
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, SEML.cs);
    int ncs = invman.getNCS();
    //int dcs = default_coordinate_system(ncs);

    //------------------------------------------------------------------------------------
    // ODE
    //------------------------------------------------------------------------------------
    //Hard precision
    Config::configManager().C_PREC_HARD();

    //Driver
    OdeStruct driver;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    int coll;
    OdeParams odeParams(&coll, &SEML);
    //Init ode structure
    init_ode_structure(&driver, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------------------------------------------------------------------------
    //Orbit
    //------------------------------------------------------------------------------------
    Orbit orbit(&invman, &SEML, &driver, OFTS_ORDER, OFS_ORDER, t0, tf);

    //Orbit IC
    orbit.update_ic(st0, t0);

    //------------------------------------------------------------------------------------
    //Plotting parameters
    //------------------------------------------------------------------------------------
    int N         = floor(fabs(tf - t0)/dt);
    double **yNCE = dmatrix(0, 5, 0, N);
    double **ySYS = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);
    double *tnSYS = dvector(0, N);
    double *dHSYS = dvector(0, N);

    //====================================================================================
    //Integration
    //====================================================================================
    int status = orbit.traj_int_grid(tf, yNCE, tNCE, N, true);

    //To SYS coordinates
    NCtoSYS_vec(yNCE, tNCE, ySYS, N, &SEML);

    //====================================================================================
    //Hamiltonian:
    //====================================================================================
    double y[6], H0 = 0.0;

    //Loop on all positions
    for(int k= 0; k <= N; k++)
    {
        for(int i = 0; i <6; i++) y[i] = ySYS[i][k];
        if(k == 0)
        {
             dHSYS[k] = 0.0;
             H0 = qbcp_H(tNCE[k], y, &SEML);
        }else{
            dHSYS[k] = qbcp_H(tNCE[k], y, &SEML) - H0;
        }

        //Corresponding normalized times
        tnSYS[k] = tNCE[k]/SEML.us.T;
    }


    //Plotting
    gnuplot_ctrl  *h3;
    h3 = gnuplot_init();
    gnuplot_cmd(h3, "set logscale y");
    gnuplot_setstyle(h3,   (char*)"lines");
    gnuplot_set_xlabel(h3, (char*)"t [x T]");
    gnuplot_set_ylabel(h3, (char*)"H [-]");
    gnuplot_plot_xy(h3, tnSYS, dHSYS, N+1, (char*)"H(t) - H(0)", "lines", "1", "1", 1);


    //====================================================================================
    //Save in file
    //====================================================================================
    writeOrbit(tNCE, yNCE, OFTS_ORDER, floor(fabs(st0[0])), N+1, TYPE_ORBIT);


    //====================================================================================
    //Plotting
    //====================================================================================
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init();

    gnuplot_setstyle(h1,   (char*)"lines");
    gnuplot_set_xlabel(h1, (char*)"x [-]");
    gnuplot_set_ylabel(h1, (char*)"y [-]");
    gnuplot_plot_xyz(h1, yNCE[0], yNCE[1], yNCE[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);


    //User check
    pressEnter(true);


    //====================================================================================
    //Correction procedure
    //====================================================================================
    gnuplot_ctrl  *h2;
    h2 = gnuplot_init();

    //------------------------------------------------------------------------------------
    // Coordinates for the correction
    //------------------------------------------------------------------------------------
    int new_ncd = VEM;
    int new_dcs = default_coordinate_system(new_ncd);

    //------------------------------------------------------------------------------------
    // Init
    //------------------------------------------------------------------------------------
    double dtd = 2*86400*2*M_PI/SEML.cs_em.cr3bp.T; //every two days
    N          = floor(fabs(tf - t0)/dtd);

    //Containers
    double **yt = dmatrix(0, 41, 0, N);
    double **yc = dmatrix(0, 41, 0, N);
    double *tc  = dvector(0, N);

    // Time grid
    for(int k = 0; k <= N; k++) tc[k] = t0 + (double) k *(tf-t0)/(1.0*N);

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    status = orbit.traj_int_grid(tf, yt, tc, N, true);
    qbcp_coc_vec(yt, tc, yc, N, ncs, new_ncd);


    //------------------------------------------------------------------------------------
    //Trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    int mPlot = 100, ode78coll;
    double** ymc        = dmatrix(0, 5, 0, mPlot);
    double* tmc         = dvector(0, mPlot);
    double** ymc_v      = dmatrix(0, 5, 0, mPlot*N);
    double* tmc_v       = dvector(0, mPlot*N);
    double yv[6];
    for(int k = 0; k < N; k++)
    {
        //--------------------------------------------------------------------------------
        //Integration segment by segment
        //--------------------------------------------------------------------------------
        for(int i = 0; i < 6; i++) yv[i] = yc[i][k];
        //cout << tc[k] << " " << yc[0][k] << endl;

        ode78(ymc, tmc, &ode78coll, tc[k], tc[k+1], yv, 6, mPlot, new_dcs, new_ncd, new_ncd);

        //--------------------------------------------------------------------------------
        //Store
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*N+1, "final trajectory", "lines", "1", "1", 2);

    //User check
    pressEnter(true);

    //------------------------------------------------------------------------------------
    //Diff Corr!
    //------------------------------------------------------------------------------------
    multiple_shooting_direct(yc, tc, yc, tc, 42, N, new_ncd, true, h2);

    //------------------------------------------------------------------------------------
    //Trajectory on lines, segment by segment
    //------------------------------------------------------------------------------------
    for(int k = 0; k < N; k++)
    {
        //--------------------------------------------------------------------------------
        //Integration segment by segment
        //--------------------------------------------------------------------------------
        for(int i = 0; i < 6; i++) yv[i] = yc[i][k];
        ode78(ymc, tmc, &ode78coll, tc[k], tc[k+1], yv, 6, mPlot, new_dcs, new_ncd, new_ncd);

        //--------------------------------------------------------------------------------
        //Store
        //--------------------------------------------------------------------------------
        for(int p = 0; p <= mPlot; p++)
        {
            for(int i = 0; i < 6; i++) ymc_v[i][k*mPlot + p] = ymc[i][p];
            tmc_v[k*mPlot + p] = tmc[p];
        }
    }

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    //Plot on h2
    gnuplot_plot_X(h2, ymc_v, mPlot*N+1, "final trajectory", "lines", "1", "1", 2);

    //User check
    pressEnter(true);
    gnuplot_close(h1);
    gnuplot_close(h2);

    //====================================================================================
    //Free
    //====================================================================================
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(ySYS, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    return status;

    return FTC_SUCCESS;
}
