#include "ode.h"

void update_ode_structure(OdeStruct *ode_s,
                          const gsl_odeiv2_step_type *T,
                          gsl_odeiv2_step *s,
                          gsl_odeiv2_control *c,
                          gsl_odeiv2_evolve *e,
                          gsl_odeiv2_system sys,
                          gsl_odeiv2_driver * d,
                          gsl_root_fsolver *s_root,
                          double eps_root,
                          double eps_diff)
{
    ode_s->dim = d->sys->dimension;
    ode_s->T = T;
    ode_s->s = s;
    ode_s->c = c;
    ode_s->e = e;
    ode_s->sys = sys;
    ode_s->d = d;
    ode_s->s_root = s_root;

    ode_s->eps_root = eps_root;
    ode_s->eps_diff = eps_diff;
    ode_s->h = d->h;

    //Integration tolerances are set to -1 because they are hidden in d and e (be careful during initialization!)
    ode_s->eps_int_abs = -1;
    ode_s->eps_int_rel = -1;
}

void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        double eps_int_abs,
                        double eps_int_rel,
                        double eps_root,
                        double eps_diff,
                        size_t dim,
                        double h,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params),
                        void *odeParams)
{
    //Precisions and initial step
    ode_s->eps_root = eps_root;
    ode_s->eps_diff = eps_diff;
    ode_s->h = h;
    ode_s->eps_int_abs = eps_int_abs;
    ode_s->eps_int_rel = eps_int_rel;
    ode_s->dim = dim;

    //Stepper
    ode_s->T = T;
    ode_s->s = gsl_odeiv2_step_alloc(T, dim);
    if (ode_s->s == NULL)
    {
        cout << "failed to allocate stepper object in init_ode_structure" << endl;
    }

    //Control
    ode_s->c = gsl_odeiv2_control_y_new(eps_int_abs, eps_int_rel);
    if (ode_s->c == NULL)
    {
        cout << "failed to allocate control object in init_ode_structure"  << endl;
    }

    //Evolution
    ode_s->e = gsl_odeiv2_evolve_alloc(dim);
    if (ode_s->e == NULL)
    {
        cout << "failed to allocate evolution object in init_ode_structure"  << endl;
    }


    //System
    ode_s->sys.function  = func;
    ode_s->sys.jacobian  = jacobian;
    ode_s->sys.dimension = dim;
    ode_s->sys.params    = (void *) odeParams;

    ode_s->d = gsl_odeiv2_driver_alloc_y_new (&ode_s->sys, ode_s->T, ode_s->h, ode_s->eps_int_abs, ode_s->eps_int_rel);
    if (ode_s->d == NULL)
    {
        cout << "failed to allocate driver object in init_ode_structure"  << endl;
    }

    ode_s->s_root = gsl_root_fsolver_alloc (T_root);
    if (ode_s->s_root == NULL)
    {
        cout << "failed to allocate root_fsolver object in init_ode_structure"  << endl;
    }

}

void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        size_t dim,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        void *odeParams)
{
    init_ode_structure(ode_s, T, T_root,
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL(),
                       Config::configManager().G_PREC_ROOT(),
                       Config::configManager().G_PREC_DIFF(),
                       dim,
                       Config::configManager().G_PREC_HSTART(),
                       func, NULL, odeParams);
    ode_s->warnings=1;
}

void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        size_t dim,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        void *odeParams,
                        int warnings)
{
    init_ode_structure(ode_s, T, T_root,
                       Config::configManager().G_PREC_ABS(),
                       Config::configManager().G_PREC_REL(),
                       Config::configManager().G_PREC_ROOT(),
                       Config::configManager().G_PREC_DIFF(),
                       dim,
                       Config::configManager().G_PREC_HSTART(),
                       func, NULL, odeParams);
    ode_s->warnings=warnings;
}

void reset_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_reset(ode_s->s);
    gsl_odeiv2_evolve_reset(ode_s->e);
    gsl_odeiv2_driver_reset(ode_s->d);
}


void free_ode_structure(OdeStruct *ode_s)
{
    gsl_odeiv2_step_free(ode_s->s);
    gsl_odeiv2_evolve_free(ode_s->e);
    gsl_odeiv2_driver_free(ode_s->d);
}


//========================================================================================
// ODEZERO functions & structures
//========================================================================================
int custom_odezero_2(double y[],
                     double** ye,
                     double *tcross,
                     double t0,
                     double t1,
                     OdeStruct *ode_s,
                     struct value_function fvalue)
{
    //------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------
    int last_indix;                     //indix of the final position to be returned
    int events = 0;                     //number of events during integration
    int dim = ode_s->d->sys->dimension; //system dimension
    int i, k;                           //loop parameter
    int status;                         //status for GSL function integer output
    double previous_yv[dim], yv[dim];   //for integration purposes, copy of initial values
    double ys[dim];                     //for integration purposes
    double previous_t, t;               //times
    //Root finding
    gsl_function F;                     //called function in the root finding routine
    struct ode_params params;           //params for F
    double t_low, t_high;               //times for root bracketing
    //Loop parameters
    double r;
    double fy;
    int iter;
    //Event condition
    double yv_mat[dim][MAX_EVENTS];
    double t_mat[2][MAX_EVENTS];
    double previous_s;
    double new_s;

    //Copy of initial values for integration purposes
    for(i=0; i<dim; i++) yv[i] = y[i];

    //------------------------------------------------------------------------------
    //Reaching y=0
    //------------------------------------------------------------------------------
    //Previous value of yv an t are kept
    t = t0;
    previous_t = t;
    for(i=0; i<dim; i++) previous_yv[i] = yv[i];


    //Stepping until y=0 is crossed
    do
    {
        //Updating the previous y and t values
        previous_t = t;
        for(i=0; i<dim; i++) previous_yv[i] = yv[i];

        //Evolve one step
        status = gsl_odeiv2_evolve_apply (ode_s->e, ode_s->c, ode_s->s, &ode_s->sys, &t, t1, &ode_s->h, yv);

        //Break if evolution has gone wrong
        if (status != GSL_SUCCESS) break;

        //event detection
        previous_s = fvalue.value(previous_t,previous_yv,fvalue.val_par);
        new_s      = fvalue.value(t,yv,fvalue.val_par);

        if(previous_s*new_s<0)  //a zero of the value function has been crossed
        {
            switch(fvalue.val_par->direction) //switch on the direction value
            {
            case -1:  //decreasing zeros
                if(new_s < previous_s)
                {
                    //Storage of values
                    for(i=0; i<dim; i++) yv_mat[i][events] = previous_yv[i];
                    t_mat[0][events] = previous_t;
                    t_mat[1][events] = t;
                    //Add one event
                    events++;
                }
                break;
            case 1:  //increasing zeros
                if(new_s > previous_s)
                {
                    //Storage of values
                    for(i=0; i<dim; i++) yv_mat[i][events] = previous_yv[i];
                    t_mat[0][events] = previous_t;
                    t_mat[1][events] = t;
                    //Add one event
                    events++;
                }
                break;
            case 0:  //all zeros
                //Storage of values
                for(i=0; i<dim; i++) yv_mat[i][events] = previous_yv[i];
                t_mat[0][events] = previous_t;
                t_mat[1][events] = t;
                //Add one event
                events++;
                break;
            default: //all zeros
                //Storage of values
                for(i=0; i<dim; i++) yv_mat[i][events] = previous_yv[i];
                t_mat[0][events] = previous_t;
                t_mat[1][events] = t;
                //Add one event
                events++;
                break;
            }
        }

        //Termination detection
        if(events>=fvalue.val_par->max_events)
        {
            //mexPrintf("Maximum number of events reached. break.\n");
            break;
        }

    }
    while(fabs(t)<fabs(t1));

    if(fabs(t)>=fabs(t1))
    {
        //mexPrintf("custom_odezero_2. Final time was reached, last state is returned.\n");
        for(i=0; i<dim; i++){
            ye[i][events] = yv[i];
        }
        last_indix = events;
        tcross[events] = t;
    }
    else last_indix = events-1;

    //------------------------------------------------------------------------------------
    //Root finding for all events
    //------------------------------------------------------------------------------------
    if(events==0)
    {
        //mexPrintf("No events was found, last state is returned. \n");
        //Storage in ye
        for(i=0; i<dim; i++) ye[i][0] = yv[i];
        *tcross = t;
        last_indix = 0;

    }
    else
    {
        for(k=0; k<events; k++)
        {

            //Copy of initial values for storage in ye after the root finding
            for(i=0; i<dim; i++)
            {
                yv[i] = yv_mat[i][k];
            }


            //New start
            for(i=0; i<dim; i++) ys[i] = yv_mat[i][k];

            //Initialization of the parameters for F
            params.t0 = t_mat[0][k];       //new initial time is previous_t - epsilon
            params.y0 = ys;
            params.d  = ode_s->d;
            params.fvalue = fvalue;
            params.nvar = dim;

            //Initialization of F
            F.function = &cr3bp_event;
            F.params   = &params;

            //Bracketing the root
            t_low  = (t_mat[0][k] < t_mat[1][k])? t_mat[0][k] : t_mat[1][k];
            t_high = (t_mat[0][k] < t_mat[1][k])? t_mat[1][k] : t_mat[0][k];

            //Setting the solver
            status = gsl_root_fsolver_set (ode_s->s_root, &F, t_low, t_high);

            //Loop
            iter = 0;
            do
            {
                status = gsl_root_fsolver_iterate (ode_s->s_root);         //updating the solver
                r = gsl_root_fsolver_root (ode_s->s_root);                 //updating the root
                previous_t = gsl_root_fsolver_x_lower (ode_s->s_root);     //updating t_low
                t = gsl_root_fsolver_x_upper (ode_s->s_root);              //updating t_high

                //Checking convergence
                fy = cr3bp_event (r, &params);
                status = gsl_root_test_residual (fy , ode_s->eps_root);

            }
            while (status == GSL_CONTINUE && (++iter)<100);

            if(iter>=100)
            {
                printf("custom_odezero_2. Warning: number of iter max exceeded in custom_odezero. Last precision is %5.5e.\n", fy);
                //return GSL_FAILURE;
            }

            //Updating the crossing time
            tcross[k] = r;
            //Updating the outputs (state + STM)
            t = t_mat[0][k];
            status = gsl_odeiv2_driver_apply(ode_s->d, &t, tcross[k], yv);  //updating y  //TO BE CHANGED: INTEGRATION FROM PREVIOUS_YV !!
            //Storage in ye
            for(i=0; i<dim; i++) ye[i][k] = yv[i];


            //fy = fvalue.value(t,yv,fvalue.val_par).val;
            //mexPrintf("Obtained precision: fy = %5.5e.\n", fy);

            //----------------------------------------------------------------------------
            //Reset
            //----------------------------------------------------------------------------
            reset_ode_structure(ode_s);
            //gsl_odeiv2_step_reset(s);
            //gsl_odeiv2_evolve_reset(e);
            //gsl_odeiv2_driver_reset(d);
        }

    }
    return last_indix;

}

double linear_intersection(double t, double yv[], void *params)
{
    struct value_params *p = (struct value_params *) params;
    double val =  yv[p->dim] - p->value;
    return val;
}


double radius_intersection(double t, double yv[], void *params)
{
    struct value_params *p = (struct value_params *) params;

    //Coordinate relative to the center
    double xl = yv[0] - p->center[0];
    double yl = yv[1] - p->center[1];
    double zl = yv[2] - p->center[2];

    double val =  sqrt(xl*xl + yl*yl + zl*zl) - p->value;
    return val;
}

double angle_intersection(double t, double yv[], void *params)
{
    struct value_params *p = (struct value_params *) params;

    //Coordinate relative to the center
    double xl = yv[0] - p->center[0];
    double yl = yv[1] - p->center[1];

    double val;
    if(fabs(p->value) == M_PI/2)
    {
        val = xl;

    }else
    {
        val = yl - tan(p->value)*xl;
    }

    return val;
}

double null_flight_path_angle(double t, double yv[], void *params)
{
    struct value_params *p = (struct value_params *) params;

    //Coordinates relative to the center:
    //Position (x) and velocity (v)
    double x[3], v[3];
    for(int i = 0; i<3; i++)
    {
      x[i] = yv[i]-p->center[i];
      v[i] = yv[i+3];
    }

    //Zero function:  value = - vrt * xrt - event.value;
    double val = 0.0;
    for(int i = 0; i<3; i++) val += -v[i]*x[i];
    val -= p->value;

    return val;
}


/**
 * \brief return the value y[1] = y for a given integration t from t0
 *        and for given initial conditions provided in the structure pointer *params.
 **/
double cr3bp_event (double t, void *params)
{
    //------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------
    int i;
    struct ode_params *p = (struct ode_params *) params;
    double ystart[p->nvar];
    for(i=0; i<p->nvar; i++) ystart[i] = (p->y0)[i];
    double t0 = p->t0;

    //------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------
    //Starting in the right direction
    p->d->h = (t>t0) ? fabs(p->d->h) : -fabs(p->d->h);
    gsl_odeiv2_driver_apply (p->d, &t0, t, ystart);

    //------------------------------------------------------------------------------
    //Reset
    //------------------------------------------------------------------------------
    gsl_odeiv2_driver_reset(p->d);

    //------------------------------------------------------------------------------
    //Output
    //------------------------------------------------------------------------------
    return p->fvalue.value(t, ystart, p->fvalue.val_par);
}
