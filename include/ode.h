#ifndef CUSTOM_ODE_H_INCLUDED
#define CUSTOM_ODE_H_INCLUDED

//C++
#include <iostream>

//GSL
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

//Custom
#include "Config.h"
#include "env.h"

#define MAX_EVENTS 1000

using namespace std;

typedef struct OdeStruct OdeStruct;
struct OdeStruct
{
    //Stepper
    const gsl_odeiv2_step_type *T;
    gsl_odeiv2_step *s;
    //Control
    gsl_odeiv2_control *c;
    //Evolution
    gsl_odeiv2_evolve *e;
    //System
    gsl_odeiv2_system sys;
    //Driver
    gsl_odeiv2_driver * d;
    //Solver for root finding
    gsl_root_fsolver *s_root;

    //Precisions
    double eps_int_rel;  //for integration (relative precision)
    double eps_int_abs;  //for integration (absolute precision)
    double eps_root;     //for root finding
    double eps_diff;     //for differential correction

    //Initial step
    double h;

    //Dimension
    int dim;

    //Warnings
    int warnings;
};


void update_ode_structure(OdeStruct *ode_s,
                          const gsl_odeiv2_step_type *T,
                          gsl_odeiv2_step *s,
                          gsl_odeiv2_control *c,
                          gsl_odeiv2_evolve *e,
                          gsl_odeiv2_system sys,
                          gsl_odeiv2_driver * d,
                          gsl_root_fsolver *s_root,
                          double eps_root,
                          double eps_diff);

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
                        void *odeParams);

void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        size_t dim,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        void *odeParams);

void init_ode_structure(OdeStruct *ode_s,
                        const gsl_odeiv2_step_type *T,
                        const gsl_root_fsolver_type *T_root,
                        size_t dim,
                        int (* func) (double t, const double y[], double dydt[], void *params),
                        void *odeParams,
                        int warnings);

void reset_ode_structure(OdeStruct *ode_s);
void free_ode_structure(OdeStruct *ode_s);




//========================================================================================
// ODEZERO functions & structures
//========================================================================================
//Set the parameters to trigger an event
typedef struct value_params value_params;
struct value_params
{
    int max_events;
    int direction;
    int dim;
    double value;
    double *center;
    int type;
};

//Structure of the output of a value function for event procedure
typedef struct value_output value_output;
struct value_output
{
    double val;
    //int max_events;
    //int direction;
};

//Structure of a value function for event procedure
typedef struct value_function value_function;
struct value_function
{
    struct value_params *val_par;
    double (*value)(double, double [], void *);
};

//Structure dedicated to ode parameters
typedef struct ode_params ode_params;
struct ode_params
{
    int nvar;
    double t0;
    double* y0;
    gsl_odeiv2_driver *d;
    struct value_function fvalue;
};

int custom_odezero_2(double y[],
                     double** ye,
                     double *tcross,
                     double t0,
                     double t1,
                     OdeStruct *ode_s,
                     struct value_function fvalue);

double linear_intersection(double t, double yv[], void *params);
double angle_intersection(double t, double yv[], void *params);
double null_flight_path_angle(double t, double yv[], void *params);
double cr3bp_event (double t, void *params);

#endif // CUSTOM_ODE_H_INCLUDED
