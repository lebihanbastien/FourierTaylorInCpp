#include "single_orbit.h"
int COMPLETION = 0;

//========================================================================================
//
//          Integration - main routines
//
//========================================================================================
/**
 * \brief Integrates a generic vector field from any input type to any output type.
 **/
int ode78(double **yv, double *tv, OdeEvent *odeEvent,
          double t0NC, double tfNC, const double *y0NC,
          int nvar, int nGrid, int dcs,
          int inputType, int outputType)
{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    string fname = "ode78";

    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
    cerr << fname << ". The number of variables (nvar) must be of size 6 or 42." << endl;

    //Type of default coordinate system (dcs) for integration
    if(dcs > I_NJ2000) cerr << fname << ". Unknown dcs." << endl;

    //Type of inputs
    if(inputType > NJ2000) cerr << fname << ". Unknown inputType" << endl;

    //Type of outputs
    if(outputType > NJ2000)
        cerr << fname << ". Unknown outputType" << endl;

    //If nvar == 42, the dcs and the outputType must match,
    //otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs != default_coordinate_system(outputType))
        {
            cerr << fname << ". If the variational equations are desired, the "  << endl;
            cerr << "outputType must match the default coordinate system (dcs)." << endl;
        }
    }

    //If the JPL integration is required, the dcs, the outputType and the inputType
    //must match! This is to simplify  the change of coordinates afterwards
    if(dcs > I_ECISEM)
    {
        if(dcs != default_coordinate_system(outputType) ||
        dcs != default_coordinate_system(inputType))
        {
            cerr << fname << ". If a JPL integration is desired, the "  << endl;
            cerr << "outputType must match the default coordinate system (dcs)." << endl;
        }
    }


    //====================================================================================
    // 2. Define the framework from the default coordinate system
    //    Define also the default variable type that will be used throughout the computation
    //====================================================================================
    int varType = default_coordinate_type(dcs);

    //====================================================================================
    // 4. Selection of the vector field
    //====================================================================================
    vfptr vf = ftc_select_vf(dcs, nvar);

    //====================================================================================
    // 5. ODE system
    //====================================================================================
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs, odeEvent->detection, odeEvent->warnings);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, nvar, vf, &odeParams);


    //====================================================================================
    // 6. to NCFWRK coordinates.
    //====================================================================================
    double y0[nvar];
    double t0 = 0, tf = 0;

    //------------------------------------------------------------------------------------
    //From inputType to varType, if we are using a QBCP integration
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If we are dealing with a JPL integration, there is no need for a COC
        //--------------------------------------------------------------------------------
        case  I_ECLI:
        case  I_J2000:
        case  I_NJ2000:
        {
            for(int i =  0; i < nvar; i++) y0[i] = y0NC[i];
            break;
        }
        //--------------------------------------------------------------------------------
        // Else, we proceed
        //--------------------------------------------------------------------------------
        default:
        qbcp_coc(t0NC, y0NC, y0, inputType, varType);
    }


    //------------------------------------------------------------------------------------
    //For the initial and final time, a switch is necessary
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If I_PSEM/I_VSEM, etc... the system is focused on the SEM system via SEML
        //--------------------------------------------------------------------------------
    case I_NCSEM:
    case I_VNCSEM:
    case I_PSEM:
    case I_VSEM:
    case I_INSEM:
    case I_ECISEM:
    {
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
        case VEM:
        case INEM:
            //Time is now in SEM units
            t0 = t0NC*SEML.us_em.ns;
            tf = tfNC*SEML.us_em.ns;
            break;
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // If I_PEM/I_VEM, etc, the system is focused on the EM system via SEML
        //--------------------------------------------------------------------------------
    case I_NCEM:
    case I_VNCEM:
    case I_PEM:
    case I_VEM:
    case I_INEM:
    {
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
        case INSEM:
        case VSEM:
        case ECISEM:
            //Time is now in EM units
            t0 = t0NC/SEML.us_em.ns;
            tf = tfNC/SEML.us_em.ns;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
        case VEM:
        case INEM:
            //Time is already in EM units
            t0 = t0NC;
            tf = tfNC;
            break;
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // The default case will apply on the JPL integration schemes, for which the time
        // is already in the right units (seconds, or normalized time)
        //--------------------------------------------------------------------------------
    default:
    {
        t0 = t0NC;
        tf = tfNC;
    }
    }

    //------------------------------------------------------------------------------------
    // Add the identity matrix if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //====================================================================================
    // 7. Integration, for 2 outputs
    //====================================================================================
    double **yvIN, *tvIN;
    yvIN  = dmatrix(0, nvar-1, 0, nGrid);
    tvIN  = dvector(0, nGrid);

    //Integration in yvIN/tvIN
    int status = ode78_grid(&odestruct, t0, tf, y0, yvIN, tvIN, nGrid);

    //====================================================================================
    // 8. To the right outputs: varType to outputType, if necessary
    //====================================================================================
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If we are dealing with a JPL integration, there is no need for a COC
        //--------------------------------------------------------------------------------
    case  I_ECLI:
    case  I_J2000:
    case  I_NJ2000:
    {
        for(int p = 0; p <= nGrid; p++)
        {
            for(int i = 0; i < nvar; i++) yv[i][p] = yvIN[i][p];
            tv[p] = tvIN[p];
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // Else, we proceed to the COC
        //--------------------------------------------------------------------------------
    default:
    {
        qbcp_coc_vec(yvIN, tvIN, yv, tv, nGrid, varType, outputType);
        // Add the STM if necessary
        if(nvar == 42)
        {
            for(int p = 0; p <= nGrid; p++)
                for(int i = 6; i < 42; i++) yv[i][p] = yvIN[i][p];
        }
    }
    }


    //====================================================================================
    // Update the collisionner
    //====================================================================================
    odeEvent->coll           = odeParams.event.coll;
    odeEvent->crossings      = odeParams.event.crossings;
    odeEvent->crossings_moon = odeParams.event.crossings_moon;

    //====================================================================================
    // 10. Free memory
    //====================================================================================
    free_dmatrix(yvIN, 0, nvar-1, 0, nGrid);
    free_dvector(tvIN, 0, nGrid);

    //====================================================================================
    // 11. Something went wrong inside the stepper?
    //====================================================================================
    return status;
}


/**
 * \brief Integrates a generic vector field from any input type to any output type.
 **/
int ode78(double **yv, double *tv, int *ode78coll,
          double t0NC, double tfNC, const double *y0NC,
          int nvar, int nGrid, int dcs,
          int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // Use of the other ode78 with OdeEvent.
    //------------------------------------------------------------------------------------
    OdeEvent odeEvent;
    int status = ode78(yv, tv, &odeEvent, t0NC, tfNC, y0NC,
                       nvar, nGrid, dcs, inputType, outputType);

    //------------------------------------------------------------------------------------
    // Update of ode78coll from the information in the OdeEvent.
    //------------------------------------------------------------------------------------
    *ode78coll = odeEvent.coll;

    return status;
}


//========================================================================================
//
//          Integration - with event
//
//========================================================================================
/**
 * \brief Integrates the QBCP vector field from any input type to any output type.
 *        The output are states and time corresponding to a particular event stored
 *        in the structure val_par.
 **/
int ode78_qbcp_event(double **ye, double *te, int *ode78coll,
                     double t0NC, double tfNC, const double *y0NC,
                     int nvar, int dcs,
                     int inputType, int outputType,
                     value_params* val_par)
{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
        perror("ode78. The number of variables (nvar) must be of size 6 or 42.");

    //Type of default coordinate system (dcs) for integration
    if(dcs > 5)
        perror("ode78. Wrong dcs.");

    //Type of inputs
    if(inputType > 7)
        perror("ode78. Unknown inputType");

    //Type of outputs
    if(outputType > 8)
        perror("ode78. Unknown outputType");

    //If nvar == 42, the dcs and the outputType must match, otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs == I_PEM && outputType != PEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_NCEM && outputType != NCEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VNCEM && outputType != VNCEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");

        if(dcs == I_PSEM && outputType != PSEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_NCSEM && outputType != NCSEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VNCSEM && outputType != VNCSEM)
            perror("ode78. If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
    }


    //====================================================================================
    // 2. Define the framework from the default coordinate system
    //    Define also the default variable type that will be used throughout the computation
    //====================================================================================
    int varType = 0;
    switch(dcs)
    {
    case I_PSEM:
        varType = PSEM;
        break;
    case I_NCSEM:
        varType = NCSEM;
        break;
    case I_VNCSEM:
        varType = VNCSEM;
        break;
    case I_PEM:
        varType = PEM;
        break;
    case I_NCEM:
        varType = NCEM;
        break;
    case I_VNCEM:
        varType = VNCEM;
        break;
    }


    //====================================================================================
    // 4. Selection of the vector field
    //====================================================================================
    vfptr vf = ftc_select_vf(dcs, nvar);


    //====================================================================================
    // 5. ODE system
    //====================================================================================
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, nvar, vf, &odeParams);


    //====================================================================================
    // 6. to NCFWRK coordinates.
    //====================================================================================
    double y0[nvar];
    double t0 = 0, tf = 0;

    //From inputType to varType
    qbcp_coc(t0NC, y0NC, y0, inputType, varType);

    //For the initial and final time, a switch is necessary
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If I_PSEM/I_VSEM, the system is focused on the SEM system via
        // SEML
        //--------------------------------------------------------------------------------
    case I_PSEM:
    case I_NCSEM:
    case I_VNCSEM:
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
            //Time is now in SEM units
            t0 = t0NC*SEML.us_em.ns;
            tf = tfNC*SEML.us_em.ns;
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // If I_PEM/I_VEM, the system is focused on the EM system via
        // SEML
        //--------------------------------------------------------------------------------
    case I_PEM:
    case I_NCEM:
    case I_VNCEM:
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            //Time is now in EM units
            t0 = t0NC/SEML.us_em.ns;
            tf = tfNC/SEML.us_em.ns;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
            //Time is already in EM units
            t0 = t0NC;
            tf = tfNC;
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Add the identity matrix if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //====================================================================================
    // 7. Integration, for 2 outputs
    //====================================================================================
    double **yvIN, *tvIN;
    int nGrid = val_par->max_events;
    yvIN  = dmatrix(0, nvar-1, 0, nGrid);
    tvIN  = dvector(0, nGrid);

    //------------------------------------------------------------------------------------
    // Integration until t = t1, and saving the events
    //------------------------------------------------------------------------------------
    struct value_function fvalue;
    fvalue.val_par = val_par;

    switch(val_par->type)
    {
        case 'X':
        case 'Y':
        case 'Z':
        {
            fvalue.value   = &linear_intersection;
            break;
        }

        case 'R':
        {
            fvalue.value   = &radius_intersection;
            break;
        }


        case 'A':
        {
            fvalue.value   = &angle_intersection;
            break;
        }

        case 'F':
        {
            fvalue.value   = &null_flight_path_angle;
            break;
        }

        default:
        {
            cerr << "Bad event type." << endl;
            return FTC_EDOM;
        }
    }

    //Integration in yvIN/tvIN, with event detection
    int nEvents = custom_odezero_2(y0, yvIN, tvIN, t0, tf, &odestruct, fvalue);
    //ode78_grid(&odestruct, t0, tf, y0, yvIN, tvIN, nGrid);

    //====================================================================================
    // 8. To the right outputs: varType to outputType
    //====================================================================================
    qbcp_coc_vec(yvIN, tvIN, ye, te, nGrid, varType, outputType);


    //------------------------------------------------------------------------------------
    // Add the STM if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        for(int p = 0; p <= nGrid; p++) for(int i = 6; i < 42; i++) ye[i][p] = yvIN[i][p];
    }

    //====================================================================================
    // Update the collisionner
    //====================================================================================
    *ode78coll = odeParams.event.coll;

    //====================================================================================
    // 10. Free memory
    //====================================================================================
    free_dmatrix(yvIN, 0, nvar-1, 0, nGrid);
    free_dvector(tvIN, 0, nGrid);

    return nEvents+1;
}


//========================================================================================
//
//          Integration - variable and given time grids (not really used)
//
//========================================================================================
/**
 * \brief Integrates the QBCP vector field from any input type to any output type. Return the last position that is filled on the grid.
 **/
int ode78_qbcp_vg(double **yv, double *tv, int *ode78coll,
                  double t0NC, double tfNC, const double *y0NC,
                  int nvar, int nGridmax, int dcs,
                  int inputType, int outputType, int dli)
{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
        perror("The number of variables (nvar) must be of size 6 or 42.");

    //Type of default coordinate system (dcs) for integration
    if(dcs > 5)
        perror("Wrong dcs.");

    //Type of inputs
    if(inputType > 7)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > 8)
        perror("Unknown outputType");

    //If nvar == 42, the dcs and the outputType must match, otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs == I_PEM && outputType != PEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_NCEM && outputType != NCEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VNCEM && outputType != VNCEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");

        if(dcs == I_PSEM && outputType != PSEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_NCSEM && outputType != NCSEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VNCSEM && outputType != VNCSEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
    }


    //====================================================================================
    // 2. Define the framework from the default coordinate system
    //    Define also the default variable type that will be used throughout the computation
    //====================================================================================
    int varType = default_coordinate_type(dcs);

    //====================================================================================
    // 4. Selection of the vector field
    //====================================================================================
    vfptr vf = ftc_select_vf(dcs, nvar);

    //====================================================================================
    // 5. ODE system
    //====================================================================================
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, nvar, vf, &odeParams);


    //====================================================================================
    // 6. to NCFWRK coordinates.
    //====================================================================================
    double y0[nvar];
    double t0 = 0, tf = 0;

    //From inputType to varType
    qbcp_coc(t0NC, y0NC, y0, inputType, varType);

    //For the initial and final time, a switch is necessary
    switch(dcs)
    {
        //-------------------------------------------------------------------------------------
        // If I_PSEM/I_VSEM, the system is focused on the SEM system via
        // SEML
        //-------------------------------------------------------------------------------------
    case I_PSEM:
    case I_NCSEM:
    case I_VNCSEM:
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            //Time is already in SEM units
            t0 = t0NC;
            tf = tfNC;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
            //Time is now in SEM units
            t0 = t0NC*SEML.us_em.ns;
            tf = tfNC*SEML.us_em.ns;
            break;
        }
        break;

        //-------------------------------------------------------------
        // If I_PEM/I_VEM, the system is focused on the EM system via
        // SEML
        //-------------------------------------------------------------
    case I_PEM:
    case I_NCEM:
    case I_VNCEM:
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
            //Time is now in EM units
            t0 = t0NC/SEML.us_em.ns;
            tf = tfNC/SEML.us_em.ns;
            break;
        case PEM:
        case NCEM:
        case VNCEM:
            //Time is already in EM units
            t0 = t0NC;
            tf = tfNC;
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Add the identity matrix if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //====================================================================================
    // 7. Integration, for 2 outputs
    //====================================================================================
    double **yvIN, *tvIN, **yvINc, *tvINc;
    yvIN   = dmatrix(0, nvar-1, 0, nGridmax);
    tvIN   = dvector(0, nGridmax);
    yvINc  = dmatrix(0, nvar-1, 0, nGridmax);
    tvINc  = dvector(0, nGridmax);

    //Integration in yvIN/tvIN
    int lastind = ode78_variable_grid(&odestruct, t0, tf, y0, yvIN, tvIN, nGridmax);

    //====================================================================================
    // 8. Get only a certain number of points of points if necessary
    //====================================================================================
    if(dli > 0)
    {
        int freq = floor(lastind/dli);
        // Select a subset
        for(int kman = 0; kman <= dli; kman++)
        {
            for(int i = 0; i < nvar; i++) yvINc[i][kman] = yvIN[i][freq*kman];
            tvINc[kman] = tvIN[freq*kman];
        }

        // Put the last point if it was not selected in the previous subset
        if(freq*dli < lastind)
        {
            for(int i = 0; i < nvar; i++) yvINc[i][dli] = yvIN[i][lastind];
            tvINc[dli] = tvIN[lastind];
        }

        //Put the subset back in yIN
        for(int kman = 0; kman <= dli; kman++)
        {
            for(int i = 0; i < nvar; i++) yvIN[i][kman] = yvINc[i][kman];
            tvIN[kman] = tvINc[kman];
        }

        //New lastindix is dli
        lastind = dli;
    }


    //====================================================================================
    // 8. To the right outputs: varType to outputType
    //====================================================================================
    qbcp_coc_vec(yvIN, tvIN, yv, tv, lastind, varType, outputType);

    //------------------------------------------------------------------------------------
    // Add the STM if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        for(int p = 0; p <= lastind; p++) for(int i = 6; i < 42; i++) yv[i][p] = yvIN[i][p];
    }


    //====================================================================================
    // Update the collisionner
    //====================================================================================
    *ode78coll = odeParams.event.coll;

    //====================================================================================
    // 10. Free memory
    //====================================================================================
    free_dmatrix(yvIN,  0, nvar-1, 0, nGridmax);
    free_dvector(tvIN,  0, nGridmax);
    free_dmatrix(yvINc, 0, nvar-1, 0, nGridmax);
    free_dvector(tvINc, 0, nGridmax);

    return lastind;
}

/**
 * \brief Integrates the QBCP vector field from any input type to any output type.
 *        The time grid is given in inputs (tvi).
 **/
int ode78_qbcp_gg(double **yv, double *tvf, int *ode78coll,
                  const double *tvi,
                  const double *y0NC,
                  int nvar, int nGrid, int dcs,
                  int inputType, int outputType)
{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    string fname = "ode78_qbcp_gg";

    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
    cerr << fname << ". The number of variables (nvar) must be of size 6 or 42." << endl;

    //Type of default coordinate system (dcs) for integration
    if(dcs > I_NJ2000) cerr << fname << ". Unknown dcs." << endl;

    //Type of inputs
    if(inputType > NJ2000) cerr << fname << ". Unknown inputType" << endl;

    //Type of outputs
    if(outputType > NJ2000)
        cerr << fname << ". Unknown outputType" << endl;

    //If nvar == 42, the dcs and the outputType must match,
    //otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs != default_coordinate_system(outputType))
        {
            cerr << fname << ". If the variational equations are desired, the "  << endl;
            cerr << "outputType must match the default coordinate system (dcs)." << endl;
        }
    }

    //If the JPL integration is required, the dcs, the outputType and the inputType
    //must match! This is to simplify  the change of coordinates afterwards
    if(dcs > I_ECISEM)
    {
        if(dcs != default_coordinate_system(outputType) ||
        dcs != default_coordinate_system(inputType))
        {
            cerr << fname << ". If a JPL integration is desired, the "  << endl;
            cerr << "outputType must match the default coordinate system (dcs)." << endl;
        }
    }


    //====================================================================================
    // 2. Define the framework from the default coordinate system
    //    Define also the default variable type that will be used throughout the computation
    //====================================================================================
    int varType = default_coordinate_type(dcs);


    //====================================================================================
    // 4. Selection of the vector field
    //====================================================================================
    vfptr vf = ftc_select_vf(dcs, nvar);

    //====================================================================================
    // 5. ODE system
    //====================================================================================
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, nvar, vf, &odeParams);


    //====================================================================================
    // 6. to NCFWRK coordinates.
    //====================================================================================
    double y0[nvar];
    double *tvIN = dvector(0, nGrid);

    //------------------------------------------------------------------------------------
    //From inputType to varType, if we are using a QBCP integration
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If we are dealing with a JPL integration, there is no need for a COC
        //--------------------------------------------------------------------------------
        case  I_ECLI:
        case  I_J2000:
        case  I_NJ2000:
        {
            for(int i =  0; i < nvar; i++) y0[i] = y0NC[i];
            break;
        }
        //--------------------------------------------------------------------------------
        // Else, we proceed
        //--------------------------------------------------------------------------------
        default:
        qbcp_coc(tvIN[0], y0NC, y0, inputType, varType);
    }

    //------------------------------------------------------------------------------------
    //For the time vector a switch is necessary
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If I_PSEM/I_VSEM, etc... the system is focused on the SEM system via SEML
        //--------------------------------------------------------------------------------
    case I_NCSEM:
    case I_VNCSEM:
    case I_PSEM:
    case I_VSEM:
    case I_INSEM:
    case I_ECISEM:
    {
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //Time is already in SEM units
            for(int i = 0; i <= nGrid; i++) tvIN[i] = tvi[i];
            break;
        case PEM:
        case NCEM:
        case VNCEM:
        case VEM:
        case INEM:
            //Time is now in SEM units
            for(int i = 0; i <= nGrid; i++) tvIN[i] = tvi[i]*SEML.us_em.ns;
            break;
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // If I_PEM/I_VEM, etc, the system is focused on the EM system via SEML
        //--------------------------------------------------------------------------------
    case I_NCEM:
    case I_VNCEM:
    case I_PEM:
    case I_VEM:
    case I_INEM:
    {
        //Time
        switch(inputType)
        {
        case PSEM:
        case NCSEM:
        case VNCSEM:
        case INSEM:
        case VSEM:
        case ECISEM:
            //Time is already in EM units
            for(int i = 0; i <= nGrid; i++) tvIN[i] = tvi[i];
            break;
        case PEM:
        case NCEM:
        case VNCEM:
        case VEM:
        case INEM:
            //Time is already in EM units
            for(int i = 0; i <= nGrid; i++) tvIN[i] = tvi[i];
            break;
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // The default case will apply on the JPL integration schemes, for which the time
        // is already in the right units (seconds, or normalized time)
        //--------------------------------------------------------------------------------
    default:
    {
        for(int i = 0; i <= nGrid; i++) tvIN[i] = tvi[i];
    }
    }

    //------------------------------------------------------------------------------------
    // Add the identity matrix if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix *Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //====================================================================================
    // 7. Integration, for 2 outputs
    //====================================================================================
    double **yvIN;
    yvIN  = dmatrix(0, nvar-1, 0, nGrid);

    //Integration in yvIN/tvIN

    vector_printf(tvIN, nGrid);
    ode78_grid_gg(&odestruct, y0, yvIN, tvIN, nGrid);

    //====================================================================================
    // 8. To the right outputs: varType to outputType, if necessary
    //====================================================================================
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If we are dealing with a JPL integration, there is no need for a COC
        //--------------------------------------------------------------------------------
    case  I_ECLI:
    case  I_J2000:
    case  I_NJ2000:
    {
        for(int p = 0; p <= nGrid; p++)
        {
            for(int i = 0; i < nvar; i++) yv[i][p] = yvIN[i][p];
            tvf[p] = tvIN[p];
        }
        break;
    }
        //--------------------------------------------------------------------------------
        // Else, we proceed to the COC
        //--------------------------------------------------------------------------------
    default:
    {
        qbcp_coc_vec(yvIN, tvIN, yv, tvf, nGrid, varType, outputType);
        // Add the STM if necessary
        if(nvar == 42)
        {
            for(int p = 0; p <= nGrid; p++)
                for(int i = 6; i < 42; i++) yv[i][p] = yvIN[i][p];
        }
    }
    }

    //====================================================================================
    // Update the collisionner
    //====================================================================================
    *ode78coll = odeParams.event.coll;


    //====================================================================================
    // 10. Free memory
    //====================================================================================
    free_dmatrix(yvIN, 0, nvar-1, 0, nGrid);
    free_dvector(tvIN, 0, nGrid);

    return GSL_SUCCESS;
}


//========================================================================================
//
//          Integration - JPL (deprecated)
//
//========================================================================================
/**
 * \brief Integrates the Ephemerides vector field. Only VECLI, VSEM and VEM coordinates are used for now.
 **/
int ode78_jpl(double **yv, double *tv, int *ode78coll,
              double t0NC, double tfNC, const double *y0NC,
              int nvar, int nGrid, int dcs,
              int inputType, int outputType)

{
    //====================================================================================
    // 1. Do some checks on the inputs
    //====================================================================================
    //Size of the variable vector
    if(nvar!=6 && nvar!=42)
        perror("ode78_jpl. The number of variables (nvar) must be of size 6 or 42.");

    //Type of default coordinate system (dcs) for integration
    if(dcs != I_ECLI && dcs != I_VSEM && dcs != I_VEM && dcs != I_J2000 && dcs != I_NJ2000)
        perror("ode78_jpl. Wrong dcs (only I_ECLI, I_VSEM, I_J2000 & I_NJ2000 are allowed).");

    //Type of inputs
    if(inputType != VECLI && inputType != VSEM && inputType != VEM && inputType != J2000 && inputType != NJ2000)
        perror("ode78_jpl. Unknown inputType");

    //Type of outputs
    if(outputType != VECLI && outputType != VSEM && outputType != VEM && outputType != J2000 && inputType != NJ2000)
        perror("ode78_jpl. Unknown outputType");

    //If nvar == 42, the dcs and the outputType must match, otherwise the variational equations make no sense with the outputs
    if(nvar == 42)
    {
        if(dcs == I_ECLI && outputType != VECLI)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VSEM && outputType != VSEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_VEM && outputType != VEM)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_J2000 && outputType != J2000)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
        if(dcs == I_NJ2000 && outputType != NJ2000)
            perror("If the variational equations are desired, the outputType must match the default coordinate system (dcs).");
    }


    //====================================================================================
    // 4. Selection of the vector field
    //====================================================================================
    int (*vf)(double, const double*, double*, void*) = jpl_vf; //by default, to avoid warning from gcc compiler
    switch(nvar)
    {
        //--------------------------------------------------------------------------------
        // 6 variables: only the state is integrated
        //--------------------------------------------------------------------------------
    case 6:
        switch(dcs)
        {
        case I_ECLI:
            vf = jpl_vf;       //vector field with a state (X, VX)
            break;
        case I_J2000:
            vf = jpl_vf_eci;      //vector field with a state (X, VX)
            break;
        case I_NJ2000:
            vf = jpl_vf_neci;
            break;
        case I_VSEM:
        case I_VEM:
            vf = jpl_vf_syn; //vector field with a state (x, vx)
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // 42 variables: the state + var. eq. are integrated
        //--------------------------------------------------------------------------------
    case 42:
        switch(dcs)
        {
        case I_ECLI:
            vf = jpl_vf_var; //vector field with a state (X, VX)
            break;
        case I_J2000:
            vf = jpl_vf_eci_var; //vector field with a state (X, VX)
            break;
        case I_NJ2000:
            vf = jpl_vf_neci_var; //vector field with a state (X, VX)
            break;
        case I_VSEM:
        case I_VEM:
            vf = jpl_vf_syn_var; //vector field with a state (x, vx)
            break;
        default:
            perror("ode78_jpl. Unknown dcs with 42 states (no corresponding vf).");
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // Default case: should NOT be reached
        //--------------------------------------------------------------------------------
    default:
        vf = jpl_vf;
        break;
    }

    //====================================================================================
    // 5. ODE system
    //====================================================================================
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, nvar, vf, &odeParams);

    //====================================================================================
    // 6. to FWRK coordinates.
    //====================================================================================
    double y0[nvar]; //y0_VSEM[nvar],
    double et0 = 0, etf = 0;

    //------------------------------------------------------------------------------------
    //For the initial and final time, a switch is necessary
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //--------------------------------------------------------------------------------
        // If I_ECLI, the time must be set in seconds
        //--------------------------------------------------------------------------------
    case I_ECLI:
        //Time
        switch(inputType)
        {
        case VECLI:
            //Time is already in seconds
            et0 = t0NC;
            etf = tfNC;
            break;
        case VSEM:
        case VEM:
        case J2000:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //Time is already in SEM units, so we only need to set it in seconds
            //            et0 = t0NC/SEML.n_sem;
            //            etf = tfNC/SEML.n_sem;
            //////////////////////////////////////////////////
            break;
        }
        break;

        //-------------------------------------------------------------------------------------
        // If I_J2000, the time must NOT be normalized
        //-------------------------------------------------------------------------------------
    case I_J2000:
        //Time
        switch(inputType)
        {
        case J2000:
            //Time is already in seconds
            et0 = t0NC;
            etf = tfNC;
            break;
        case VSEM:
        case VEM:
        case VECLI:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //Time is already in SEM units, so we only need to set it in seconds
            //            et0 = t0NC/SEML.n_sem;
            //            etf = tfNC/SEML.n_sem;
            //////////////////////////////////////////////////
            break;
        }
        break;

        //-------------------------------------------------------------------------------------
        // If I_NJ2000, the time must be normalized
        //-------------------------------------------------------------------------------------
    case I_NJ2000:
        //Time
        switch(inputType)
        {
        case NJ2000:
            //Time is already normalized
            et0 = t0NC;
            etf = tfNC;
            break;
        case VSEM:
        case VEM:
        case VECLI:
        case J2000:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //Time is already in SEM units, so we only need to set it in seconds
            //            et0 = t0NC/SEML.n_sem;
            //            etf = tfNC/SEML.n_sem;
            //////////////////////////////////////////////////
            break;
        }
        break;
        //-------------------------------------------------------------------------------------
        // If I_VSEM, the time must be normalized
        //-------------------------------------------------------------------------------------
    case I_VSEM:
        //Time
        switch(inputType)
        {
        case VECLI:
        case VEM:
        case J2000:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //Time is in seconds, we need it normalized
            //            et0 = t0NC*SEML.n_sem;
            //            etf = tfNC*SEML.n_sem;
            //////////////////////////////////////////////////
            break;
        case VSEM:
            //Time is already in SEM units
            et0 = t0NC;
            etf = tfNC;
            break;
        }
        break;
        //-------------------------------------------------------------------------------------
        // If I_PSEM, the time must be normalized
        //-------------------------------------------------------------------------------------
    case I_VEM:
        //Time
        switch(inputType)
        {
        case VECLI:
        case VSEM:
        case J2000:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //Time is in seconds, we need it normalized
            //            et0 = t0NC*SEML.n_sem;
            //            etf = tfNC*SEML.n_sem;
            //////////////////////////////////////////////////
            break;
        case VEM:
            //Time is already in EM units
            et0 = t0NC;
            etf = tfNC;
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    //From inputType to dcs coordinates
    //------------------------------------------------------------------------------------
    switch(dcs)
    {
        //-------------------------------------------------------------------------------------
        // If I_ECLI
        //-------------------------------------------------------------------------------------
    case I_ECLI:
        //Time
        switch(inputType)
        {
        case VECLI:
            //Simple copy
            for(int i = 0; i < 6; i++) y0[i] = y0NC[i];
            break;
        case VSEM:
        case VEM:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //If inputType is not VECLI, then the initial state is -> eph_coord(inputType)
            //            qbcp_coc(t0NC, y0NC, y0_VSEM, inputType, eph_coord(inputType));
            //            //Then, to inputType
            //            syn2eclstate(y0_VSEM, y0, et0, eph_coord(inputType));
            //////////////////////////////////////////////////
            break;
        }
        break;

        //-------------------------------------------------------------------------------------
        // If I_J2000
        //-------------------------------------------------------------------------------------
    case I_J2000:
        //Time
        switch(inputType)
        {
        case J2000:
            //Simple copy
            for(int i = 0; i < 6; i++) y0[i] = y0NC[i];
            break;
        case VSEM:
        case VEM:
        case VECLI:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //If inputType is not VECLI, then the initial state is -> eph_coord(inputType)
            //            qbcp_coc(t0NC, y0NC, y0_VSEM, inputType, eph_coord(inputType));
            //            //Then, to inputType
            //            syn2eclstate(y0_VSEM, y0, et0, eph_coord(inputType));
            //////////////////////////////////////////////////
            break;
        }
        break;

        //-------------------------------------------------------------------------------------
        // If I_NJ2000
        //-------------------------------------------------------------------------------------
    case I_NJ2000:
        //Time
        switch(inputType)
        {
        case NJ2000:
            //Simple copy
            for(int i = 0; i < 6; i++) y0[i] = y0NC[i];
            break;
        case VSEM:
        case VEM:
        case VECLI:
        case J2000:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            //If inputType is not VECLI, then the initial state is -> eph_coord(inputType)
            //            qbcp_coc(t0NC, y0NC, y0_VSEM, inputType, eph_coord(inputType));
            //            //Then, to inputType
            //            syn2eclstate(y0_VSEM, y0, et0, eph_coord(inputType));
            //////////////////////////////////////////////////
            break;
        }
        break;


        //-------------------------------------------------------------------------------------
        // If I_VSEM
        //-------------------------------------------------------------------------------------
    case I_VSEM:
        //Time
        switch(inputType)
        {
        case VECLI:
        case VEM:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            ecl2synstate(y0NC, t0NC y0, eph_coord(inputType));
            //            qbcp_coc(et0, y0NC, y0_VSEM, inputType, VSEM);
            //////////////////////////////////////////////////
            break;
        case VSEM:
            //Simple copy
            for(int i = 0; i < 6; i++) y0[i] = y0NC[i];
            break;
        }
        break;
        //-------------------------------------------------------------------------------------
        // If I_VEM
        //-------------------------------------------------------------------------------------
    case I_VEM:
        //Time
        switch(inputType)
        {
        case VECLI:
        case VSEM:
            //////////////////////////////////////////////////
            //            //€€TODO: probably BAD
            //            ecl2synstate(y0NC, t0NC y0, eph_coord(inputType));
            //            qbcp_coc(et0, y0NC, y0_VSEM, inputType, VSEM);
            //////////////////////////////////////////////////
            break;
        case VEM:
            //Simple copy
            for(int i = 0; i < 6; i++) y0[i] = y0NC[i];
            break;
        }
        break;
    }

    //------------------------------------------------------------------------------------
    // Add the identity matrix if necessary
    //------------------------------------------------------------------------------------
    if(nvar == 42)
    {
        //Identity matrix eye(6)
        gsl_matrix* Id = gsl_matrix_alloc(6,6);
        gsl_matrix_set_identity (Id);
        //Storing eye(6) into the initial vector
        gslc_matrixToVector(y0, Id, 6, 6, 6);
        gsl_matrix_free(Id);
    }

    //====================================================================================
    // 7. Integration, for 2 outputs
    //====================================================================================
    double** yvIN, *tvIN;
    yvIN  = dmatrix(0, nvar-1, 0, nGrid);
    tvIN  = dvector(0, nGrid);
    //Integration in yvIN/tvIN
    ode78_grid(&odestruct, et0, etf, y0, yvIN, tvIN, nGrid);

    //====================================================================================
    // 8. To the right outputs: varType to outputType. Fo now, simple copy
    //====================================================================================
    //qbcp_coc_vec(yvIN, tvIN, yv, tv, nGrid, varType, outputType);
    for(int p = 0; p <= nGrid; p++)
    {
        for(int i = 0; i < nvar; i++) yv[i][p] = yvIN[i][p];
        tv[p] = tvIN[p];
    }


    //====================================================================================
    // Update the collisionner
    //====================================================================================
    *ode78coll = odeParams.event.coll;


    //====================================================================================
    // 10. Free memory
    //====================================================================================
    free_dmatrix(yvIN, 0, nvar-1, 0, nGrid);
    free_dvector(tvIN, 0, nGrid);

    return GSL_SUCCESS;
}


//========================================================================================
//
//          Integration - subroutines
//
//========================================================================================
/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int ode78_grid(OdeStruct *odestruct, double t0, double tf, double *y0, double **yNCE, double *tNCE, int N)
{
    //====================================================================================
    //Initialization
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Parameters
    //------------------------------------------------------------------------------------
    OdeParams* odP = (OdeParams*) odestruct->sys.params;

    //------------------------------------------------------------------------------------
    //Initialization of the parameters
    //------------------------------------------------------------------------------------
    int nvar = odestruct->dim;
    int status;            //current status
    double yv[nvar], t;       //current state and time
    int nt;
    double ti;

    //------------------------------------------------------------------------------------
    //Initialization of the odestruct
    //------------------------------------------------------------------------------------
    //Change sign of step if necessary
    if((tf < t0 && odestruct->d->h>0) || (tf > t0 && odestruct->d->h<0)) odestruct->d->h *= -1;
    if((tf < t0 && odestruct->h>0)    || (tf > t0 && odestruct->h<0))    odestruct->h    *= -1;

    //Reset ode structure.
    reset_ode_structure(odestruct);

    //------------------------------------------------------------------------------------
    //Initialization of the IC
    //------------------------------------------------------------------------------------
    //Init the state & time
    for(int k = 0; k < nvar; k++) yv[k] = y0[k];
    t = t0;
    nt = 1;

    //====================================================================================
    //Creation of the time grid
    //====================================================================================
    for(int k = 0; k <= N; k++)
    {
        tNCE[k] = t0 + (double) k *(tf-t0)/N;
    }

    //====================================================================================
    //Loop
    //====================================================================================
    //First position
    for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];

    //Limit the max step
    //gsl_odeiv2_driver_set_hmax (odestruct->d, 1.0);

    //Other positions
    do
    {
        //Evolve up to ti
        ti = tNCE[nt];

        if(odP->event.detection)
        {
            // New version: gslc_event_evolve, including event detection via odP
            status = gslc_event_evolve(odestruct, odP, yv, &t, ti);
        }
        else
        {
            // Old version: gsl_odeiv2_driver_apply from GSL
            status = gsl_odeiv2_driver_apply(odestruct->d, &t , ti, yv);
        }

        //Update the storage
        for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];

        //Advance one step
        nt++;
    }
    while((nt<=N) && (status == GSL_SUCCESS));

    //------------------------------------------------------------------------------------
    //Something went wrong inside the stepper?
    //------------------------------------------------------------------------------------
    if(status)
    {
        if(odP->event.warnings) cout << "Warning in ode78_grid: the stepper went wrong. No output is produced.";
        return FTC_FAILURE;
    }

    return FTC_SUCCESS;
}

/**
 *  \brief Stepper that includes detections of event as described in OdeParams.
 **/
int gslc_event_step(OdeStruct* odestruct, OdeParams* odP, double yv[], double* t, double t1)
{
    //------------------------------------------------------------------------------------
    //Initialize
    //------------------------------------------------------------------------------------
    double x1, x2;
    double dist_to_Bem;

    //------------------------------------------------------------------------------------
    //Evolve one step
    //------------------------------------------------------------------------------------
    int status = gsl_odeiv2_evolve_apply (odestruct->e, odestruct->c, odestruct->s, &odestruct->sys, t, t1, &odestruct->h, yv);

    //------------------------------------------------------------------------------------
    //Check crossings of x = -1, in EML2 NCSEM coordinates
    //------------------------------------------------------------------------------------
    x1 = odP->event.x1;
    if(x1 == -1.0) odP->event.x1 = yv[0] + 1.0; //first use
    else
    {
        x2 = yv[0] + 1.0;

        if(x1 > 0 && x2 < 0)
        {
            if(yv[1] > 0) odP->event.crossings += 1.0;    //clockwise
            else          odP->event.crossings += 0.1;    //counterclockwise
        }
        else if(x1 < 0 && x2 > 0)
        {

            if(yv[1] < 0) odP->event.crossings += 1.0;    //clockwise
            else          odP->event.crossings += 0.1;    //counterclockwise
        }
        odP->event.x1 = x2;
    }

    //------------------------------------------------------------------------------------
    // Check crossings with the radius of the Moon, in EML2 NCSEM coordinates
    // Mean radius of the lunar orbit:
    //          odP->event.mr_moon = 2.518747349676265e-01 (see env.h)
    // being centered on (-1.0, 0, 0), the position of Bem, the Earth-Moon center of mass.
    //------------------------------------------------------------------------------------
    x1 = odP->event.x1_moon;
    if(x1 == -1.0)
    {
        //We compare the distance to (-1.0, 0) wrt to the radius of the lunar orbit
        dist_to_Bem = sqrt((-1.0 - yv[0])*(-1.0 - yv[0]) + yv[1]*yv[1] + yv[2]*yv[2]);
        odP->event.x1_moon = dist_to_Bem - odP->event.mr_moon; //first use
    }
    else
    {
        //We compare the distance to (-1.0, 0) wrt to the radius of the lunar orbit
        dist_to_Bem = sqrt((-1.0 - yv[0])*(-1.0 - yv[0]) + yv[1]*yv[1] + yv[2]*yv[2]);
        x2 = dist_to_Bem - odP->event.mr_moon; //first use

        //Detecting a cross of the lunar orbit
        if(x1*x2 < 0)
        {
            odP->event.crossings_moon += 1.0;
        }

        odP->event.x1_moon = x2;
    }

    return status;
}

/**
 *   \brief Integrates the current state yv up to t = t1, includes detections of event as described in OdeParams.
 **/
int gslc_event_evolve(OdeStruct* odestruct, OdeParams* odP, double yv[], double* t, double t1)
{
    int status;
    do
    {
            //Evolve one step
            status = gslc_event_step(odestruct, odP, yv, t, t1);

    }while((status == GSL_SUCCESS) && fabs(*t - t1) != 0);

    return status;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid. Return the last position that is filled on the grid.
 **/
int ode78_variable_grid(OdeStruct *odestruct, double t0, double tf, double *y0, double **yNCE, double *tNCE, int N)
{
    //====================================================================================
    //Initialization
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Initialization of the parameters
    //------------------------------------------------------------------------------------
    int nvar = odestruct->dim;
    int status;               //current status
    double yv[nvar], t;       //current state and time
    int nt;

    //------------------------------------------------------------------------------------
    //Initialization of the odestruct
    //------------------------------------------------------------------------------------
    //Change sign of step if necessary
    if((tf < t0 && odestruct->d->h>0) || (tf > t0 && odestruct->d->h<0)) odestruct->d->h *= -1;
    if((tf < t0 && odestruct->h>0) || (tf > t0 && odestruct->h<0)) odestruct->h *= -1;

    //Reset ode structure.
    reset_ode_structure(odestruct);

    //------------------------------------------------------------------------------------
    //Initialization of the IC
    //------------------------------------------------------------------------------------
    //Init the state & time
    for(int k = 0; k < nvar; k++) yv[k] = y0[k];
    t = t0;
    nt = 1;

    //====================================================================================
    //Loop
    //====================================================================================
    //First position and time
    for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];
    tNCE[0] = t;

    //Other positions
    do
    {
        status = gsl_odeiv2_evolve_apply (odestruct->e, odestruct->c, odestruct->s, &odestruct->sys, &t, tf, &odestruct->h, yv);
        //Update the storage
        for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
        //Evolve up to ti
        tNCE[nt] = t;
        //Advance one step
        nt++;
    }
    while((t < tf) && (nt<=N) && (status == GSL_SUCCESS));

    //------------------------------------------------------------------------------------
    //Something went wrong inside the stepper?
    //------------------------------------------------------------------------------------
    if(status)
    {
        cout << "Warning in ode78_variable_grid: the stepper went wrong. No output is produced.";
        return FTC_FAILURE;
    }

    //------------------------------------------------------------------------------------
    //The maximum number of points is reached
    //------------------------------------------------------------------------------------
    if(nt == N)
    {
        cout << "Warning in ode78_variable_grid: the final time was not reached because the maximum number of points is reached." << endl;
    }


    return nt-1;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a given grid. The time grid is given as input
 **/
int ode78_grid_gg(OdeStruct *odestruct, double *y0, double **yNCE, const double *tNCEi, int N)
{
    //====================================================================================
    //Initialization
    //====================================================================================

    //------------------------------------------------------------------------------------
    //Initialization of the parameters
    //------------------------------------------------------------------------------------
    int nvar = odestruct->dim;
    int status;               //current status
    double yv[nvar], t;       //current state and time
    int nt;
    double ti;

    //------------------------------------------------------------------------------------
    //Initialization of the odestruct
    //------------------------------------------------------------------------------------
    //Change sign of step if necessary
    if((tNCEi[N] < tNCEi[0] && odestruct->d->h>0) || (tNCEi[N] > tNCEi[0] && odestruct->d->h<0)) odestruct->d->h *= -1;
    if((tNCEi[N] < tNCEi[0] && odestruct->h>0)    || (tNCEi[N] > tNCEi[0] && odestruct->h<0)) odestruct->h *= -1;

    //Reset ode structure.
    reset_ode_structure(odestruct);

    //------------------------------------------------------------------------------------
    //Initialization of the IC
    //------------------------------------------------------------------------------------
    //Init the state & time
    for(int k = 0; k < nvar; k++) yv[k] = y0[k];
    t = tNCEi[0];
    nt = 1;


    //====================================================================================
    //Loop
    //====================================================================================
    //First position
    for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];

    //Other positions
    do
    {
        //Evolve up to ti
        ti = tNCEi[nt];
        status = gsl_odeiv2_driver_apply (odestruct->d, &t , ti, yv);

        //Update the storage
        for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];

        //Advance one step
        nt++;
    }
    while((nt<=N) && (status == GSL_SUCCESS));

    //------------------------------------------------------------------------------------
    //Something went wrong inside the stepper?
    //------------------------------------------------------------------------------------
    if(status)
    {
        cout << "Warning in ode78_grid_gg: the stepper went wrong. No output is produced.";
        return FTC_FAILURE;
    }


    return FTC_SUCCESS;
}

//========================================================================================
//
//          Integration - SingleOrbit (deprecated)
//
//========================================================================================
/**
 *   \brief Integrates a given trajectory up to tf, on a given grid
 **/
int trajectory_integration_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------------
    int nvar = orbit.driver->dim;
    int status;               //current status
    double yv[nvar], t;       //current state and time

    //Projection tools
    double proj_dist_SEM;
    int nreset, nt;

    //Plot
    double ti;

    //------------------------------------------------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------------------------------------------------
    do
    {
        //Reset ode structure.
        reset_ode_structure(orbit.driver);

        //Change sign of step if necessary
        if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;
        if((tf < t0 && orbit.driver->d->h>0) || (tf > t0 && orbit.driver->d->h<0)) orbit.driver->d->h *= -1;

        //Init the state & time
        for(int i = 0; i < nvar; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 1;
        nreset = 1;

        //First position
        for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];
        tNCE[0] = t0;

        //Loop
        do
        {

            ti = t0 + (double) nt *(tf-t0)/(1.0*N);
            if(isResetOn) status = gslc_proj_evolve(orbit, yv, &t, t0, ti, &proj_dist_SEM, &nreset, isResetOn);
            else  status = gsl_odeiv2_driver_apply (orbit.driver->d, &t, ti, yv);

            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = ti;

            //Advance one step
            nt++;
        }
        while((nt<=N) && (status == FTC_SUCCESS) && (orbit.tproj > orbit.tprojmin));

        //If a new reset is necessary
        if (status == ORBIT_EPROJ && isResetOn)
        {
            cout << "Warning in trajectory_integration_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

        //Something went wrong inside the stepper
        if(status == ORBIT_EINT)
        {
            cout << "Warning in trajectory_integration_grid: the stepper went wrong. No output is produced.";
            return FTC_FAILURE;
        }

        //Something unknown went wrong
        if(status == FTC_FAILURE)
        {
            cout << "Warning in trajectory_integration_grid: something went wrong. No output is produced.";
            return FTC_FAILURE;
        }

    }
    while(status!= FTC_SUCCESS && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return FTC_FAILURE;
    }

    return FTC_SUCCESS;
}

/**
 *   \brief Integrates a given trajectory up to tf, on a variable grid of maximum size N. Return the last position that is filled on the grid.
 **/
int trajectory_integration_variable_grid(SingleOrbit &orbit, double t0, double tf, double **yNCE, double *tNCE, int N, int isResetOn)
{
    //------------------------------------------------------------------------------------
    //Initialization
    //------------------------------------------------------------------------------------
    int nvar = orbit.driver->dim;   //number of variables
    int status;                     //current status
    double yv[nvar], t;             //current state and time

    //Projection tools
    double proj_dist;
    int nreset, nt;

    //------------------------------------------------------------------------------------
    //Evolving yv(t) up to tf
    //------------------------------------------------------------------------------------
    do
    {
        //Change sign of step if necessary
        if((tf < t0 && orbit.driver->h>0) || (tf > t0 && orbit.driver->h<0)) orbit.driver->h *= -1;

        //Reset ode structure.
        reset_ode_structure(orbit.driver);

        //Init the state & time
        for(int i = 0; i < nvar; i++) yv[i] = orbit.z0[i];
        t = t0;
        nt = 1;
        nreset = 1;

        //First step
        for(int k = 0; k < nvar; k++) yNCE[k][0] = yv[k];
        tNCE[0] = t0;

        //Loop
        do
        {
            status = gslc_proj_step(orbit, yv, &t, t0, tf, &proj_dist, &nreset, isResetOn);
            for(int k = 0; k < nvar; k++) yNCE[k][nt] = yv[k];
            tNCE[nt] = t;
            //Advance one step
            nt++;
        }
        while((t != tf) && (nt <= N) && (status == 0) && (orbit.tproj > orbit.tprojmin));

        if (status == ORBIT_EPROJ && isResetOn) //something went wrong during projection
        {
            cout << "Warning in trajectory_integration_variable_grid: the interval of projection has to be reduced: ";
            cout << setprecision(3) << "orbit.tproj : " << orbit.tproj << " -> ";
            orbit.tproj *= 0.5;
            cout << orbit.tproj << setprecision(15) << endl;
        }

        if(status == ORBIT_EINT) //something went wrong inside the stepper
        {
            cout << "Warning in trajectory_integration_variable_grid: the stepper went wrong. No output is produced.";
            return FTC_FAILURE;
        }

        //Something unknown went wrong
        if(status == FTC_FAILURE)
        {
            cout << "Warning in trajectory_integration_variable_grid: something went wrong. No output is produced.";
            return FTC_FAILURE;
        }

        //The maximum number of points is reached
        if(nt == N)
        {
            cout << "Warning in trajectory_integration_variable_grid: the final time was not reached because the maximum number of points is reached." << endl;
        }

    }
    while(status!= FTC_SUCCESS && orbit.tproj > orbit.tprojmin);

    if(orbit.tproj < orbit.tprojmin)
    {
        cout << "Error in trajectory_integration_grid: the interval of projection is too small." << endl;
        return FTC_FAILURE;
    }

    return nt-1;
}

/**
 *   \brief Integrates one step the current state yv using projection on the CM if necessary
 **/
int gslc_proj_step(SingleOrbit &orbit,
                   double yv[],
                   double *t,
                   double t0,
                   double t1,
                   double *proj_dist_SEM,
                   int *nreset,
                   int isResetOn)
{
    int status;
    double yvp[6], yvi[6];
    cdouble scp[orbit.reduced_nv];

    //----------------------
    //Evolve one step of z(t)
    //----------------------
    orbit.driver->h = (t1 >= *t)? fabs(orbit.driver->h):-fabs(orbit.driver->h);
    status = gsl_odeiv2_evolve_apply (orbit.driver->e, orbit.driver->c, orbit.driver->s, &orbit.driver->sys, t, t1, &orbit.driver->h, yv);
    if (status)
    {
        cout << "error in gslc_proj_step: integration of z(t) has gone wrong. break." << endl;
        return ORBIT_EINT;
    }

    //----------------------
    //Projection if necessary
    //----------------------
    if(isResetOn && fabs(*t-t0) > fabs(*nreset*orbit.tproj))
    {
        //----------------------
        //Projection tools
        //----------------------
        double omega1 = cimag(orbit.Wh->at(0).getCA(1,0)->getCoef(0));
        double omega3 = cimag(orbit.Wh->at(2).getCA(1,1)->getCoef(0));

        //----------------------
        // Projection on the center manifold
        //----------------------
        //Get the closest point on the center manifold
        NCprojCCM(yv, *t, orbit.n, OFS_ORDER, *orbit.CQ, *orbit.V, omega1, omega3, scp, orbit.reduced_nv);
        //Update the state
        CCMtoNCbyTFC(scp, *t, orbit.n, orbit.order,  orbit.ofs_order,  *orbit.Wh,  *orbit.PC, *orbit.V,  yvp,  orbit.isGS);
        //For comparison
        for(int i = 0; i <6; i++) yvi[i] = yv[i];
        // Copy of yvp in current state
        for(int i=0; i<6; i++) yv[i]  = yvp[i];

        //-----------------
        // Get the current projection error
        //-----------------
        //Get the current error
        *proj_dist_SEM = (yvi[0] - yv[0])*(yvi[0] - yv[0]);
        for(int i = 1; i < 6 ; i++)
        {
            *proj_dist_SEM += (yvi[i] - yv[i])*(yvi[i] - yv[i]);
        }
        *proj_dist_SEM  = sqrt(*proj_dist_SEM);

        if(*proj_dist_SEM > orbit.ePmax)
        {
            cout << "Warning: Reset n° " << *nreset << ". Error (NC) = " << *proj_dist_SEM << endl;
            cout << "Error (SYS) = " << *proj_dist_SEM*orbit.qbcp_l->cs->gamma << endl;
            cout << "Error (km) = " << *proj_dist_SEM*orbit.qbcp_l->cs->gamma*orbit.qbcp_l->cs->cr3bp.L << endl;
            return ORBIT_EPROJ;
        }

        //-----------------
        //Reset ode structure for next step
        //-----------------
        //reset_ode_structure(orbit.driver);

        //-----------------
        //One additional reset
        //-----------------
        *nreset = *nreset +1;
    }


    return FTC_SUCCESS;
}


/**
 *   \brief Integrates the current state yv up to t = t1, using projection on the CM if necessary
 **/
int gslc_proj_evolve(SingleOrbit &orbit,
                     double yv[],
                     double *t,
                     double t0,
                     double t1,
                     double *proj_dist_SEM,
                     int *nreset,
                     int isResetOn)
{
    //reset_ode_structure(orbit.driver);
    int status;
    do
    {
        status = gslc_proj_step(orbit, yv, t, t0, t1, proj_dist_SEM, nreset, isResetOn);
    }
    while(status == FTC_SUCCESS && fabs(*t - t1) > 1e-16);

    return status;
}


//========================================================================================
//
//          Projection on (un)stable manifold - deprecated
//
//========================================================================================
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CM_TFC, MIcoc, Vcoc). Then the CCM state is extended by adding a non-null direction
 *        along the hyperbolic direction (sti[4]).
 **/
void NCprojCCMtoCUS(double *yv, double tv, double n, double sti[5], vector<Oftsc> &CM_TFC, double epsilon, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc)
{
    //Projection tools
    double omega1 = cimag(CM_TFC[0].getCA(1,0)->getCoef(0));
    double omega3 = cimag(CM_TFC[2].getCA(1,1)->getCoef(0));
    cdouble scp[5];
    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, 5);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, 5);
    //Add a given quantity on the hyperbolic direction
    sti[4] = epsilon;
}

/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given by (CM_TFC, MIcoc, Vcoc).
 **/
void NCprojCCMtoCM(double *yv, double tv, double n, double sti[5], vector<Oftsc> &CM_TFC, matrix<Ofsc>  &MIcoc, vector<Ofsc>  &Vcoc)
{
    //Projection tools
    double omega1 = cimag(CM_TFC[0].getCA(1,0)->getCoef(0));
    double omega3 = cimag(CM_TFC[2].getCA(1,1)->getCoef(0));
    cdouble scp[4];
    //Get the closest point on the center manifold, in scp[4]
    NCprojCCM(yv, tv, n, OFS_ORDER, MIcoc, Vcoc, omega1, omega3, scp, 4);
    //Get the correspondance in RCM coordinates
    CCMtoRCM(scp, sti, 4);
}


//========================================================================================
//
//                  Energy on vectors
//
//========================================================================================
/**
 *  \brief Hamiltonian along one orbit in SEM units and SEM coordinates
 **/
void HSEM_vec(double *tSEM, double **ySEM, double *Hvec, int N, QBCP_L *qbcp_l)
{
    double yv[6];
    for(int i = 0; i <= N; i++)
    {
        for(int p = 0; p < 6; p++) yv[p] = ySEM[p][i];
        Hvec[i] = qbcp_H_SEM(tSEM[i], yv, qbcp_l); //careful: time in SEM units!
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
        NCSEMmtoSEMm(tSEM[i], yv, yvSEM, qbcp_l);
        Hvec[i] = qbcp_H_SEM(tSEM[i], yvSEM, qbcp_l);
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
        Hvec[i] = qbcp_H_SEM(tEM[i]*qbcp_l->us_em.ns, yvSEM, qbcp_l); //careful, time must be in SEM units at this step!
    }
}



//========================================================================================
//
//         Update points
//
//========================================================================================
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

    //SEMLi position in EM coordinates
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
    //SEMLi position
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

/**
 *  \brief Update some key positions of notable points in the EM system (Normalized version)
 **/
void emNCPoints(double t, double **emP)
{
    double temp1[6], temp2[6];
    //----------------------------------
    //All in EM coordinates, at time t in EM units
    //----------------------------------
    //Earth position
    evaluateCoef(emP[0], t, SEML.us_em.n, SEML.nf, SEML.cs_em.pe, 3);
    //Moon position
    evaluateCoef(emP[1], t, SEML.us_em.n, SEML.nf, SEML.cs_em.pm, 3);
    //Sun position
    evaluateCoef(emP[2], t, SEML.us_em.n, SEML.nf, SEML.cs_em.ps, 3);

    //EML1 position
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l1.position[i];
    temp1[0] *= -1.0;
    EMmtoNCEMm(t, temp1, temp2, &SEML);
    for(int i = 0; i < 3; i++) emP[3][i] = temp2[i];

    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l2.position[i];
    temp1[0] *= -1.0;
    EMmtoNCEMm(t, temp1, temp2, &SEML);
    for(int i = 0; i < 3; i++) emP[4][i] = temp2[i];


    //SEML1 position in EM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l1.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in EM coordinates, starting with a time in SEM units !
    SEMmtoNCEMm(t*SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) emP[5][i] = temp2[i];

    //SEMLi position in EM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l2.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in EM coordinates, starting with a time in SEM units !
    SEMmtoNCEMm(t*SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) emP[6][i] = temp2[i];
}

/**
 *  \brief Update some key positions of notable points in the SEM system (Normalized version)
 **/
void semNCPoints(double t, double **semP)
{
    double temp1[6], temp2[6];
    //----------------------------------
    //All in SEM coordinates, at time t in SEM units
    //----------------------------------
    //Earth position
    evaluateCoef(semP[0], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.pe, 3);
    //Moon position
    evaluateCoef(semP[1], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.pm, 3);
    //Sun position
    evaluateCoef(semP[2], t, SEML.us_sem.n, SEML.nf, SEML.cs_sem.ps, 3);

    //SEML1 position
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l1.position[i];
    temp1[0] *= -1.0;
    qbcp_coc(t, temp1, temp2, PSEM, NCSEM);
    for(int i = 0; i < 3; i++) semP[5][i] = temp2[i];

    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_sem.cr3bp.l2.position[i];
    temp1[0] *= -1.0;
    qbcp_coc(t, temp1, temp2, PSEM, NCSEM);
    for(int i = 0; i < 3; i++) semP[6][i] = temp2[i];


    //EML1 position in SEM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l1.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in SEM coordinates, starting with a time in EM units !
    EMmtoNCSEMm(t/SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) semP[3][i] = temp2[i];

    //EML2 position in SEM coordinates
    for(int i = 0; i < 3; i++) temp1[i] = SEML.cs_em.cr3bp.l2.position[i];
    for(int i = 4; i < 6; i++) temp1[i] = 0.0;  //arbitrary momenta
    temp1[0] *= -1.0;
    //Back in SEM coordinates, starting with a time in EM units !
    EMmtoNCSEMm(t/SEML.us_em.ns, temp1, temp2, &SEML);
    //Store in emP
    for(int i = 0; i < 3; i++) semP[4][i] = temp2[i];
}

//========================================================================================
//
//        Plots
//
//========================================================================================
/**
 *  \brief Sets notable points in EM system on the gnuplot window ctrl h1
 **/
void emPlot(gnuplot_ctrl *h1, double **emP, int isPlot)
{
    //Misc
    gnuplot_cmd(isPlot, h1, "set grid");
    gnuplot_set_xlabel(isPlot, h1, (char*)"Xem [-]");
    gnuplot_set_ylabel(isPlot, h1, (char*)"Yem [-]");
    gnuplot_set_zlabel(isPlot, h1, (char*)"Zem [-]");
    //EML point
    switch(SEML.li_EM)
    {
    case 1:
        gnuplot_plot_xyz(isPlot, h1, emP[3], emP[3]+1,  emP[3]+2, 1, (char*)"EML1", "points", "1", "2", 8);
        break;
    case 2:
        gnuplot_plot_xyz(isPlot, h1, emP[4], emP[4]+1,  emP[4]+2, 1, (char*)"EML2", "points", "1", "2", 8);
        break;
    }
    //Primaries
    gnuplot_plot_xyz(isPlot, h1, emP[1], emP[1]+1,  emP[1]+2, 1, (char*)"Moon", "points", "7", "2", 9);
}

/**
 *  \brief Sets notable points in SEM system on the gnuplot window ctrl h2
 **/
void semPlot(gnuplot_ctrl *h2, double **semP, int isPlot)
{
    //Misc
    gnuplot_cmd(isPlot, h2, "set grid");
    gnuplot_set_xlabel(isPlot, h2, (char*)"Xsem [-]");
    gnuplot_set_ylabel(isPlot, h2, (char*)"Ysem [-]");
    gnuplot_set_zlabel(isPlot, h2, (char*)"Zsem [-]");

    // SEML points
    switch(SEML.li_SEM)
    {
    case 1:
        gnuplot_plot_xyz(isPlot, h2, semP[5], semP[5]+1,  semP[5]+2, 1, (char*)"SEML1", "points", "1", "2", 9);
        break;
    case 2:
        gnuplot_plot_xyz(isPlot, h2, semP[6], semP[6]+1,  semP[6]+2, 1, (char*)"SEMLi", "points", "2", "2", 9);
        break;
    }
    //EML points
    gnuplot_plot_xyz(isPlot, h2, semP[3], semP[3]+1,  semP[3]+2, 1, (char*)"EML1", "points", "1", "2", 8);
    gnuplot_plot_xyz(isPlot, h2, semP[4], semP[4]+1,  semP[4]+2, 1, (char*)"EML2", "points", "2", "2", 8);
    //Primaries
    gnuplot_plot_xyz(isPlot, h2, semP[0], semP[0]+1,  semP[0]+2, 1, (char*)"Earth", "points", "7", "2", 6);
}

//========================================================================================
//
//          I/O orbit
//
//========================================================================================
/**
 *  \brief Computes a data filename as a string, depending on the ofts_order, size, and type of the data.
 **/
string filenameOrbit(int ofts_order, int sizeOrbit, int type)
{
    switch(type)
    {
    case TYPE_ORBIT:
        return SEML.cs->F_PLOT+"orbit_order_"+numTostring(ofts_order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_CU:
        return SEML.cs->F_PLOT+"cu_order_"+numTostring(ofts_order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_CS:
        return SEML.cs->F_PLOT+"cs_order_"+numTostring(ofts_order)+"_size_"+numTostring(sizeOrbit)+".txt";
    case TYPE_MAN_PROJ:
        return SEML.cs->F_PLOT+"man_proj_"+numTostring(ofts_order)+"_size_"+numTostring(sizeOrbit)+".bin";
    default:
        cout << "filenameOrbit: unknown type." << endl;
        return "";
    }
}

/**
 * \brief Store the orbit (tNCE, yNCE) in the  data file filenameOrbit(ofts_order, sizeOrbit, type)
 **/
int writeOrbit(double *tNCE, double **yNCE, int ofts_order, int sizeOrbit, int N, int type)
{
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);
    gnuplot_fplot_txp(tNCE, yNCE[0], yNCE[1], yNCE[2], yNCE[3], yNCE[4], yNCE[5], N, filename.c_str());
    return FTC_SUCCESS;
}


/**
 * \brief Get the length of the data file filenameOrbit(ofts_order, sizeOrbit, type).
 **/
int getLineNumber(int ofts_order, int sizeOrbit, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    string line;
    ifstream filestream (filename.c_str());
    int Nfile = 0;
    if (filestream.is_open())
    {
        while ( getline (filestream,line) ) Nfile++;
        filestream.close();
    }
    else return FTC_FAILURE;

    return Nfile;
}

/**
 * \brief Read the orbit (tNCE, yNCE) in the  data file filenameOrbit(ofts_order, size, type)
 **/
int readOrbit(double *tNCE, double **yNCE, int ofts_order, int sizeOrbit, int N, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    int Nfile = getLineNumber(ofts_order, sizeOrbit, type);
    //If N is greater than Nfile, the end of the vectors will not be updated
    if(Nfile < N)
    {
        cout << "readOrbit: wrong input, N = " << N << ", but Nfile = " << Nfile << endl;
        return FTC_FAILURE;
    }
    else Nfile = N;

    //------------------------
    //Updating the vectors
    //------------------------
    ifstream filestream (filename.c_str());
    if (filestream.is_open())
    {
        for(int i = 0; i < Nfile; i++)
        {
            filestream >> tNCE[i];
            for(int k = 0; k < 6; k++) filestream >> yNCE[k][i];
        }
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

/**
 * \brief Read the orbit (tNCE, yNCE, sNCE) in the  data file filenameOrbit(ofts_order, sizeOrbit, type)
 **/
int readOrbit(double *tNCE, double **yNCE, double **sNCE, int ofts_order, int sizeOrbit, int N, int type)
{
    //------------------------
    //Name of the file
    //------------------------
    string filename = filenameOrbit(ofts_order, sizeOrbit, type);

    //------------------------
    //Getting the sizeOrbit of the file
    //------------------------
    int Nfile = getLineNumber(ofts_order, sizeOrbit, type);
    //If N is greater than Nfile, the end of the vectors will not be updated
    if(Nfile < N)
    {
        cout << "readOrbit: wrong input, N = " << N << ", but Nfile = " << Nfile << endl;
        return FTC_FAILURE;
    }
    else Nfile = N;

    //------------------------
    //Updating the vectors
    //------------------------
    ifstream filestream (filename.c_str());
    if (filestream.is_open())
    {
        for(int i = 0; i < Nfile; i++)
        {
            filestream >> tNCE[i];
            for(int k = 0; k < 6; k++) filestream >> yNCE[k][i];
            for(int k = 0; k < 5; k++) filestream >> sNCE[k][i];
        }
    }
    else return FTC_FAILURE;

    return FTC_SUCCESS;
}

//========================================================================================
//
//          QBCP test
//
//========================================================================================
/**
 *  \brief Derivatives of the QBCP in EM inertial coordinates (primaries + fourth body motion)
 *          - First 8 variable is the primaries' motion (z, Z in real form).
 *          - Last  6 variables is the state (Xin, Vin).
 */
int qbcp_derivatives_em_in(double t, const double y[], double f[], void *params)
{
    //=======================================================================
    //Retrieving the parameters
    //=======================================================================
    QBCP_L* qbp = (QBCP_L *) params;
    double ms   = qbp->us_em.ms;
    double me   = qbp->us_em.me;
    double mm   = qbp->us_em.mm;
    double mu   = qbp->us_em.mu_EM;

    //=======================================================================
    // 1. Primaries' motion
    //=======================================================================
    //Reconstruction of z and Z
    cdouble z = y[0] + I*y[1];
    cdouble Z = y[2] + I*y[3];

    cdouble temp1 = Z-mu*z;
    temp1 = temp1/pow(cabs(temp1), 3.0);

    cdouble temp2 = Z+(1-mu)*z;
    temp2 = temp2/pow(cabs(temp2), 3.0);

    cdouble zdd = +0.0*I-z/pow(cabs(z), 3.0) + ms*(temp1-temp2);
    cdouble Zdd = -(1+ms)*(mu*temp2 + (1-mu)*temp1);

    //-----------------------
    //Phase space derivatives
    //-----------------------
    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];
    f[4] = creal(zdd);
    f[5] = cimag(zdd);
    f[6] = creal(Zdd);
    f[7] = cimag(Zdd);


    //=======================================================================
    // 1. S/C's motion
    //=======================================================================
    int shift = 8;
    double yIE[3], yIM[3], yIS[3];
    //Earth position
    yIE[0] = -ms/(1.0+ms)*creal(Z) + mu*creal(z);
    yIE[1] = -ms/(1.0+ms)*cimag(Z) + mu*cimag(z);
    yIE[2] = 0.0;

    //Moon position
    yIM[0] = -ms/(1.0+ms)*creal(Z) - (1-mu)*creal(z);
    yIM[1] = -ms/(1.0+ms)*cimag(Z) - (1-mu)*cimag(z);
    yIM[2] = 0.0;

    //Sun position
    yIS[0] = 1.0/(1.0+ms)*creal(Z);
    yIS[1] = 1.0/(1.0+ms)*cimag(Z);
    yIS[2] = 0.0;

    //Distance to Earth, Moon and Sun
    double dIE2 = 0.0, dIM2 = 0.0, dIS2 = 0.0;
    for(int i = 0; i < 3; i++)
    {
        dIE2 += (y[i+shift] - yIE[i])*(y[i+shift] - yIE[i]);
        dIM2 += (y[i+shift] - yIM[i])*(y[i+shift] - yIM[i]);
        dIS2 += (y[i+shift] - yIS[i])*(y[i+shift] - yIS[i]);
    }

    //-----------------------
    //Phase space derivatives
    //-----------------------
    f[0+shift] = y[3+shift];
    f[1+shift] = y[4+shift];
    f[2+shift] = y[5+shift];

    f[3+shift] = - me/pow(dIE2, 3.0/2)*(y[0 + shift] - yIE[0])
                 - mm/pow(dIM2, 3.0/2)*(y[0 + shift] - yIM[0])
                 - ms/pow(dIS2, 3.0/2)*(y[0 + shift] - yIS[0]);

    f[4+shift] = - me/pow(dIE2, 3.0/2)*(y[1 + shift] - yIE[1])
                 - mm/pow(dIM2, 3.0/2)*(y[1 + shift] - yIM[1])
                 - ms/pow(dIS2, 3.0/2)*(y[1 + shift] - yIS[1]);

    f[5+shift] = - me/pow(dIE2, 3.0/2)*(y[2 + shift] - yIE[2])
                 - mm/pow(dIM2, 3.0/2)*(y[2 + shift] - yIM[2])
                 - ms/pow(dIS2, 3.0/2)*(y[2 + shift] - yIS[2]);

    return GSL_SUCCESS;
}


/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbcp_test(QBCP_L &qbcp_l)
{
    double yEx[6];
    //=======================================================================
    // 1. Initialization
    //=======================================================================
    //IC in NCEM coordinates
    double yvNCEM[6];
    yvNCEM[0] = -3.190475469637594e-02;
    yvNCEM[1] = -3.459923857900106e-01;
    yvNCEM[2] = +0.000000000000000e+00;
    yvNCEM[3] = +1.137240172435051e-01;
    yvNCEM[4] = -7.499822805182552e-02;
    yvNCEM[5] = +0.000000000000000e+00;

    //Initial an final times in EM units
    double t0EM = +6.016997770510415e+00;
    double tfEM = t0EM + 3*SEML.us_em.T;

    //Size of the vectors
    int man_grid_size = 10000;

    //Collisionner (not used here)
    int ode78coll;

    //=======================================================================
    // 2. Integration in NCEM system
    //=======================================================================
    double **ymNCEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmNCEM_IN   = dvector(0, man_grid_size);
    ode78(ymNCEM_IN, tmNCEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, I_NCEM, NCEM, INEM);

    //=======================================================================
    // 2. Integration in EM system
    //=======================================================================
    double **ymEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmEM_IN   = dvector(0, man_grid_size);
    ode78(ymEM_IN, tmEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, I_PEM, NCEM, INEM);

    //=======================================================================
    // 3. Integration in VNCEM system
    //=======================================================================
    double **ymVNCEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmVNCEM_IN   = dvector(0, man_grid_size);
    ode78(ymVNCEM_IN, tmVNCEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, I_VNCEM, NCEM, INEM);

    //=======================================================================
    // 4. Integration in NCSEM system
    //=======================================================================
    double **ymNCSEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmNCSEM_IN   = dvector(0, man_grid_size);
    ode78(ymNCSEM_IN, tmNCSEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, I_NCSEM, NCEM, INEM);

    //=======================================================================
    // 5. Integration in SEM system
    //=======================================================================
    double **ymSEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmSEM_IN   = dvector(0, man_grid_size);
    ode78(ymSEM_IN, tmSEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, F_SEM, NCEM, INEM);

    //=======================================================================
    // 6. Integration in VNCSEM system
    //=======================================================================
    double **ymVNCSEM_IN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmVNCSEM_IN   = dvector(0, man_grid_size);
    ode78(ymVNCSEM_IN, tmVNCSEM_IN, &ode78coll, t0EM, tfEM, yvNCEM, 6, man_grid_size, I_VNCSEM, NCEM, INEM);

    //=======================================================================
    // 7. Integration in "true" system
    //=======================================================================
    double **ymIN  = dmatrix(0, 5, 0, man_grid_size);
    double *tmIN   = dvector(0, man_grid_size);

    double n   = qbcp_l.us_em.n;
    double ns  = qbcp_l.us_em.ns;
    double as  = qbcp_l.us_em.as;
    double ni  = qbcp_l.us_em.ni;
    double ai  = qbcp_l.us_em.ai;

    //--------------------------------------
    // 7.1. Init the integration tools
    //--------------------------------------
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    OdeStruct ode_s;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method
    //Parameters
    OdeParams odeParams(&SEML, I_INEM);
    //General structures
    init_ode_structure(&ode_s, T, T_root, 14, qbcp_derivatives_em_in, &odeParams);

    //--------------------------------------
    // 7.2. Initital conditions
    //--------------------------------------
    double t = t0EM;
    //z(0) and Z(0)
    cdouble z0    = evz(qbcp_l.cs_em.zt, t0EM, n, ni, ai);
    cdouble Z0    = evz(qbcp_l.cs_em.Zt, t0EM, n, ns, as);
    cdouble zdot0 = evzdot(qbcp_l.cs_em.zt, qbcp_l.cs_em.ztdot, t0EM, n, ni, ai);
    cdouble Zdot0 = evzdot(qbcp_l.cs_em.Zt, qbcp_l.cs_em.Ztdot, t0EM, n, ns, as);

    //Put the IC in real form
    double yv[14];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    //Initial state of the S/C
    int shift = 8;
    double yINEM[6];
    qbcp_coc(t0EM, yvNCEM, yINEM, NCEM, INEM);
    for(int i = 0; i <6; i++) yv[i+shift] = yINEM[i];


    //--------------------------------------
    // 7.3. Integration on a loop of size man_grid_size
    //--------------------------------------
    double ti;
    for(int k = 0; k <= man_grid_size; k++)
    {
        reset_ode_structure(&ode_s);
        //--------------------------------------
        //Integration from t to ti
        //--------------------------------------
        ti = t0EM + 1.0*(tfEM - t0EM)*k/man_grid_size;
        if(k > 0) gsl_odeiv2_driver_apply (ode_s.d, &t , ti, yv);

        //--------------------------------------
        //Store
        //--------------------------------------
        for(int i = 0; i < 6; i++) ymIN[i][k] = yv[i+shift];
        tmIN[k] = (ti - t0EM)/qbcp_l.us_em.T;
    }


    //=======================================================================
    // 8. Postprocess
    //=======================================================================
    int kep = 6;
    double *dNCEM_IN   = dvector(0, man_grid_size);
    double *dEM_IN     = dvector(0, man_grid_size);
    double *dVNCEM_IN  = dvector(0, man_grid_size);
    double *dNCSEM_IN  = dvector(0, man_grid_size);
    double *dSEM_IN    = dvector(0, man_grid_size);
    double *dVNCSEM_IN = dvector(0, man_grid_size);
    double *dSEM_EM    = dvector(0, man_grid_size);
    //Loop of postprocess
    for(int k = 0; k <= man_grid_size; k++)
    {
        //NCEM
        for(int i = 0; i < kep; i++) yEx[i] = ymNCEM_IN[i][k] - ymIN[i][k];
        dNCEM_IN[k] = ENorm(yEx, kep);

        //EM
        for(int i = 0; i < kep; i++) yEx[i] = ymEM_IN[i][k] - ymIN[i][k];
        dEM_IN[k] = ENorm(yEx, kep);

        //VNCEM
        for(int i = 0; i < kep; i++) yEx[i] = ymVNCEM_IN[i][k] - ymIN[i][k];
        dVNCEM_IN[k] = ENorm(yEx, kep);

        //NCSEM
        for(int i = 0; i < kep; i++) yEx[i] = ymNCSEM_IN[i][k] - ymIN[i][k];
        dNCSEM_IN[k] = ENorm(yEx, kep);

        //SEM
        for(int i = 0; i < kep; i++) yEx[i] = ymSEM_IN[i][k] - ymIN[i][k];
        dSEM_IN[k] = ENorm(yEx, kep);

        //VNCSEM
        for(int i = 0; i < kep; i++) yEx[i] = ymVNCSEM_IN[i][k] - ymIN[i][k];
        dVNCSEM_IN[k] = ENorm(yEx, kep);

        //NCSEM vs NCEM
        for(int i = 0; i < kep; i++) yEx[i] = ymNCSEM_IN[i][k] - ymNCEM_IN[i][k];
        dSEM_EM[k] = ENorm(yEx, kep);
    }


    //=======================================================================
    // 8. Plot
    //=======================================================================
    gnuplot_ctrl *h1;
    h1 = gnuplot_init(true);
    gnuplot_cmd(true, h1, "set title \"Variations of the error\" ");
    gnuplot_cmd(true, h1, "set grid");
    gnuplot_set_xlabel(true, h1, (char*)"t [x Tsem]");
    gnuplot_set_ylabel(true, h1, (char*)"Error [EM units]");
    gnuplot_cmd(true, h1, "set logscale y");
    gnuplot_cmd(true, h1, "set format y \"1e\%%L\"");

    int color = 1;
    gnuplot_plot_xy(true, h1, tmIN, dNCEM_IN, man_grid_size+1,   (char*)"NCEM",  "lines", "3", "2", color++);
    gnuplot_plot_xy(true, h1, tmIN, dVNCEM_IN, man_grid_size+1,  (char*)"VNCEM", "lines", "3", "2", color++);
    gnuplot_plot_xy(true, h1, tmIN, dEM_IN, man_grid_size+1,     (char*)"EM",    "lines", "3", "2", color++);
    gnuplot_plot_xy(true, h1, tmIN, dNCSEM_IN, man_grid_size+1,  (char*)"NCSEM",  "lines", "dashed", "2", color++);
    gnuplot_plot_xy(true, h1, tmIN, dVNCSEM_IN, man_grid_size+1, (char*)"VNCSEM", "lines", "dashed", "2", color++);
    gnuplot_plot_xy(true, h1, tmIN, dSEM_IN, man_grid_size+1,    (char*)"SEM",    "lines", "dashed", "2", color++);

    //=======================================================================
    // 9. fPlot
    //=======================================================================
    gnuplot_fplot_xy(tmIN, dNCEM_IN, man_grid_size+1,  (char*) (qbcp_l.cs->F_PLOT+"QBCP_NCEM_vs_IN.txt").c_str());
    gnuplot_fplot_xy(tmIN, dEM_IN, man_grid_size+1,    (char*) (qbcp_l.cs->F_PLOT+"QBCP_EM_vs_IN.txt").c_str());
    gnuplot_fplot_xy(tmIN, dVNCEM_IN, man_grid_size+1, (char*) (qbcp_l.cs->F_PLOT+"QBCP_VNCEM_vs_IN.txt").c_str());
    gnuplot_fplot_xy(tmIN, dNCSEM_IN, man_grid_size+1,  (char*) (qbcp_l.cs->F_PLOT+"QBCP_NCSEM_vs_IN.txt").c_str());
    gnuplot_fplot_xy(tmIN, dSEM_IN, man_grid_size+1,    (char*) (qbcp_l.cs->F_PLOT+"QBCP_SEM_vs_IN.txt").c_str());
    gnuplot_fplot_xy(tmIN, dVNCSEM_IN, man_grid_size+1, (char*) (qbcp_l.cs->F_PLOT+"QBCP_VNCSEM_vs_IN.txt").c_str());

    char ch;
    gnuplot_cmd(true, h1, "set logscale y");
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(true, h1);

    //=======================================================================
    // 10. Free
    //=======================================================================
    free_dmatrix(ymNCEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmNCEM_IN, 0, man_grid_size);
    free_dmatrix(ymEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmEM_IN, 0, man_grid_size);
    free_dmatrix(ymVNCEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmVNCEM_IN, 0, man_grid_size);
    free_dmatrix(ymNCSEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmNCSEM_IN, 0, man_grid_size);
    free_dmatrix(ymSEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmSEM_IN, 0, man_grid_size);
    free_dmatrix(ymVNCSEM_IN, 0, 5, 0, man_grid_size);
    free_dvector(tmVNCSEM_IN, 0, man_grid_size);
    free_dmatrix(ymIN, 0, 5, 0, man_grid_size);
    free_dvector(tmIN, 0, man_grid_size);
    free_dvector(dNCEM_IN, 0, man_grid_size);
    free_dvector(dEM_IN, 0, man_grid_size);
    free_dvector(dVNCEM_IN, 0, man_grid_size);
    free_dvector(dNCSEM_IN, 0, man_grid_size);
    free_dvector(dSEM_IN, 0, man_grid_size);
    free_dvector(dVNCSEM_IN, 0, man_grid_size);
    free_dvector(dSEM_EM, 0, man_grid_size);

}


//========================================================================================
//
//          Orbit on a grid - deprecated
//
//========================================================================================
/**
 *  \brief Computes an orbit on a given time grid [t0 t1], with initial conditions st0, on the center manifold CM.
 **/
int gridOrbit(double st0[], double t0, double tf, double dt, vector<Oftsc> &CM,
              vector<Oftsc> &CM_TFC, matrix<Ofsc>  &Mcoc, matrix<Ofsc>  &MIcoc,
              vector<Ofsc>  &Vcoc)
{
    //------------------------------------------------------------------------------------
    // Computation
    //------------------------------------------------------------------------------------
    int dcs = SEML.fwrk == F_EM ? I_NCEM:I_NCSEM;

    //------------------
    // ODE
    //------------------
    OdeStruct odestruct;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, dcs);
    //Init ode structure
    init_ode_structure(&odestruct, T, T_root, 6, qbcp_vfn, &odeParams);

    //------------------
    //Orbit structure
    //------------------
    Ofsc orbit_ofs(OFS_ORDER);
    SingleOrbit orbit;

    //The default interval of projection is set to Tproj = T/5
    double tproj = SEML.us->T/5;

    //Init routine
    init_orbit(orbit, &CM, &CM_TFC, &Mcoc, &MIcoc, &Vcoc, &orbit_ofs, OFTS_ORDER, OFS_ORDER, 1, t0, tf, tproj, &odestruct, &SEML);

    //------------------
    //Orbit IC
    //------------------
    orbit_update_ic(orbit, st0, orbit.t0);

    //------------------------------------------------------------------------------------
    //Plotting parameters
    //------------------------------------------------------------------------------------
    int N         = floor(fabs(orbit.tf - orbit.t0)/dt);
    double **yNCE = dmatrix(0, 5, 0, N);
    double **ySYS = dmatrix(0, 5, 0, N);
    double *tNCE  = dvector(0, N);

    //------------------------------------------------------------------------------------
    //Integration
    //------------------------------------------------------------------------------------
    int status = trajectory_integration_grid(orbit, orbit.t0, orbit.tf, yNCE, tNCE, N, 1);

    //------------------------------------------------------------------------------------
    //To SYS coordinates
    //------------------------------------------------------------------------------------
    NCtoSYS_vec(yNCE, tNCE, ySYS,N, &SEML);

    //------------------------------------------------------------------------------------
    //Plotting
    //------------------------------------------------------------------------------------
    gnuplot_ctrl  *h1;
    h1 = gnuplot_init(true);

    gnuplot_setstyle(true, h1, (char*)"lines");
    gnuplot_set_xlabel(true, h1, (char*)"x [-]");
    gnuplot_set_ylabel(true, h1, (char*)"y [-]");
    gnuplot_plot_xyz(true, h1, ySYS[0], ySYS[1], ySYS[2], N+1, (char*)"NC coordinates", "lines", "1", "1", 1);


    //------------------------------------------------------------------------------------
    //Save in file
    //------------------------------------------------------------------------------------
    writeOrbit(tNCE, yNCE, OFTS_ORDER, floor(fabs(st0[0])), N+1, TYPE_ORBIT);

    //------------------------------------------------------------------------------------
    //User check
    //------------------------------------------------------------------------------------
    char ch;
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c",&ch);
    gnuplot_close(true, h1);

    //------------------------------------------------------------------------------------
    //Free
    //------------------------------------------------------------------------------------
    free_orbit(&orbit);
    //Free
    free_dmatrix(yNCE, 0, 5, 0, N);
    free_dmatrix(ySYS, 0, 5, 0, N);
    free_dvector(tNCE, 0, N);
    return status;
}

//========================================================================================
//
//          SingleOrbit structure - deprecated
//
//========================================================================================
/**
 *   \brief Initialize one SingleOrbit structure
 **/
void init_orbit(SingleOrbit &orbit, vector<Oftsc>*  W, vector<Oftsc>*  Wh,
                matrix<Ofsc>*  PC, matrix<Ofsc>*  CQ, vector<Ofsc>*  V,
                Ofsc* orbit_ofs, int ofts_order, int ofs_order, int isGS, double t0,
                double tf, double tproj, OdeStruct *driver, QBCP_L *qbcp_l)
{
    //-----------
    //Parameterization (common to all orbits)
    //-----------
    orbit.W   =  W;            //z(t) = W(s(t), t)
    orbit.Wh  =  Wh;           //zh(t) = Wh(s(t), t)

    //-----------
    //COC (common to all orbits)
    //-----------
    orbit.PC  = PC;            //COC matrix
    orbit.CQ  = CQ;            //inv COC matrix
    orbit.V   = V;             //COC vector

    //-----------
    //Orders
    //-----------
    orbit.reduced_nv =  (*Wh)[0].getNV();   //reduced number of variables, taken from the parameterization directly
    orbit.ofs        =  orbit_ofs;          //Auxiliary Ofs object
    orbit.order      =  ofts_order;         //ofts_order of the expansions
    orbit.ofs_order  =  ofs_order;          //ofts_order of the Fourier coefficients
    orbit.isGS       =  isGS;               //Was the pm obtained through graph style?

    //-----------
    //Characteristics
    //-----------
    orbit.z0  = dvector(0, 5);                        //Initial position in NC coordinates dim = 6
    orbit.si  = dvector(0, orbit.reduced_nv-1);       //Initial RCM configuration dim = 4
    orbit.s0d = dvector(0, 2*orbit.reduced_nv-1);     //Initial position in CCM8 coordinates (real+imag part) dim = 8
    orbit.xf  = dvector(0, 5);                        //Final position NC dim = 6
    orbit.s0  = dcvector(0,orbit.reduced_nv-1);       //Initial position in CCM4 coordinates (real+imag part) dim = 4
    orbit.t0  = t0;                                   //Initial time
    orbit.tf  = tf;                                   //Final time after computation
    orbit.tproj  = tproj;                             //default time between each projection
    orbit.tprojmin  = 1e-4*qbcp_l->us->T;              //minimum time between each projection
    orbit.ePmax = (qbcp_l->fwrk == F_SEM)? 5e-3:5e-5; //maximum projection distance allowed

    //-----------
    //OLD VALUE
    //-----------
    //  orbit.ePmax = (qbcp_l->fwrk == F_SEM)? 5e-4:1e-5;


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
    orbit.n = qbcp_l->us->n;
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

//========================================================================================
//
//          Initial conditions - deprecated
//
//========================================================================================
/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], double t0)
{
    //------------------------------------------------------------------------------------
    // 1. Update si
    //------------------------------------------------------------------------------------
    for(int p = 0; p < orbit.reduced_nv; p++) orbit.si[p] = si[p];

    //------------------------------------------------------------------------------------
    // 2. Update s0
    //------------------------------------------------------------------------------------
    RCMtoCCM(si, orbit.s0, orbit.reduced_nv);

    //------------------------------------------------------------------------------------
    // 2. Update s0d
    //------------------------------------------------------------------------------------
    RCMtoCCM8(si, orbit.s0d, orbit.reduced_nv);

    //------------------------------------------------------------------------------------
    // 4. Update z0
    //------------------------------------------------------------------------------------
    //z0 = W(si, t0)
    RCMtoNCbyTFC(si,
                 t0,
                 orbit.n,
                 orbit.order,
                 orbit.ofs_order,
                 orbit.reduced_nv,
                *orbit.Wh,
                *orbit.PC,
                *orbit.V,
                 orbit.z0,
                 orbit.isGS);
}

/**
 *   \brief Update the initial conditions (si, s0, z0 and s0d) of the orbit given an array of initial TFC conditions si
 *          and an array of initial NC conditions z0.
 **/
void orbit_update_ic(SingleOrbit &orbit, const double si[], const double z0[])
{
    //------------------------------------------------------------------------------------
    // 1. Update si
    //------------------------------------------------------------------------------------
    for(int p = 0; p < orbit.reduced_nv; p++) orbit.si[p] = si[p];

    //------------------------------------------------------------------------------------
    // 2. Update s0
    //------------------------------------------------------------------------------------
    RCMtoCCM(si, orbit.s0, orbit.reduced_nv);

    //------------------------------------------------------------------------------------
    // 2. Update s0d
    //------------------------------------------------------------------------------------
    RCMtoCCM8(si, orbit.s0d, orbit.reduced_nv);

    //------------------------------------------------------------------------------------
    // 4. Update z0
    //------------------------------------------------------------------------------------
    //z0 = W(si, 0.0)
    for(int p = 0; p < NV; p++) orbit.z0[p] = z0[p];
}
