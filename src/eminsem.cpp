/**
 * \file  eminsem.cpp
 * \brief Contains all the routines to perform changes of coordinates between the EM and SEM frameworks. Including
 * \author BLB.
 * \date   2016
 * \version 1.0
 */

#include "eminsem.h"

//========================================================================================
//
// Super COCs
//
//========================================================================================
/**
 *  \brief COC: from inputType to outputType. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc(double t, const double y0[], double yout[], int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // 1. Do some checks on the inputs
    //------------------------------------------------------------------------------------
    //Type of inputs
    if(inputType > ECISEM)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > ECISEM)
        perror("Unknown outputType");

    //------------------------------------------------------------------------------------
    // 2. Define the default framework wrt the inputType
    //------------------------------------------------------------------------------------
    int fwrk = 0;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        fwrk = F_EM;
        break;
    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        fwrk = F_SEM;
        break;
    }


    //------------------------------------------------------------------------------------
    // 2. Check that the focus in SEML is
    // in accordance with the inputType.
    //------------------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);


    //------------------------------------------------------------------------------------
    // 3. Updating output
    //------------------------------------------------------------------------------------
    qbcp_coc_fwrk(t, y0, yout, inputType, outputType);


    //------------------------------------------------------------------------------------
    // 4. Reset the focus in SEML, if necessary
    //------------------------------------------------------------------------------------
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);
}

/**
 *  \brief COC: from inputType to outputType. With a time update in tout. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc(double t, const double y0[], double yout[], double *tout, int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // 1. Do some checks on the inputs
    //------------------------------------------------------------------------------------
    //Type of inputs
    if(inputType > ECISEM)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > ECISEM)
        perror("Unknown outputType");

    //------------------------------------------------------------------------------------
    // 2. Define the default framework wrt the inputType
    //    Define the factor to apply to the time vector
    //------------------------------------------------------------------------------------
    int fwrk = 0;
    double tfactor = 1.0;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        //EM framework
        fwrk = F_EM;
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //EM to EM units
            tfactor = 1.0;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //EM to SEM units
            tfactor = SEML.us_em.ns;
            break;
        }
        break;

    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        //SEM framework
        fwrk = F_SEM;
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //SEM to EM units
            tfactor = 1.0/SEML.us_em.ns;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //SEM to SEM units
            tfactor = 1.0;
            break;
        }
        break;
    }


    //------------------------------------------------------------------------------------
    // 2. Check that the focus in SEML is
    // in accordance with the inputType.
    //------------------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);


    //------------------------------------------------------------------------------------
    // 3. Updating output
    //------------------------------------------------------------------------------------
    qbcp_coc_fwrk(t, y0, yout, inputType, outputType);
    *tout = tfactor*t;


    //------------------------------------------------------------------------------------
    // 4. Reset the focus in SEML, if necessary
    //------------------------------------------------------------------------------------
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);
}

/**
 *  \brief COC: from inputType to outputType, in vector form. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc_vec(double **y0, double *t0, double **yout, double *tout, int N, int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // 1. Do some checks on the inputs
    //------------------------------------------------------------------------------------
    //Type of inputs
    if(inputType > ECISEM)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > ECISEM)
        perror("Unknown outputType");

    //------------------------------------------------------------------------------------
    // 2. Define the default framework wrt the inputType
    //    Define the factor to apply to the time vector
    //------------------------------------------------------------------------------------
    int fwrk = 0;
    double tfactor = 1.0;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        //EM framework
        fwrk = F_EM;
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //EM to EM units
            tfactor = 1.0;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //EM to SEM units
            tfactor = SEML.us_em.ns;
            break;
        }
        break;

    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        //SEM framework
        fwrk = F_SEM;
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //SEM to EM units
            tfactor = 1.0/SEML.us_em.ns;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //SEM to SEM units
            tfactor = 1.0;
            break;
        }
        break;
    }


    //------------------------------------------------------------------------------------
    // 2. Check that the focus in SEML is
    // in accordance with the inputType.
    //------------------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);


    //------------------------------------------------------------------------------------
    // 3. Updating output
    //------------------------------------------------------------------------------------
    double yt[6], yt2[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yt
        for(int k = 0; k < 6; k++) yt[k] = y0[k][p];
        //Perform COC
        qbcp_coc_fwrk(t0[p], yt, yt2, inputType, outputType);
        //Copy result in yout
        for(int k = 0; k < 6; k++) yout[k][p] = yt2[k];
        //time
        tout[p] = t0[p]*tfactor;
    }


    //------------------------------------------------------------------------------------
    // 4. Reset the focus in SEML, if necessary
    //------------------------------------------------------------------------------------
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);
}

/**
 *  \brief COC: from inputType to outputType, in vector form. Version with no time vector output. Some specific checks are made.
 *         In particular, this routine makes sure that SEML is focused on the right system (either SEM or EM, depending on the inputType).
 *         The routine is able to make the COC between 8 different types of outputs: NCEM, VNCEM, PEM, VEM, and their equivalents in SEM coordinates.
 *         All 64 possibilities are available.
 **/
void qbcp_coc_vec(double **y0, double *t0, double **yout, int N, int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // 1. Do some checks on the inputs
    //------------------------------------------------------------------------------------
    //Type of inputs
    if(inputType > ECISEM)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > ECISEM)
        perror("Unknown outputType");

    //------------------------------------------------------------------------------------
    // 2. Define the default framework wrt the inputType
    //    Define the factor to apply to the time vector
    //------------------------------------------------------------------------------------
    int fwrk = 0;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        //EM framework
        fwrk = F_EM;
        break;

    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        //SEM framework
        fwrk = F_SEM;
        break;
    }


    //------------------------------------------------------------------------------------
    // 2. Check that the focus in SEML is
    // in accordance with the inputType.
    //------------------------------------------------------------------------------------
    int fwrk0 = SEML.fwrk;
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk);


    //------------------------------------------------------------------------------------
    // 3. Updating output
    //------------------------------------------------------------------------------------
    double yt[6], yt2[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yt
        for(int k = 0; k < 6; k++) yt[k] = y0[k][p];
        //Perform COC
        qbcp_coc_fwrk(t0[p], yt, yt2, inputType, outputType);
        //Copy result in yout
        for(int k = 0; k < 6; k++) yout[k][p] = yt2[k];
    }


    //------------------------------------------------------------------------------------
    // 4. Reset the focus in SEML, if necessary
    //------------------------------------------------------------------------------------
    if(fwrk0 != fwrk) changeDCS(SEML, fwrk0);
}

/**
 *  \brief COC: from inputType to outputType.
 *         Used ONLY in qbcp_coc and qbcp_coc_vec, since specific checks
 *         are made in these routines prior to any computations.
 *         €€TODO: simplify the COCs. Done for now: VNCEM.
 **/
void qbcp_coc_fwrk(double t, const double y0[], double yout[], int inputType, int outputType)
{
    double ytemp[6], ytemp2[6], ytemp3[6], ttemp;
    switch(inputType)
    {
        //--------------------------------------------------------------------------------
        // For inputType = VNCEM, NCEM, PEM, VEM, INEM, the focus of SEML is
        // on the EM SYSTEM
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        // VNCEM
        //--------------------------------------------------------------------------------
    case VNCEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //VNCEM -> INEM
            //----------------------------------------------------------------------------
            //VNCEM -> VEM
            NCvtoSYSv(y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //VNCEM -> VNCEM nothing to do
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //VNCEM -> NCEM nothing to do
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //VNCEM -> PEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //VNCEM -> VEM
            //----------------------------------------------------------------------------
            NCEMvtoEMv(y0, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> NCSEM
            //----------------------------------------------------------------------------
            NCEMvtoNCSEMm(t, y0, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> PSEM
            //----------------------------------------------------------------------------
            NCEMvtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> VSEM
            //----------------------------------------------------------------------------
            NCEMvtoSEMv(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> VNCSEM
            //----------------------------------------------------------------------------
            NCEMvtoNCSEMv(t, y0, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> INSEM
            //----------------------------------------------------------------------------
            //VNCEM -> VEM
            NCvtoSYSv(y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, yout, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //VNCEM -> ECISEM
            //----------------------------------------------------------------------------
            //VNCEM -> VEM
            NCvtoSYSv(y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, ytemp2, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, ytemp2, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(ttemp, ytemp2, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // NCEM
        //--------------------------------------------------------------------------------
    case NCEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //NCEM -> INEM
            //----------------------------------------------------------------------------
            //NCEM -> VNCEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VNCEM -> VEM
            NCvtoSYSv(ytemp, ytemp2, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp2, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //NCEM -> VNCEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, yout, &SEML);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //NCEM -> NCEM nothing to do
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //NCEM -> PEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //NCEM -> VEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //NCEM -> NCSEM
            //----------------------------------------------------------------------------
            NCEMmtoNCSEMm(t, y0, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //NCEM -> PSEM
            //----------------------------------------------------------------------------
            NCEMmtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //NCEM -> VSEM
            //----------------------------------------------------------------------------
            NCEMmtoSEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            MOMtoVEL(t*SEML.us_em.ns, ytemp, yout, &SEML_SEM);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //NCEM -> VNCSEM
            //----------------------------------------------------------------------------
            NCEMmtoNCSEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            MOMtoVEL(t*SEML.us_em.ns, ytemp, yout, &SEML_SEM);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //NCEM -> INSEM
            //----------------------------------------------------------------------------
            //NCEM -> VNCEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VNCEM -> VEM
            NCvtoSYSv(ytemp, ytemp2, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp2, yout, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //NCEM -> ECISEM
            //----------------------------------------------------------------------------
            //NCEM -> VNCEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VNCEM -> VEM
            NCvtoSYSv(ytemp, ytemp2, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp2, ytemp3, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, ytemp3, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(ttemp, ytemp3, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // PEM
        //--------------------------------------------------------------------------------
    case PEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //PEM -> INEM
            //----------------------------------------------------------------------------
            //PEM -> VEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, yout, &SEML);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //PEM -> NCEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //PEM -> PEM
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case VEM:
            //----------------------------------------------------------------------------
            // PEM -> VEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            // PEM -> VNCEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //PEM -> NCSEM
            //----------------------------------------------------------------------------
            EMmtoSEMm(t, y0, ytemp, &SEML);
            //Careful here: the time must be in SEM units.
            SEMmtoNCSEMm(t*SEML.us_em.ns, ytemp, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //PEM -> PSEM
            //----------------------------------------------------------------------------
            EMmtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //PEM -> VSEM
            //----------------------------------------------------------------------------
            EMmtoSEMv(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //PEM -> VNCSEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCEMmtoNCSEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in SEM units.
            // 2. The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            MOMtoVEL(t*SEML.us_em.ns, ytemp2, yout, &SEML_SEM);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //PEM -> INSEM
            //----------------------------------------------------------------------------
            //PEM -> VEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, yout, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //PEM -> ECISEM
            //----------------------------------------------------------------------------
            //PEM -> VEM
            MOMtoVEL(t, y0, ytemp, &SEML);
            //VEM -> INEM
            EMvtoIN(t, ytemp, ytemp2, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, ytemp2, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(ttemp, ytemp2, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // VEM
        //--------------------------------------------------------------------------------
    case VEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //VEM -> INEM
            //----------------------------------------------------------------------------
            EMvtoIN(t, y0, yout, &SEML);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            // VEM -> NCEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            // VEM -> PEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            // VEM -> VEM
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //VEM -> VNCEM
            //----------------------------------------------------------------------------
            EMvtoNCEMv(y0, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            // VEM -> NCSEM
            //----------------------------------------------------------------------------
            EMvtoSEMm(t, y0, ytemp, &SEML);
            //Careful here: the time must be in SEM units.
            SEMmtoNCSEMm(t*SEML.us_em.ns, ytemp, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            // VEM -> PSEM
            //----------------------------------------------------------------------------
            EMvtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            // VEM -> VSEM
            //----------------------------------------------------------------------------
            EMvtoSEMv(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            // VEM -> VNCSEM
            //----------------------------------------------------------------------------
            EMvtoSEMv(t, y0, ytemp, &SEML);
            //Careful here:
            // The SEML_SEM structure must be used, in order to
            // get the SEM coefficients
            SYSvtoNCv(ytemp, yout, &SEML_SEM);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //VEM -> INSEM
            //----------------------------------------------------------------------------
            EMvtoIN(t, y0, yout, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //VEM -> ECISEM
            //----------------------------------------------------------------------------
            EMvtoIN(t, y0, ytemp, &SEML);
            //INEM -> INSEM
            ttemp = t;
            usem2ussem(&ttemp, ytemp, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(ttemp, ytemp, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // INEM
        //--------------------------------------------------------------------------------
    case INEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //INEM -> INEM
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            // INEM -> NCEM
            //----------------------------------------------------------------------------
            INtoEM(t, y0, ytemp, &SEML);
            VELtoMOM(t, ytemp, ytemp2, &SEML);
            SYStoNC(t, ytemp2, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            // INEM -> PEM
            //----------------------------------------------------------------------------
            INtoEM(t, y0, ytemp, &SEML);
            VELtoMOM(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            // INEM -> VEM
            //----------------------------------------------------------------------------
            INtoEM(t, y0, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //INEM -> VNCEM
            //----------------------------------------------------------------------------
            INtoEM(t, y0, ytemp, &SEML);
            EMvtoNCEMv(ytemp, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            // INEM -> NCSEM
            //----------------------------------------------------------------------------
            //INEM -> VSEM
            INEMvtoSEMv(t, y0, ytemp, &SEML);
            //VSEM to NCSEM. Careful here: the time must be in SEM units.
            SEMvtoNCSEMm(t*SEML.us_em.ns, ytemp, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            // INEM -> PSEM
            //----------------------------------------------------------------------------
            INEMvtoSEMm(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            // INEM -> VSEM
            //----------------------------------------------------------------------------
            INEMvtoSEMv(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            // INEM -> VNCSEM
            //----------------------------------------------------------------------------
            //INEM -> VSEM
            INEMvtoSEMv(t, y0, ytemp, &SEML);
            //VSEM -> VNCSEM.
            SEMvtoNCSEMv(ytemp, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //INEM -> INSEM
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            ttemp = t;
            usem2ussem(&ttemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //INEM -> ECISEM
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) ytemp[i] = y0[i];
            ttemp = t;
            usem2ussem(&ttemp, ytemp, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(ttemp, ytemp, yout, &SEML);
            break;
        }
        break;


        //--------------------------------------------------------------------------------
        // For inputType = VNCSEM, NCSEM, PSEM, VSEM, INSEM, the focus of SEML is
        // on the SEM SYSTEM
        //--------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------
        // VNCSEM
        //--------------------------------------------------------------------------------
    case VNCSEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> INEM
            //----------------------------------------------------------------------------
            //VNCSEM to VSEM
            NCSEMvtoSEMv(y0, ytemp, &SEML);
            //VSEM to VEM
            SEMvtoEMv(t, ytemp, ytemp2, &SEML);
            //VEM to INEM. Careful here:
            // 1. The time must be in EM units.
            EMvtoIN(t/SEML.us_em.ns, ytemp2, yout, &SEML);
            break;
        case VNCEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> VNCEM
            //----------------------------------------------------------------------------
            //VNCSEM to VSEM
            NCSEMvtoSEMv(y0, ytemp, &SEML);
            //VSEM to VEM
            SEMvtoEMv(t, ytemp, ytemp2, &SEML);
            //VEM to VNCEM.
            EMvtoNCEMv(ytemp2, yout, &SEML);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> NCEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> PEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            NCSEMmtoEMm(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> VEM
            //----------------------------------------------------------------------------
            //VNCSEM to VSEM
            NCSEMvtoSEMv(y0, ytemp, &SEML);
            //VSEM to VEM
            SEMvtoEMv(t, ytemp, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> NCSEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> PSEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            NCtoSYS(t, ytemp, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> VSEM
            //----------------------------------------------------------------------------
            NCSEMvtoSEMv(y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //VNCSEM -> VNCSEM nothing to do
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //VNCEM -> INSEM
            //----------------------------------------------------------------------------
            NCSEMvtoSEMv(y0, ytemp, &SEML);
            SEMtoIN(t, ytemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //VNCEM -> ECISEM
            //----------------------------------------------------------------------------
            NCSEMvtoSEMv(y0, ytemp, &SEML);
            SEMtoIN(t, ytemp, ytemp2, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(t, ytemp2, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // NCSEM
        //--------------------------------------------------------------------------------
    case NCSEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //NCSEM -> INEM
            //----------------------------------------------------------------------------
            //NCSEM -> EM
            NCSEMmtoEMm(t, y0, ytemp, &SEML);
            //EM -> VEM
            MOMtoVEL(t/SEML.us_em.ns, ytemp, ytemp2, &SEML_EM);
            //VEM to INEM
            EMvtoIN(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //NCSEM -> VNCEM
            //----------------------------------------------------------------------------
            NCSEMmtoNCEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            MOMtoVEL(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //NCSEM -> NCEM
            //----------------------------------------------------------------------------
            NCSEMmtoNCEMm(t, y0, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //NCSEM -> PEM
            //----------------------------------------------------------------------------
            NCSEMmtoEMm(t, y0, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //NCSEM -> VEM
            //----------------------------------------------------------------------------
            NCSEMmtoEMm(t, y0, ytemp, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            MOMtoVEL(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //NCSEM -> NCSEM nothing to do here
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++) yout[i] = y0[i];
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //NCSEM -> PSEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //NCSEM -> VSEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //NCSEM -> VNCSEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //NCSEM -> INSEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, ytemp2, &SEML);
            SEMtoIN(t, ytemp2, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //NCSEM -> ECISEM
            //----------------------------------------------------------------------------
            NCtoSYS(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, ytemp2, &SEML);
            SEMtoIN(t, ytemp2, ytemp3, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(t, ytemp3, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // PSEM
        //--------------------------------------------------------------------------------
    case PSEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //SEM -> INEM
            //----------------------------------------------------------------------------
            //SEM -> VEM
            SEMmtoEMv(t, y0, ytemp, &SEML);
            //VEM to INEM
            EMvtoIN(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //PSEM -> VNCEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, ytemp2, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            MOMtoVEL(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //PSEM -> NCEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            NCSEMmtoNCEMm(t, ytemp, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //PSEM -> PEM
            //----------------------------------------------------------------------------
            SEMmtoEMm(t, y0, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //PSEM -> VEM
            //----------------------------------------------------------------------------
            SEMmtoEMv(t, y0, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //PSEM -> NCSEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //PSEM -> PSEM nothing to do here
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //PSEM -> VSEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //PSEM -> VNCSEM
            //----------------------------------------------------------------------------
            SYStoNC(t, y0, ytemp, &SEML);
            MOMtoVEL(t, ytemp, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //PSEM -> INSEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, ytemp, &SEML);
            SEMtoIN(t, ytemp, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //PSEM -> ECISEM
            //----------------------------------------------------------------------------
            MOMtoVEL(t, y0, ytemp, &SEML);
            SEMtoIN(t, ytemp, ytemp2, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(t, ytemp2, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // VSEM
        //--------------------------------------------------------------------------------
    case VSEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //VSEM -> INEM
            //----------------------------------------------------------------------------
            //VSEM -> VEM
            SEMvtoEMv(t, y0, ytemp, &SEML);
            //VEM to INEM
            EMvtoIN(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //VSEM -> VNCEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            NCSEMmtoNCEMm(t, ytemp2, ytemp3, &SEML);
            //Careful here:
            // 1. The time must be in EM units.
            // 2. The SEML_EM structure must be used, in order to
            // get the EM coefficients
            MOMtoVEL(t/SEML.us_em.ns, ytemp3, yout, &SEML_EM);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //VSEM -> NCEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, ytemp2, &SEML);
            NCSEMmtoNCEMm(t, ytemp2, yout, &SEML);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //VSEM -> PEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            SEMmtoEMm(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //VSEM -> VEM
            //----------------------------------------------------------------------------
            SEMvtoEMv(t, y0, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //VSEM -> NCSEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, ytemp, &SEML);
            SYStoNC(t, ytemp, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //VSEM -> PSEM
            //----------------------------------------------------------------------------
            VELtoMOM(t, y0, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //VSEM -> VSEM
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //VSEM -> VNCSEM
            //----------------------------------------------------------------------------
            SEMvtoNCSEMv(y0, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //VSEM -> INSEM
            //----------------------------------------------------------------------------
            SEMtoIN(t, y0, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //VSEM -> ECISEM
            //----------------------------------------------------------------------------
            SEMtoIN(t, y0, ytemp, &SEML);
            //INSEM -> ECISEM
            INSEMtoECISEM(t, ytemp, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // INSEM
        //--------------------------------------------------------------------------------
    case INSEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //INSEM -> INEM
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            ttemp = t;
            ussem2usem(&ttemp, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //INSEM -> VNCEM
            //----------------------------------------------------------------------------
            //INSEM -> VEM
            INSEMvtoEMv(t, y0, ytemp, &SEML);
            //VEM -> VNCEM
            EMvtoNCEMv(ytemp, yout, &SEML_EM);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //INSEM -> NCEM
            //----------------------------------------------------------------------------
            //INSEM -> VEM
            INSEMvtoEMv(t, y0, ytemp, &SEML);
            //VEM -> NCEM. The time must be set in EM units
            EMvtoNCEMm(t/SEML.us_em.ns, ytemp, yout, &SEML_EM);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //INSEM -> PEM
            //----------------------------------------------------------------------------
            INSEMvtoEMm(t, y0, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //INSEM -> VEM
            //----------------------------------------------------------------------------
            INSEMvtoEMv(t, y0, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //INSEM -> NCSEM
            //----------------------------------------------------------------------------
            INtoSEM(t, y0, ytemp, &SEML);
            VELtoMOM(t, ytemp, ytemp2, &SEML);
            SYStoNC(t, ytemp2, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //INSEM -> PSEM
            //----------------------------------------------------------------------------
            INtoSEM(t, y0, ytemp, &SEML);
            VELtoMOM(t, ytemp, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //INSEM -> VSEM
            //----------------------------------------------------------------------------
            INtoSEM(t, y0, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //INSEM -> VNCSEM
            //----------------------------------------------------------------------------
            INtoSEM(t, y0, ytemp, &SEML);
            SEMvtoNCSEMv(ytemp, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //INSEM -> INSEM
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //INSEM -> ECISEM
            //----------------------------------------------------------------------------
            INSEMtoECISEM(t, y0, yout, &SEML);
            break;
        }
        break;

        //--------------------------------------------------------------------------------
        // ECISEM
        //--------------------------------------------------------------------------------
    case ECISEM:
        switch(outputType)
        {
        case INEM:
            //----------------------------------------------------------------------------
            //ECISEM -> INEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, yout, &SEML);
            ttemp = t;
            ussem2usem(&ttemp, yout, &SEML);
            break;

        case VNCEM:
            //----------------------------------------------------------------------------
            //ECISEM -> VNCEM
            //----------------------------------------------------------------------------
            //ECISEM -> INSEM
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            //INSEM -> VEM
            INSEMvtoEMv(t, ytemp, ytemp2, &SEML);
            //VEM -> VNCEM
            EMvtoNCEMv(ytemp2, yout, &SEML_EM);
            break;

        case NCEM:
            //----------------------------------------------------------------------------
            //ECISEM -> NCEM
            //----------------------------------------------------------------------------
            //ECISEM -> INSEM
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            //INSEM -> VEM
            INSEMvtoEMv(t, ytemp, ytemp2, &SEML);
            //VEM -> NCEM. The time must be set in EM units
            EMvtoNCEMm(t/SEML.us_em.ns, ytemp2, yout, &SEML_EM);
            break;

        case PEM:
            //----------------------------------------------------------------------------
            //ECISEM -> PEM
            //----------------------------------------------------------------------------
            //ECISEM -> INSEM
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            //INSEM -> PEM
            INSEMvtoEMm(t, ytemp, yout, &SEML);
            break;

        case VEM:
            //----------------------------------------------------------------------------
            //ECISEM -> VEM
            //----------------------------------------------------------------------------
            //ECISEM -> INSEM
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            //INSEM -> VEM
            INSEMvtoEMv(t, ytemp, yout, &SEML);
            break;

        case NCSEM:
            //----------------------------------------------------------------------------
            //ECISEM -> NCSEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            INtoSEM(t, ytemp, ytemp2, &SEML);
            VELtoMOM(t, ytemp2, ytemp3, &SEML);
            SYStoNC(t, ytemp3, yout, &SEML);
            break;

        case PSEM:
            //----------------------------------------------------------------------------
            //ECISEM -> PSEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            INtoSEM(t, ytemp, ytemp2, &SEML);
            VELtoMOM(t, ytemp2, yout, &SEML);
            break;

        case VSEM:
            //----------------------------------------------------------------------------
            //ECISEM -> VSEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            INtoSEM(t, ytemp, yout, &SEML);
            break;

        case VNCSEM:
            //----------------------------------------------------------------------------
            //ECISEM -> VNCSEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, ytemp, &SEML);
            INtoSEM(t, ytemp, ytemp2, &SEML);
            SEMvtoNCSEMv(ytemp2, yout, &SEML);
            break;

        case INSEM:
            //----------------------------------------------------------------------------
            //ECISEM -> INSEM
            //----------------------------------------------------------------------------
            ECISEMtoINSEM(t, y0, yout, &SEML);
            break;

        case ECISEM:
            //----------------------------------------------------------------------------
            //ECISEM -> ECISEM
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++) yout[i] = y0[i];
            break;
        }
        break;
    }
}

/**
 *  \brief COC: from inputType to outputType, in vector form, for time vectors only.
 **/
void qbcp_coc_time(double *t0, double *tout, int N, int inputType, int outputType)
{
    //------------------------------------------------------------------------------------
    // 1. Do some checks on the inputs
    //------------------------------------------------------------------------------------
    //Type of inputs
    if(inputType > ECISEM)
        perror("Unknown inputType");

    //Type of outputs
    if(outputType > ECISEM)
        perror("Unknown outputType");

    //------------------------------------------------------------------------------------
    // 2. Define the factor to apply to the time vector
    //------------------------------------------------------------------------------------
    double tfactor = 1.0;
    switch(inputType)
    {
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //EM to EM units
            tfactor = 1.0;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //EM to SEM units
            tfactor = SEML.us_em.ns;
            break;
        }
        break;

    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        //Time factor
        switch(outputType)
        {
        case VNCEM:
        case NCEM:
        case PEM:
        case VEM:
        case INEM:
            //SEM to EM units
            tfactor = 1.0/SEML.us_em.ns;
            break;
        case VNCSEM:
        case NCSEM:
        case PSEM:
        case VSEM:
        case INSEM:
        case ECISEM:
            //SEM to SEM units
            tfactor = 1.0;
            break;
        }
        break;
    }


    //------------------------------------------------------------------------------------
    // 3. Updating output
    //------------------------------------------------------------------------------------
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //time
        tout[p] = t0[p]*tfactor;
    }
}

//========================================================================================
//
// Change of coordinates: NC <-> SYS
//
//========================================================================================
//----------------------------------------------------------------------------------------
// COC: NC <--> SYS
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: NC coordinates to SYSTEM coordinates. Use in priority instead of NCEMmtoEMm or NCSEMmtoSEMm.
 **/
void NCtoSYS(double t, const double yNC[], double yEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from SYSTEM coordinates to NC coordinates. Use in priority instead of EMtoN or SEMmtoNCSEMm.
 **/
void SYStoNC(double t, const double yEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us.n;
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

//----------------------------------------------------------------------------------------
// COC: NC <--> SEM
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from SEM coordinates to NC coordinates
 **/
void SEMmtoNCSEMm(double t, const double ySEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -ySEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -ySEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +ySEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -ySEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -ySEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +ySEM[5]/gamma;
}

/**
 *  \brief COC: from NC coordinates to SEM coordinates
 **/
void NCSEMmtoSEMm(double t, const double yNC[], double ySEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_sem.n;
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    ySEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    ySEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    ySEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    ySEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    ySEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    ySEM[5] = +gamma*yNC[5];

}

//----------------------------------------------------------------------------------------
// COC: NC <--> EM
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: from NC coordinates to EM coordinates
 **/
void NCEMmtoEMm(double t, const double yNC[], double yEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma(px+a2/a1*c1)
    yEM[3] = -gamma*(yNC[3] + alpha[1]/alpha[0]*c1);
    //PY = -gamma(py-a3/a1*c1)
    yEM[4] = -gamma*(yNC[4] - alpha[2]/alpha[0]*c1);
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];

}

/**
 *  \brief COC: from EM coordinates to NC coordinates
 **/
void EMmtoNCEMm(double t, const double yEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double n     =  qbp->us_em.n;
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[8];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 8);

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma - a2/a1*c1
    yNC[3] = -yEM[3]/gamma - alpha[1]/alpha[0]*c1;
    //py = -PY/gamma + a3/a1*c1
    yNC[4] = -yEM[4]/gamma + alpha[2]/alpha[0]*c1;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}


//----------------------------------------------------------------------------------------
// COC: VNC <--> NC
//----------------------------------------------------------------------------------------
/**
 *  \brief COC: NC coordinates to SYSTEM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCvtoSYSv(const double yNC[], double yEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma*px
    yEM[3] = -gamma*yNC[3];
    //PY = -gamma*py
    yEM[4] = -gamma*yNC[4];
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from SYSTEM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void SYSvtoNCv(const double yEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs.gamma;
    double c1    =  qbp->cs.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma
    yNC[3] = -yEM[3]/gamma;
    //py = -PY/gamma
    yNC[4] = -yEM[4]/gamma;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

/**
 *  \brief COC: NC coordinates to EM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCEMvtoEMv(const double yNC[], double yEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma*px
    yEM[3] = -gamma*yNC[3];
    //PY = -gamma*py
    yEM[4] = -gamma*yNC[4];
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from EM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void EMvtoNCEMv(const double yEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs_em.gamma;
    double c1    =  qbp->cs_em.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma
    yNC[3] = -yEM[3]/gamma;
    //py = -PY/gamma
    yNC[4] = -yEM[4]/gamma;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

/**
 *  \brief COC: NC coordinates to SEM coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void NCSEMvtoSEMv(const double yNC[], double yEM[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //X = -gamma*(x-c1)
    yEM[0] = -gamma*(yNC[0] - c1);
    //Y = -gamma*y
    yEM[1] = -gamma*yNC[1];
    //Z = +gamma*z
    yEM[2] = +gamma*yNC[2];

    //PX = -gamma*px
    yEM[3] = -gamma*yNC[3];
    //PY = -gamma*py
    yEM[4] = -gamma*yNC[4];
    //PZ = +gamma*pz
    yEM[5] = +gamma*yNC[5];
}

/**
 *  \brief COC: from SEM coordinates to NC coordinates, with a state of the form (x, v), instead of (x, p).
 **/
void SEMvtoNCSEMv(const double yEM[], double yNC[], QBCP_L *qbp)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    double gamma =  qbp->cs_sem.gamma;
    double c1    =  qbp->cs_sem.c1;

    //------------------------------------------------------------------------------------
    //CoC
    //------------------------------------------------------------------------------------
    //x = -X/gamma + c1
    yNC[0] = -yEM[0]/gamma  + c1;
    //y = -Y/gamma
    yNC[1] = -yEM[1]/gamma;
    //z = +Z/gamma
    yNC[2] = +yEM[2]/gamma;
    //px = -PX/gamma
    yNC[3] = -yEM[3]/gamma;
    //py = -PY/gamma
    yNC[4] = -yEM[4]/gamma;
    //pz = +PZ/gamma
    yNC[5] = +yEM[5]/gamma;
}

//========================================================================================
//
// COC: Velocities <--> Momenta
//
//========================================================================================

//----------------------------------------------------------------------------------------
// SEM
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the SEM velocities into SEM momenta
 **/
void SEMvtoSEMm(double t, const double ySEv[], double ySEm[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp = (QBCP_L*) params_void;
    double n    =  qbp->us_sem.n;

    //------------------------------------------------------------------------------------
    //evaluate the deltas
    //------------------------------------------------------------------------------------
    double delta[3];
    evaluateCoef(delta, t, n, qbp->nf, qbp->cs_sem.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySEm[i] = ySEv[i];

    //------------------------------------------------------------------------------------
    //Velocity to Momenta
    //------------------------------------------------------------------------------------
    //px
    ySEm[3] = 1.0/delta[0] * (ySEv[3] - delta[1]*ySEv[0] - delta[2]*ySEv[1]);
    //py
    ySEm[4] = 1.0/delta[0] * (ySEv[4] - delta[1]*ySEv[1] + delta[2]*ySEv[0]);
    //pz
    ySEm[5] = 1.0/delta[0] * (ySEv[5] - delta[1]*ySEv[2] );
}

/**
 *  \brief Change the SEM momenta into SEM velocities
 **/
void SEMmtoSEMv(double t, const double ySEm[], double ySEv[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp = (QBCP_L*) params_void;
    double n    =  qbp->us_sem.n;

    //------------------------------------------------------------------------------------
    //evaluate the deltas
    //------------------------------------------------------------------------------------
    double delta[3];
    evaluateCoef(delta, t, n, qbp->nf, qbp->cs_sem.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySEv[i] = ySEm[i];

    //------------------------------------------------------------------------------------
    //Momenta to velocities
    //------------------------------------------------------------------------------------
    //vx
    ySEv[3] = delta[0]*ySEm[3] + delta[1]*ySEm[0] + delta[2]*ySEm[1];
    //vy
    ySEv[4] = delta[0]*ySEm[4] + delta[1]*ySEm[1] - delta[2]*ySEm[0];
    //vz
    ySEv[5] = delta[0]*ySEm[5] + delta[1]*ySEm[2];
}

//----------------------------------------------------------------------------------------
// EM
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the EM velocities into EM momenta
 **/
void EMvtoEMm(double t, const double yEMv[], double yEMm[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us_em.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas
    //------------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) yEMm[i] = yEMv[i];


    //------------------------------------------------------------------------------------
    //Velocity to Momenta
    //------------------------------------------------------------------------------------
    //px
    yEMm[3] = 1.0/alpha[0] * (yEMv[3] - alpha[1]*yEMv[0] - alpha[2]*yEMv[1]);
    //py
    yEMm[4] = 1.0/alpha[0] * (yEMv[4] - alpha[1]*yEMv[1] + alpha[2]*yEMv[0]);
    //pz
    yEMm[5] = 1.0/alpha[0] * (yEMv[5] - alpha[1]*yEMv[2] );
}

/**
 *  \brief Change the EM momenta into EM velocities
 **/
void EMmtoEMv(double t, const double yEMm[], double yEMv[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us_em.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas
    //------------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_em.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) yEMv[i] = yEMm[i];

    //------------------------------------------------------------------------------------
    //Momenta to velocities
    //------------------------------------------------------------------------------------
    //vx
    yEMv[3] = alpha[0]*yEMm[3] + alpha[1]*yEMm[0] + alpha[2]*yEMm[1];
    //vy
    yEMv[4] = alpha[0]*yEMm[4] + alpha[1]*yEMm[1] - alpha[2]*yEMm[0];
    //vz
    yEMv[5] = alpha[0]*yEMm[5] + alpha[1]*yEMm[2];
}


//----------------------------------------------------------------------------------------
// SYS & NC
//----------------------------------------------------------------------------------------
/**
 *  \brief Change the velocities into momenta. Works for both SYS (SEM, EM) and NC (NCEM, NCSEM) coordinates.
 *         For this reason, the routine is not called SYSvSYSm, but rather VELtoMOM.
 **/
void VELtoMOM(double t, const double ySYSv[], double ySYSm[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySYSm[i] = ySYSv[i];

    //------------------------------------------------------------------------------------
    //Velocity to Momenta
    //------------------------------------------------------------------------------------
    //px
    ySYSm[3] = 1.0/alpha[0] * (ySYSv[3] - alpha[1]*ySYSv[0] - alpha[2]*ySYSv[1]);
    //py
    ySYSm[4] = 1.0/alpha[0] * (ySYSv[4] - alpha[1]*ySYSv[1] + alpha[2]*ySYSv[0]);
    //pz
    ySYSm[5] = 1.0/alpha[0] * (ySYSv[5] - alpha[1]*ySYSv[2] );
}

/**
 *  \brief Change the momenta into velocities. Works for both SYS (SEM, EM) and NC (NCEM, NCSEM) coordinates.
 *         For this reason, the routine is not called SYSmSYSv, but rather MOMtoVEL.
 **/
void MOMtoVEL(double t, const double ySYSm[], double ySYSv[], void *params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    QBCP_L* qbp  = (QBCP_L*) params_void;
    double n     =  qbp->us.n;

    //------------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //------------------------------------------------------------------------------------
    double alpha[3];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //------------------------------------------------------------------------------------
    //Position to position
    //------------------------------------------------------------------------------------
    for(int i =0; i<3; i++) ySYSv[i] = ySYSm[i];

    //------------------------------------------------------------------------------------
    //Momenta to velocities
    //------------------------------------------------------------------------------------
    //vx
    ySYSv[3] = alpha[0]*ySYSm[3] + alpha[1]*ySYSm[0] + alpha[2]*ySYSm[1];
    //vy
    ySYSv[4] = alpha[0]*ySYSm[4] + alpha[1]*ySYSm[1] - alpha[2]*ySYSm[0];
    //vz
    ySYSv[5] = alpha[0]*ySYSm[5] + alpha[1]*ySYSm[2];
}


//========================================================================================
//
// COC: Velocities <--> Momenta + NC <-> SYS
//
//========================================================================================
/**
 *  \brief From VSEM to NCSEMm
 **/
void SEMvtoNCSEMm(double t, const double ySEMv[], double yNCSEMm[], QBCP_L *qbcp_l)
{
    double ySEMm[6];
    //VSEM to SEM
    SEMvtoSEMm(t, ySEMv, ySEMm, qbcp_l);
    //SEM to NCSEM
    SEMmtoNCSEMm(t, ySEMm, yNCSEMm, qbcp_l);
}

/**
 *  \brief From VEM to NCEMm
 **/
void EMvtoNCEMm(double t, const double yEMv[], double yNCEMm[], QBCP_L *qbcp_l)
{
    double yEMm[6];
    //VEM to EM
    EMvtoEMm(t, yEMv, yEMm, qbcp_l);
    //EM to NCEM
    EMmtoNCEMm(t, yEMm, yNCEMm, qbcp_l);
}

//========================================================================================
//
// Change of unit SYSTEM
//
//========================================================================================
/**
 *   \brief From SEM units to EM units for a Position/Velocity and time vector in IN coordinates and time
 **/
void ussem2usem(double *tc, double yINv[], QBCP_L *qbcp_l)
{
    //IN[SEM] to IN[EM]
    yINv[0] *= qbcp_l->us_em.as;
    yINv[1] *= qbcp_l->us_em.as;
    yINv[2] *= qbcp_l->us_em.as;
    yINv[3] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[4] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[5] *= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    *tc     /= qbcp_l->us_em.ns;
}

/**
 *   \brief From EM units to SEM units for a Position/Velocity and time vector in IN coordinates and time
 **/
void usem2ussem(double *tc, double yINv[], QBCP_L *qbcp_l)
{
    //IN[EM] to IN[SEM]
    yINv[0] /= qbcp_l->us_em.as;
    yINv[1] /= qbcp_l->us_em.as;
    yINv[2] /= qbcp_l->us_em.as;
    yINv[3] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[4] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    yINv[5] /= qbcp_l->us_em.as*qbcp_l->us_em.ns;
    *tc     *= qbcp_l->us_em.ns;
}

//========================================================================================
//
// INSEM to ECISEM
//
//========================================================================================
/**
 *   \brief Earth state in INSEM coordinates
 **/
void zEarth_INSEM(double t, double zEarth[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Parameters
    //------------------------------------------------------------------------------------
    double n  = qbcp_l->us_sem.n;
    double ms = qbcp_l->us_sem.ms;
    double mm = qbcp_l->us_sem.mm;
    double me = qbcp_l->us_sem.me;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;
    double M = mm + ms + me;

    //------------------------------------------------------------------------------------
    //Retrieve Earth position and velocity in INSEM coordinates
    //------------------------------------------------------------------------------------
    //r
    //double r1 = creal(evz(qbcp_l->cs_sem.zt, t, n, ni, ai));
    //double r2 = cimag(evz(qbcp_l->cs_sem.zt, t, n, ni, ai));
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    //rdot
    //cdouble rdot = evzdot(qbcp_l->cs_sem.zt, qbcp_l->cs_sem.ztdot, t, n, ni, ai);
    //Rdot
    cdouble Rdot = evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as);

    //------------------------------------------------------------------------------------
    //Store results
    //------------------------------------------------------------------------------------
    //    zEarth[0] = - ms / M * R1 + mm/(mm + me) * r1;
    //    zEarth[1] = - ms / M * R2 + mm/(mm + me) * r2;
    //    zEarth[2] = 0.0; //+ 0.0;
    //
    //    zEarth[3] = - ms / M * creal(Rdot) + mm/(mm + me) * creal(rdot);
    //    zEarth[4] = - ms / M * cimag(Rdot) + mm/(mm + me) * cimag(rdot);
    //    zEarth[5] = + 0.0;

    zEarth[0] = - ms / M * R1;
    zEarth[1] = - ms / M * R2;
    zEarth[2] = 0.0; //+ 0.0;

    zEarth[3] = - ms / M * creal(Rdot);
    zEarth[4] = - ms / M * cimag(Rdot);
    zEarth[5] = + 0.0;

    //Sun case!
    //    zEarth[0] = (mm + me) * R1; //- ms / M * R1 + mm/(mm + me) * r1;
    //    zEarth[1] = (mm + me) * R2; //- ms / M * R2 + mm/(mm + me) * r2;
    //    zEarth[2] = 0.0; //+ 0.0;
    //
    //    zEarth[3] = (mm + me) * creal(Rdot); //- ms / M * creal(Rdot) + mm/(mm + me) * creal(rdot);
    //    zEarth[4] = (mm + me) * cimag(Rdot); //- ms / M * cimag(Rdot) + mm/(mm + me) * cimag(rdot);
    //    zEarth[5] = + 0.0;
}

/**
 *   \brief Earth acceleration in INSEM coordinates
 **/
void aEarth_INSEM(double t, double aEarth[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Parameters
    //------------------------------------------------------------------------------------
    double n  = qbcp_l->us_sem.n;
    double ms = qbcp_l->us_sem.ms;
    double mm = qbcp_l->us_sem.mm;
    double me = qbcp_l->us_sem.me;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;
    //double ni = qbcp_l->us_sem.ni;
    //double ai = qbcp_l->us_sem.ai;
    double M  = mm + ms + me;

    //------------------------------------------------------------------------------------
    //Retrieve Earth position and velocity in INSEM coordinates
    //------------------------------------------------------------------------------------
    //rddot
    //cdouble rddot = evzddot(qbcp_l->cs_sem.zt, qbcp_l->cs_sem.ztdot, qbcp_l->cs_sem.ztddot, t, n, ni, ai);
    //double r1ddot = creal(rddot);
    //double r2ddot = cimag(rddot);
    //Rddot
    cdouble Rddot = evzddot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, qbcp_l->cs_sem.Ztddot, t, n, ns, as);
    double R1ddot = creal(Rddot);
    double R2ddot = cimag(Rddot);

    //------------------------------------------------------------------------------------
    //Store results
    //------------------------------------------------------------------------------------
    //    aEarth[0] = - ms / M * R1ddot + mm/(mm + me) * r1ddot;
    //    aEarth[1] = - ms / M * R2ddot + mm/(mm + me) * r2ddot;
    //    aEarth[2] = + 0.0;

    aEarth[0] = - ms / M * R1ddot;
    aEarth[1] = - ms / M * R2ddot;
    aEarth[2] = + 0.0;

    //Sun case!
    //    aEarth[0] = (mm + me) * R1ddot;
    //    aEarth[1] = (mm + me) * R2ddot;
    //    aEarth[2] = + 0.0;
}


/**
 *   \brief From INSEM to ECISEM (Earth-Centered) coordinates
 **/
void INSEMtoECISEM(double t, const double yINSEM[], double yECISEM[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Retrieve Earth position and velocity in INSEM coordinates
    //------------------------------------------------------------------------------------
    double zEarth[6];
    zEarth_INSEM(t, zEarth, qbcp_l);

    //------------------------------------------------------------------------------------
    //Substract from INSEM coordinates to get ECISEM
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yECISEM[i] = yINSEM[i] - zEarth[i];
}


/**
 *   \brief From ECISEM to INSEM coordinates
 **/
void ECISEMtoINSEM(double t, const double yECISEM[], double yINSEM[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Retrieve Earth position and velocity in INSEM coordinates
    //------------------------------------------------------------------------------------
    double zEarth[6];
    zEarth_INSEM(t, zEarth, qbcp_l);

    //------------------------------------------------------------------------------------
    //Add from ECISEM coordinates to get INSEM
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yINSEM[i] = yECISEM[i] + zEarth[i];
}


/**
 *   \brief From INSEM to ECISEM (Earth-Centered) coordinates
 **/
void INSEMtoECISEM(double t, const double yINSEM[], double yECISEM[], double zEarth[])
{
    //------------------------------------------------------------------------------------
    //Substract from INSEM coordinates to get ECISEM
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yECISEM[i] = yINSEM[i] - zEarth[i];
}


/**
 *   \brief From ECISEM to INSEM coordinates
 **/
void ECISEMtoINSEM(double t, const double yECISEM[], double yINSEM[], double zEarth[])
{

    //------------------------------------------------------------------------------------
    //Add from ECISEM coordinates to get INSEM
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) yINSEM[i] = yECISEM[i] + zEarth[i];
}


//========================================================================================
//
// COC: SEM <--> IN <--> EM
//
//========================================================================================

//----------------------------------------------------------------------------------------
// COC: IN <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From EM to IN (in EM units)
 **/
void EMvtoIN(double t, const double yEM[], double yIN[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_em.n;
    double ms = qbcp_l->us_em.ms;
    double ns = qbcp_l->us_em.ns;
    double as = qbcp_l->us_em.as;
    double ni = qbcp_l->us_em.ni;
    double ai = qbcp_l->us_em.ai;
    //r
    double r1 = creal(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r  = sqrt(r1*r1 + r2*r2);

    //R
    double R1 = creal(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double rdot  = 1.0/r*(r1*r1dot + r2*r2dot);
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));

    //------------------------------------------------------------------------------------
    // EM to IN: Position
    //------------------------------------------------------------------------------------
    yIN[0] = r1*yEM[0] - r2*yEM[1] - ms/(1.0+ms)*R1;
    yIN[1] = r2*yEM[0] + r1*yEM[1] - ms/(1.0+ms)*R2;
    yIN[2] = r *yEM[2];

    //------------------------------------------------------------------------------------
    // EM to IN: Velocity
    //------------------------------------------------------------------------------------
    yIN[3] = r1dot*yEM[0] - r2dot*yEM[1] + r1*yEM[3] - r2*yEM[4] - ms/(1.0+ms)*R1dot;
    yIN[4] = r2dot*yEM[0] + r1dot*yEM[1] + r2*yEM[3] + r1*yEM[4] - ms/(1.0+ms)*R2dot;
    yIN[5] = rdot *yEM[2] + r *yEM[5];
}

/**
 * \brief From IN to EM (in EM units)
 **/
void INtoEM(double t, const double yIN[], double yEM[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_em.n;
    double ms = qbcp_l->us_em.ms;
    double ns = qbcp_l->us_em.ns;
    double as = qbcp_l->us_em.as;
    double ni = qbcp_l->us_em.ni;
    double ai = qbcp_l->us_em.ai;

    //r
    double r1 = creal(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r2 = cimag(evz(qbcp_l->cs_em.zt, t, n, ni, ai));
    double r = sqrt(r1*r1 + r2*r2);
    //R
    double R1 = creal(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_em.Zt, t, n, ns, as));
    //rdot
    double r1dot = creal(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    double r2dot = cimag(evzdot(qbcp_l->cs_em.zt, qbcp_l->cs_em.ztdot, t, n, ni, ai));
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_em.Zt, qbcp_l->cs_em.Ztdot, t, n, ns, as));

    //Additional parameters
    double a = +pow(r, -4.0)*(r1dot*r*r - 2*r1*(r1*r1dot + r2*r2dot)); //dot(r1/(r*r))
    double b = +pow(r, -4.0)*(r2dot*r*r - 2*r2*(r1*r1dot + r2*r2dot)); //dot(r2/(r*r))
    double c = -pow(r, -3.0)*(r1*r1dot + r2*r2dot);                    //dot(1/r)

    //------------------------------------------------------------------------------------
    // EM to IN: Position
    //------------------------------------------------------------------------------------
    yEM[0] = 1.0/(r*r) * ( +r1*(yIN[0] + ms/(1.0+ms)*R1) +  r2*(yIN[1] + ms/(1.0+ms)*R2) );
    yEM[1] = 1.0/(r*r) * ( -r2*(yIN[0] + ms/(1.0+ms)*R1) +  r1*(yIN[1] + ms/(1.0+ms)*R2) );
    yEM[2] = 1.0/r * yIN[2];

    //------------------------------------------------------------------------------------
    // EM to IN: Velocity
    //------------------------------------------------------------------------------------
    yEM[3] = +a*(yIN[0] + ms/(1.0+ms)*R1) + b*(yIN[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( +r1*(yIN[3] + ms/(1.0+ms)*R1dot) +  r2*(yIN[4] + ms/(1.0+ms)*R2dot) );
    yEM[4] = -b*(yIN[0] + ms/(1.0+ms)*R1) + a*(yIN[1] + ms/(1.0+ms)*R2)
             + 1.0/(r*r) * ( -r2*(yIN[3] + ms/(1.0+ms)*R1dot) +  r1*(yIN[4] + ms/(1.0+ms)*R2dot) );
    yEM[5] = c*yIN[2] + 1.0/r * yIN[5];
}

//----------------------------------------------------------------------------------------
// COC: IN <--> SEM
//----------------------------------------------------------------------------------------
/**
 * \brief From SEM to IN (in SEM units)
 **/
void SEMtoIN(double t, const double ySE[], double yIN[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_sem.n;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;

    //------------------------------------------------------------------------------------
    //r & R
    //------------------------------------------------------------------------------------
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R = sqrt(R1*R1 + R2*R2);

    //------------------------------------------------------------------------------------
    //Derivatives
    //------------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double Rdot  = 1.0/R*(R1*R1dot + R2*R2dot);

    //------------------------------------------------------------------------------------
    // SE to IN: Position
    //------------------------------------------------------------------------------------
    yIN[0] = R1*ySE[0] - R2*ySE[1];
    yIN[1] = R2*ySE[0] + R1*ySE[1];
    yIN[2] = R *ySE[2];

    //------------------------------------------------------------------------------------
    // SE to IN: Velocity
    //------------------------------------------------------------------------------------
    yIN[3] = R1dot*ySE[0] - R2dot*ySE[1] + R1*ySE[3] - R2*ySE[4];
    yIN[4] = R2dot*ySE[0] + R1dot*ySE[1] + R2*ySE[3] + R1*ySE[4];
    yIN[5] = Rdot *ySE[2] + R *ySE[5];
}

/**
 * \brief From IN to SEM (in SEM units)
 **/
void INtoSEM(double t, const double yIN[], double ySE[],
             QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_sem.n;
    double ns = qbcp_l->us_sem.ns;
    double as = qbcp_l->us_sem.as;

    //------------------------------------------------------------------------------------
    //r & R
    //------------------------------------------------------------------------------------
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    //h
    double h1 = R1;
    double h2 = R2;
    double h  = sqrt(h1*h1 + h2*h2);

    //------------------------------------------------------------------------------------
    //Derivatives
    //------------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    //hdot
    double h1dot = R1dot;
    double h2dot = R2dot;

    //------------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates and SE units
    //------------------------------------------------------------------------------------
    double yesIN[6];
    yesIN[0] = 0.0;
    yesIN[1] = 0.0;
    yesIN[2] = 0.0;
    yesIN[3] = 0.0;
    yesIN[4] = 0.0;
    yesIN[5] = 0.0;

    //Additional parameters
    double a = +pow(h, -4.0)*(h1dot*h*h - 2*h1*(h1*h1dot + h2*h2dot)); //dot(h1/(h*h))
    double b = +pow(h, -4.0)*(h2dot*h*h - 2*h2*(h1*h1dot + h2*h2dot)); //dot(h2/(h*h))
    double c = -pow(h, -3.0)*(h1*h1dot  + h2*h2dot);                   //dot(1/h)

    //------------------------------------------------------------------------------------
    // SE to IN: Position
    //------------------------------------------------------------------------------------
    ySE[0] = 1.0/(h*h) * ( +h1*(yIN[0] - yesIN[0]) +  h2*(yIN[1] - yesIN[1]) );
    ySE[1] = 1.0/(h*h) * ( -h2*(yIN[0] - yesIN[0]) +  h1*(yIN[1] - yesIN[1]) );
    ySE[2] = 1.0/h * yIN[2];

    //------------------------------------------------------------------------------------
    // SE to IN: Velocity
    //------------------------------------------------------------------------------------
    ySE[3] = +a*(yIN[0] - yesIN[0]) + b*(yIN[1] - yesIN[1])
             + 1.0/(h*h) * ( +h1*(yIN[3] - yesIN[3]) +  h2*(yIN[4] - yesIN[4]) );
    ySE[4] = -b*(yIN[0] - yesIN[0]) + a*(yIN[1] - yesIN[1])
             + 1.0/(h*h) * ( -h2*(yIN[3] - yesIN[3]) +  h1*(yIN[4] - yesIN[4]) );
    ySE[5] = +c*yIN[2] + 1.0/h * yIN[5];
}

/**
 * \brief From SEM to ECIN (in SEM units)
 **/
void SEMtoECISEM(double t, const double ySE[], double yIN[], QBCP_L *qbcp_l)
{
    //------------------------------------------------------------------------------------
    //Init
    //------------------------------------------------------------------------------------
    //Param
    double n  = qbcp_l->us_sem.n;
    double ns = qbcp_l->us_sem.ns;
    //double ni = qbcp_l->us_sem.ni;
    double as = qbcp_l->us_sem.as;
    //double ai = qbcp_l->us_sem.ai;
    //double ms = qbcp_l->us_sem.ms;
    //double mm = qbcp_l->us_sem.mm;
    //double me = qbcp_l->us_sem.me;

    //------------------------------------------------------------------------------------
    //r & R
    //------------------------------------------------------------------------------------
    //R
    double R1 = creal(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R2 = cimag(evz(qbcp_l->cs_sem.Zt, t, n, ns, as));
    double R  = sqrt(R1*R1 + R2*R2);
    //r
    //double r1 = creal(evz(qbcp_l->cs_sem.zt, t, n, ni, ai));
    //double r2 = cimag(evz(qbcp_l->cs_sem.zt, t, n, ni, ai));

    //------------------------------------------------------------------------------------
    //Derivatives
    //------------------------------------------------------------------------------------
    //Rdot
    double R1dot = creal(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double R2dot = cimag(evzdot(qbcp_l->cs_sem.Zt, qbcp_l->cs_sem.Ztdot, t, n, ns, as));
    double Rdot  = 1.0/R*(R1*R1dot + R2*R2dot);
    //rdot
    //cdouble rdot = evzdot(qbcp_l->cs_sem.zt, qbcp_l->cs_sem.ztdot, t, n, ni, ai);

    //------------------------------------------------------------------------------------
    //Position & Velocity of the SE barycenter in inertial coordinates (in SE units)
    //------------------------------------------------------------------------------------
    double yesIN[6];
    yesIN[0] = 0.0;
    yesIN[1] = 0.0;
    yesIN[2] = 0.0;
    yesIN[3] = 0.0;
    yesIN[4] = 0.0;
    yesIN[5] = 0.0;

    //------------------------------------------------------------------------------------
    // Earth state
    //------------------------------------------------------------------------------------
    double zEarth[6];
    zEarth_INSEM(t, zEarth, qbcp_l);

    //------------------------------------------------------------------------------------
    // SE to IN: Position
    //------------------------------------------------------------------------------------
    yIN[0] = R1*ySE[0] - R2*ySE[1] + yesIN[0] - zEarth[0];
    yIN[1] = R2*ySE[0] + R1*ySE[1] + yesIN[1] - zEarth[1];
    yIN[2] = R *ySE[2]                        - zEarth[2];

    //------------------------------------------------------------------------------------
    // SE to IN: Velocity
    //------------------------------------------------------------------------------------
    yIN[3] = R1dot*ySE[0] - R2dot*ySE[1] + R1*ySE[3] - R2*ySE[4] + yesIN[3] - zEarth[3];
    yIN[4] = R2dot*ySE[0] + R1dot*ySE[1] + R2*ySE[3] + R1*ySE[4] + yesIN[4] - zEarth[4];
    yIN[5] = Rdot *ySE[2] + R *ySE[5]                                       - zEarth[5];
}


//========================================================================================
//
// COC: INSEM <--> EM, INEM <--> SEM
//
//========================================================================================
/**
 * \brief From INEM to SEM
 **/
void INEMvtoSEMm(double t, const double yINEMv[], double ySEMm[], QBCP_L *qbcp_l)
{
    double tc = t;
    double yINv[6];
    //Copy of yINEMv
    for(int i = 0; i <6; i++) yINv[i] = yINEMv[i];

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SE
    double ySEMv[6];
    INtoSEM(tc, yINv, ySEMv, qbcp_l);

    //Velocities to momenta
    SEMvtoSEMm(tc, ySEMv, ySEMm, qbcp_l);
}

/**
 * \brief From INEM to SEMv
 **/
void INEMvtoSEMv(double t, const double yINEMv[], double ySEMv[], QBCP_L *qbcp_l)
{
    double tc = t;
    double yINv[6];
    //Copy of yINEMv
    for(int i = 0; i <6; i++) yINv[i] = yINEMv[i];

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SE
    INtoSEM(tc, yINv, ySEMv, qbcp_l);
}

/**
 * \brief From INSEM to VEM
 **/
void INSEMvtoEMv(double t, const double yINSEMv[], double yEMv[], QBCP_L *qbcp_l)
{
    double tc = t;
    //SEMv --> IN
    double yINv[6];
    //Copy of yINEMv
    for(int i = 0; i <6; i++) yINv[i] = yINSEMv[i];

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EMv
    INtoEM(tc, yINv, yEMv, qbcp_l);
}

/**
 * \brief From INSEM to EM
 **/
void INSEMvtoEMm(double t, const double yINSEMv[], double yEMm[], QBCP_L *qbcp_l)
{
    double tc = t;
    //SEMv --> IN
    double yINv[6], yEMv[6];
    //Copy of yINEMv
    for(int i = 0; i <6; i++) yINv[i] = yINSEMv[i];

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EMv
    INtoEM(tc, yINv, yEMv, qbcp_l);

    //Velocities to momenta
    EMvtoEMm(tc, yEMv, yEMm, qbcp_l);
}

//========================================================================================
//
// COC: SEM <--> EM
//
//========================================================================================

//----------------------------------------------------------------------------------------
// COC: SEMv <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMvtoEMm(double t, const double ySEMv[], double yEMm[], QBCP_L *qbcp_l)
{
    double tc = t;

    //SEMv --> IN
    double yINv[6];
    SEMtoIN(tc, ySEMv, yINv, qbcp_l);

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EM
    double yEMv[6];
    INtoEM(tc, yINv, yEMv, qbcp_l);

    //Velocities to momenta
    EMvtoEMm(tc, yEMv, yEMm, qbcp_l);
}

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMv(double t, const double yEMm[], double ySEMv[], QBCP_L *qbcp_l)
{
    double tc = t;
    //Momenta to velocities
    double yEMv[6];
    EMmtoEMv(tc, yEMm, yEMv, qbcp_l);

    //EM-->IN
    double yINv[6];
    EMvtoIN(tc, yEMv, yINv, qbcp_l);

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SEMv
    INtoSEM(tc, yINv, ySEMv, qbcp_l);
}



//----------------------------------------------------------------------------------------
// COC: SEMm <--> EMv
//----------------------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMv(double t, const double ySEMm[], double yEMv[], QBCP_L *qbcp_l)
{
    double tc = t;

    //Momenta to velocities
    double ySEv[6];
    SEMmtoSEMv(tc, ySEMm, ySEv, qbcp_l);

    //SE --> IN
    double yINv[6];
    SEMtoIN(tc, ySEv, yINv, qbcp_l);

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EM
    INtoEM(tc, yINv, yEMv, qbcp_l);
}

/**
 * \brief From EMv to SEM
 **/
void EMvtoSEMm(double t, const double yEMv[], double ySEMm[], QBCP_L *qbcp_l)
{
    double tc = t;

    //EM-->IN
    double yINv[6];
    EMvtoIN(tc, yEMv, yINv, qbcp_l);

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SE
    double ySEMv[6];
    INtoSEM(tc, yINv, ySEMv, qbcp_l);

    //Velocities to momenta
    SEMvtoSEMm(tc, ySEMv, ySEMm, qbcp_l);
}

//----------------------------------------------------------------------------------------
// COC: SEMv <--> EMv
//----------------------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMvtoEMv(double t, const double ySEMv[], double yEMv[], QBCP_L *qbcp_l)
{
    double tc = t;
    //SEv --> IN
    double yINv[6];
    SEMtoIN(tc, ySEMv, yINv, qbcp_l);

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EMv
    INtoEM(tc, yINv, yEMv, qbcp_l);
}

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMvtoSEMv(double t, const double yEMv[], double ySEMv[], QBCP_L *qbcp_l)
{
    double tc = t;

    //EMv-->IN
    double yINv[6];
    EMvtoIN(tc, yEMv, yINv, qbcp_l);

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SEMv
    INtoSEM(tc, yINv, ySEMv, qbcp_l);
}



//----------------------------------------------------------------------------------------
// COC: SEM <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From SEM to EM (both in position/momenta form)
 **/
void SEMmtoEMm(double t, const double ySEm[], double yEMm[], QBCP_L *qbcp_l)
{
    double tc = t;

    //Momenta to velocities
    double ySEv[6];
    SEMmtoSEMv(tc, ySEm, ySEv, qbcp_l);

    //SE --> IN
    double yINv[6];
    SEMtoIN(tc, ySEv, yINv, qbcp_l);

    //IN[SE] to IN[SEM]
    ussem2usem(&tc, yINv, qbcp_l);

    //IN --> EM
    double yEMv[6];
    INtoEM(tc, yINv, yEMv, qbcp_l);

    //Velocities to momenta
    EMvtoEMm(tc, yEMv, yEMm, qbcp_l);
}

/**
 * \brief From EM to SEM (both in position/momenta form)
 **/
void EMmtoSEMm(double t, const double yEMm[], double ySEMm[], QBCP_L *qbcp_l)
{
    double tc = t;
    //Momenta to velocities
    double yEMv[6];
    EMmtoEMv(tc, yEMm, yEMv, qbcp_l);

    //EM-->IN
    double yINv[6];
    EMvtoIN(tc, yEMv, yINv, qbcp_l);

    //IN[EM] to IN[SEM]
    usem2ussem(&tc, yINv, qbcp_l);

    //IN-->SE
    double ySEMv[6];
    INtoSEM(tc, yINv, ySEMv, qbcp_l);

    //Velocities to momenta
    SEMvtoSEMm(tc, ySEMv, ySEMm, qbcp_l);
}

//----------------------------------------------------------------------------------------
// COC: SEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From VNCEM to SEMm
 **/
void NCEMvtoSEMm(double tEM, const double yNCEMv[], double ySEMm[], QBCP_L *qbcp_l)
{
    double yEMv[6];
    //VNCEM to VEM
    NCEMvtoEMv(yNCEMv, yEMv, qbcp_l);
    //VEM to SEM
    EMvtoSEMm(tEM, yEMv, ySEMm, qbcp_l);
}

/**
 * \brief From VNCEM to SEMv
 **/
void NCEMvtoSEMv(double tEM, const double yNCEMv[], double ySEMv[], QBCP_L *qbcp_l)
{
    double yEMv[6];
    //VNCEM to VEM
    NCEMvtoEMv(yNCEMv, yEMv, qbcp_l);
    //VEM to SEM
    EMvtoSEMv(tEM, yEMv, ySEMv, qbcp_l);
}


//----------------------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From VNCEM to  NCSEMm
 **/
void NCEMvtoNCSEMm(double tEM, const double yNCEMv[], double yNCSEMm[], QBCP_L *qbcp_l)
{
    double ySEMm[6];
    //VNCEM to SEM
    NCEMvtoSEMm(tEM, yNCEMv, ySEMm, qbcp_l);
    //SEM to NCSEM
    SEMmtoNCSEMm(tEM*qbcp_l->us_em.ns, ySEMm, yNCSEMm, qbcp_l);
}

/**
 * \brief From NC EM to  NC SEM (both in position/momenta form)
 **/
void NCEMmtoNCSEMm(double tEM, const double yNCEMm[], double yNCSEM[], QBCP_L *qbcp_l)
{
    double yEMm[6], ySEMm[6];
    //NC to EM
    NCEMmtoEMm(tEM, yNCEMm, yEMm, qbcp_l);
    //EM to SEM
    EMmtoSEMm(tEM, yEMm, ySEMm, qbcp_l);
    //SEM to NC (careful, the time should be set into SEM units!)
    SEMmtoNCSEMm(tEM*qbcp_l->us_em.ns, ySEMm, yNCSEM, qbcp_l);
}

/**
 * \brief From NC SEM to  NC EM (both in position/momenta form)
 **/
void NCSEMmtoNCEMm(double tSEM, const double yNCSEMm[], double yNCEM[], QBCP_L *qbcp_l)
{
    double ySEMm[6], yEMm[6];
    //NC to EM
    NCSEMmtoSEMm(tSEM, yNCSEMm, ySEMm, qbcp_l);
    //SEM to EM
    SEMmtoEMm(tSEM, ySEMm, yEMm, qbcp_l);
    //EM to NC (careful, the time should be set into EM units!)
    EMmtoNCEMm(tSEM/qbcp_l->us_em.ns, yEMm, yNCEM, qbcp_l);
}

/**
 * \brief From VNCEM to VNCSEM
 **/
void NCEMvtoNCSEMv(double tEM, const double yNCEMv[], double yNCSEMv[], QBCP_L *qbcp_l)
{
    double yEMv[6], ySEMv[6];
    //VNCEM to VEM
    NCEMvtoEMv(yNCEMv, yEMv, qbcp_l);
    //VEM to VSEM
    EMvtoSEMv(tEM, yEMv, ySEMv, qbcp_l);
    //VSEM to VNCSEM
    SEMvtoNCSEMv(ySEMv, yNCSEMv, qbcp_l);
}


//----------------------------------------------------------------------------------------
// COC: SEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM (both in position/momenta form)
 **/
void NCEMmtoSEMm(double t, const double yNCEMm[], double ySEMm[], QBCP_L *qbcp_l)
{
    double yEMm[6];
    //NC to EM
    NCEMmtoEMm(t, yNCEMm, yEMm, qbcp_l);
    //EM to SEM
    EMmtoSEMm(t, yEMm, ySEMm, qbcp_l);
}

/**
 * \brief From SEM to NCEM (both in position/momenta form)
 **/
void SEMmtoNCEMm(double t, const double ySEMm[], double yNCEMm[], QBCP_L *qbcp_l)
{
    double yEMm[6];
    //SEM to EM
    SEMmtoEMm(t, ySEMm, yEMm, qbcp_l);
    //EM to NC (careful, the time should be set into EM units!)
    EMmtoNCEMm(t/qbcp_l->us_em.ns, yEMm, yNCEMm, qbcp_l);
}


//----------------------------------------------------------------------------------------
// COC: NCSEM <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From NCSEM to EM (both in position/momenta form)
 **/
void NCSEMmtoEMm(double t, const double yNCSEMm[], double yEMm[], QBCP_L *qbcp_l)
{
    double ySEMm[6];
    //NC to SEM
    NCSEMmtoSEMm(t, yNCSEMm, ySEMm, qbcp_l);
    //SEM to EM
    SEMmtoEMm(t, ySEMm, yEMm, qbcp_l);
}

/**
 * \brief From EM to  NCSEM (both in position/momenta form)
 **/
void EMmtoNCSEMm(double tEM, const double yEMm[], double yNCSEM[], QBCP_L *qbcp_l)
{
    double ySEMm[6];
    //EM to SEM
    EMmtoSEMm(tEM, yEMm, ySEMm, qbcp_l);
    //SEM to NC (careful, the time should be set into SEM units!)
    SEMmtoNCSEMm(tEM*qbcp_l->us_em.ns, ySEMm, yNCSEM, qbcp_l);
}


//=============================================================================================
//
//          Change of coordinates on vectors
//
//=============================================================================================

//----------------------------------------------------------------------------------------
// COC: SEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NC EM to SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
    }
}

/**
 * \brief From NC EM to SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMm_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMd[k];
        //time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC EM to SEMv for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoSEMv_vec(double **yNCEM, double *tNCEM, double **ySEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySEMm[6], ySEMv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to SEM
        NCEMmtoSEMm(tNCEM[p], yNCd, ySEMm, qbcp_l);
        //SEM to SEMv
        SEMmtoSEMv(tSEM[p], ySEMm, ySEMv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySEM[k][p] = ySEMv[k];
    }
}

//----------------------------------------------------------------------------------------
// COC: NCSEM <--> NCEM
//----------------------------------------------------------------------------------------
/**
 * \brief From NC EM to NC SEM for vectors. The time vector tNCEM is transformed into SEM units and stored into tSEM.
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, double *tSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCSEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to NC SEM
        NCEMmtoNCSEMm(tNCEM[p], yNCd, yNCSEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCSEM[k][p] = yNCSEMd[k];
        //time
        tSEM[p] = tNCEM[p]*qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC EM to NC SEM for vectors. Note that the time vector is left unchanged (still in EM units).
 **/
void NCEMmtoNCSEMm_vec(double **yNCEM, double *tNCEM, double **yNCSEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCSEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCEM[k][p];
        //NC EM to NC SEM
        NCEMmtoNCSEMm(tNCEM[p], yNCd, yNCSEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCSEM[k][p] = yNCSEMd[k];
    }
}

/**
 * \brief From NC SEM to NC EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoNCEMm(tNCSEM[p], yNCd, yNCEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCEM[k][p] = yNCEMd[k];
    }
}

/**
 * \brief From NC SEM to NC EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoNCEMm_vec(double **yNCSEM, double *tNCSEM, double **yNCEM, double *tEM,  int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoNCEMm(tNCSEM[p], yNCd, yNCEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yNCEM[k][p] = yNCEMd[k];
        //time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
    }
}


//----------------------------------------------------------------------------------------
// COC: NCSEM <--> EM
//----------------------------------------------------------------------------------------
/**
 * \brief From NC SEM to EM for vectors. Note that the time vector is left unchanged (still in SEM units).
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMd[k];
    }
}

/**
 * \brief From NC SEM to EM for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMm_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC EM to NC SEM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMd, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMd[k];
        //time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
    }
}

/**
 * \brief From NC SEM to EMv for vectors. The time vector tNCSEM is transformed into EM units and stored into tEM.
 **/
void NCSEMmtoEMv_vec(double **yNCSEM, double *tNCSEM, double **yEM, double *tEM, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yEMm[6], yEMv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Time
        tEM[p] = tNCSEM[p]/qbcp_l->us_em.ns;
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNCSEM[k][p];
        //NC SEM to EM
        NCSEMmtoEMm(tNCSEM[p], yNCd, yEMm, qbcp_l);
        //EM to EMv
        EMmtoEMv(tEM[p], yEMm, yEMv, qbcp_l);
        //Copy result in yNCSEM
        for(int k = 0; k <6; k++) yEM[k][p] = yEMv[k];
    }
}


//----------------------------------------------------------------------------------------
// COC: NC <--> SYS
//----------------------------------------------------------------------------------------
/**
 *  \brief From NC to SYS, in SYS units (either EM or SEM).
 **/
void NCtoSYS_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSd[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC EM to SEM
        NCtoSYS(tNC[p], yNCd, ySYSd, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySYS[k][p] = ySYSd[k];
    }
}

//----------------------------------------------------------------------------------------
// COC: NC <--> SYS (NOT USED, HERE FOR BACKUP)
//----------------------------------------------------------------------------------------
/**
 *  \brief From NC to SYSv, in SYS units (either EM or SEM). NOT USED
 **/
void NCtoSYSv_vec(double **yNC, double *tNC, double **ySYS, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSm[6], ySYSv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC to SYS
        NCtoSYS(tNC[p], yNCd, ySYSm, qbcp_l);
        //SYSm to SYSv
        MOMtoVEL(tNC[p], ySYSm, ySYSv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) ySYS[k][p] = ySYSv[k];
    }
}

/**
 *  \brief From SYS to NC, in SYS units (either EM or SEM). NOT USED.
 **/
void SYStoNC_vec(double **ySYS, double *tSYS, double **yNC, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], ySYSd[6];
    //Loop on all elements in ySYS
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in ySYS
        for(int k = 0; k <6; k++) ySYSd[k] = ySYS[k][p];
        //SYS to NC
        SYStoNC(tSYS[p], ySYSd, yNCd, qbcp_l);
        //Copy result in yNC
        for(int k = 0; k <6; k++) yNC[k][p] = yNCd[k];
    }
}


//----------------------------------------------------------------------------------------
// COC: NC <--> VNC (NOT USED, HERE FOR BACKUP)
//----------------------------------------------------------------------------------------
/**
 *  \brief From NC to VNC, in NC units (either EM or SEM). NOT USED
 **/
void NCtoVNC_vec(double **yNC, double *tNC, double **yVNC, int N, QBCP_L *qbcp_l)
{
    double yNCd[6], yNCv[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCd[k] = yNC[k][p];
        //NC to VNC
        MOMtoVEL(tNC[p], yNCd, yNCv, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) yVNC[k][p] = yNCv[k];
    }
}

/**
 *  \brief From VNC to NC, in NC units (either EM or SEM). NOT USED
 **/
void VNCtoNC_vec(double **yVNC, double *tNC, double **yNC, int N, QBCP_L *qbcp_l)
{
    double yNCv[6], yNCm[6];
    //Loop on all elements in yNCEM
    for(int p = 0; p <= N; p++)
    {
        //Copy step p in yNCd
        for(int k = 0; k <6; k++) yNCv[k] = yVNC[k][p];
        //VNC to NC
        VELtoMOM(tNC[p], yNCv, yNCm, qbcp_l);
        //Copy result in ySEM
        for(int k = 0; k <6; k++) yNC[k][p] = yNCm[k];
    }
}

//------------------------------------------------------------------------------------
//   Precisions
//------------------------------------------------------------------------------------
/**
 *  \brief Sets an average precision in cout.
 **/
void coutmp()
{
    cout <<  setprecision(2) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a big precision in cout.
 **/
void coutlp()
{
    cout <<  setprecision(15) << std::showpos  <<  setiosflags(ios::scientific);
}

/**
 *  \brief Sets a small precision in cout.
 **/
void coutsp()
{
    cout <<   setprecision(3) << resetiosflags(ios::scientific);
}
//------------------------------------------------------------------------------------
//   Print
//------------------------------------------------------------------------------------
/**
 *  \brief Prints an array of double using cout. More precise than vector_printf.
 **/
void vector_printf_prec(const double *y, int n)
{
    coutlp();
    for(int i = 0; i < n; i ++) cout << i << "  " << y[i] << endl;
}

/**
 *  \brief Prints an array of double using cout.
 **/
void vector_printf(double const *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << y[i] << endl;
}

/**
 *  \brief Prints an array of complex double using cout.
 **/
void vector_complex_printf(const cdouble *y, int n)
{
    coutmp();
    for(int i = 0; i < n; i ++) cout << i << "  " << creal(y[i]) << "  " << cimag(y[i]) << endl;
}

/**
 *  \brief Prints an array of complex double using cout. More precise than vector_complex_printf.
 **/
void vector_complex_printf_prec(const cdouble *y, int n)
{
    coutlp();
    for(int i = 0; i < n; i ++) cout << i << "  " << creal(y[i]) << "  " << cimag(y[i]) << endl;
}

/**
 *  \brief Prints a matrix of complex double using cout.
 **/
void matrix_complex_printf(const cdouble **y, int n, int m)
{
    coutlp();
    for(int i = 0; i < n; i ++)
    {
        cout << i;
        for(int j = 0; j < m; j++) printf(" %+1.2e%+1.2ei", creal(y[i][j]), cimag(y[i][j]));
        cout << endl;
    }
}

//------------------------------------------------------------------------------------
//   Norm
//------------------------------------------------------------------------------------
/**
 *  Euclidian norm computed on the first k components of a complex vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(const cdouble z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

/**
 *  Euclidian norm computed on the first k components of a double vector: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  z_0[p] ^2 \right)^{-1/2} \f$
 **/
double ENorm(const double z0[], int k)
{
    double res = cabs(z0[0])*cabs(z0[0]);
    for(int p = 1; p < k; p++) res += cabs(z0[p])*cabs(z0[p]);
    return(sqrt(res));
}

/**
 *  Euclidian norm of the difference of two double vectors, computed on the first k: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  (z_1[p] - z_2[p]) ^2 \right)^{-1/2} \f$
 **/
double DENorm(const double z1[], const double z2[], int k)
{
    double diff[k];
    for(int p = 0; p < k; p++) diff[p] = z1[p] - z2[p];
    return(ENorm(diff, k));
}

/**
 *  Euclidian norm of the difference of two cdouble vectors, computed on the first k: \f$ ENorm(z_0, k) = \left( \sum_{p = 0}^{k-1}  (z_1[p] - z_2[p]) ^2 \right)^{-1/2} \f$
 **/
double DENorm(const cdouble z1[], const cdouble z2[], int k)
{
    cdouble diff[k];
    for(int p = 0; p < k; p++) diff[p] = z1[p] - z2[p];
    return(ENorm(diff, k));
}

/**
 *  Euclidian norm of the difference of two cdouble matrices.
 **/
double DENorm(const cdouble **z1, const cdouble **z2, int k, int l)
{
    cdouble diff[k*l];
    for(int i = 0; i < k; i++)
    for(int j = 0; j < l; j++) diff[i*l+j] = z1[i][j] - z2[i][j];
    return(ENorm(diff, k*l));
}

/**
 *  Normalize the vector z0, containing k components
 **/
void ENorm(const double z0[], double z0n[], int k)
{
    double nz0 = ENorm(z0, k);
    for(int p = 0; p < k; p++) z0n[p] = z0[p] / nz0;
}

//------------------------------------------------------------------------------------
//   Default parameterization
//------------------------------------------------------------------------------------
/**
 * \brief Get the default coordinates system for variational equations, from the coord_type
 **/
int default_coordinate_system(int coord_type)
{
    switch(coord_type)
    {
    case PSEM:
        return I_PSEM;
    case NCSEM:
        return I_NCSEM;
    case VNCSEM:
        return I_VNCSEM;
    case PEM:
        return I_PEM;
    case NCEM:
        return I_NCEM;
    case VNCEM:
        return I_VNCEM;
    case VECLI:
        return I_ECLI;
    case J2000:
        return I_J2000;
    case NJ2000:
        return I_NJ2000;
    case VSEM:
        return I_VSEM;
    case VEM:
        return I_VEM;
    case INSEM:
        return I_INSEM;
    case ECISEM:
        return I_ECISEM;
    case INEM:
        return I_INEM;
    default:
        cerr << "default_coordinate_system. Unknown type of coordinates." << endl;
        return FTC_FAILURE;
    }
}

/**
 * \brief Get the default type of coordinates, from the dcs (inverse of default_coordinate_system).
 **/
int default_coordinate_type(int coord_type)
{
    switch(coord_type)
    {
    case I_PSEM:
        return PSEM;
    case I_NCSEM:
        return NCSEM;
    case I_VNCSEM:
        return VNCSEM;
    case I_PEM:
        return PEM;
    case I_NCEM:
        return NCEM;
    case I_VNCEM:
        return VNCEM;
    case I_ECLI:
        return VECLI;
    case I_J2000:
        return J2000;
    case I_NJ2000:
        return NJ2000;
    case I_VSEM:
        return VSEM;
    case I_VEM:
        return VEM;
    case I_INSEM:
        return INSEM;
    case I_ECISEM:
        return ECISEM;
    case I_INEM:
        return INEM;
    default:
        cerr << "default_coordinate_type. Unknown dcs." << endl;
        return FTC_FAILURE;
    }
}

/**
 * \brief Get the default framework, from the coord_type
 **/
int default_framework(int coord_type)
{
    switch(coord_type)
    {
    case VNCSEM:
    case NCSEM:
    case PSEM:
    case VSEM:
    case INSEM:
    case ECISEM:
        return F_SEM;
    case VNCEM:
    case NCEM:
    case PEM:
    case VEM:
    case INEM:
        return F_EM;
    case VECLI:
    case J2000:
    case NJ2000:
        return SEML.fwrk;
    default:
        cerr << "default_framework. Unknown type of coordinates." << endl;
        return FTC_FAILURE;
    }
}


/**
 * \brief Get the default framework, from the dcs
 **/
int default_framework_dcs(int dcs)
{
    switch(dcs)
    {
    case I_VNCSEM:
    case I_NCSEM:
    case I_PSEM:
    case I_VSEM:
    case I_INSEM:
    case I_ECISEM:
        return F_SEM;
    case I_VNCEM:
    case I_NCEM:
    case I_PEM:
    case I_VEM:
    case I_INEM:
        return F_EM;
    case I_ECLI:
    case I_J2000:
    case I_NJ2000:
        return SEML.fwrk;
    default:
        cerr << "default_framework_dcs. Unknown tdcs." << endl;
        return FTC_FAILURE;
    }
}
