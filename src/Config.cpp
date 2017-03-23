#include "Config.h"

//----------------------------------------------------------------------------------------
//Constructor
//----------------------------------------------------------------------------------------
/**
 *  \brief Initialize all necessary configuration parameters
 **/
Config::Config()
{
    PREC_ABS    = 1e-10;
    PREC_REL    = 1e-10;
    PREC_ROOT   = 1e-13;
    PREC_DIFF   = 1e-12;
    PREC_HSTART = 1e-8;
    DC_ITERMAX  = 20;

    PREC_ABS_BACKUP = PREC_ABS;
    PREC_REL_BACKUP = PREC_REL;

    //C_PREC_HARD(); force hard precision?
}


//----------------------------------------------------------------------------------------
//Setters
//----------------------------------------------------------------------------------------
/**
 *  \brief Hard precision for integration purposes
 **/
void Config::C_PREC_HARD()
{
    //Backup
    PREC_ABS_BACKUP = PREC_ABS;
    PREC_REL_BACKUP = PREC_REL;

    //New values
    PREC_ABS = 1e-14;
    PREC_REL = 1e-14;
    cout << "Config manager: the integration precisions (absolute and relative)";
    cout << " have been set to 1e-14" << endl;
}

/**
 *  \brief Soft precision for integration purposes, to fasten correction loops
 *         It is the default case.
 **/
void Config::C_PREC_SOFT()
{
    //Backup
    PREC_ABS_BACKUP = PREC_ABS;
    PREC_REL_BACKUP = PREC_REL;

    //New values
    PREC_ABS = 1e-10;
    PREC_REL = 1e-10;
    cout << "Config manager: the integration precisions (absolute and relative)";
    cout << " have been set to 1e-10" << endl;
}

/**
 *  \brief Use backup values for integration precision.
 **/
void Config::C_PREC_BACK()
{
    //New values
    PREC_ABS = PREC_ABS_BACKUP;
    PREC_REL = PREC_REL_BACKUP;
    cout << "Config manager: the integration precisions (absolute and relative)";
    cout << " have been set to their backup values" << endl;
}

//----------------------------------------------------------------------------------------
//Destructor
//----------------------------------------------------------------------------------------
Config::~Config()
{
    //dtor
}
