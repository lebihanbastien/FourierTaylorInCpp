#include "ftc_errno.h"


const char*
ref_strerror (const int ref_errno)
{
    switch (ref_errno)
    {
    //------------------------------------------------------------------------------------
    // Generic codes
    //------------------------------------------------------------------------------------
    case FTC_SUCCESS:
        return "success" ;
    case FTC_FAILURE:
        return "failure" ;
    case FTC_EDOM:
        return "input domain error" ;
    case FTC_ENOENT:
        return "no such file or directory" ;

    //------------------------------------------------------------------------------------
    //For refinement procedures
    //------------------------------------------------------------------------------------
    case REF_EOUTOFDPC:
        return "out of the domain of practical convergence (DPC) of the semi-analytical tools";
    case REF_EMAXITREACH:
        return "maximum iteration number reached" ;

    //------------------------------------------------------------------------------------
    //For orbit procedures
    //------------------------------------------------------------------------------------
    case ORBIT_EINT:
        return "integration failed";
    case ORBIT_EPROJ:
        return "projection procedure failed" ;
    default:
        return "unknown error code" ;
    }
}



