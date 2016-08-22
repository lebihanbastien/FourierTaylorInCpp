#include "timec.h"

/**
 * \file timec.cpp
 * \brief Time manipulation similar to MATLAB routines.
 * \author BLB, Zubin Olikara
 * \date May 2015
 * \version 1.0
 */

struct timespec TIC_TIME;

//Equivalent to the tic function in Matlab
void tic(void)
{
    clock_gettime(CLOCK_REALTIME, &TIC_TIME);
}

//Equivalent to the toc function in Matlab
double toc(void)
{
    struct timespec toc_time;
    clock_gettime(CLOCK_REALTIME, &toc_time);
    return (double)(toc_time.tv_sec - TIC_TIME.tv_sec) +
           (double)(toc_time.tv_nsec - TIC_TIME.tv_nsec) / 1e9;
}
