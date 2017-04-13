#include "ephemerides.h"


//========================================================================================
//                         Subroutines to choose between SEM and EM representation
//========================================================================================
/**
 * \brief Get the ephemerides coordinate system associated
 *        to the current coordinate system.
 *        For now, only the only possible outputs are VSEM/VEM.
 *        In the future, VNCSEM and VNCEM should be made available.
 **/
int eph_coord(int coord_type)
{
    switch(coord_type)
    {
    case NCSEM  :
    case VNCSEM :
    case PSEM   :
    case VSEM   :
    case ECISEM :
        return VSEM;

    case NCEM   :
    case VNCEM  :
    case PEM    :
    case VEM    :
        return VEM;

    default:
        cout << "eph_coord. Unknown type of coordinates. GSL_FAILURE is returned." << endl;
        return GSL_FAILURE;
    }
}

/**
 * \brief Get the ephemerides fwrk associated to the current coordinate system.
 *        For now, only the only possible outputs are I_VSEM/I_VEM.
 *        In the future, I_VNCSEM and I_VNCEM should be made available.
 **/
int eph_fwrk(int coord_type)
{
    switch(coord_type)
    {
    case NCSEM  :
    case VNCSEM :
    case PSEM   :
    case VSEM   :
        return I_VSEM;

    case NCEM   :
    case VNCEM  :
    case PEM    :
    case VEM    :
        return I_VEM;

    default:
        cout << "eph_coord. Unknown type of coordinates. GSL_FAILURE is returned." << endl;
        return GSL_FAILURE;
    }
}

/**
 *    \brief Mean motion, in rad/s, assciated to the ephemerides coordinates coord_eph
 **/
double mean_motion(int coord_eph)
{
    switch(coord_eph)
    {
    case VNCSEM :
    case VSEM   :
        return SEML.n_sem;

    case VNCEM  :
    case VEM    :
        return SEML.n_em;
    case J2000  :
    case NJ2000 :
        return SEML.ss->n;
    case VECLI  : //no normalization!
        return 1.0;
    default:
        cout << "mean_motion. Unknown type of coordinates. GSL_FAILURE is returned." << endl;
        return GSL_FAILURE;
    }
}


//========================================================================================
//                         Looking for the best fit wrt to a given
//                              QBCP configuration
//========================================================================================
/**
 *  \brief Find the epoch that gives the best correspondance between the Sun-Earth-Moon QBCP configuration at time t, given in SEM/EM coordinates, and the Sun-Earth-Moon true motion
 *         given by the JPL ephemerides. The input tSYS must be in SEM/EM normalized coordinates. The output epoch et is given in seconds.
 **/
void qbcp2jpl(double tSYS, double *et, int coord_type)
{
    //====================================================================================
    //Precheck on the coord_type
    //====================================================================================
    int coord_eph = eph_coord(coord_type);
    if(coord_eph == GSL_FAILURE)
    {
        cout << "qbcp2jpl. There is no ephemerides coordinates" << endl;
        cout << " associated to the current coordinate system. return." << endl;
        return;
    }


    //====================================================================================
    //Compute the positions of primaries in system coordinates
    //====================================================================================
    double Pe[3], Pm[3], Ps[3];
    switch(coord_eph)
    {
        case VEM:
            evaluateCoef(Pe, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Pe, 3);
            evaluateCoef(Pm, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Pm, 3);
            evaluateCoef(Ps, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Ps, 3);
            break;
        case VSEM:
            evaluateCoef(Pe, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
            evaluateCoef(Pm, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
            evaluateCoef(Ps, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Ps, 3);
            break;
    }


    //====================================================================================
    // Init
    //====================================================================================
    //============================================
    // Starting and end time
    //============================================
    //Covered time span in 432s.bsp is:
    //Start:  1949 DEC 14 00:00:00.000 (TDB)
    //Stop:   2050 JAN 02 00:00:00.000 (TDB)
    string epoch_start, epoch_end;
    SpiceDouble  et_start, et_end;

    epoch_start = "2046 MAR 01 00:00:00.000"; //We start in 1950
    epoch_end   = "2046 JUN 01 00:00:00.000"; //We end in 1950

    //In ephemeris time
    str2et_c(epoch_start.c_str(), &et_start);
    str2et_c(epoch_end.c_str(), &et_end);

    //============================================
    // Initialize the change of coordinates:
    // Ecliptic coordinates  <-> SEM/EM synodic in position.
    //============================================
    double B[3];
    double C[3][3], k;
    SpiceDouble lt, RE[6], RM[6], RS[6], re[3], rm[3], rs[3];

    //============================================
    //Initialize the minimum distance
    //============================================
    double minDist    = 0.0;
    double argMinDist = 0.0;
    double dist       = 0.0;

    //============================================
    //Frequence of check
    //============================================
    int freq = 60;//spd_c();

    //====================================================================================
    // Loop in the Kernel, each day
    //====================================================================================
    *et = et_start;
    do
    {
        //Initialize the change of coordinates
        init_ecl2synpos(*et, B, C, &k, coord_eph);

        //Compute the ecliptic positions of the Earth & the Moon
        spkezr_c ("EARTH", *et, DEFFRAME,  "NONE",  DEFOBS, RE, &lt);
        spkezr_c ("MOON",  *et, DEFFRAME,  "NONE",  DEFOBS, RM, &lt);
        spkezr_c ("SUN",   *et, DEFFRAME,  "NONE",  DEFOBS, RS, &lt);

        //Back in synodical coordinates
        ecl2synpos(RE, re, B, C, k);
        ecl2synpos(RM, rm, B, C, k);
        ecl2synpos(RS, rs, B, C, k);

        //Compute the minimum distance
        dist = 0.0;
        for(int i = 0; i < 3; i++) dist += (Pe[i] - re[i])*(Pe[i] - re[i]);
        for(int i = 0; i < 3; i++) dist += (Pm[i] - rm[i])*(Pm[i] - rm[i]);
        for(int i = 0; i < 3; i++) dist += (Ps[i] - rs[i])*(Ps[i] - rs[i]);

        if(*et == et_start)
        {
            minDist    = dist;
            argMinDist = *et;
        }
        else
        {
            if(dist < minDist)
            {
                minDist    = dist;
                argMinDist = *et;
            }
        }

        //Add one freq
        *et += freq;//spd_c();
    }while(*et <= et_end);


    //============================================
    // Update results
    //============================================
    *et = argMinDist;

    //============================================
    // Display results
    //============================================
    //        //Initialize the change of coordinates
    //        init_ecl2synpos(*et, B, C, &k, coord_eph);
    //
    //        //Compute the ecliptic positions of the Earth & the Moon
    //        spkezr_c ("EARTH", *et, DEFFRAME,  "NONE", DEFOBS, RE, &lt);
    //        spkezr_c ("MOON",  *et, DEFFRAME,  "NONE", DEFOBS, RM, &lt);
    //        spkezr_c ("SUN",   *et, DEFFRAME,  "NONE", DEFOBS, RS, &lt);
    //
    //        //Back in synodical coordinates
    //        ecl2synpos(RE, re, B, C, k);
    //        ecl2synpos(RM, rm, B, C, k);
    //        ecl2synpos(RS, rs, B, C, k);
    //
    //        cout << "Earth position at best fit: " << endl;
    //        vector_printf_prec(re, 3);
    //
    //        cout << "Moon position at best fit: " << endl;
    //        vector_printf_prec(rm, 3);
    //
    //        cout << "Sun position at best fit: " << endl;
    //        vector_printf_prec(rs, 3);
    //    //============================================

}

/**
 *  \brief Find the epoch that gives the best correspondance between the Sun-Earth-Moon QBCP configuration at time t, given in SEM/EM coordinates, and the Sun-Earth-Moon true motion
 *         given by the JPL ephemerides. The input tSYS must be in SEM/EM normalized coordinates. The output epoch et is given in seconds. Display routine (test).
 **/
void qbcp2jpl_disp(double tSYS, double *et, int coord_type)
{
    //====================================================================================
    //Precheck on the coord_type
    //====================================================================================
    int coord_eph = eph_coord(coord_type);
    if(coord_eph == GSL_FAILURE)
    {
        cout << "qbcp2jpl. There is no ephemerides coordinates" << endl;
        cout << " associated to the current coordinate system. return." << endl;
        return;
    }


    //====================================================================================
    //Compute the positions of primaries in system coordinates
    //====================================================================================
    double Pe[3], Pm[3], Ps[3];
    switch(coord_eph)
    {
        case VEM:
            evaluateCoef(Pe, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Pe, 3);
            evaluateCoef(Pm, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Pm, 3);
            evaluateCoef(Ps, tSYS, SEML.us_em.n, SEML.nf, SEML.cs_em.Ps, 3);
            break;
        case VSEM:
            evaluateCoef(Pe, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pe, 3);
            evaluateCoef(Pm, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Pm, 3);
            evaluateCoef(Ps, tSYS, SEML.us_sem.n, SEML.nf, SEML.cs_sem.Ps, 3);
            break;
    }

    //====================================================================================
    //Display
    //====================================================================================
    cout << "-----------------------------"<< endl;
    cout << "Pe = " << endl;
    vector_printf_prec(Pe, 3);
    cout << "-----------------------------"<< endl;
    cout << "Pm = " << endl;
    vector_printf_prec(Pm, 3);
    cout << "-----------------------------"<< endl;
    cout << "Ps = " << endl;
    vector_printf_prec(Ps, 3);

    //====================================================================================
    //Gnuplot
    //====================================================================================
    gnuplot_ctrl *h1, *h2;
    h1 = gnuplot_init(true);
    h2 = gnuplot_init(true);

    //Best fit
    gnuplot_cmd(true, h1, "set title \"Best fit (SEM coordinates)\"");
    gnuplot_set_xlabel(true, h1, (char*) "Xsem");
    gnuplot_set_ylabel(true, h1, (char*) "Ysem");
    gnuplot_set_zlabel(true, h1, (char*) "Zsem");
    gnuplot_plot_xyz(true, h1, &Pe[0], &Pe[1], &Pe[2], 1, (char*) "Earth", "point", "1", "2", 1);
    gnuplot_plot_xyz(true, h1, &Pm[0], &Pm[1], &Pm[2], 1, (char*) "Moon", "point",  "1", "2", 2);
    gnuplot_plot_xyz(true, h1, &Ps[0], &Ps[1], &Ps[2], 1, (char*) "Sun", "point",  "1", "2",  3);

    //Best fit
    gnuplot_cmd(true, h2, "set title \"Fit error\"");
    gnuplot_set_xlabel(true, h2, (char*) "t (s)");
    gnuplot_set_ylabel(true, h2, (char*) "Fit error");


    //====================================================================================
    // Init
    //====================================================================================
    //============================================
    // Starting and end time
    //============================================
    //Covered time span in 432s.bsp is:
    //Start:  1949 DEC 14 00:00:00.000 (TDB)
    //Stop:   2050 JAN 02 00:00:00.000 (TDB)
    string epoch_start, epoch_end;
    SpiceDouble  et_start, et_end;

    epoch_start = "2000 JAN 01 00:00:00.000"; //We start in 2000
    epoch_end   = "2001 JAN 01 00:00:00.000"; //We end in 2001

    //In ephemeris time
    str2et_c(epoch_start.c_str(), &et_start);
    str2et_c(epoch_end.c_str(), &et_end);

    //============================================
    // Initialize the change of coordinates:
    // Ecliptic coordinates  <-> SEM/EM synodic in position.
    //============================================
    double B[3];
    double C[3][3], k;
    SpiceDouble lt, RE[6], RM[6], RS[6], re[3], rm[3], rs[3];

    //============================================
    //Initialize the minimum distance
    //============================================
    double minDist    = 0.0;
    double argMinDist = 0.0;
    double dist       = 0.0;

    //============================================
    //Frequence of check
    //============================================
    int freq = 60;//spd_c();

    //============================================
    //Store the results
    //============================================
    int n_sol  = floor((et_end - et_start)/freq) + 1;
    double *fit_error   = dvector(0, n_sol-1);
    double *fit_error_t = dvector(0, n_sol-1);
    int id = 0;

    //====================================================================================
    // Loop in the Kernel, each day
    //====================================================================================
    *et = et_start;
    do
    {
        //Initialize the change of coordinates
        init_ecl2synpos(*et, B, C, &k, coord_eph);

        //Compute the ecliptic positions of the Earth & the Moon
        spkezr_c ("EARTH", *et, DEFFRAME,  "NONE",  DEFOBS, RE, &lt);
        spkezr_c ("MOON",  *et, DEFFRAME,  "NONE",  DEFOBS, RM, &lt);
        spkezr_c ("SUN",   *et, DEFFRAME,  "NONE",  DEFOBS, RS, &lt);

        //Back in synodical coordinates
        ecl2synpos(RE, re, B, C, k);
        ecl2synpos(RM, rm, B, C, k);
        ecl2synpos(RS, rs, B, C, k);

        //Compute the minimum distance
        dist = 0.0;
        for(int i = 0; i < 3; i++) dist += (Pe[i] - re[i])*(Pe[i] - re[i]);
        for(int i = 0; i < 3; i++) dist += (Pm[i] - rm[i])*(Pm[i] - rm[i]);
        for(int i = 0; i < 3; i++) dist += (Ps[i] - rs[i])*(Ps[i] - rs[i]);

        if(*et == et_start)
        {
            minDist    = dist;
            argMinDist = *et;
        }
        else
        {
            if(dist < minDist)
            {
                minDist    = dist;
                argMinDist = *et;
            }
        }

        //Store result
        fit_error[id]   = dist;
        fit_error_t[id] = *et;

        //Add one freq
        *et += freq;//spd_c();
        id++;

    }while(*et <= et_end);


    //============================================
    // Update results
    //============================================
    *et = argMinDist;

    //====================================================================================
    // Display results
    //====================================================================================
    //Initialize the change of coordinates
    init_ecl2synpos(*et, B, C, &k, coord_eph);

    //Compute the ecliptic positions of the Earth & the Moon
    spkezr_c ("EARTH", *et, DEFFRAME,  "NONE", DEFOBS, RE, &lt);
    spkezr_c ("MOON",  *et, DEFFRAME,  "NONE", DEFOBS, RM, &lt);
    spkezr_c ("SUN",  *et, DEFFRAME,  "NONE", DEFOBS, RS, &lt);

    //Back in synodical coordinates
    ecl2synpos(RE, re, B, C, k);
    ecl2synpos(RM, rm, B, C, k);
    ecl2synpos(RS, rs, B, C, k);

    //Best fit
    gnuplot_plot_xyz(true, h1, &re[0], &re[1], &re[2], 1, (char*) "Earth (JPL)", "point",  "2", "2", 1);
    gnuplot_plot_xyz(true, h1, &rm[0], &rm[1], &rm[2], 1, (char*) "Moon (JPL)",  "point",  "2", "2", 2);
    gnuplot_plot_xyz(true, h1, &rs[0], &rs[1], &rs[2], 1, (char*) "Sun (JPL)",   "point",  "2", "2", 3);

    //All the results
    string title = "Best fit = " + numTostring(minDist);
    gnuplot_plot_xy(true, h2, fit_error_t, fit_error, n_sol, (char*) "",  "lines",  "2", "2", 2);
    gnuplot_plot_xy(true, h2, &argMinDist, &minDist, 1, (char*) title.c_str(),  "point",  "2", "2", 2);

    cout << "--------------------------------------------------------" << endl;
    cout << "qbcp2jpl. Closest Positions of the primaries from JPL ephemerides" << endl;
    cout << "Earth:" << endl;
    vector_printf_prec(re, 3);
    cout << "Moon:" << endl;
    vector_printf_prec(rm, 3);
    cout << "--------------------------------------------------------" << endl;

    //====================================================================================
    // Separation angle between BemS and ME, in QBCP
    //====================================================================================
    double Bem[3], BemS[3], ME[3];
    for(int i = 0; i <3; i++)
    {
        //Barycenter of the Earth-Moon system
        Bem[i] = (SEML.us_sem.me*Pe[i] + SEML.us_sem.mm*Pm[i])/(SEML.us_sem.me + SEML.us_sem.mm);
        //BemS = BS - BBem
        BemS[i] = Ps[i] - Bem[i];
        //ME = BE - BM
        ME[i] = Pe[i] - Pm[i];
    }
    double gamma = vsep_c(BemS, ME)*180/M_PI;


    //====================================================================================
    // Separation angle between BemS and ME, in JPL
    //====================================================================================
    double me[3], r[3], rp[3], mag, e3[3];
    //Initialize the change of coordinates for EM state
    init_ecl2synpos(*et, B, C, &k, VEM);

    //Barycenter of the Earth-Moon system, in JPL ephemerides
    double RB[6], RBS[3], RME[3];
    spkezr_c ("EARTH MOON BARYCENTER", *et, DEFFRAME,  "NONE", DEFOBS, RB, &lt);
    for(int i = 0; i <3;i++) RBS[i] = RS[i] - RB[i];
    //ME line
    for(int i = 0; i <3;i++) RME[i] = RE[i] - RM[i];

    //r = RBS/|RBS|
    unorm_c(RBS, r, &mag);
    //em = RME/|RME|
    unorm_c(RME, me, &mag);
    //e3 = normal to the Earth-Moon plane
    for(int i = 0; i <3; i++) e3[i] = C[i][2];
    //rp = r - n*(r.n)
    for(int i = 0; i <3; i++) rp[i] = r[i] - vdot_c(r, e3)*e3[i];
    //gammap = angle between rp and ME
    double gammap = vsep_c(rp, me)*180/M_PI;

    //Angle between BemS and ME
    cout << "Separation angle between BemS and ME (QBCP)= " << gamma  << endl;
    cout << "Separation angle between BemS and ME (JPL) = " << gammap << endl;


    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
}

/**
 *  \brief Find the epoch that gives the best correspondance between the Sun-Earth-Moon QBCP configuration at time t, given in SEM/EM coordinates, and the Sun-Earth-Moon true motion
 *         given by the JPL ephemerides. The input tSYS must be in SEM/EM normalized coordinates. The output epoch et is given in seconds.
 **/
void qbcp2jpl_inertial(double tSYS, double *et, int coord_type)
{
    //====================================================================================
    //Precheck on the coord_type
    //====================================================================================
    int coord_eph = eph_coord(coord_type);
    if(coord_eph == GSL_FAILURE)
    {
        cout << "qbcp2jpl. There is no ephemerides coordinates" << endl;
        cout << " associated to the current coordinate system. return." << endl;
        return;
    }

    //====================================================================================
    //Compute the positions of primaries in ECISEM coordinates
    //====================================================================================
    double Pe[3], Pm[3], Ps[3];
    //double Pnorm = 0.0;
    switch(coord_eph)
    {
        case VEM:
            cout << "VEM cannot be used here!" << endl;
            break;
        case VSEM:
            double n  = SEML.us_sem.n;
            double ni = SEML.us_sem.ni;
            double ns = SEML.us_sem.ns;
            double ai = SEML.us_sem.ai;
            double as = SEML.us_sem.as;
            double me = SEML.us_sem.me;
            double mm = SEML.us_sem.mm;

            //----------------------------------------------------------------------------
            //r & R
            //----------------------------------------------------------------------------
            //R
            double R1 = creal(evz(SEML.cs_sem.Zt, tSYS, n, ns, as));
            double R2 = cimag(evz(SEML.cs_sem.Zt, tSYS, n, ns, as));
            //r
            double r1 = creal(evz(SEML.cs_sem.zt, tSYS, n, ni, ai));
            double r2 = cimag(evz(SEML.cs_sem.zt, tSYS, n, ni, ai));

            //----------------------------------------------------------------------------
            //Position of the primaries
            //----------------------------------------------------------------------------
            Ps[0] = R1;
            Ps[1] = R2;
            Ps[2] = 0.0;

            Pm[0] = - me/(mm + me)* r1;
            Pm[1] = - me/(mm + me)* r2;
            Pm[2] = 0.0;

            Pe[0] = + mm/(mm + me)* r1;
            Pe[1] = + mm/(mm + me)* r2;
            Pe[2] = 0.0;

            //----------------------------------------------------------------------------
            //Normalization
            //----------------------------------------------------------------------------
            //            Pnorm = vnorm_c(Pm);
            //            for(int i = 0; i < 3; i++) Pm[i] /= Pnorm;
            //
            //            Pnorm = vnorm_c(Ps);
            //            for(int i = 0; i < 3; i++) Ps[i] /= Pnorm;
            //
            //            Pnorm = vnorm_c(Pe);
            //            for(int i = 0; i < 3; i++) Pe[i] /= Pnorm;

        break;
    }


    //====================================================================================
    // Init
    //====================================================================================
    //============================================
    // Starting and end time
    //============================================
    //Covered time span in 432s.bsp is:
    //Start:  1949 DEC 14 00:00:00.000 (TDB)
    //Stop:   2050 JAN 02 00:00:00.000 (TDB)
    string epoch_start, epoch_end;
    SpiceDouble  et_start, et_end;

    epoch_start = "2047 JAN 01 00:00:00.000"; //We start in 2047
    epoch_end   = "2048 JAN 01 00:00:00.000"; //We end in 2048

    //In ephemeris time
    str2et_c(epoch_start.c_str(), &et_start);
    str2et_c(epoch_end.c_str(), &et_end);

    //============================================
    // Initialize the change of coordinates:
    // Ecliptic coordinates  <-> SEM/EM synodic in position.
    //============================================
    double B[3];
    double C[3][3], k;
    SpiceDouble lt, RE[6], RM[6], RS[6], re[3], rm[3], rs[3], REM[6];

    //============================================
    //Initialize the minimum distance
    //============================================
    double minDist    = 0.0;
    double argMinDist = 0.0;
    double dist       = 0.0;

    //============================================
    //Frequence of check
    //============================================
    int freq = 60;//spd_c();

    //====================================================================================
    // Loop in the Kernel, each day
    //====================================================================================
    *et = et_start;
    double dSEM= 0.0;
    double ddist = 0.0;
    do
    {
        //Initialize the change of coordinates
        init_ecl2synpos(*et, B, C, &k, coord_eph);

        //Compute the ecliptic positions of the Earth & the Moon
        spkezr_c ("EARTH", *et, DEFFRAME,  "NONE",  DEFOBS, RE, &lt);
        spkezr_c ("MOON",  *et, DEFFRAME,  "NONE",  DEFOBS, RM, &lt);
        spkezr_c ("SUN",   *et, DEFFRAME,  "NONE",  DEFOBS, RS, &lt);
        spkezr_c ("EARTH MOON BARYCENTER", *et, DEFFRAME,  "NONE",  DEFOBS, REM, &lt);

        //Back in ECI coordinates
        //        ecl2neci(RE, REM, re, SEML.ss);
        //        ecl2neci(RM, REM, rm, SEML.ss);
        //        ecl2neci(RS, REM, rs, SEML.ss);

        //Instead, we compute the Sun-Bem distance each time
        dSEM= 0.0;
        for(int i = 0; i < 3; i++) dSEM += (RS[i] - REM[i])*(RS[i] - REM[i]);
        dSEM = sqrt(dSEM);

        for(int i = 0; i < 3; i++)
        {
            re[i] = (RE[i] - REM[i])/dSEM;
            rm[i] = (RM[i] - REM[i])/dSEM;
            rs[i] = (RS[i] - REM[i])/dSEM;
        }

        //Normalisation
        //        Pnorm = vnorm_c(re);
        //        for(int i = 0; i < 3; i++) re[i] /= Pnorm;
        //        Pnorm = vnorm_c(rs);
        //        for(int i = 0; i < 3; i++) rs[i] /= Pnorm;
        //        Pnorm = vnorm_c(rm);
        //        for(int i = 0; i < 3; i++) rm[i] /= Pnorm;

        //Compute the minimum distance
        dist = 0.0;
        ddist = 0.0;
        for(int i = 0; i < 3; i++) ddist += (Pe[i] - re[i])*(Pe[i] - re[i]);
        dist += sqrt(ddist);
        ddist = 0.0;
        for(int i = 0; i < 3; i++) ddist += (Pm[i] - rm[i])*(Pm[i] - rm[i]);
        dist += sqrt(ddist);
        ddist = 0.0;
        for(int i = 0; i < 3; i++) ddist += (Ps[i] - rs[i])*(Ps[i] - rs[i]);
        dist += sqrt(ddist);

        if(*et == et_start)
        {
            minDist    = dist;
            argMinDist = *et;
        }
        else
        {
            if(dist < minDist)
            {
                minDist    = dist;
                argMinDist = *et;
            }
        }

        //Add one freq
        *et += freq;//spd_c();
    }while(*et <= et_end);


    //============================================
    // Update results
    //============================================
    *et = argMinDist;
}



//========================================================================================
//                         Name of the primaries
//========================================================================================
/**
 *      \brief Return the name of the first primary associated to the coord_type.
 *             Not included (since it makes no sense):
 *             case INEM   case INSEM  case VECLI
 **/
int m1name(int coord_type)
{
    switch(coord_type)
    {
    case NCSEM  :
    case VNCSEM :
    case PSEM   :
    case VSEM   :
    case ECISEM :
            return SUN;
    case NCEM   :
    case VNCEM  :
    case PEM    :
    case VEM    :
            return EARTH;
    default:
        cout << "m1name. Unknown coord_type. SUN is used by default." << endl;
        return SUN;
    }
}

/**
 *      \brief Return the name of the second primary associated to the coord_type.
 *             Not included (since it makes no sense):
 *             case INEM   case INSEM  case VECLI
 **/
int m2name(int coord_type)
{
    switch(coord_type)
    {
    case NCSEM  :
    case VNCSEM :
    case PSEM   :
    case VSEM   :
    case ECISEM :
            return EARTH_MOON_BARYCENTER;
    case NCEM   :
    case VNCEM  :
    case PEM    :
    case VEM    :
            return MOON;
    default:
        cout << "m2name. Unknown coord_type. EARTH_MOON_BARYCENTER is used by default." << endl;
        return EARTH_MOON_BARYCENTER;
    }
}

/**
 *      \brief Return the name of the second primary associated to the coord_type.
 *             Not included (since it makes no sense):
 *             case INEM   case INSEM  case VECLI
 **/
int m2name2(int coord_type)
{
    switch(coord_type)
    {
    case NCSEM  :
    case VNCSEM :
    case PSEM   :
    case VSEM   :
            return EARTH;
    case NCEM   :
    case VNCEM  :
    case PEM    :
    case VEM    :
            return MOON;
    default:
        cout << "m2name2. Unknown coord_type. EARTH is used by default." << endl;
        return EARTH;
    }
}

//========================================================================================
//                          Change of coordinates:
//              Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic
//                            Only in position
//========================================================================================
/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + kC r        and    r = 1/k inv(C) (R - B)
 *
 *        The time is given as a string, recognized by SPICE
 **/
void init_ecl2synpos(string epoch, double B[3], double C[3][3], double *k, int coord_eph)
{
    //============================================
    //To ephemeris time
    //============================================
    SpiceDouble  et;
    str2et_c (epoch.c_str(), &et);

    //============================================
    //Initialize
    //============================================
    init_ecl2synpos(et, B, C, k, coord_eph);
}

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + kC r        and    r = 1/k inv(C) (R - B)
 *
 *        The time is given as a double, in seconds (SPICE ephemeris time)s
 **/
void init_ecl2synpos(double et, double B[3], double C[3][3], double *k, int coord_eph)
{
    //----------------------------
    //Retrieve the positions of m1 and m2
    //----------------------------
    SpiceDouble lt;
    SpiceDouble RS[6], RB[6];

    spkez_c (m1name(coord_eph), et, DEFFRAME, "NONE", SSB, RS, &lt);
    spkez_c (m2name(coord_eph), et, DEFFRAME, "NONE", SSB, RB, &lt);

    //============================================
    //Build the necessary objects
    //for the COC
    //============================================
    //----------------------------
    //Init
    //----------------------------
    double RBS[3], dRBS[3];

    //Masses. Note that the system of units is not important here
    //since only mass ratios are used.
    double m1 = 1.0, m2 = 1.0;
    switch(coord_eph)
    {
        case VEM:
            m1 = SEML.us_em.me; //Earth mass
            m2 = SEML.us_em.mm; //Moon mass
            break;
        case VSEM:
            m1 = SEML.us_sem.ms; //Sun mass
            m2 = SEML.us_sem.me + SEML.us_sem.mm; //Earth+Moon masses
            break;
        default:
            cout << "init_ecl2synpos. Unknown coord. system." << endl;
    }


    //----------------------------
    //Position B of the barycenter
    //----------------------------
    for(int i = 0; i < 3; i++) B[i] = (m1*RS[i] + m2*RB[i])/(m1 + m2);

    //----------------------------
    //Building the orthogonal matrix C & k factor
    //----------------------------
    double e1[3], e2[3], e3[3];

    //Difference between the two position vectors
    for(int i = 0; i < 3; i++) RBS[i] = RS[i] - RB[i];

    //Derivative of the difference wrt to time
    for(int i = 0; i < 3; i++) dRBS[i] = RS[i+3] - RB[i+3];

    //Normalized RBS in e1, norm in k
    unorm_c(RBS, e1, k);
    //Normalized crossed vector in e3
    ucrss_c(RBS, dRBS, e3);
    //Normalized crossed vector in e2
    ucrss_c(e3, e1, e2);

    //Building C
    for(int i = 0; i <3; i++) C[i][0] = e1[i];
    for(int i = 0; i <3; i++) C[i][1] = e2[i];
    for(int i = 0; i <3; i++) C[i][2] = e3[i];

    //----------------------------
    //Test
    //----------------------------
    //    double vout[3], vout2[3];
    //    cout << "Sun position in ecliptic coordinates:" << endl;
    //    vector_printf_prec(RS, 3);
    //    ecl2synpos(RS, vout, B, C, *k);
    //    cout << "Sun position in synodical coordinates:" << endl;
    //    vector_printf_prec(vout, 3);
    //    syn2eclpos(vout, vout2, B, C, *k);
    //    cout << "Sun position back in ecliptic coordinates:" << endl;
    //    vector_printf_prec(vout2, 3);
    //
    //    double err;
    //    int dim = 1;
    //    double dXS = dfridr(xBSecl, et, dim, 1.0, &err);
    //
    //    cout << "Estimated RS[dim] - RB[dim] = " << dXS << ", with error " << err << endl;
    //    cout << "Real value = " << RS[dim+3] - RB[dim+3] << endl;

}

/**
 * \brief Change of coordinates: Ecliptic coordinates  <- Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *      R = B + kC r
 **/
void syn2eclpos(double vin[3], double vout[3], double B[3], double C[3][3], double k)
{
    for(int i = 0; i < 3; i++)
    {
        vout[i] = 0.0;
        for(int j = 0; j < 3; j++) vout[i] += k*C[i][j]*vin[j];
        vout[i] += B[i];
    }
}

/**
 * \brief Change of coordinates: Ecliptic coordinates  -> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *                  r = 1/k inv(C) (R - B)
 **/
void ecl2synpos(double vin[3], double vout[3], double B[3], double C[3][3], double k)
{
    //============================================
    //K1 =  vin - B
    //============================================
    gsl_vector *K1 = gsl_vector_alloc(3);
    gsl_vector *K2 = gsl_vector_alloc(3);
    for(int i = 0; i < 3; i++) gsl_vector_set(K1, i, vin[i] - B[i]);

    //============================================
    //Inverse the system  Cinv*K2 = K1 with LU decomposition
    //============================================
    gsl_matrix *Cinv = gsl_matrix_alloc(3,3);
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) gsl_matrix_set(Cinv, i, j, C[i][j]);
    int s;
    gsl_permutation * p = gsl_permutation_alloc(3);
    gsl_linalg_LU_decomp (Cinv, p, &s);
    gsl_linalg_LU_solve(Cinv, p, K1, K2);


    //============================================
    // Normalization
    //============================================
    for(int i = 0; i < 3; i++) vout[i] = gsl_vector_get(K2, i)/k;

    //============================================
    // Free
    //============================================
    gsl_vector_free(K1);
    gsl_vector_free(K2);
    gsl_matrix_free(Cinv);
    gsl_permutation_free (p);
}


/**
 * \brief Change of coordinates: Ecliptic coordinates  -> Sun-(Earth+Moon) synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relation:
 *
 *                  r = inv(C) (R - B)
 **/
void ecl2syndpos(double vin[3], double vout[3], double B[3], double C[3][3])
{
    //============================================
    //K1 =  vin - B
    //============================================
    gsl_vector *K1 = gsl_vector_alloc(3);
    gsl_vector *K2 = gsl_vector_alloc(3);
    for(int i = 0; i < 3; i++) gsl_vector_set(K1, i, vin[i] - B[i]);

    //============================================
    //Inverse the system  Cinv*K2 = K1 with LU decomposition
    //============================================
    gsl_matrix *Cinv = gsl_matrix_alloc(3,3);
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) gsl_matrix_set(Cinv, i, j, C[i][j]);
    int s;
    gsl_permutation * p = gsl_permutation_alloc(3);
    gsl_linalg_LU_decomp (Cinv, p, &s);
    gsl_linalg_LU_solve(Cinv, p, K1, K2);

    //============================================
    // Storage
    //============================================
    for(int i = 0; i < 3; i++) vout[i] = gsl_vector_get(K2, i);

    //============================================
    // Free
    //============================================
    gsl_vector_free(K1);
    gsl_vector_free(K2);
    gsl_matrix_free(Cinv);
    gsl_permutation_free (p);
}

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) dimensionalized synodic in position.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        R = B + C r        and    r = inv(C) (R - B)
 *
 *        The time is given as a double, in seconds (SPICE ephemeris time)s
 **/
void init_ecl2syndpos(double et, double B[3], double C[3][3], int coord_eph)
{
    //----------------------------
    //Retrieve the positions of m1 and m2
    //----------------------------
    SpiceDouble lt;
    SpiceDouble R1[6], R2[6];

    spkez_c (m1name(coord_eph),  et, DEFFRAME, "NONE", SSB, R1, &lt);
    spkez_c (m2name2(coord_eph), et, DEFFRAME, "NONE", SSB, R2, &lt);

    //============================================
    //Build the necessary objects
    //for the COC
    //============================================
    double R12[3], dR12[3];

    //----------------------------
    //Position B is the position of the biggest primary
    //----------------------------
    for(int i = 0; i < 3; i++) B[i] = R1[i];

    //----------------------------
    //Building the orthogonal matrix C
    //----------------------------
    double e1[3], e2[3], e3[3];

    //Difference between the two position vectors
    for(int i = 0; i < 3; i++) R12[i] = - R1[i] + R2[i];

    //Derivative of the difference wrt to time
    for(int i = 0; i < 3; i++) dR12[i] = - R1[i+3] + R2[i+3];

    //Normalized RBS in e1, norm in k
    double k;
    unorm_c(R12, e1, &k);
    //Normalized crossed vector in e3
    ucrss_c(R12, dR12, e3);
    //Normalized crossed vector in e2
    ucrss_c(e3, e1, e2);

    //Building C
    for(int i = 0; i <3; i++) C[i][0] = e1[i];
    for(int i = 0; i <3; i++) C[i][1] = e2[i];
    for(int i = 0; i <3; i++) C[i][2] = e3[i];
}

//========================================================================================
//                          Building the change of coordinates
//              between the ecliptic coordinates and the synodic ones
//              mainly from Dei Tos 2014 and Gomez et al. 2002.
//========================================================================================
//----------------------------------------------------------------------------------------
// Necessary parameters
// From Dei Tos (2014)
//----------------------------------------------------------------------------------------
/**
 * \brief Init k, k', and k" coefficients, with k = ||r21||, r21
 *        being the relative position of the primaries.
 *        See equations (3.4) of Dei Tos (2014)
 **/
void init_k(const double r21[3], const double v21[3], const double a21[3],
            double *k, double *kp, double *kpp)
{
    //k = ||r21|
    *k = vnorm_c(r21);
    //kp = r21 . v21 / k
    *kp = vdot_c(r21, v21)/ (*k);
    //kpp = (v21 . v21 + r21.a21 - kp^2)/k
    *kpp = (vdot_c(v21, v21) + vdot_c(r21, a21) - (*kp)* (*kp))/(*k);
}

/**
 *  \brief Init the cross products r21 x v21, r21 x a21, etc.
 **/
void init_cp(const double r21[3], const double v21[3],
             const double a21[3], const double j21[3],
             double r21_v21[3], double r21_a21[3],
             double r21_j21[3], double v21_a21[3])
{
    vcrss_c(r21, v21, r21_v21);
    vcrss_c(r21, a21, r21_a21);
    vcrss_c(r21, j21, r21_j21);
    vcrss_c(v21, a21, v21_a21);
}

/**
 * \brief Init h, h', and h" coefficients, with h = ||r21 x v21||, r21/v21
 *        being the relative position/velocity of the primaries.
 *        See equations (3.5) of Dei Tos (2014)
 **/
void init_h(const double r21_v21[3], const double r21_a21[3],
            const double r21_j21[3], const double v21_a21[3],
            double *h, double *hp, double *hpp, double etemp[3])
{
    //h = ||r21 x v21||
    *h = vnorm_c(r21_v21);

    //hp = (r21 x v21) . (r21 x a21) / h
    *hp = vdot_c(r21_v21, r21_a21)/(*h);

    //et = v21 x a21 + r21 x j21
    vadd_c(v21_a21, r21_j21, etemp);

    //hpp = hp'
    *hpp  = + ( vdot_c(r21_a21, r21_a21) + vdot_c(r21_v21, etemp)  - (*hp)*(*hp) )/(*h);
}

/**
 *  \brief Init the orthonormal basis (e1, e2, e3)
 **/
void init_ei(const double r21[3], const double v21[3],
             double e1[3],  double e2[3], double e3[3])
{
    double k;
    //First, k = ||r21||, e1 = (r21)/k
    unorm_c(r21, e1, &k);
    //e3 = (r21 x v21)/||r21 x v21||
    ucrss_c(r21, v21, e3);
    //e2 = e3 x e1
    ucrss_c(e3, e1, e2);
}

/**
 *  \brief Init the derivative of the orthonormal basis (e1', e2', e3')
 **/
void init_dei(const double r21[3], const double v21[3], const double a21[3],
              const double k,      const double kp,
              const double h,      const double hp,
              const double  e1[3], const double e3[3],
              double de1[3], double de2[3], double de3[3], double etemp[3])
{
    //de1 = (k * v21 - kp * r21)/k^2 = (v21 - kp * e1)/k
    for(int i = 0; i < 3; i++) de1[i] = (v21[i] - kp*e1[i])/k;

    //de3 = (r21 x a21 - hp * e3)/h
    vcrss_c(r21, a21, de3);
    for(int i = 0; i < 3; i++) de3[i] = (de3[i] - hp * e3[i])/h;

    //de2 = e3' x e1 + e3 x e1'
    vcrss_c(de3, e1,  de2); //e3' x e1
    vcrss_c(e3,  de1, etemp);  //e3 x e1'
    for(int i = 0; i < 3; i++) de2[i] = de2[i] + etemp[i];
}


/**
 *  \brief Init the double derivative of the orthonormal basis (e1", e2", e3")
 **/
void init_ddei(const double r21[3], const double v21[3], const double a21[3],
               const double r21_v21[3], const double r21_a21[3],
               const double r21_j21[3], const double v21_a21[3],
               const double k, const double kp, const double kpp,
               const double h, const double hp, const double hpp,
               const double  e1[3],   const double  e3[3],
               const double  de1[3],  const double  de3[3],
               double dde1[3], double dde2[3], double dde3[3],
               double  etemp[3])
{
    //dde1
    for(int i = 0; i < 3; i++)
    {
        dde1[i]  = +a21[i]/k;
        dde1[i] += -2*kp*v21[i]/(k*k);
        dde1[i] += +(2*kp*kp - k*kpp)*r21[i]/(k*k*k);
    }

    //dde3
    for(int i = 0; i < 3; i++)
    {
        dde3[i]  = +(v21_a21[i] + r21_j21[i])/h;
        dde3[i] += -2*hp*r21_a21[i]/(h*h);
        dde3[i] += +(2*hp*hp - h*hpp)*r21_v21[i]/(h*h*h);
    }

    //dde2 = e3" x e1 + 2 e3' x e1' + e3 x e1"
    vcrss_c(dde3, e1, dde2);    //e3" x e1
    vcrss_c(de3, de1, etemp);   //e3'x e1'
    for(int i = 0; i < 3; i++) dde2[i] += 2*etemp[i];
    vcrss_c(e3, dde1, etemp);   //e3x e1"
    for(int i = 0; i < 3; i++) dde2[i] += etemp[i];
}

/**
 *  \brief Init the derivative of the orthonormal basis times k:
 *           (kC)' = ((ke1)', (ke2)', (ke3)')
 **/
void init_kdei(const double  e1[3], const double  e2[3], const double  e3[3],
               const double de1[3], const double de2[3], const double de3[3],
               const double k,      const double kp,
               double dke1[3], double dke2[3], double dke3[3])
{
    for(int i = 0; i < 3; i++)
    {
        dke1[i] = k*de1[i] + kp*e1[i];
        dke2[i] = k*de2[i] + kp*e2[i];
        dke3[i] = k*de3[i] + kp*e3[i];
    }
}

/**
 *  \brief Init the matrices C, kCprim and the scalar k necessary to compute
 *         the cocs: ecliptic -> synodic and synodic -> ecliptic.
 **/
void init_coc(const double r21[3],  const double v21[3],
              const double a21[3],  const double j21[3],
              double C[3][3], double *k, double kCprim[3][3])
{
    //------------------------------------------------------------------------------------
    //Orthonormal basis (e1, e2, e3) + h and k factors and their derivatives
    //------------------------------------------------------------------------------------
    double e1[3], e2[3], e3[3], et[3], kp, kpp, h, hp, hpp;
    double r21_v21[3], r21_a21[3], r21_j21[3], v21_a21[3];
    //Building the cross products
    init_cp(r21, v21, a21, j21, r21_v21, r21_a21, r21_j21, v21_a21);
    //Building h, h', h"
    init_h(r21_v21, r21_a21, r21_j21, v21_a21, &h, &hp, &hpp, et);
    //Building k, k', k''
    init_k(r21, v21, a21, k, &kp, &kpp);
    //Building e1, e2, e3
    init_ei(r21, v21, e1, e2, e3);

    //------------------------------------------------------------------------------------
    //Orthonormal basis derivatives (e1', e2', e3')
    //------------------------------------------------------------------------------------
    double de1[3], de2[3], de3[3];
    init_dei(r21, v21, a21, *k, kp, h, hp, e1, e3, de1, de2, de3, et);

    //------------------------------------------------------------------------------------
    //Orthonormal basis derivatives (kC)' = ((ke1)', (ke2)', (ke3)')
    //------------------------------------------------------------------------------------
    double dke1[3], dke2[3], dke3[3];
    init_kdei(e1, e2, e3, de1, de2, de3, *k, kp, dke1, dke2, dke3);

    //------------------------------------------------------------------------------------
    //Building C
    //------------------------------------------------------------------------------------
    for(int i = 0; i <3; i++) C[i][0] = e1[i];
    for(int i = 0; i <3; i++) C[i][1] = e2[i];
    for(int i = 0; i <3; i++) C[i][2] = e3[i];

    //------------------------------------------------------------------------------------
    //Building (kC)'
    //------------------------------------------------------------------------------------
    for(int i = 0; i <3; i++) kCprim[i][0] = dke1[i];
    for(int i = 0; i <3; i++) kCprim[i][1] = dke2[i];
    for(int i = 0; i <3; i++) kCprim[i][2] = dke3[i];



}



//========================================================================================
//                          Change of coordinates:
//              Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic
//                            Position + Velocity
//========================================================================================
//----------------------------------------------------------------------------------------
// Init
//----------------------------------------------------------------------------------------
/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate(string epoch, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph)
{
    //=================================================================================
    //To ephemeris time
    //=================================================================================
    SpiceDouble  et;
    str2et_c (epoch.c_str(), &et);

    //=================================================================================
    //Initialize
    //=================================================================================
    init_ecl2synstate(et, B, C, k, Bprim, kCprim, coord_eph);
}

/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate(double et, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph)
{
    //=================================================================================
    //Build the necessary objects for the COC
    //=================================================================================
    //------------------------------------------------------------------------------------
    //Position, velocity, acceleration, and jerk of the primaries
    //------------------------------------------------------------------------------------
    double lt;
    double R1[6];  //state of m1
    double R2[6];  //state of m2
    //double A1[3];  //acceleration of m1
    //double A2[3];  //acceleration of m2
    double J1[3];  //jerk of m1: set to ZERO because not needed here
    double J2[3];  //jerk of m2: set to ZERO because not needed here

    //--------------------------------------------
    //Position, velocity of the primaries, in dimensionnalized units
    //--------------------------------------------
    spkez_c (m1name(coord_eph), et, DEFFRAME, "NONE", SSB, R1, &lt);
    spkez_c (m2name(coord_eph), et, DEFFRAME, "NONE", SSB, R2, &lt);

    //--------------------------------------------
    //Acceleration of the primaries, in dimensionnalized units
    //--------------------------------------------
    //First possibility: from ephemerides
    double err;
    //for(int i = 0; i < 3; i++)  A1[i] = dfridr(xm1ecl, et, i+3, 1.0, &err, coord_eph);
    //for(int i = 0; i < 3; i++)  A2[i] = dfridr(xm2ecl, et, i+3, 1.0, &err, coord_eph);

    //--------------------------------------------
    //Jerk of the primaries, set to ZERO because not needed here
    //--------------------------------------------
    for(int i = 0; i < 3; i++) J1[i] = 0.0;
    for(int i = 0; i < 3; i++) J2[i] = 0.0;

    //----------------------------
    //Init
    //----------------------------
    //Masses. Note that the system of units is not important here
    //since only mass ratios are used.
    double m1 = 1.0, m2 = 1.0;
    switch(coord_eph)
    {
        case VEM:
            m1 = SEML.us_em.me; //Earth mass
            m2 = SEML.us_em.mm; //Moon mass
            break;
        case VSEM:
            m1 = SEML.us_sem.ms; //Sun mass
            m2 = SEML.us_sem.me + SEML.us_sem.mm; //Earth+Moon masses
            break;
        default:
            cout << "init_ecl2synstate. Unknown coord. system." << endl;
    }

    //====================================================================================
    // Barycenter B, B'
    //====================================================================================
    //Position B of the barycenter
    for(int i = 0; i < 3; i++) B[i] = (m1*R1[i] + m2*R2[i])/(m1 + m2);

    //Position Bprim of the barycenter = B'
    for(int i = 0; i < 3; i++) Bprim[i] = (m1*R1[i+3] + m2*R2[i+3])/(m1 + m2);

    //------------------------------------------------------------------------------------
    //Relative position, velocity, acceleration, and jerk
    //------------------------------------------------------------------------------------
    double r21[3]; //relative position
    double v21[3]; //relative velocity
    double a21[3]; //relative acceleration
    double j21[3]; //relative jerk
    for(int i = 0; i < 3; i++)
    {
        //Difference between the two position vectors
        r21[i]  = R1[i] - R2[i];
        //Derivative of the difference wrt to dimentionalized time
        v21[i]  = R1[i+3] - R2[i+3];
        //Difference between the two accelerations
        a21[i]  = dfridr(xBSecl, et, i+3, 1.0, &err, coord_eph);//A1[i] - A2[i];
        //Difference between the two jerks
        j21[i]  = J1[i] - J2[i];
    }

    //------------------------------------------------------------------------------------
    //Orthonormal basis (e1, e2, e3) + h and k factors and their derivatives
    //------------------------------------------------------------------------------------
    init_coc(r21, v21, a21, j21, C, k, kCprim);
}


/**
 * \brief Initialize the change of coordinates Ecliptic coordinates  <-> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r                      and    r  = 1/k inv(C) (R - B)
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'       and    r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot   and    rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 **/
void init_ecl2synstate_2(double et, double B[3], double C[3][3], double *k, double Bprim[3], double kCprim[3][3], int coord_eph)
{
    //=================================================================================
    //Build the necessary objects for the COC
    //=================================================================================
    //============================================
    //Retrieve the positions of the Sun and B
    //============================================
    SpiceDouble lt;
    SpiceDouble RS[6], RB[6];
    spkez_c (m1name(coord_eph), et, DEFFRAME, "NONE", SSB, RS, &lt);
    spkez_c (m2name(coord_eph), et, DEFFRAME, "NONE", SSB, RB, &lt);

    //----------------------------
    //Init
    //----------------------------
    double RBS[3], dRBS[3], ddRBS[3];

    //Masses. Note that the system of units is not important here
    //since only mass ratios are used.
    double m1 = 1.0, m2 = 1.0;
    switch(coord_eph)
    {
        case VEM:
            m1 = SEML.us_em.me; //Earth mass
            m2 = SEML.us_em.mm; //Moon mass
            break;
        case VSEM:
            m1 = SEML.us_sem.ms; //Sun mass
            m2 = SEML.us_sem.me + SEML.us_sem.mm; //Earth+Moon masses
            break;
        default:
            cout << "init_ecl2synstate. Unknown coord. system." << endl;
    }


    //----------------------------
    //Position B of the barycenter
    //----------------------------
    for(int i = 0; i < 3; i++) B[i] = (m1*RS[i] + m2*RB[i])/(m1 + m2);

    //============================================
    //Building the orthogonal matrix C & k factor
    //============================================
    double e1[3], e2[3], e3[3];

    //Difference between the two position vectors
    for(int i = 0; i < 3; i++) RBS[i]  = RS[i] - RB[i];

    //Derivative of the difference wrt to time
    for(int i = 0; i < 3; i++) dRBS[i] = RS[i+3] - RB[i+3];


    //Normalized RBS in e1, norm in k
    unorm_c(RBS, e1, k);
    //Normalized crossed vector in e3
    ucrss_c(RBS, dRBS, e3);
    //Normalized crossed vector in e2
    ucrss_c(e3, e1, e2);

    //Building C
    for(int i = 0; i <3; i++) C[i][0] = e1[i];
    for(int i = 0; i <3; i++) C[i][1] = e2[i];
    for(int i = 0; i <3; i++) C[i][2] = e3[i];

    //=================================================================================
    //Build the derivatives of the necessary objects for the COC
    //=================================================================================

    //============================================
    //Position Bprim of the barycenter = B'
    //============================================
    for(int i = 0; i < 3; i++) Bprim[i] = (m1*RS[i+3] + m2*RB[i+3])/(m1 + m2);

    //============================================
    //Building the orthogonal matrix kCprim = (k C)'
    //============================================
    double dke1[3], dke2[3], dke3[3], dk, de1[3], de2[3], de3[3], A[3], An, err;

    //dk = k' = 1/k * RBS.RBS'
    dk = 0.0;
    for(int i = 0; i < 3; i++) dk += RBS[i]*dRBS[i]/(*k);

    //----------------------------
    //dke1 = RBS'
    //----------------------------
    for(int i = 0; i < 3; i++) dke1[i] = dRBS[i];


    //----------------------------
    //dke3 = k'*e3 + k * e3'
    //----------------------------
    //Double derivatives of RBS = RBS", using Ridders’ method of polynomial extrapolation;
    for(int i = 0; i < 3; i++)  ddRBS[i] = dfridr(xBSecl, et, i+3, 1.0, &err, coord_eph);

    cout << "dRBS = " << endl;
    vector_printf_prec(dRBS, 3);
    cout << "------------------" << endl;

    cout << "ddRBS = " << endl;
    vector_printf_prec(ddRBS, 3);
    cout << "------------------" << endl;

    //Cross products
    double  RBS_dRBS[3];
    double  dRBS_dRBS[3];
    double  RBS_ddRBS[3];
    vcrss_c(RBS,  dRBS,  RBS_dRBS);  //RBS  x  RBS'
    vcrss_c(dRBS, dRBS,  dRBS_dRBS); //RBS' x  RBS'
    vcrss_c(RBS,  ddRBS, RBS_ddRBS); //RBS  x  RBS"

    //Magnitudes of RBS  x  RBS'
    double RBS_dRBS_n = vnorm_c(RBS_dRBS);

    //A = RBS' x RBS' + RBS x RBS"
    for(int i = 0; i < 3; i++) A[i] = dRBS_dRBS[i] + RBS_ddRBS[i];

    //An = (RBS x  RBS') . A
    An = 0.0;
    for(int i = 0; i < 3; i++) An += RBS_dRBS[i] * A[i];

    //de3 = 1/|RBS  x  RBS'| * A -  1/|RBS  x  RBS'|^3 * An * RBS  x  RBS'
    for(int i = 0; i < 3; i++)
    {
        de3[i]  = A[i]/RBS_dRBS_n;
        de3[i] -= An/pow(RBS_dRBS_n,3.0)*RBS_dRBS[i];
    }

    //dke3 = k'*e3 + k * e3'
    for(int i = 0; i < 3; i++) dke3[i] = dk * e3[i] + (*k) * de3[i];

    //----------------------------
    //dke2 = k'*e2 + k * e2'
    //----------------------------
    //de1 = 1/|RBS| * RBS' -  1/|RBS|^3 * RBS * (RBS .  RBS')
    for(int i = 0; i < 3; i++)
    {
        de1[i]  = dRBS[i]/(*k);
        de1[i] -= (RBS[0]*dRBS[0] + (RBS[1]*dRBS[1]) + (RBS[2]*dRBS[2]))/pow((*k),3.0)*RBS[i];
    }

    //Cross products
    double  de3_e1[3];
    double  e3_de1[3];
    vcrss_c(de3, e1,  de3_e1); //e3' x e1
    vcrss_c(e3,  de1, e3_de1); //e3 x e1'

    //de2 = e3' x e1 + e3 x e1'
    for(int i = 0; i < 3; i++) de2[i] = de3_e1[i] + e3_de1[i];

    //dke2 = k'*e2 + k * e2'
    for(int i = 0; i < 3; i++) dke2[i] = dk * e2[i] + (*k) * de2[i];

    //Building C
    for(int i = 0; i <3; i++) kCprim[i][0] = dke1[i];
    for(int i = 0; i <3; i++) kCprim[i][1] = dke2[i];
    for(int i = 0; i <3; i++) kCprim[i][2] = dke3[i];


    cout << "C = " << endl;
    vector_printf_prec(e1, 3);
    vector_printf_prec(e2, 3);
    vector_printf_prec(e3, 3);

    cout << "de1 = " << endl;
    vector_printf_prec(de1, 3);
    cout << "de2 = " << endl;
    vector_printf_prec(de2, 3);
    cout << "de3 = " << endl;
    vector_printf_prec(de3, 3);

    cout << "(kC)' = " << endl;
    vector_printf_prec(dke1, 3);
    vector_printf_prec(dke2, 3);
    vector_printf_prec(dke3, 3);
}


//----------------------------------------------------------------------------------------
// Ecliptic coordinates -> Sun-(Earth+Moon) synodic
//----------------------------------------------------------------------------------------
/**
 * \brief Change of coordinates Ecliptic coordinates  -> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  r  = 1/k inv(C) (R - B)
 *        - For the velocities: r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 *
 *         Note: n_sem is the Sun-Bem mean motion, in rad/s.
 **/
void ecl2synstate(double vin[6], double vout[6], double et, int coord_eph)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3];
    init_ecl2synstate(et, B, C, &k, Bprim, kCprim, coord_eph);

    //=================================================================================
    //COC
    //=================================================================================
    ecl2synstate(vin, vout, B, C, k, Bprim, kCprim, mean_motion(coord_eph));
}

/**
 * \brief Change of coordinates Ecliptic coordinates  -> Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  r  = 1/k inv(C) (R - B)
 *        - For the velocities: r' = 1/k inv(C) (R' - B' - (kC)'r)
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: rdot =  = 1/(nk) inv(C) (R' - B' - (kC)'r)
 *
 *         Note: n_sem is the Sun-Bem mean motion, in rad/s.
 **/
void ecl2synstate(double vin[6], double vout[6], double B[3], double C[3][3], double k, double Bprim[3], double kCprim[3][3], double n_sys)
{
    //=================================================================================
    //For the position
    //=================================================================================
    ecl2synpos(vin, vout, B, C, k);

    //=================================================================================
    //For the velocity
    //=================================================================================
    gsl_vector *K1 = gsl_vector_alloc(3);
    gsl_vector *K2 = gsl_vector_alloc(3);
    gsl_matrix *Cinv = gsl_matrix_alloc(3,3);
    int s; gsl_permutation * p = gsl_permutation_alloc(3);

    //============================================
    //K1 =  vin' - B' - (kC)'r
    //============================================
    for(int i = 0; i < 3; i++) gsl_vector_set(K1, i, vin[i+3] - Bprim[i] - kCprim[i][0]*vout[0] - kCprim[i][1]*vout[1] - kCprim[i][2]*vout[2]);

    //============================================
    //Inverse the system  Cinv*K2 = K1 with LU decomposition
    //============================================
    for(int i = 0; i < 3; i++) for(int j = 0; j < 3; j++) gsl_matrix_set(Cinv, i, j, C[i][j]);
    gsl_linalg_LU_decomp (Cinv, p, &s);
    gsl_linalg_LU_solve(Cinv, p, K1, K2);


    //============================================
    // Normalization
    // CAREFUL: the time normalization is also performed here:
    //       t = nt*,
    // where t* is the time in seconds, and t is the normalized SEM/EM time
    //============================================
    for(int i = 0; i < 3; i++) vout[i+3] = gsl_vector_get(K2, i)/(k*n_sys);

    //============================================
    // Free
    //============================================
    gsl_vector_free(K1);
    gsl_vector_free(K2);
    gsl_matrix_free(Cinv);
    gsl_permutation_free (p);

}

//----------------------------------------------------------------------------------------
// Ecliptic coordinates  <- Sun-(Earth+Moon) synodic
//----------------------------------------------------------------------------------------
/**
 * \brief Change of coordinates Ecliptic coordinates  <- Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot
 **/
void syn2eclstate(double vin[6], double vout[6], double et, int coord_eph)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3];
    init_ecl2synstate(et, B, C, &k, Bprim, kCprim, coord_eph);

    //=================================================================================
    //COC
    //=================================================================================
    syn2eclstate(vin, vout, B, C, k, Bprim, kCprim, mean_motion(coord_eph));
}

/**
 * \brief Change of coordinates Ecliptic coordinates  <- Sun-(Earth+Moon) synodic coordinates in position+velocity state.
 *        For a position vector R in ecliptic coordinates, and a position vector r in Sun-(Earth+Moon) synodical coordinates, we have the following relations:
 *
 *        - For the positions:  R = B + kC r
 *        - For the velocities: R' = B' + (kC)' r + (kC) r'
 *
 *        However, the SEM time is normalized: t = nt*, where t is the normalized time and t* is the dimensionalized time (usually in seconds).
 *        For the SEM case, we have n ~ 1.9909866091e-7 rad/s.
 *
 *        Therefore, denoting rdot = dr/dt = 1/n r' the velocity in normalized units, we have:
 *
 *        - For the velocities: R' = B' + (kC)' r + n (kC) rdot
 **/
void syn2eclstate(double vin[6], double vout[6], double B[3], double C[3][3], double k, double Bprim[3], double kCprim[3][3], double n_sys)
{
    //=================================================================================
    //For the position
    //=================================================================================
    syn2eclpos(vin, vout, B, C, k);

    //=================================================================================
    //For the velocity
    // CAREFUL: the time normalization is also performed here:
    //       t = nt*,
    // where t* is the time in seconds, and t is the normalized (Sun-Earth+Moon) time
    //=================================================================================
    //vout' = B' + (kC)'r + kC r' = B' + (kC)'r + n kC rdot
    for(int i = 0; i < 3; i++)
    {
        vout[i+3]  = Bprim[i] + kCprim[i][0]*vin[0] + kCprim[i][1]*vin[1] + kCprim[i][2]*vin[2];
        vout[i+3] += n_sys*k*C[i][0]*vin[3] + n_sys*k*C[i][1]*vin[4] + n_sys*k*C[i][2]*vin[5];
    }

}


//----------------------------------------------------------------------------------------
// On vectors
//----------------------------------------------------------------------------------------
//---------------------------------
// Ecliptic coordinates <-> any other native type (NCEM, NCSEM...)
//---------------------------------
/**
 *  From Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void ecl2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //Ecliptic coordinates to VSEM/VEM coordinates
    //=================================================================================
    ecl2synstate_vec(yecl, etecl, ytemp, ttemp, N, et0, tsys0, coord_eph);

    //=================================================================================
    //VSEM/VEM coordinates to coord_type coordinates
    //=================================================================================
    qbcp_coc_vec(ytemp, ttemp, yout, tout, N, coord_eph, coord_type);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}

/**
 *  To Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2eclstate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //coord_type to VSEM/VEM coordinates
    //=================================================================================
    qbcp_coc_vec(yin, tin, ytemp, ttemp, N, coord_type, coord_eph);

    //=================================================================================
    //VSEM/VEM to ecliptic coordinates
    //=================================================================================
    syn2eclstate_vec(ytemp, ttemp, yecl, etecl, N, et0, tsys0, coord_eph);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}

//---------------------------------
// Ecliptic coordinates <-> any other synodical type (VEM, VSEM...)
//---------------------------------
/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void ecl2synstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3];
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //Init the COC
        init_ecl2synstate(etecl[p], B, C, &k, Bprim, kCprim, coord_eph);

        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yecl[i][p];

        //COC
        ecl2synstate(vin, vout, B, C, k, Bprim, kCprim, n_sys);

        //Update yout
        for(int i = 0; i <6; i++) yout[i][p] = vout[i];

        //Update tout (normalized time)
        tout[p] = (etecl[p] - et0)*n_sys + tsys0;
    }
}

/**
 * \brief From synodic to ecliptic coordinates, vector format
 **/
void syn2eclstate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, double et0, double tsys0, int coord_eph)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3];
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //Update etecl (non-normalized time)
        etecl[p] = et0 + (tin[p]-tsys0)/n_sys;

        //Init the COC
        init_ecl2synstate(etecl[p], B, C, &k, Bprim, kCprim, coord_eph);

        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yin[i][p];

        //COC
        syn2eclstate(vin, vout, B, C, k, Bprim, kCprim, n_sys);

        //Update yecl
        for(int i = 0; i <6; i++) yecl[i][p] = vout[i];
    }
}


//---------------------------------
// Earth-centered inertial (ECI) coordinates <-> any other native type (NCEM, NCSEM...)
//---------------------------------
/**
 *  From Normalized-Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void eci2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph)//, SS &ss)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //Ecliptic coordinates to VSEM/VEM coordinates
    //=================================================================================
    eci2synstate_vec(yecl, etecl, ytemp, ttemp, N, et0, tsys0, coord_eph);//, ss);

    //=================================================================================
    //VSEM/VEM coordinates to coord_type coordinates
    //=================================================================================
    qbcp_coc_vec(ytemp, ttemp, yout, tout, N, coord_eph, coord_type);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}

/**
 *  To Normalized-Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2ecistate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph)//, SS &ss)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //coord_type to VSEM/VEM coordinates
    //=================================================================================
    qbcp_coc_vec(yin, tin, ytemp, ttemp, N, coord_type, coord_eph);

    //=================================================================================
    //VSEM/VEM to ecliptic coordinates
    //=================================================================================
    syn2ecistate_vec(ytemp, ttemp, yecl, etecl, N, et0, tsys0, coord_eph);//, ss);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}

//---------------------------------
// Earth-centered inertial (ECI)  <-> any other synodical type (VEM, VSEM...)
//---------------------------------
/**
 * \brief From synodic to earth-centered normalized coordinates, vector format.
 *        Note that coord_eph = VEM/VSEM can be chosen independantly from the ss structure
 *        that normalizes the state.
 *        Note: for now the normalization is NOT taken into account because the stepper is going wrong when the state is normalized.
 **/
void syn2ecistate_vec(double **yin, double *tin, double **yeci, double *teci, int N, double et0, double tsys0, int coord_eph)//, SS &ss)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3], etecl;
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //----------------------------
        //Update etecl (non-normalized shifted time)
        //----------------------------
        etecl = et0 + (tin[p]-tsys0)/n_sys;

        //----------------------------
        // COC
        //----------------------------
        //Init the COC
        init_ecl2synstate(etecl, B, C, &k, Bprim, kCprim, coord_eph);
        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yin[i][p];
        //COC
        syn2eclstate(vin, vout, B, C, k, Bprim, kCprim, n_sys);

        //Position of the Earth
        double Re[6], lt;
        spkez_c (EARTH, etecl, DEFFRAME, "NONE", SSB, Re, &lt);

        //----------------------------
        //Store data
        //----------------------------
        //Normalization
        ecl2eci(vout, Re, vin);
        //Update yecl
        for(int i = 0; i <6; i++) yeci[i][p] = vin[i];
        //Update tecl, which is also the non-normalized shifted time
        teci[p] = etecl;
    }
}

/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void eci2synstate_vec(double **yeci, double *teci, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph)//, SS & ss)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3], etecl;
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //----------------------------
        //Update etecl (non-normalized shifted time):
        //etecl = tecln[p]*ss.n;
        //----------------------------
        etecl = teci[p];

        //Position of the Earth
        double Re[6], lt;
        spkez_c (EARTH, etecl, DEFFRAME, "NONE", SSB, Re, &lt);

        //----------------------------
        // COC
        //----------------------------
        //Init the COC
        init_ecl2synstate(etecl, B, C, &k, Bprim, kCprim, coord_eph);
        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yeci[i][p];
        //Denormalization
        eci2ecl(vin, Re, vout);
        //COC
        ecl2synstate(vout, vin, B, C, k, Bprim, kCprim, n_sys);

        //----------------------------
        //Store data
        //----------------------------
        //Update yout
        for(int i = 0; i <6; i++) yout[i][p] = vin[i];
        //Update tout (normalized time):
        tout[p] = (teci[p] - et0)*n_sys + tsys0;
    }
}


/**
 *  \brief From full ecliptic coordinates (directly from SPICE) to Earth-centered inertial coordinates.
 **/
void ecl2eci(double YECL[6], double YEARTH[6], double yeci[6])//, SS &ss)
{
    //Position
    for(int i = 0; i < 3; i++) yeci[i] = (YECL[i]-YEARTH[i]); // /ss.a;
    //Velocity
    for(int i = 3; i < 6; i++) yeci[i] = (YECL[i]-YEARTH[i]); // /(ss.a*ss.n);
}

/**
 *  \brief To full ecliptic coordinates (directly from SPICE) from to earth-centered inertial coordinates.
 **/
void eci2ecl(double yeci[6], double YEARTH[6], double YECL[6])//, SS &ss)
{
    //Position
    //for(int i = 0; i < 3; i++) YECL[i] = yecln[i]*ss.a + YEARTH[i];
    for(int i = 0; i < 3; i++) YECL[i] = yeci[i] + YEARTH[i];
    //Velocity
    //for(int i = 3; i < 6; i++) YECL[i] = yecln[i]*(ss.a*ss.n) + YEARTH[i];
    for(int i = 3; i < 6; i++) YECL[i] = yeci[i] + YEARTH[i];
}


/**
 * \brief From ecliptic to synodic coordinates, vector format.
 *        The final time is in Julian date.
 **/
void eci2syndpos_vec(double **yeci, double *teci, double **yout, double *tout, int N, int coord_eph)//, SS & ss)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], etecl;
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    for(int p = 0; p <= N; p++)
    {
        //----------------------------
        //Update etecl (non-normalized shifted time):
        //----------------------------
        etecl = teci[p];

        //Position of the Earth
        double Re[6], lt;
        spkez_c (EARTH, etecl, DEFFRAME, "NONE", SSB, Re, &lt);

        //----------------------------
        // COC
        //----------------------------
        //Init the COC
        init_ecl2syndpos(etecl, B, C, coord_eph);
        //Update vin
        for(int i = 0; i < 6; i++) vin[i] = yeci[i][p];
        //Denormalization
        eci2ecl(vin, Re, vout);
        //COC, ONLY IN POSITION FOR NOW
        ecl2syndpos(vout, vin, B, C);


        //----------------------------
        //Store data
        //----------------------------
        //Update yout
        for(int i = 0; i <6; i++) yout[i][p] = vin[i];
        //Update tout in Julian Date
        tout[p] = unitim_c(teci[p], "TDB", "JED");
    }
}


//----------------------------------------------------------------------------------------
// Normalized Earth-centered inertial (NECI) coordinates <-> any other native type (NCEM, NCSEM...)
//----------------------------------------------------------------------------------------
/**
 *  From Normalized-Ecliptic coordinates to a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void neci2coordstate_vec(double **yecl, double *etecl, double **yout, double *tout, int N, int coord_type, double et0,  double tsys0, int coord_eph, SS &ss)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //NECI coordinates to VSEM/VEM coordinates
    //=================================================================================
    neci2synstate_vec(yecl, etecl, ytemp, ttemp, N, et0, tsys0, coord_eph, ss);

    //=================================================================================
    //VSEM/VEM coordinates to coord_type coordinates
    //=================================================================================
    qbcp_coc_vec(ytemp, ttemp, yout, tout, N, coord_eph, coord_type);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}

/**
 *  To Normalized-Ecliptic coordinates from a generic other coordinate system. The initial time tsys0 is given in SEM/EM coordinates, and yields the correspondance
 *  between the ephemerides epoch and the normalized QBCP time.
 **/
void coord2necistate_vec(double **yin, double *tin, double **yecl, double *etecl, int N, int coord_type, double et0,  double tsys0, int coord_eph, SS &ss)
{
    //=================================================================================
    //Allocate memory
    //=================================================================================
    double **ytemp = dmatrix(0, 6, 0, N);
    double *ttemp  = dvector(0, N);

    //=================================================================================
    //coord_type to VSEM/VEM coordinates
    //=================================================================================
    qbcp_coc_vec(yin, tin, ytemp, ttemp, N, coord_type, coord_eph);

    //=================================================================================
    //VSEM/VEM to NECI coordinates
    //=================================================================================
    syn2necistate_vec(ytemp, ttemp, yecl, etecl, N, et0, tsys0, coord_eph, ss);

    //=================================================================================
    //Free memory
    //=================================================================================
    free_dmatrix(ytemp, 0, 6, 0, N);
    free_dvector(ttemp, 0, N);
}



//---------------------------------
// Normalized Earth-centered inertial (NECI)  <-> any other synodical type (VEM, VSEM...)
//---------------------------------
/**
 *  \brief From full ecliptic coordinates (directly from SPICE) to Earth-centered inertial coordinates.
 **/
void ecl2neci(double YECL[6], double YEARTH[6], double yneci[6], SS &ss)
{
    //Position
    for(int i = 0; i < 3; i++) yneci[i] = (YECL[i]-YEARTH[i])/ss.a;
    //Velocity
    for(int i = 3; i < 6; i++) yneci[i] = (YECL[i]-YEARTH[i])/(ss.a*ss.n);
}

/**
 *  \brief To full ecliptic coordinates (directly from SPICE) from to earth-centered inertial coordinates.
 **/
void neci2ecl(double yneci[6], double YEARTH[6], double YECL[6], SS &ss)
{
    //Position
    for(int i = 0; i < 3; i++) YECL[i] = yneci[i]*ss.a + YEARTH[i];
    //Velocity
    for(int i = 3; i < 6; i++) YECL[i] = yneci[i]*(ss.a*ss.n) + YEARTH[i];
}

/**
 * \brief From synodic to earth-centered normalized coordinates, vector format.
 *        Note that coord_eph = VEM/VSEM can be chosen independantly from the ss structure
 *        that normalizes the state.
 *        Note: for now the normalization is NOT taken into account because the stepper is going wrong when the state is normalized.
 **/
void syn2necistate_vec(double **yin, double *tin, double **yneci, double *tneci, int N, double et0, double tsys0, int coord_eph, SS &ss)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3], etecl;
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //----------------------------
        //Update etecl (non-normalized shifted time)
        //----------------------------
        etecl = et0 + (tin[p]-tsys0)/n_sys;

        //----------------------------
        // COC
        //----------------------------
        //Init the COC
        init_ecl2synstate(etecl, B, C, &k, Bprim, kCprim, coord_eph);
        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yin[i][p];
        //COC
        syn2eclstate(vin, vout, B, C, k, Bprim, kCprim, n_sys);

        //Position of the Center
        double Re[6], lt;
        spkez_c (ss.center, etecl, DEFFRAME, "NONE", SSB, Re, &lt);

        //----------------------------
        //Store data
        //----------------------------
        //Normalization
        ecl2neci(vout, Re, vin, ss);
        //Update yecl
        for(int i = 0; i <6; i++) yneci[i][p] = vin[i];
        //Update tecl, which is the normalized shifted time
        tneci[p] = etecl*ss.n;
    }
}

/**
 * \brief From ecliptic to synodic coordinates, vector format
 **/
void neci2synstate_vec(double **yneci, double *tneci, double **yout, double *tout, int N, double et0, double tsys0, int coord_eph, SS & ss)
{
    //=================================================================================
    //Init
    //=================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3], etecl;
    double vin[6], vout[6];

    //=================================================================================
    //Loop
    //=================================================================================
    double n_sys = mean_motion(coord_eph);
    for(int p = 0; p <= N; p++)
    {
        //----------------------------
        //Update etecl (non-normalized shifted time):
        //----------------------------
        etecl = tneci[p]/ss.n;

        //Position of the Center
        double Re[6], lt;
        spkez_c (ss.center, etecl, DEFFRAME, "NONE", SSB, Re, &lt);

        //----------------------------
        // COC
        //----------------------------
        //Init the COC
        init_ecl2synstate(etecl, B, C, &k, Bprim, kCprim, coord_eph);
        //Update vin
        for(int i = 0; i <6; i++) vin[i] = yneci[i][p];
        //Denormalization
        neci2ecl(vin, Re, vout, ss);
        //COC
        ecl2synstate(vout, vin, B, C, k, Bprim, kCprim, n_sys);

        //----------------------------
        //Store data
        //----------------------------
        //Update yout
        for(int i = 0; i <6; i++) yout[i][p] = vin[i];
        //Update tout (normalized time):
        tout[p] = (etecl - et0)*n_sys + tsys0;
    }
}




//========================================================================================
//                          Compute the acceleration of the
//                     m1-m2 line, where Bem is the Earth-Moon barycenter
//                           (in ecliptic coordinates)
//========================================================================================
/**
 *  \brief Return a component (along the dimension dim) of the m1-m2 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xBSecl(double et, int dim, int coord_type)
{
    SpiceDouble lt;
    SpiceDouble RS[6], RB[6];
    spkez_c (m1name(coord_type), et, DEFFRAME, "NONE", SSB, RS, &lt);
    spkez_c (m2name(coord_type), et, DEFFRAME, "NONE", SSB, RB, &lt);
    return RS[dim] - RB[dim];
}

/**
 *  \brief Return a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xm1ecl(double et, int dim, int coord_type)
{
    SpiceDouble lt;
    SpiceDouble RS[6];
    spkez_c (m1name(coord_type), et, DEFFRAME, "NONE", SSB, RS, &lt);
    return RS[dim];
}

/**
 *  \brief Return a component (along the dimension dim) of the m2 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double xm2ecl(double et, int dim, int coord_type)
{
    SpiceDouble lt;
    SpiceDouble RS[6];
    spkez_c (m2name(coord_type), et, DEFFRAME, "NONE", SSB, RS, &lt);
    return RS[dim];
}


/**
 *  \brief Return the derivative of a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double dxm1ecl(double et, int dim, int coord_type)
{
    double err;
    return dfridr(xm1ecl, et, dim, 1.0, &err, coord_type);
}

/**
 *  \brief Return the derivative of a component (along the dimension dim) of the m1 vector in ecliptic coordinates, at epoch et given in ephemeris time.
 **/
double dxm2ecl(double et, int dim, int coord_type)
{
    double err;
    return dfridr(xm2ecl, et, dim, 1.0, &err, coord_type);
}


/**
 *   \brief Returns the derivative of a function func at a point x by Ridders’ method of polynomial
 *          extrapolation. The value h is input as an estimated initial stepsize; it need not be small, but
 *          rather should be an increment in x over which func changes substantially. An estimate of the
 *          error in the derivative is returned as err. (from Numerical Recipe in C).
 **/
double dfridr(double (*func)(double, int, int), double x, int dim, double h, double *err, int coord_type)
{
    #define CON 1.4
    #define CON2 (CON*CON)
    #define BIG 1.0e30
    #define NTAB 10
    #define SAFE 2.0

    int i,j;
    double errt,fac,hh,**a,ans = 0.0;

    if (h == 0.0) nrerror((char*)"h must be nonzero in dfridr.");
    a=dmatrix(1,NTAB,1,NTAB);
    hh=h;
    a[1][1]=((*func)(x+hh, dim, coord_type)-(*func)(x-hh, dim, coord_type))/(2.0*hh);
    *err=BIG;

    //Successive columns in the Neville tableau will go to smaller stepsizes and higher orders of
    //extrapolation.
    for (i=2; i<=NTAB; i++)
    {
        hh /= CON;
        a[1][i]=((*func)(x+hh, dim, coord_type)-(*func)(x-hh, dim, coord_type))/(2.0*hh); //Try new, smaller stepsize
        fac=CON2;

        //Compute extrapolations of various orders, requiring no new function evaluations.
        for (j=2; j<=i; j++)
        {

            a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
            fac=CON2*fac;
            errt=max(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));

            //The error strategy is to compare each new extrapolation to one order lower, both
            //at the present stepsize and the previous one.
            if (errt <= *err)
            {
                //If error is decreased, save the improved answer.
                *err=errt;
                ans=a[j][i];
            }
        }
        if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break; //If higher order is worse by a significant factor SAFE, then quit early.
    }
    free_dmatrix(a,1,NTAB,1,NTAB);
    return ans;
}


//========================================================================================
//                         Test of the vector field
//========================================================================================
/**
 *  \brief Test of the JPL vector field using the LUTETIA asteroid
 **/
void test_asteroid()
{
    string asteroid = "SIDING SPRING";
    cout << "------------------------------------------------------------------------" << endl;
    cout << "Test of the JPL vector field using the  asteroid/comet " << asteroid      << endl;
    cout << "------------------------------------------------------------------------" << endl;

    //=========================================================================
    // Retrieve the initial position of the asteroid/comet
    //=========================================================================
    //    Summary for: C_G_1000012_2012_2017.bsp
    //
    //    Body: CHURYUMOV-GERASIMENKO (1000012) w.r.t. SUN (10)
    //          Start of Interval (ET)              End of Interval (ET)
    //          -----------------------------       -----------------------------
    //          2012 JAN 01 00:00:00.000            2017 JAN 01 00:00:00.000
    //
    //    Summary for: siding_spring_8-19-14.bsp
    //
    //    Body: SIDING SPRING (1003228) w.r.t. SUN (10)
    //          Start of Interval (ET)              End of Interval (ET)
    //          -----------------------------       -----------------------------
    //          2014 AUG 19 00:00:00.000            2014 DEC 30 00:00:00.000
    //
    //        Summary for: ceres-2003-2016.bsp
    //
    //        Bodies                       Start of Interval (ET)          End of Interval (ET)
    //        -------                      -----------------------------   -----------------------------
    //        2000001 CERES w.r.t. 10 SUN  2003 JAN 01 00:00:00.000        2016 JAN 01 00:00:00.000
    //
    //        Summary for: vesta-2003-2013.bsp
    //
    //        Bodies                       Start of Interval (ET)          End of Interval (ET)
    //        -------                      -----------------------------   -----------------------------
    //        2000004 VESTA w.r.t. 10 SUN  2003 JAN 01 00:00:00.000        2013 JAN 01 00:00:00.000
    //=========================================================================
    //Time
    SpiceDouble et_start, et_end, et, eti;

    //For "SIDING SPRING"
    string epoch_start  = "2014 AUG 20 00:00:00.000";
    string epoch_end    = "2014 DEC 29 00:00:00.000";


    str2et_c(epoch_start.c_str(), &et_start);
    str2et_c(epoch_end.c_str(),   &et_end);

    //Initial state
    SpiceDouble Rjpl[6], Ri[6],  lt;
    spkezr_c (asteroid.c_str(), et_start, DEFFRAME,  "NONE", DEFOBS, Rjpl, &lt);


    cout << asteroid << " initial state: " << endl;
    vector_printf_prec(Rjpl, 6);
    cout << "-----------------------------------------------------" << endl;


    //====================================================================================
    // Initialize the vector field
    //====================================================================================
    OdeStruct odestruct_JPL;
    //Root-finding
    const gsl_root_fsolver_type *T_root = gsl_root_fsolver_brent;
    //Stepper
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //Parameters
    OdeParams odeParams(&SEML, I_ECLI);
    //Init ode structure
    init_ode_structure(&odestruct_JPL, T, T_root, 6, jpl_vf, &odeParams);

    //====================================================================================
    // Save values
    //====================================================================================
    int Npoints = 500;
    double **y_from_int = dmatrix(0, 5, 0, Npoints);
    double *t_from_int  = dvector(0, Npoints);

    //====================================================================================
    // Integration
    //====================================================================================
    //Initialise the state
    //-----------------------
    for(int i = 0; i <6;i++) Ri[i] = Rjpl[i];

    //-----------------------
    //Loop
    //-----------------------
    et = et_start;
    int nt = 0, status = 0;
    do
    {
            eti = et_start + (double) nt *(et_end-et_start)/Npoints;
            status = gsl_odeiv2_driver_apply (odestruct_JPL.d, &et, eti, Ri);

            for(int k = 0; k < 6; k++) y_from_int[k][nt] = Ri[k];
            t_from_int[nt] = eti;

            //Advance one step
            nt++;
    }while((nt<=Npoints) && status == GSL_SUCCESS);


    //=========================================================================
    // Earth & Moon & Asteroid trajectory
    //=========================================================================
    double **y_earth_spice = dmatrix(0, 5, 0, Npoints);
    double **y_moon_spice = dmatrix(0, 5, 0, Npoints);
    double **y_lutetia_spice = dmatrix(0, 5, 0, Npoints);

    for(nt = 0; nt <= Npoints; nt++)
    {
        et = t_from_int[nt];

        //EARTH
        spkezr_c ("EARTH", et, DEFFRAME,  "NONE", DEFOBS, Rjpl, &lt);
        for(int k = 0; k < 6; k++) y_earth_spice[k][nt] = Rjpl[k];

        //MOON
        spkezr_c ("MOON", et, DEFFRAME,  "NONE", DEFOBS, Rjpl, &lt);
        for(int k = 0; k < 6; k++) y_moon_spice[k][nt] = Rjpl[k];

        //LUTETIA
        spkezr_c (asteroid.c_str(), et, DEFFRAME,  "NONE", DEFOBS, Rjpl, &lt);
        for(int k = 0; k < 6; k++) y_lutetia_spice[k][nt] = Rjpl[k];
    }

    //=========================================================================
    // Plot
    //=========================================================================
    gnuplot_ctrl *h1;
    h1 = gnuplot_init(true);

    gnuplot_cmd(true, h1, (char*) ("set title \"" + asteroid + " trajectory\"").c_str());
    gnuplot_set_xlabel(true, h1, (char*) "X [km]");
    gnuplot_set_ylabel(true, h1, (char*) "Y [km]");
    gnuplot_set_zlabel(true, h1, (char*) "Z [km]");

    gnuplot_plot_xyz(true, h1, y_from_int[0], y_from_int[1], y_from_int[2], Npoints+1, (char*) (asteroid + " (int)").c_str(), "lines", "2", "2", 1);
    gnuplot_plot_xyz(true, h1, y_lutetia_spice[0], y_lutetia_spice[1], y_lutetia_spice[2], Npoints+1, (char*) (asteroid + " (JPL)").c_str(), "lines", "2", "2", 2);
    gnuplot_plot_xyz(true, h1, y_earth_spice[0], y_earth_spice[1], y_earth_spice[2], Npoints+1, (char*) "EARTH", "lines", "2", "2", 3);
    gnuplot_plot_xyz(true, h1, y_moon_spice[0], y_moon_spice[1], y_moon_spice[2], Npoints+1, (char*) "MOON", "lines", "2", "2", 4);


    //=========================================================================
    // Final error
    //=========================================================================
    cout << "Final (relative) error = " << DENorm(Ri, Rjpl, 6)/ENorm(Ri, 6) << endl;
    cout << "-----------------------------------------------------" << endl;
    char ch;
    printf("Press ENTER to go on\n");
    scanf("%c",&ch);
}
/**
 *  \brief Generic test to display some feature of the COC
 **/
void test_coc(int coord_eph)
{
    //====================================================================================
    //Init the time
    //====================================================================================
    string epoch_start, epoch_end;
    SpiceDouble  et_start, et_end;
    epoch_start = "1950 JAN 01 00:00:00.000"; //We start in 1950
    epoch_end   = "1950 FEB 01 00:00:00.000"; //We end in 1950

    //In ephemeris time
    str2et_c(epoch_start.c_str(), &et_start);
    str2et_c(epoch_end.c_str(), &et_end);

    //====================================================================================
    //Init the COC, first case
    //====================================================================================
    double B[3], C[3][3], k, Bprim[3], kCprim[3][3];
    init_ecl2synstate(et_start, B, C, &k, Bprim, kCprim, coord_eph);

    //====================================================================================
    //Init the COC, second case
    //====================================================================================
    init_ecl2synstate_2(et_start, B, C, &k, Bprim, kCprim, coord_eph);
}


//========================================================================================
//                         Benchmark of the numerical constants
//========================================================================================
/**
 *  \brief Benchmark of the numerical constants that appear in the use of JPL eph.
 **/
void comp_num_const()
{
    //============================================
    // From QBCP, in SEM coordinates
    //============================================
    // Masses in QBCP
    double me_qbcp = SEML.us_sem.me;
    double mm_qbcp = SEML.us_sem.mm;
    double ms_qbcp = SEML.us_sem.ms;

    //============================================
    // From NASA, in SEM coordinates
    //============================================
    //The GM constants
    double GMe = SEML.cs_em.cr3bp.m1.GM;
    double GMm = SEML.cs_em.cr3bp.m2.GM;
    double GMs = SEML.cs_sem.cr3bp.m1.GM;

    //Back in SEM coordinates
    double me_jpl = GMe/(GMs + GMe + GMm);
    double mm_jpl = GMm/(GMs + GMe + GMm);
    double ms_jpl = GMs/(GMs + GMe + GMm);

    //============================================
    // Comparison
    //============================================
    cout << "--------------------------------------------------------" << endl;
    cout << "comp_num_const. Comparison of the values of the masses  " << endl;
    cout << "used throughout this software:                          " << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Earth mass, in SEM coordinates:" << endl;
    cout << "me_qbcp = "  << me_qbcp << endl;
    cout << "me_jpl = "   << me_jpl  << endl;
    cout << "relative error = " << (me_qbcp - me_jpl)/me_jpl << endl;
    cout << "absolute error = " << (me_qbcp - me_jpl) << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Moon mass, in SEM coordinates:" << endl;
    cout << "mm_qbcp = "  << mm_qbcp << endl;
    cout << "mm_jpl = "   << mm_jpl  << endl;
    cout << "relative error = " << (mm_qbcp - mm_jpl)/mm_jpl << endl;
    cout << "absolute error = " << (mm_qbcp - mm_jpl) << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Sun mass, in SEM coordinates:" << endl;
    cout << "ms_qbcp = "  << ms_qbcp << endl;
    cout << "ms_jpl = "   << ms_jpl  << endl;
    cout << "relative error = " << (ms_qbcp - ms_jpl)/ms_jpl << endl;
    cout << "absolute error = " << (ms_qbcp - ms_jpl) << endl;
    cout << "--------------------------------------------------------" << endl;

    cout << "--------------------------------------------------------" << endl;
    cout << "The uncertainty on the GMsun is (http://ssd.jpl.nasa.gov/?constants): " << 8 << " km^3/s^2" << endl;
    cout << "Back in SEM coordinates, this value is: " << 8/(GMm + GMe + GMs) << endl;
    cout << "--------------------------------------------------------" << endl;


    //============================================
    // Comparison for the mean Sun-Bem mean motion in EM coordinates
    //============================================
    double ns_qbcp   = SEML.us_em.ns;//from QBCP
    double ns_jpl    = SEML.cs_sem.cr3bp.T/SEML.cs_em.cr3bp.T; //from JPL data sheet
    double ns_gomez  = 0.01720209883844/86400/(2*M_PI)*SEML.cs_em.cr3bp.T; //from Gomez et al. 2002

    cout << "--------------------------------------------------------" << endl;
    cout << "Sun-Bem mean motion, in rad/s:" << endl;
    cout << "ns_qbcp = "    << ns_qbcp  << endl;
    cout << "ns_jpl = "     << ns_jpl   << endl;
    cout << "ns_gomez = "   << ns_gomez << endl;
    cout << "relative error (with jpl) = " << (ns_qbcp - ns_jpl)/ns_jpl << endl;
    cout << "absolute error (with jpl) = " << (ns_qbcp - ns_jpl) << endl;
    cout << "--------------------------------------------------------" << endl;
    cout << "Earth-Moon mean motion use for conversion = " << 2*M_PI/(SEML.cs_em.cr3bp.T) << "rad/s" << endl;
    cout << "Conclusion:" << endl;
    cout << "It seems like QBCP value has a great uncertainty. In fact, there is" << endl;
    cout << "no justification whatsoever for this value, neither in Andreu (1998)," << endl;
    cout << "nor in Gomez et al. Book IV." << endl;
    cout << "--------------------------------------------------------" << endl;
}


//========================================================================================
//                          Display information on SPICE Kernels.
//                          Routine are based on SPICE tutorials.
//========================================================================================
/**
 *  \brief Display information on a given kernel (currently de432s). The SRC code may be changed to include a user input for the name of the kernel
 **/
void displayKernelFeatures()
{
    /*
    Local parameters
    */
    #define  FILSIZ         256
    #define  MAXIV          1000
    #define  WINSIZ         ( 2 * MAXIV )
    #define  TIMLEN         51
    #define  MAXOBJ         1000

    /*
    Local variables
    */
    SPICEDOUBLE_CELL        ( cover, WINSIZ );
    SPICEINT_CELL           ( ids,   MAXOBJ );

    SpiceChar               *lsk;
    SpiceChar               *spk;
    SpiceChar               timstr  [ TIMLEN ];

    SpiceDouble             b;
    SpiceDouble             e;

    SpiceInt                i;
    SpiceInt                j;
    SpiceInt                niv;
    SpiceInt                obj;


    /*
    Load a leapseconds kernel for output time conversion.
    SPKCOV itself does not require a leapseconds kernel.
    */
    //prompt_c ( "Name of leapseconds kernel > ", FILSIZ, lsk );
    lsk = (char*) "spice/kernels/naif0008.tls";
    furnsh_c (lsk);

    /*
    Get name of SPK file.
    */
    //prompt_c ( "Name of SPK file           > ", FILSIZ, spk    );

    /*
    Find the set of objects in the SPK file.
    */
    spk = ((char*) "spice/kernels/de432s.bsp");
    spkobj_c ( spk, &ids );

    /*
    We want to display the coverage for each object. Loop over
    the contents of the ID code set, find the coverage for
    deach item in the set, and display the coverage.
    */
    for (i = 0;  i < card_c( &ids );  i++  )
    {
        /*
        Find the coverage window for the current object.
        Empty the coverage window each time so we don't
        include data for the previous object.
        */
        obj  =  SPICE_CELL_ELEM_I( &ids, i );

        scard_c  ( 0,        &cover );
        spkcov_c ( spk, obj, &cover );

        /*
        Get the number of intervals in the coverage window.
        */
        niv = wncard_c ( &cover );

        /*
        Display a simple banner.
        */
        printf ( "%s\n", "========================================" );

        printf ( "Coverage for object %i\n", obj );

        /*
        Convert the coverage interval start and stop times to TDB
        calendar strings.
        */
        for ( j = 0;  j < niv;  j++  )
        {
            /*
            Get the endpoints of the jth interval.
            */
            wnfetd_c ( &cover, j, &b, &e );

            /*
            Convert the endpoints to TDB calendar
            format time strings and display them.
            */
            timout_c ( b,
                       "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                       TIMLEN,
                       timstr                                  );

            printf ( "\n"
                     "Interval:  %i\n"
                     "Start:     %s\n",
                     j,
                     timstr            );

            timout_c ( e,
                       "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                       TIMLEN,
                       timstr                                  );
            printf ( "Stop:      %s\n", timstr );

        }

    }

}

/**
 *  \brief Display information on a given body accross multiple kernels loaded via a meta kernel (eg: spice/kernels/metakernel.furnsh).
 **/
void displayKernelFeaturesOneBody()
{

    /*
    Local parameters
    */
    #define  FILSIZ         256
    #define  LNSIZE         81
    #define  MAXCOV         100000
    #define  WINSIZS        ( 2 * MAXCOV )
    #define  TIMLEN         51

    /*
    Local variables
    */
    SPICEDOUBLE_CELL        ( cover, WINSIZS );

    SpiceBoolean            found;

    SpiceChar               file    [ FILSIZ ];
    SpiceChar               idch    [ LNSIZE ];
    SpiceChar               meta    [ FILSIZ ];
    SpiceChar               source  [ FILSIZ ];
    SpiceChar               timstr  [ TIMLEN ];
    SpiceChar               type    [ LNSIZE ];

    SpiceDouble             b;
    SpiceDouble             e;

    SpiceInt                count;
    SpiceInt                handle;
    SpiceInt                i;
    SpiceInt                idcode;
    SpiceInt                niv;


    /*
    Prompt for the metakernel name; load the metakernel.
    The metakernel lists the SPK files whose coverage
    for `idcode' we'd like to determine.  The metakernel
    must also specify a leapseconds kernel.
    */
    prompt_c ( "Name of metakernel > ", FILSIZ, meta );
    furnsh_c ( meta );


    /*
    Get the ID code of interest.
    */
    prompt_c ( "Enter ID code      > ", LNSIZE, idch );
    prsint_c ( idch,  &idcode );

    /*
    Find out how many kernels are loaded.  Loop over the
    kernels:  for each loaded SPK file, add its coverage
    for `idcode', if any, to the coverage window.
    */
    ktotal_c ( "SPK", &count );

    for ( i = 0;  i < count;  i++  )
    {
        kdata_c  ( i,     "SPK",   FILSIZ,  LNSIZE,   FILSIZ,
                   file,  type,    source,  &handle,  &found );

        spkcov_c ( file,  idcode,  &cover );
    }

    /*
    Display results.

    Get the number of intervals in the coverage window.
    */
    niv = wncard_c ( &cover );

    /*
    Display a simple banner.
    */
    printf ( "\nCoverage for object %i\n", idcode );

    /*
    Convert the coverage interval start and stop times to TDB
    calendar strings.
    */
    for ( i = 0;  i < niv;  i++  )
    {
        /*
        Get the endpoints of the ith interval.
        */
        wnfetd_c ( &cover, i, &b, &e );

        /*
        Convert the endpoints to TDB calendar
        format time strings and display them.
        */
        timout_c ( b,
                   "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                   TIMLEN,
                   timstr                                  );

        printf ( "\n"
                 "Interval:  %d\n"
                 "Start:     %s\n",
                 i,
                 timstr            );

        timout_c ( e,
                   "YYYY MON DD HR:MN:SC.### (TDB) ::TDB",
                   TIMLEN,
                   timstr                                  );
        printf ( "Stop:      %s\n", timstr );

    }
}
