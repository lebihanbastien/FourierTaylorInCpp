#include "vf.h"

//=============================================================================================
//
// Integration without STM
//
//=============================================================================================
//----------------------------------------------------------------
// Ecliptic coordinates for JPL ephemerides (X, X')
//----------------------------------------------------------------
/**
 * \brief Vector field of the solar system (JPL ephemerides) in ecliptic coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO
 **/
int jpl_vf(double et, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;

    double Rb[6], lt;

    //Loop on all the primaries
    for(int i = 0; i < qbp->ss.maxBodies; i++)
    {
        //Retrieve posiiton from SPICE kernel
        spkez_c(qbp->ss.id[i], et, DEFFRAME,  "NONE", SSB, Rb, &lt);
        //Build the vector field
        f[3] += - qbp->ss.Gmi[i] / pow((y[0] - Rb[0]) * (y[0] - Rb[0]) + (y[1] - Rb[1]) * (y[1] - Rb[1]) + (y[2] - Rb[2]) * (y[2] - Rb[2]), 3.0 / 2) * (y[0] - Rb[0]);
        f[4] += - qbp->ss.Gmi[i] / pow((y[0] - Rb[0]) * (y[0] - Rb[0]) + (y[1] - Rb[1]) * (y[1] - Rb[1]) + (y[2] - Rb[2]) * (y[2] - Rb[2]), 3.0 / 2) * (y[1] - Rb[1]);
        f[5] += - qbp->ss.Gmi[i] / pow((y[0] - Rb[0]) * (y[0] - Rb[0]) + (y[1] - Rb[1]) * (y[1] - Rb[1]) + (y[2] - Rb[2]) * (y[2] - Rb[2]), 3.0 / 2) * (y[2] - Rb[2]);
    }


    //-------------------------------------------------------------------------------
    //Display
    //-------------------------------------------------------------------------------
                cout << "y = " << endl;
                vector_printf_prec((double*) y, 6);
                cout << "et = " << et << endl;
                cout << "f[3] = " << f[3] << endl;
                cout << "f[4] = " << f[4] << endl;
                cout << "f[5] = " << f[5] << endl;
                char ch;
                printf("Press ENTER to plot the Initial Guess\n");
                scanf("%c",&ch);


    return 0;
}

/**
 *  \brief Computes the acceleration vector of the Earth.
 **/
void acc_earth_from_vf(double et, double Ae[3],  double Rj[11][6], QBCP_L* qbp)
{
    //------------------------------------------------------------------------------------
    //Reinit
    //------------------------------------------------------------------------------------
    double DR[3], DR2, Re[6], lt;
    for(int i = 0; i < 3; i++) Ae[i] = 0;

    //------------------------------------------------------------------------------------
    //Position of the Earth
    //------------------------------------------------------------------------------------
    spkez_c (EARTH, et, DEFFRAME, "NONE", SSB, Re, &lt);

    //------------------------------------------------------------------------------------
    //Acceleration of the Earth
    //------------------------------------------------------------------------------------
    for(int p = 0; p < qbp->ss.maxBodies; p++)
    {
        if(p != 3)
        {
            //DR = R1 - Rj
            for(int i = 0; i < 3; i++) DR[i]  = Re[i] - Rj[p][i];

            //DR2 = DR . DR
            DR2 = vdot_c(DR, DR);

            //A1 += Gmj/|DR|^3*DR
            for(int i = 0; i < 3; i++) Ae[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
        }
    }
}

/**
 *   \brief Update the positions of the primaries of the Solar system at epoch et.
 **/
void usspos(double et, double Rj[11][6], QBCP_L* qbp)
{
    double Rb[6], lt;

    for(int i = 0; i < qbp->ss.maxBodies; i++)
    {
        //Retrieve position from SPICE kernel
        spkez_c(qbp->ss.id[i], et, DEFFRAME,  "NONE", SSB, Rb, &lt);

        //Update Rj
        for(int j = 0; j < 6; j++) Rj[i][j] = Rb[j];
    }
}

/**
 * \brief Vector field of the solar system (JPL ephemerides) in earth-centered normalized. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 *         Note: for now the normalization is NOT taken into account because
 *         the stepper is going wrong when the state is normalized.
 **/
int jpl_vfn(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    double rb[6], et, rbe2;
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;

    //------------------------------------------------------------------------------------
    //Epoch in seconds
    //------------------------------------------------------------------------------------
    //et = t*qbp->ss.n;
    et = t;

    //------------------------------------------------------------------------------------
    // Update the positions of all the bodies
    //------------------------------------------------------------------------------------
    double Rj[11][6];
    usspos(et, Rj, qbp);

    //------------------------------------------------------------------------------------
    //Acceleration of the Earth
    //------------------------------------------------------------------------------------
    double Ae[3], ae[3];
    acc_earth_from_vf(et, Ae, Rj, qbp);
    //for(int i = 0; i <3; i++) ae[i] = Ae[i]/(qbp->ss.a*qbp->ss.n*qbp->ss.n);
    for(int i = 0; i <3; i++) ae[i] = Ae[i];

    //Position of the Earth
    double Re[6], lt;
    spkez_c (EARTH, et, DEFFRAME, "NONE", SSB, Re, &lt);

    //------------------------------------------------------------------------------------
    //Loop on all the primaries
    //------------------------------------------------------------------------------------
    for(int i = 0; i < qbp->ss.maxBodies; i++)
    {
        //To normalized coordinates
        ecl2ecln(Rj[i], Re, rb);//, qbp->ss);

        //Distance
        rbe2 = (y[0] - rb[0]) * (y[0] - rb[0]) + (y[1] - rb[1]) * (y[1] - rb[1]) + (y[2] - rb[2]) * (y[2] - rb[2]);

        //Build the vector field
        f[3] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[0] - rb[0]);
        f[4] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[1] - rb[1]);
        f[5] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[2] - rb[2]);
        //        f[3] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[0] - rb[0]);
        //        f[4] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[1] - rb[1]);
        //        f[5] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[2] - rb[2]);
    }

    //------------------------------------------------------------------------------------
    //Last contribution: acceleration of the Earth
    //------------------------------------------------------------------------------------
    f[3] += -ae[0];
    f[4] += -ae[1];
    f[5] += -ae[2];

    //-------------------------------------------------------------------------------
    //Display
    //-------------------------------------------------------------------------------
    //            double y0[6];
    //            ecln2ecl((double*)y, Re, y0, qbp->ss);
    //            cout << "y = " << endl;
    //            vector_printf_prec((double*) y0, 6);
    //            cout << "et = " << et << endl;
    //            cout << "f[3] = " << f[3]* pow(qbp->ss.n, 2.0)* pow(qbp->ss.a, 1.0) + Ae[0] << endl;
    //            cout << "f[4] = " << f[4]* pow(qbp->ss.n, 2.0)* pow(qbp->ss.a, 1.0) + Ae[1] << endl;
    //            cout << "f[5] = " << f[5]* pow(qbp->ss.n, 2.0)* pow(qbp->ss.a, 1.0) + Ae[2] << endl;
    //            cout << "-----------------------" << endl;
    //            cout << "y = " << endl;
    //            vector_printf_prec((double*) y, 6);
    //            cout << "f[3] = " << f[3] << endl;
    //            cout << "f[4] = " << f[4] << endl;
    //            cout << "f[5] = " << f[5] << endl;
    //            cout << "-----------------------" << endl;
    //            cout << "-----------------------" << endl;

    //            char ch;
    //            printf("Press ENTER to plot the Initial Guess\n");
    //            scanf("%c",&ch);


    return 0;
}

//----------------------------------------------------------------
// Synodic coordinates for JPL ephemerides (x, x')
//----------------------------------------------------------------

/**
 *  \brief Computes the acceleration vectors A1 and A2 of the two primaries that defines the synodic coordinates.
 *         They are computed from the newtonian vector field, i.e. the gravitational influence of the massive bodies of the solar system.
 *         If the smaller primary is the Earth-Moon Barycenter, the acceleration is computed as the barycenter of the acceleration of the Earth and the Moon.
 **/
void acc_from_vf(double et, int pos1, int pos2, double A1[3], double A2[3], double R1[3], double R2[3], double Rj[11][6], QBCP_L* qbp)
{
    //------------------------------------------------------------------------------------
    //Reinit
    //------------------------------------------------------------------------------------
    double DR[3], DR2;

    for(int i = 0; i < 3; i++) A1[i] = 0;

    for(int i = 0; i < 3; i++) A2[i] = 0;

    //------------------------------------------------------------------------------------
    //Switch
    //------------------------------------------------------------------------------------
    if(pos1 == 0 && pos2 == 11)
    {
        //--------------------------------------------------------------------------------
        //Sun-(Earth+Moon) case
        //--------------------------------------------------------------------------------
        double Re[6], Rm[6], lt;

        //------------------------------------------------
        // Acceleration of the Earth-Moon barycenter
        //------------------------------------------------
        spkez_c (EARTH, et, DEFFRAME, "NONE", SSB, Re, &lt);
        spkez_c (MOON,  et, DEFFRAME, "NONE", SSB, Rm, &lt);

        //Loop on all primaries
        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            //Earth
            if(p != 3)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = Re[i] - Rj[p][i];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);

                //A1 += Gmj/|DR|^3*DR
                for(int i = 0; i < 3; i++) A1[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
            }

            //Moon
            if(p != 4)
            {
                //DR = R2 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = Rm[i] - Rj[p][i];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);

                //A2 += Gmj/|DR|^3*DR
                for(int i = 0; i < 3; i++) A2[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
            }
        }

        //Barycenter of Earth & Moon
        double Gme = qbp->ss.Gmi[3]; //Earth mass
        double Gmm = qbp->ss.Gmi[4]; //Moon mass

        for(int i = 0; i < 3; i++) A2[i] = (Gme * A1[i] + Gmm * A2[i]) / (Gme + Gmm);

        //Sun
        for(int i = 0; i < 3; i++) A1[i] = 0;

        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            if(p != pos1)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R1[i] - Rj[p][i];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);

                //A1 += Gmj/|DR|^3*DR
                for(int i = 0; i < 3; i++) A1[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
            }
        }
    }
    else
    {
        //--------------------------------------------------------------------------------
        //Classic case
        //--------------------------------------------------------------------------------
        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            //First primary
            if(p != pos1)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R1[i] - Rj[p][i];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);

                //A1 += Gmj/|DR|^3*DR
                for(int i = 0; i < 3; i++) A1[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
            }

            //Second primary
            if(p != pos2)
            {
                //DR = R2 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R2[i] - Rj[p][i];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);

                //A2 += Gmj/|DR|^3*DR
                for(int i = 0; i < 3; i++) A2[i] += - qbp->ss.Gmi[p] / pow(DR2, 3.0 / 2) * DR[i];
            }

        }
    }

}

/**
 *  \brief Computes the jerk vectors A1 and A2 of the two primaries that defines the synodic coordinates.
 *         They are computed from the derivative of the newtonian vector field, i.e. the gravitational influence of the massive bodies of the solar system.
 *         If the smaller primary is the Earth-Moon Barycenter, the jerk is computed as the barycenter of the jerk of the Earth and the Moon.
 **/
void jerk_from_vf(double et, int pos1, int pos2, double J1[3], double J2[3], double R1[3], double R2[3], double Rj[11][6], QBCP_L* qbp)
{
    //------------------------------------------------------------------------------------
    //Reinit
    //------------------------------------------------------------------------------------
    double DR[3], DV[3], DR2, DRV;

    for(int i = 0; i < 3; i++) J1[i] = 0;

    for(int i = 0; i < 3; i++) J2[i] = 0;

    //------------------------------------------------------------------------------------
    //Switch
    //------------------------------------------------------------------------------------
    if(pos1 == 0 && pos2 == 11)
    {
        //--------------------------------------------------------------------------------
        //Sun-(Earth+Moon) case
        //--------------------------------------------------------------------------------
        double Re[6], Rm[6], lt;

        //------------------------------------------------
        // Jerk of the Earth-Moon barycenter
        //------------------------------------------------
        spkez_c (EARTH, et, DEFFRAME, "NONE", SSB, Re, &lt);
        spkez_c (MOON,  et, DEFFRAME, "NONE", SSB, Rm, &lt);

        //Loop on all primaries
        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            //Earth
            if(p != 3)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = Re[i]   - Rj[p][i];

                //DV = V1 - Vj
                for(int i = 0; i < 3; i++) DV[i]  = Re[i + 3] - Rj[p][i + 3];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);
                //DRV = DR . DV
                DRV = vdot_c(DR, DV);

                //J1 += Gmj*( 3 * DR * (DR . DV) / |DR|^5 - DV / |DR|^3)
                for(int i = 0; i < 3; i++) J1[i] += qbp->ss.Gmi[p] * (3.0 * DR[i] * DRV / pow(DR2, 5.0 / 2) - DV[i] / pow(DR2, 3.0 / 2));
            }

            //Moon
            if(p != 4)
            {
                //DR = R2 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = Rm[i]   - Rj[p][i];

                //DV = V2 - Vj
                for(int i = 0; i < 3; i++) DV[i]  = Rm[i + 3] - Rj[p][i + 3];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);
                //DRV = DR . DV
                DRV = vdot_c(DR, DV);

                //J2 += Gmj*( 3 * DR * (DR . DV) / |DR|^5 - DV / |DR|^3)
                for(int i = 0; i < 3; i++) J2[i] += qbp->ss.Gmi[p] * (3.0 * DR[i] * DRV / pow(DR2, 5.0 / 2) - DV[i] / pow(DR2, 3.0 / 2));
            }

        }

        //Barycenter of Earth & Moon
        double Gme = qbp->ss.Gmi[3]; //Earth mass
        double Gmm = qbp->ss.Gmi[4]; //Moon mass

        for(int i = 0; i < 3; i++) J2[i] = (Gme * J1[i] + Gmm * J2[i]) / (Gme + Gmm);


        //------------------------------------------------
        // Jerk of the Sun
        //------------------------------------------------
        for(int i = 0; i < 3; i++) J1[i] = 0;

        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            if(p != pos1)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R1[i]   - Rj[p][i];

                //DV = V1 - Vj
                for(int i = 0; i < 3; i++) DV[i]  = R1[i + 3] - Rj[p][i + 3];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);
                //DRV = DR . DV
                DRV = vdot_c(DR, DV);

                //J1 += Gmj*( 3 * DR * (DR . DV) / |DR|^5 - DV / |DR|^3)
                for(int i = 0; i < 3; i++) J1[i] += qbp->ss.Gmi[p] * (3.0 * DR[i] * DRV / pow(DR2, 5.0 / 2) - DV[i] / pow(DR2, 3.0 / 2));
            }

        }

    }
    else
    {
        //--------------------------------------------------------------------------------
        //Classic case
        //--------------------------------------------------------------------------------
        for(int p = 0; p < qbp->ss.maxBodies; p++)
        {
            if(p != pos1)
            {
                //DR = R1 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R1[i]   - Rj[p][i];

                //DV = V1 - Vj
                for(int i = 0; i < 3; i++) DV[i]  = R1[i + 3] - Rj[p][i + 3];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);
                //DRV = DR . DV
                DRV = vdot_c(DR, DV);

                //J1 += Gmj*( 3 * DR * (DR . DV) / |DR|^5 - DV / |DR|^3)
                for(int i = 0; i < 3; i++) J1[i] += qbp->ss.Gmi[p] * (3.0 * DR[i] * DRV / pow(DR2, 5.0 / 2) - DV[i] / pow(DR2, 3.0 / 2));
            }

            if(p != pos2)
            {
                //DR = R2 - Rj
                for(int i = 0; i < 3; i++) DR[i]  = R2[i]   - Rj[p][i];

                //DV = V2 - Vj
                for(int i = 0; i < 3; i++) DV[i]  = R2[i + 3] - Rj[p][i + 3];

                //DR2 = DR . DR
                DR2 = vdot_c(DR, DR);
                //DRV = DR . DV
                DRV = vdot_c(DR, DV);

                //J2 += Gmj*( 3 * DR * (DR . DV) / |DR|^5 - DV / |DR|^3)
                for(int i = 0; i < 3; i++) J2[i] += qbp->ss.Gmi[p] * (3.0 * DR[i] * DRV / pow(DR2, 5.0 / 2) - DV[i] / pow(DR2, 3.0 / 2));
            }

        }
    }

}

/**
 * \brief Vector field of the solar system (JPL ephemerides) in synodical coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO (11 bodies)
 **/
int jpl_vf_syn(double t, const double y[], double f[], void* params_void)
{
    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp   = (QBCP_L*) params_void;
    int coord_eph = qbp->ss.coord_eph;        //Current coord_eph
    double a      = qbp->ss.a;                //mean semi-major axis (in km)
    double et0    = qbp->ss.et0;              //starting epoch
    double t0     = qbp->ss.t0;               //starting time
    double mu1    = qbp->ss.mu1;              //mass ratio
    double mu2    = qbp->ss.mu2;              //1-mass ratio
    double nsys   = qbp->ss.n;                //mean motion (in s)
    double pos1   = qbp->ss.pos1;             //indix of m1
    double pos2   = qbp->ss.pos2;             //indix of m2
    double et     = et0 + (t-t0)/nsys;        //Current time in seconds

    //------------------------------------------------------------------------------------
    //Display:
    //------------------------------------------------------------------------------------
    //    cout << " a      = " << qbp->ss.a << endl;                //mean semi-major axis (in km)
    //    cout << " et0    = " << qbp->ss.et0 << endl;              //starting epoch
    //    cout << " t0     = " << qbp->ss.t0 << endl;               //starting time
    //    cout << " mu1    = " << qbp->ss.mu1 << endl;              //mass ratio
    //    cout << " mu2    = " << qbp->ss.mu2 << endl;              //1-mass ratio
    //    cout << " nsys   = " << qbp->ss.n << endl;                //mean motion (in s)
    //    cout << " pos1   = " << qbp->ss.pos1 << endl;             //indix of m1
    //    cout << " pos2   = " << qbp->ss.pos2 << endl;             //indix of m2
    //    cout << " et     = " << et0 + (t-t0)/nsys << endl;        //Current time in seconds

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z'
    //------------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //------------------------------------------------------------------------------------
    // Update the positions of all the bodies
    //------------------------------------------------------------------------------------
    double Rj[11][6];
    double rj[11][3];
    usspos(et, Rj, qbp);

    //------------------------------------------------------------------------------------
    //Position, velocity, acceleration, and jerk of the primaries
    //------------------------------------------------------------------------------------
    double lt;
    double R1[6];  //state of m1
    double R2[6];  //state of m2
    double A1[3];  //acceleration of m1
    double A2[3];  //acceleration of m2
    double J1[3];  //jerk of m1
    double J2[3];  //jerk of m2

    //--------------------------------------------
    //Position, velocity of the primaries, in dimensionnalized units
    //--------------------------------------------
    spkez_c (m1name(coord_eph), et, DEFFRAME, "NONE", SSB, R1, &lt);
    spkez_c (m2name(coord_eph), et, DEFFRAME, "NONE", SSB, R2, &lt);

    //--------------------------------------------
    //Acceleration of the primaries, in dimensionnalized units
    //--------------------------------------------
    double DR[3], DR2;
    //First possibility: from ephemerides
    //    double err;
    //    for(int i = 0; i < 3; i++)  A1[i] = dfridr(xm1ecl, et, i+3, 1.0, &err, coord_eph);
    //    for(int i = 0; i < 3; i++)  A2[i] = dfridr(xm2ecl, et, i+3, 1.0, &err, coord_eph);
    //Second possibility: from the vector field
    acc_from_vf(et, pos1, pos2, A1, A2, R1, R2, Rj, qbp);

    //--------------------------------------------
    //Jerk of the primaries, in dimensionnalized units
    //--------------------------------------------
    //First possibility: from ephemerides
    //    for(int i = 0; i < 3; i++)  J1[i] = dfridr(dxm1ecl, et, i+3, 1.0, &err, coord_eph);
    //    for(int i = 0; i < 3; i++)  J2[i] = dfridr(dxm2ecl, et, i+3, 1.0, &err, coord_eph);
    //Second possibility: from vector field
    jerk_from_vf(et, pos1, pos2, J1, J2, R1, R2, Rj, qbp);


    //------------------------------------------------------------------------------------
    //Relative position, velocity, acceleration, and jerk
    //------------------------------------------------------------------------------------
    double r21[3]; //relative position
    double v21[3]; //relative velocity
    double a21[3]; //relative acceleration
    double j21[3]; //relative jerk
    double B[3];   //position of the barycenter vector
    double Bpp[3]; //acceleration of the barycenter vector

    for(int i = 0; i < 3; i++)
    {
        //Difference between the two position vectors
        r21[i]  = R1[i] - R2[i];
        //Derivative of the difference wrt to dimentionalized time
        v21[i]  = R1[i + 3] - R2[i + 3];
        //Difference between the two accelerations
        a21[i]  = A1[i] - A2[i];
        //Difference between the two jerks
        j21[i]  = J1[i] - J2[i];
        //Position of the barycenter
        B[i]    = (mu1 * R1[i] + mu2 * R2[i]) / (mu1 + mu2);
        //Acceleration of the barycenter vector
        Bpp[i]  = (mu1 * A1[i] + mu2 * A2[i]) / (mu1 + mu2);
    }

    //------------------------------------------------------------------------------------
    //Orthonormal basis (e1, e2, e3) + h and k factors and their derivatives
    //------------------------------------------------------------------------------------
    double e1[3], e2[3], e3[3], e4[3], k, kp, kpp, h, hp, hpp;
    double r21_v21[3], r21_a21[3], r21_j21[3], v21_a21[3];
    //Building the cross products
    init_cp(r21, v21, a21, j21, r21_v21, r21_a21, r21_j21, v21_a21);
    //Building h, h', h"
    init_h(r21_v21, r21_a21, r21_j21, v21_a21, &h, &hp, &hpp, e4);
    //Building k, k', k''
    init_k(r21, v21, a21, &k, &kp, &kpp);
    //Building e1, e2, e3
    init_ei(r21, v21, e1, e2, e3);

    //------------------------------------------------------------------------------------
    //Orthonormal basis derivatives (e1', e2', e3')
    //------------------------------------------------------------------------------------
    double de1[3], de2[3], de3[3];
    init_dei(r21, v21, a21, k, kp, h, hp, e1, e3, de1, de2, de3, e4);

    //------------------------------------------------------------------------------------
    //Orthonormal basis double derivatives (e1", e2", e3")
    //------------------------------------------------------------------------------------
    double dde1[3], dde2[3], dde3[3];
    init_ddei(r21, v21, a21, r21_v21, r21_a21, r21_j21, v21_a21,
              k, kp, kpp, h, hp, hpp, e1, e3, de1, de3,
              dde1, dde2, dde3, e4);

    //------------------------------------------------------------------------------------
    // b coefficients. Note that a factor 1/n^2 is added to account for
    // the time normalization
    //------------------------------------------------------------------------------------
    double b[13];

    //b[0] = - B" . e1 /(k * n^2);
    b[0] = -vdot_c(Bpp, e1) / (k * nsys * nsys);

    //b[1] = - B" . e2 /(k * n^2);
    b[1] = -vdot_c(Bpp, e2) / (k * nsys * nsys);

    //b[2] = - B" . e3 /(k * n^2);
    b[2] = -vdot_c(Bpp, e3) / (k * nsys * nsys);

    //b[3] = - 2*k'/(n*k)
    b[3] = -2 * kp / (k * nsys);

    //b[4] = 2/n e2 . e1'
    b[4] = +2 * vdot_c(e2, de1) / nsys;

    //b[5] = 2/n e3 . e2'
    b[5] = +2 * vdot_c(e3, de2) / nsys;

    //b[6] = -1/n^2*(k"/k - e1' . e1')
    b[6] = - (kpp / k - vdot_c(de1, de1)) / (nsys * nsys);

    //b[7] = +1/n^2*(e1'.e3')
    b[7] = vdot_c(de1, de3) / (nsys * nsys);

    //b[8] = +1/n^2*(2*k'/k e2 . e1' + e2 . e1")
    b[8] = (2 * kp / k * vdot_c(e2, de1) + vdot_c(e2, dde1)) / (nsys * nsys);

    //b[9] = -1/n^2*(k"/k - e2' . e2')
    b[9] = - (kpp / k - vdot_c(de2, de2)) / (nsys * nsys);

    //b[10] = +1/n^2*(2*k'/k e3 . e3' + e3 . e2")
    b[10] = (2 * kp / k * vdot_c(e3, de2) + vdot_c(e3, dde2)) / (nsys * nsys);

    //b[11] = -1/n^2*(k"/k - e3' . e3')
    b[11] = - (kpp / k - vdot_c(de3, de3)) / (nsys * nsys);

    //b[12] = a^3/k^3
    b[12] = a * a * a / (k * k * k);

    //////////////////////////////////////////////////////////////////////////////////////
    //
    //€€TODO:
    //read the double derivatives of the ei + hpp and kpp
    //set a latex document that sum up these calculus (see p44 of Dei Tos).
    //
    //////////////////////////////////////////////////////////////////////////////////////


    //------------------------------------------------------------------------------------
    //Reduced positions of the primaries
    //------------------------------------------------------------------------------------
    double C[3][3];

    for(int i = 0; i < 3; i++) C[i][0] = e1[i];

    for(int i = 0; i < 3; i++) C[i][1] = e2[i];

    for(int i = 0; i < 3; i++) C[i][2] = e3[i];

    double vin[3], vout[3];

    for(int p = 0; p < qbp->ss.maxBodies; p++)
    {
        for(int i = 0; i < 3; i++) vin[i] = Rj[p][i];
        ecl2synpos(vin, vout, B, C, k);
        for(int i = 0; i < 3; i++) rj[p][i] = vout[i];
    }


    //------------------------------------------------------------------------------------
    //Phase space derivatives: x", y", z"
    //------------------------------------------------------------------------------------
    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = b[0] + b[3] * y[3] + b[4] * y[4]               + b[6] * y[0] + b[8] * y[1] + b[7]  * y[2];
    f[4] = b[1] - b[4] * y[3] + b[3] * y[4] + b[5] * y[5] - b[8] * y[0] + b[9] * y[1] + b[10] * y[2];
    f[5] = b[2]               - b[5] * y[4] + b[3] * y[5] + b[7] * y[0] - b[10]* y[1] + b[11] * y[2];

    //Loop on all the primaries
    for(int p = 0; p < qbp->ss.maxBodies; p++)
    {
        //DR = (x, y, z)^T - rj
        for(int i = 0; i < 3; i++) DR[i]  = y[i] - rj[p][i];

        //DR2 = DR . DR
        DR2 = vdot_c(DR, DR);
        //Build the vector field
        f[3] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[0];
        f[4] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[1];
        f[5] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[2];
    }

    return 0;


}

//----------------------------------------------------------------
// Inertial SEM coordinates, (X, V)
//----------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in SEM inertial coordinates
 **/
int qbcp_vf_insem(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    double ms = qbp->us_sem.ms;
    double me = qbp->us_sem.me;
    double mm = qbp->us_sem.mm;
    double n  = qbp->us_sem.n;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;


    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @t, in SEM coordinates
    //-------------------------------------------------------------------------------
    double Ps_syn[3];
    evaluateCoef(Ps_syn, t, n, qbp->nf, qbp->cs_sem.Ps, 3);
    double Pe_syn[3];
    evaluateCoef(Pe_syn, t, n, qbp->nf, qbp->cs_sem.Pe, 3);
    double Pm_syn[3];
    evaluateCoef(Pm_syn, t, n, qbp->nf, qbp->cs_sem.Pm, 3);

    //-------------------------------------------------------------------------------
    //Back to  INSEM coordinates
    //-------------------------------------------------------------------------------
    double Ps[6];
    SEMtoIN(t, Ps_syn, Ps, qbp);
    double Pe[6];
    SEMtoIN(t, Pe_syn, Pe, qbp);
    double Pm[6];
    SEMtoIN(t, Pm_syn, Pm, qbp);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);


    //-------------------------------------------------------------------------------
    // Update the vector field
    //-------------------------------------------------------------------------------
    if(me != 0) f[3] += - me / pow(qPe2, 3.0 / 2) * (y[0] - Pe[0]);
    if(mm != 0) f[3] += - mm / pow(qPm2, 3.0 / 2) * (y[0] - Pm[0]);
    if(ms != 0) f[3] += - ms / pow(qPs2, 3.0 / 2) * (y[0] - Ps[0]);

    if(me != 0) f[4] += - me / pow(qPe2, 3.0 / 2) * (y[1] - Pe[1]);
    if(mm != 0) f[4] += - mm / pow(qPm2, 3.0 / 2) * (y[1] - Pm[1]);
    if(ms != 0) f[4] += - ms / pow(qPs2, 3.0 / 2) * (y[1] - Ps[1]);

    if(me != 0) f[5] += - me / pow(qPe2, 3.0 / 2) * (y[2] - Pe[2]);
    if(mm != 0) f[5] += - mm / pow(qPm2, 3.0 / 2) * (y[2] - Pm[2]);
    if(ms != 0) f[5] += - ms / pow(qPs2, 3.0 / 2) * (y[2] - Ps[2]);

    return 0;
}

//----------------------------------------------------------------
// Inertial EM coordinates, (X, V)
//----------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in EM inertial coordinates
 **/
int qbcp_vf_inem(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    double ms = qbp->us_em.ms;
    double me = qbp->us_em.me;
    double mm = qbp->us_em.mm;
    double n  = qbp->us_em.n;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;


    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @t, in SEM coordinates
    //-------------------------------------------------------------------------------
    double Ps_syn[3];
    evaluateCoef(Ps_syn, t, n, qbp->nf, qbp->cs_em.Ps, 3);
    double Pe_syn[3];
    evaluateCoef(Pe_syn, t, n, qbp->nf, qbp->cs_em.Pe, 3);
    double Pm_syn[3];
    evaluateCoef(Pm_syn, t, n, qbp->nf, qbp->cs_em.Pm, 3);

    //-------------------------------------------------------------------------------
    //Back to  INSEM coordinates
    //-------------------------------------------------------------------------------
    double Ps[6];
    EMvtoIN(t, Ps_syn, Ps, qbp);
    double Pe[6];
    EMvtoIN(t, Pe_syn, Pe, qbp);
    double Pm[6];
    EMvtoIN(t, Pm_syn, Pm, qbp);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);


    //-------------------------------------------------------------------------------
    // Update the vector field
    //-------------------------------------------------------------------------------
    if(me != 0) f[3] += - me / pow(qPe2, 3.0 / 2) * (y[0] - Pe[0]);
    if(mm != 0) f[3] += - mm / pow(qPm2, 3.0 / 2) * (y[0] - Pm[0]);
    if(ms != 0) f[3] += - ms / pow(qPs2, 3.0 / 2) * (y[0] - Ps[0]);

    if(me != 0) f[4] += - me / pow(qPe2, 3.0 / 2) * (y[1] - Pe[1]);
    if(mm != 0) f[4] += - mm / pow(qPm2, 3.0 / 2) * (y[1] - Pm[1]);
    if(ms != 0) f[4] += - ms / pow(qPs2, 3.0 / 2) * (y[1] - Ps[1]);

    if(me != 0) f[5] += - me / pow(qPe2, 3.0 / 2) * (y[2] - Pe[2]);
    if(mm != 0) f[5] += - mm / pow(qPm2, 3.0 / 2) * (y[2] - Pm[2]);
    if(ms != 0) f[5] += - ms / pow(qPs2, 3.0 / 2) * (y[2] - Ps[2]);

    return 0;
}


//----------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------
//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations)
int qbfbp_vfn_novar(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    int noc      = qbp->numberOfCoefs;
    int nf       = qbp->nf;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @ t
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0] - pe[0]) * (y[0] - pe[0]) + (y[1] - pe[1]) * (y[1] - pe[1]) + (y[2] - pe[2]) * (y[2] - pe[2]);
    double qps2 = (y[0] - ps[0]) * (y[0] - ps[0]) + (y[1] - ps[1]) * (y[1] - ps[1]) + (y[2] - ps[2]) * (y[2] - ps[2]);
    double qpm2 = (y[0] - pm[0]) * (y[0] - pm[0]) + (y[1] - pm[1]) * (y[1] - pm[1]) + (y[2] - pm[2]) * (y[2] - pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return 0;
}

//Update the normalized vector field of the state. Note that alpha[14] (alpha15) is zero for the QBCP
int vfn_state(const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm,
              double gamma)
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0] * y[3] + alpha[1] * y[0] + alpha[2] * y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0] * y[4] + alpha[1] * y[1] - alpha[2] * y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0] * y[5] + alpha[1] * y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1] * y[3] + alpha[2] * y[4] + alpha[14] * y[0]
           + alpha[12];

    if(me != 0) f[3] += - alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[0] - pe[0]);

    if(mm != 0) f[3] += - alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[0] - pm[0]);

    if(ms != 0) f[3] += - alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[0] - ps[0]);

    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1] * y[4] - alpha[2] * y[3] + alpha[14] * y[1]
           + alpha[13];

    if(me != 0) f[4] += - alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[1] - pe[1]);

    if(mm != 0) f[4] += - alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[1] - pm[1]);

    if(ms != 0) f[4] += - alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[1] - ps[1]);


    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1] * y[5] + alpha[14] * y[2];

    if(me != 0) f[5] += - alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[2] - pe[2]);

    if(mm != 0) f[5] += - alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[2] - pm[2]);

    if(ms != 0) f[5] += - alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[2] - ps[2]);

    return 0;
}


//----------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------
//Update the normalized vector field of the state. The state is here (x, xdot), and not (x, px)
int vfn_state_xv(const double y[], double f[], double alpha[], double alphad[],
                 double ps[], double pe[], double pm[],
                 double qps2, double qpe2, double qpm2,
                 double ms, double me, double mm,
                 double gamma)
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = y[3];
    //f[1] = y'
    //------------------------------------
    f[1] = y[4];
    //f[2] = z'
    //------------------------------------
    f[2] = y[5];

    //-------------------------------------------------------------------------------
    //Intermediate coefficients
    //-------------------------------------------------------------------------------
    double alpha16 = alphad[1] - alphad[0] * alpha[1] / alpha[0] + alpha[1] * alpha[1] + alpha[2] * alpha[2]; //
    double alpha17 = alphad[2] - alphad[0] * alpha[2] / alpha[0]; //
    double alpha18 = alphad[0] / alpha[0]; //

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x'', y'', z''
    //-------------------------------------------------------------------------------
    //f[3] = x''
    //------------------------------------
    f[3] = alpha16 * y[0] + alpha17 * y[1] + alpha18 * y[3] + 2 * alpha[2] * y[4] + alpha[0] * alpha[12];

    if(me != 0) f[3] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[0] - pe[0]);

    if(mm != 0) f[3] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[0] - pm[0]);

    if(ms != 0) f[3] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[0] - ps[0]);

    //f[4] = y''
    //------------------------------------
    f[4] = - alpha17 * y[0] + alpha16 * y[1] - 2 * alpha[2] * y[3] + alpha18 * y[4] + alpha[0] * alpha[13];

    if(me != 0) f[4] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[1] - pe[1]);

    if(mm != 0) f[4] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[1] - pm[1]);

    if(ms != 0) f[4] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[1] - ps[1]);


    //f[5] = z''
    //------------------------------------
    f[5] = (alpha16 - alpha[2] * alpha[2]) * y[2] + alpha18 * y[5];

    if(me != 0) f[5] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * me / pow(qpe2, 3.0 / 2) * (y[2] - pe[2]);

    if(mm != 0) f[5] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * mm / pow(qpm2, 3.0 / 2) * (y[2] - pm[2]);

    if(ms != 0) f[5] -=  alpha[0] * alpha[5] / pow(gamma, 3.0) * ms / pow(qps2, 3.0 / 2) * (y[2] - ps[2]);

    return 0;
}

//Vector field of the QBCP with EM units and Normalized-Li centered coordinates (no variationnal equations). The state is here (x, xdot), and not (x, px)
int qbfbp_vfn_novar_xv(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    int noc      = qbp->numberOfCoefs;
    int nf       = qbp->nf;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @ t, as well as the derivatives of the first 3 alphas
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, nf, qbp->cs.coeffs, noc);
    double alphad[3];
    evaluateCoefDerivatives(alphad, t, n, nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, nf, qbp->cs.pm, 3);


    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0] - pe[0]) * (y[0] - pe[0]) + (y[1] - pe[1]) * (y[1] - pe[1]) + (y[2] - pe[2]) * (y[2] - pe[2]);
    double qps2 = (y[0] - ps[0]) * (y[0] - ps[0]) + (y[1] - ps[1]) * (y[1] - ps[1]) + (y[2] - ps[2]) * (y[2] - ps[2]);
    double qpm2 = (y[0] - pm[0]) * (y[0] - pm[0]) + (y[1] - pm[1]) * (y[1] - pm[1]) + (y[2] - pm[2]) * (y[2] - pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state_xv(y, f, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    return 0;
}

//----------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame.
 **/
int qbfbp_vf(double t, const double y[], double f[], void* params_void)
{

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vf_state(y, f, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    return GSL_SUCCESS;
}

/**
 *  \brief Update the vector field of the state in system coordinates (non-normalized). Note that alpha[14] (alpha15) is zero for the QBCP
 **/
int vf_state( const double y[], double f[], double alpha[],
              double ps[], double pe[], double pm[],
              double qps2, double qpe2, double qpm2,
              double ms, double me, double mm )
{
    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    //f[0] = x'
    //------------------------------------
    f[0] = alpha[0] * y[3] + alpha[1] * y[0] + alpha[2] * y[1];
    //f[1] = y'
    //------------------------------------
    f[1] = alpha[0] * y[4] + alpha[1] * y[1] - alpha[2] * y[0];
    //f[2] = z'
    //------------------------------------
    f[2] = alpha[0] * y[5] + alpha[1] * y[2];
    //f[3] = px'
    //------------------------------------
    f[3] = - alpha[1] * y[3] + alpha[2] * y[4] + alpha[14] * y[0]
           - alpha[3];

    if(me != 0) f[3] += - alpha[5] * me / pow(qpe2, 3.0 / 2) * (y[0] - pe[0]);

    if(mm != 0) f[3] += - alpha[5] * mm / pow(qpm2, 3.0 / 2) * (y[0] - pm[0]);

    if(ms != 0) f[3] += - alpha[5] * ms / pow(qps2, 3.0 / 2) * (y[0] - ps[0]);


    //f[4] = py'
    //------------------------------------
    f[4] = - alpha[1] * y[4] - alpha[2] * y[3] + alpha[14] * y[1]
           - alpha[4];

    if(me != 0) f[4] += - alpha[5] * me / pow(qpe2, 3.0 / 2) * (y[1] - pe[1]);

    if(mm != 0) f[4] += - alpha[5] * mm / pow(qpm2, 3.0 / 2) * (y[1] - pm[1]);

    if(ms != 0) f[4] += - alpha[5] * ms / pow(qps2, 3.0 / 2) * (y[1] - ps[1]);

    //f[5] = pz'
    //------------------------------------
    f[5] = - alpha[1] * y[5] + alpha[14] * y[2];

    if(me != 0) f[5] += - alpha[5] * me / pow(qpe2, 3.0 / 2) * (y[2] - pe[2]);

    if(mm != 0) f[5] += - alpha[5] * mm / pow(qpm2, 3.0 / 2) * (y[2] - pm[2]);

    if(ms != 0) f[5] += - alpha[5] * ms / pow(qps2, 3.0 / 2) * (y[2] - ps[2]);


    return GSL_SUCCESS;
}



//=============================================================================================
//
// Integration with STM
//
//=============================================================================================
/**
 *  \brief Computes the second-order derivatives of the potential of one given primary
 **/
double Uij(const double y[], double pc[], double qpc2, double mc, double factor, int i, int j)
{
    if(mc != 0.0)
    {
        if(i == j) return factor * (3 * mc / pow(qpc2, 5.0 / 2) * pow(y[i] - pc[i], 2.0) - mc / pow(qpc2, 3.0 / 2));
        else return factor * 3 * mc / pow(qpc2, 5.0 / 2) * (y[i] - pc[i]) * (y[j] - pc[j]);
    }
    else return 0.0;
}


//----------------------------------------------------------------
// Ecliptic coordinates for JPL ephemerides (X, X')
//----------------------------------------------------------------
/**
 * \brief Vector field of the solar system (JPL ephemerides) in ecliptic coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO
 **/
int jpl_vf_var(double et, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);
    double Rb[6], lt, Rbe2;
    double b[6];

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];
    //f[3] = x", f[4] = y", f[5] = z", just init, the rest is done in the sequel
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;

    //For differential equations
    for(int i = 0; i < 6; i++) b[i] = 0.0;

    //-------------------------------------------------------------------------------
    //Loop on all the primaries
    //-------------------------------------------------------------------------------
    for(int i = 0; i < qbp->ss.maxBodies; i++)
    {
        //Retrieve posiiton from SPICE kernel
        spkez_c(qbp->ss.id[i], et, DEFFRAME,  "NONE", SSB, Rb, &lt);

        //Distance
        Rbe2 = (y[0] - Rb[0]) * (y[0] - Rb[0]) + (y[1] - Rb[1]) * (y[1] - Rb[1]) + (y[2] - Rb[2]) * (y[2] - Rb[2]);

        //Build the vector field
        f[3] += - qbp->ss.Gmi[i] / pow(Rbe2, 3.0 / 2) * (y[0] - Rb[0]);
        f[4] += - qbp->ss.Gmi[i] / pow(Rbe2, 3.0 / 2) * (y[1] - Rb[1]);
        f[5] += - qbp->ss.Gmi[i] / pow(Rbe2, 3.0 / 2) * (y[2] - Rb[2]);

        //For differential equations
        b[0] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 0, 0);
        b[1] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 0, 1);
        b[2] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 1, 1);
        b[3] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 0, 2);
        b[4] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 1, 2);
        b[5] += Uij(y, Rb, Rbe2, qbp->ss.Gmi[i], 1.0, 2, 2);
    }


    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);

    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);
    return 0;
}
/**
 * \brief Vector field of the solar system (JPL ephemerides) in earth-centered normalized
 *        coordinates with variationnal equations. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO.
 *        Note: for now the normalization is NOT taken into account because
 *        the stepper is going wrong when the state is normalized.
 **/
int jpl_vfn_var(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);
    double rb[6], rbe2, et, b[6];

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z", just init, the rest is done in the sequel
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;


    //------------------------------------------------------------------------------------
    //Epoch in seconds
    //------------------------------------------------------------------------------------
    //et = t*qbp->ss.n;
    et = t;

    //------------------------------------------------------------------------------------
    // Update the positions of all the bodies
    //------------------------------------------------------------------------------------
    double Rj[11][6];
    usspos(et, Rj, qbp);

    //------------------------------------------------------------------------------------
    //Acceleration of the Earth
    //------------------------------------------------------------------------------------
    double Ae[3], ae[3];
    acc_earth_from_vf(et, Ae, Rj, qbp);
    //for(int i = 0; i <3; i++) ae[i] = Ae[i]/(qbp->ss.a*qbp->ss.n*qbp->ss.n);
    for(int i = 0; i <3; i++) ae[i] = Ae[i];
    //Position of the Earth
    double Re[6], lt;
    spkez_c (EARTH, et, DEFFRAME, "NONE", SSB, Re, &lt);

    //------------------------------------------------------------------------------------
    //For differential equations
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) b[i] = 0.0;

    //-------------------------------------------------------------------------------
    //Loop on all the primaries
    //-------------------------------------------------------------------------------
    for(int i = 0; i < qbp->ss.maxBodies; i++)
    {
        //To normalized coordinates
        ecl2ecln(Rj[i], Re, rb);//, qbp->ss);

        //Distance
        rbe2 = (y[0] - rb[0]) * (y[0] - rb[0]) + (y[1] - rb[1]) * (y[1] - rb[1]) + (y[2] - rb[2]) * (y[2] - rb[2]);

        //Build the vector field
        f[3] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[0] - rb[0]);
        f[4] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[1] - rb[1]);
        f[5] += - qbp->ss.Gmi[i] / pow(rbe2, 3.0 / 2) * (y[2] - rb[2]);

        //For differential equations
        b[0] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 0, 0);
        b[1] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 0, 1);
        b[2] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 1, 1);
        b[3] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 0, 2);
        b[4] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 1, 2);
        b[5] += Uij(y, rb, rbe2, qbp->ss.Gmi[i], 1.0, 2, 2);

        //        //Build the vector field
        //        f[3] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[0] - rb[0]);
        //        f[4] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[1] - rb[1]);
        //        f[5] += - qbp->ss.mui[i] / pow(rbe2, 3.0 / 2) * (y[2] - rb[2]);
        //
        //        //For differential equations
        //        b[0] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 0, 0);
        //        b[1] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 0, 1);
        //        b[2] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 1, 1);
        //        b[3] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 0, 2);
        //        b[4] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 1, 2);
        //        b[5] += Uij(y, rb, rbe2, qbp->ss.mui[i], 1.0, 2, 2);
    }

    //------------------------------------------------------------------------------------
    //Last contribution: acceleration of the Earth
    //------------------------------------------------------------------------------------
    f[3] += -ae[0];
    f[4] += -ae[1];
    f[5] += -ae[2];


    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);

    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);
    return 0;
}


//----------------------------------------------------------------
// Synodic coordinates for JPL ephemerides (x, x')
//----------------------------------------------------------------
/**
 * \brief Vector field of the solar system (JPL ephemerides) in synodical coordinates. Included:
 *        SUN, MERCURY, VENUS, EARTH, MOON, MARS, JUPITER, SATURN, URANUS, NEPTUNE, PLUTO (11 bodies)
 **/
int jpl_vf_syn_var(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);
    double q[6];


    //------------------------------------------------------------------------------------
    // Misc parameters
    //------------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp   = (QBCP_L*) params_void;
    int coord_eph = qbp->ss.coord_eph;        //Current coord_eph
    double a      = qbp->ss.a;                //mean semi-major axis (in km)
    double et0    = qbp->ss.et0;              //starting epoch
    double t0     = qbp->ss.t0;               //starting time
    double mu1    = qbp->ss.mu1;              //mass ratio
    double mu2    = qbp->ss.mu2;              //1-mass ratio
    double nsys   = qbp->ss.n;                //mean motion (in s)
    double pos1   = qbp->ss.pos1;             //indix of m1
    double pos2   = qbp->ss.pos2;             //indix of m2
    double et     = et0 + (t - t0) / qbp->ss.n; //Current time in seconds

    //------------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z'
    //------------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //------------------------------------------------------------------------------------
    // Update the positions of all the bodies
    //------------------------------------------------------------------------------------
    double Rj[11][6];
    usspos(et, Rj, qbp);

    //------------------------------------------------------------------------------------
    //Position, velocity, acceleration, and jerk of the primaries
    //------------------------------------------------------------------------------------
    double lt;
    double R1[6];  //state of m1
    double R2[6];  //state of m2
    double A1[3];  //acceleration of m1
    double A2[3];  //acceleration of m2
    double J1[3];  //jerk of m1
    double J2[3];  //jerk of m2

    //--------------------------------------------
    //Position, velocity of the primaries, in dimensionnalized units
    //--------------------------------------------
    spkez_c (m1name(coord_eph), et, DEFFRAME, "NONE", SSB, R1, &lt);
    spkez_c (m2name(coord_eph), et, DEFFRAME, "NONE", SSB, R2, &lt);

    //--------------------------------------------
    //Acceleration of the primaries, in dimensionnalized units
    //--------------------------------------------
    double DR[3], DR2;
    //First possibility: from ephemerides
    //    double err;
    //    for(int i = 0; i < 3; i++)  A1[i] = dfridr(xm1ecl, et, i+3, 1.0, &err, coord_eph);
    //    for(int i = 0; i < 3; i++)  A2[i] = dfridr(xm2ecl, et, i+3, 1.0, &err, coord_eph);
    //Second possibility: from the vector field
    acc_from_vf(et, pos1, pos2, A1, A2, R1, R2, Rj, qbp);

    //--------------------------------------------
    //Jerk of the primaries, in dimensionnalized units
    //--------------------------------------------
    //First possibility: from ephemerides
    //    for(int i = 0; i < 3; i++)  J1[i] = dfridr(dxm1ecl, et, i+3, 1.0, &err, coord_eph);
    //    for(int i = 0; i < 3; i++)  J2[i] = dfridr(dxm2ecl, et, i+3, 1.0, &err, coord_eph);
    //Second possibility: from vector field
    jerk_from_vf(et, pos1, pos2, J1, J2, R1, R2, Rj, qbp);


    //------------------------------------------------------------------------------------
    //Relative position, velocity, acceleration, and jerk
    //------------------------------------------------------------------------------------
    double r21[3]; //relative position
    double v21[3]; //relative velocity
    double a21[3]; //relative acceleration
    double j21[3]; //relative jerk
    double B[3];   //position of the barycenter vector
    double Bpp[3]; //acceleration of the barycenter vector

    for(int i = 0; i < 3; i++)
    {
        //Difference between the two position vectors
        r21[i]  = R1[i] - R2[i];
        //Derivative of the difference wrt to dimentionalized time
        v21[i]  = R1[i + 3] - R2[i + 3];
        //Difference between the two accelerations
        a21[i]  = A1[i] - A2[i];
        //Difference between the two jerks
        j21[i]  = J1[i] - J2[i];
        //Position of the barycenter
        B[i]    = (mu1 * R1[i] + mu2 * R2[i]) / (mu1 + mu2);
        //Acceleration of the barycenter vector
        Bpp[i]  = (mu1 * A1[i] + mu2 * A2[i]) / (mu1 + mu2);
    }

    //------------------------------------------------------------------------------------
    //Orthonormal basis (e1, e2, e3) + h and k factors and their derivatives
    //------------------------------------------------------------------------------------
    double e1[3], e2[3], e3[3], e4[3], k, kp, kpp, h, hp, hpp;
    double r21_v21[3], r21_a21[3], r21_j21[3], v21_a21[3];
    //Building the cross products
    init_cp(r21, v21, a21, j21, r21_v21, r21_a21, r21_j21, v21_a21);
    //Building h, h', h"
    init_h(r21_v21, r21_a21, r21_j21, v21_a21, &h, &hp, &hpp, e4);
    //Building k, k', k''
    init_k(r21, v21, a21, &k, &kp, &kpp);
    //Building e1, e2, e3
    init_ei(r21, v21, e1, e2, e3);

    //------------------------------------------------------------------------------------
    //Orthonormal basis derivatives (e1', e2', e3')
    //------------------------------------------------------------------------------------
    double de1[3], de2[3], de3[3];
    init_dei(r21, v21, a21, k, kp, h, hp, e1, e3, de1, de2, de3, e4);

    //------------------------------------------------------------------------------------
    //Orthonormal basis double derivatives (e1", e2", e3")
    //------------------------------------------------------------------------------------
    double dde1[3], dde2[3], dde3[3];
    init_ddei(r21, v21, a21, r21_v21, r21_a21, r21_j21, v21_a21,
              k, kp, kpp, h, hp, hpp, e1, e3, de1, de3,
              dde1, dde2, dde3, e4);

    //------------------------------------------------------------------------------------
    // b coefficients. Note that a factor 1/n^2 is added to account for
    // the time normalization
    //------------------------------------------------------------------------------------
    double b[13];
    //b[0] = - B" . e1 /(k * n^2);
    b[0] = -vdot_c(Bpp, e1) / (k * nsys * nsys);
    //b[1] = - B" . e2 /(k * n^2);
    b[1] = -vdot_c(Bpp, e2) / (k * nsys * nsys);
    //b[2] = - B" . e3 /(k * n^2);
    b[2] = -vdot_c(Bpp, e3) / (k * nsys * nsys);
    //b[3] = - 2*k'/(n*k)
    b[3] = -2 * kp / (k * nsys);
    //b[4] = 2/n e2 . e1'
    b[4] = +2 * vdot_c(e2, de1) / nsys;
    //b[5] = 2/n e3 . e2'
    b[5] = +2 * vdot_c(e3, de2) / nsys;
    //b[6] = -1/n^2*(k"/k - e1' . e1')
    b[6] = - (kpp / k - vdot_c(de1, de1)) / (nsys * nsys);
    //b[7] = +1/n^2*(e1'.e3')
    b[7] = vdot_c(de1, de3) / (nsys * nsys);
    //b[8] = +1/n^2*(2*k'/k e2 . e1' + e2 . e1")
    b[8] = (2 * kp / k * vdot_c(e2, de1) + vdot_c(e2, dde1)) / (nsys * nsys);
    //b[9] = -1/n^2*(k"/k - e2' . e2')
    b[9] = - (kpp / k - vdot_c(de2, de2)) / (nsys * nsys);
    //b[10] = +1/n^2*(2*k'/k e3 . e3' + e3 . e2")
    b[10] = (2 * kp / k * vdot_c(e3, de2) + vdot_c(e3, dde2)) / (nsys * nsys);
    //b[11] = -1/n^2*(k"/k - e3' . e3')
    b[11] = - (kpp / k - vdot_c(de3, de3)) / (nsys * nsys);
    //b[12] = a^3/k^3
    b[12] = a * a * a / (k * k * k);

    //////////////////////////////////////////////////////////////////////////////////////
    //
    //€€TODO:
    //read the double derivatives of the ei + hpp and kpp
    //set a latex document that sum up these calculus (see p44 of Dei Tos).
    //
    //////////////////////////////////////////////////////////////////////////////////////


    //------------------------------------------------------------------------------------
    //For Reduced positions of the primaries
    //------------------------------------------------------------------------------------
    double C[3][3], vin[3], rb[3];
    for(int i = 0; i < 3; i++) C[i][0] = e1[i];
    for(int i = 0; i < 3; i++) C[i][1] = e2[i];
    for(int i = 0; i < 3; i++) C[i][2] = e3[i];



    //------------------------------------------------------------------------------------
    //Phase space derivatives: x", y", z"
    //------------------------------------------------------------------------------------
    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = b[0] + b[3] * y[3] + b[4] * y[4]               + b[6] * y[0] + b[8] * y[1] + b[7]  * y[2];
    f[4] = b[1] - b[4] * y[3] + b[3] * y[4] + b[5] * y[5] - b[8] * y[0] + b[9] * y[1] + b[10] * y[2];
    f[5] = b[2]               - b[5] * y[4] + b[3] * y[5] + b[7] * y[0] - b[10]* y[1] + b[11] * y[2];

    //For differential equations
    for(int i = 0; i < 6; i++) q[i] = 0.0;


    //Loop on all the primaries
    for(int p = 0; p < qbp->ss.maxBodies; p++)
    {
        //Reduced position of the primaries
        for(int i = 0; i < 3; i++) vin[i] = Rj[p][i];
        ecl2synpos(vin, rb, B, C, k);

        //DR = (x, y, z)^T - rj
        for(int i = 0; i < 3; i++) DR[i]  = y[i] - rb[i];

        //DR2 = DR . DR
        DR2 = vdot_c(DR, DR);

        //Build the vector field
        f[3] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[0];
        f[4] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[1];
        f[5] += - b[12] * qbp->ss.mui[p] / pow(DR2, 3.0 / 2) * DR[2];


        //For differential equations
        q[0] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 0, 0);
        q[1] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 0, 1);
        q[2] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 1, 1);
        q[3] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 0, 2);
        q[4] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 1, 2);
        q[5] += Uij(y, rb, DR2, qbp->ss.mui[p], b[12], 2, 2);
    }

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);

    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  q[0]+b[6]);
    gsl_matrix_set(Q, 3, 1,  q[1]+b[8]);
    gsl_matrix_set(Q, 3, 2,  q[3]+b[7]);
    gsl_matrix_set(Q, 3, 3,  +b[3]);
    gsl_matrix_set(Q, 3, 4,  +b[4]);
    gsl_matrix_set(Q, 3, 5,  +0.0);


    gsl_matrix_set(Q, 4, 0,  q[1]-b[8]);
    gsl_matrix_set(Q, 4, 1,  q[2]+b[9]);
    gsl_matrix_set(Q, 4, 2,  q[4]+b[10]);
    gsl_matrix_set(Q, 4, 3,  -b[4]);
    gsl_matrix_set(Q, 4, 4,  +b[3]);
    gsl_matrix_set(Q, 4, 5,  +b[5]);

    gsl_matrix_set(Q, 5, 0,  q[3]+b[7]);
    gsl_matrix_set(Q, 5, 1,  q[4]-b[10]);
    gsl_matrix_set(Q, 5, 2,  q[5]+b[11]);
    gsl_matrix_set(Q, 5, 3,  +0.0);
    gsl_matrix_set(Q, 5, 4,  -b[5]);
    gsl_matrix_set(Q, 5, 5,  +b[3]);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);
    return 0;
}


//----------------------------------------------------------------
// Inertial SEM coordinates, (X, V)
//----------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in SEM inertial coordinates
 **/
int qbcp_vf_insem_var(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);
    double b[6];


    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    double ms = qbp->us_sem.ms;
    double me = qbp->us_sem.me;
    double mm = qbp->us_sem.mm;
    double n  = qbp->us_sem.n;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;


    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @t, in SEM coordinates
    //-------------------------------------------------------------------------------
    double Ps_syn[3];
    evaluateCoef(Ps_syn, t, n, qbp->nf, qbp->cs_sem.Ps, 3);
    double Pe_syn[3];
    evaluateCoef(Pe_syn, t, n, qbp->nf, qbp->cs_sem.Pe, 3);
    double Pm_syn[3];
    evaluateCoef(Pm_syn, t, n, qbp->nf, qbp->cs_sem.Pm, 3);

    //-------------------------------------------------------------------------------
    //Back to  INSEM coordinates
    //-------------------------------------------------------------------------------
    double Ps[6];
    SEMtoIN(t, Ps_syn, Ps, qbp);
    double Pe[6];
    SEMtoIN(t, Pe_syn, Pe, qbp);
    double Pm[6];
    SEMtoIN(t, Pm_syn, Pm, qbp);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);


    //-------------------------------------------------------------------------------
    // Update the vector field
    //-------------------------------------------------------------------------------
    if(me != 0) f[3] += - me / pow(qPe2, 3.0 / 2) * (y[0] - Pe[0]);
    if(mm != 0) f[3] += - mm / pow(qPm2, 3.0 / 2) * (y[0] - Pm[0]);
    if(ms != 0) f[3] += - ms / pow(qPs2, 3.0 / 2) * (y[0] - Ps[0]);

    if(me != 0) f[4] += - me / pow(qPe2, 3.0 / 2) * (y[1] - Pe[1]);
    if(mm != 0) f[4] += - mm / pow(qPm2, 3.0 / 2) * (y[1] - Pm[1]);
    if(ms != 0) f[4] += - ms / pow(qPs2, 3.0 / 2) * (y[1] - Ps[1]);

    if(me != 0) f[5] += - me / pow(qPe2, 3.0 / 2) * (y[2] - Pe[2]);
    if(mm != 0) f[5] += - mm / pow(qPm2, 3.0 / 2) * (y[2] - Pm[2]);
    if(ms != 0) f[5] += - ms / pow(qPs2, 3.0 / 2) * (y[2] - Ps[2]);

    //-------------------------------------------------------------------------------
    //For differential equations
    //-------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) b[i] = 0.0;


    //For differential equations
    b[0] += Uij(y, Pe, qPe2, me, 1.0, 0, 0) + Uij(y, Pm, qPm2, mm, 1.0, 0, 0) + Uij(y, Ps, qPs2, ms, 1.0, 0, 0);
    b[1] += Uij(y, Pe, qPe2, me, 1.0, 0, 1) + Uij(y, Pm, qPm2, mm, 1.0, 0, 1) + Uij(y, Ps, qPs2, ms, 1.0, 0, 1);
    b[2] += Uij(y, Pe, qPe2, me, 1.0, 1, 1) + Uij(y, Pm, qPm2, mm, 1.0, 1, 1) + Uij(y, Ps, qPs2, ms, 1.0, 1, 1);
    b[3] += Uij(y, Pe, qPe2, me, 1.0, 0, 2) + Uij(y, Pm, qPm2, mm, 1.0, 0, 2) + Uij(y, Ps, qPs2, ms, 1.0, 0, 2);
    b[4] += Uij(y, Pe, qPe2, me, 1.0, 1, 2) + Uij(y, Pm, qPm2, mm, 1.0, 1, 2) + Uij(y, Ps, qPs2, ms, 1.0, 1, 2);
    b[5] += Uij(y, Pe, qPe2, me, 1.0, 2, 2) + Uij(y, Pm, qPm2, mm, 1.0, 2, 2) + Uij(y, Ps, qPs2, ms, 1.0, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);

    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);
    return 0;
}

//----------------------------------------------------------------
// Inertial EM coordinates, (X, V)
//----------------------------------------------------------------
/**
 * \brief Vector field of the QBCP in EM inertial coordinates
 **/
int qbcp_vf_inem_var(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);
    double b[6];


    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    double ms = qbp->us_em.ms;
    double me = qbp->us_em.me;
    double mm = qbp->us_em.mm;
    double n  = qbp->us_em.n;

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', x", y", z"
    //-------------------------------------------------------------------------------
    //f[0] = x', f[1] = y', f[2] = z'
    //------------------------------------
    f[0] = y[3];
    f[1] = y[4];
    f[2] = y[5];

    //f[3] = x", f[4] = y", f[5] = z"
    //------------------------------------
    f[3] = 0.0;
    f[4] = 0.0;
    f[5] = 0.0;


    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @t, in SEM coordinates
    //-------------------------------------------------------------------------------
    double Ps_syn[3];
    evaluateCoef(Ps_syn, t, n, qbp->nf, qbp->cs_em.Ps, 3);
    double Pe_syn[3];
    evaluateCoef(Pe_syn, t, n, qbp->nf, qbp->cs_em.Pe, 3);
    double Pm_syn[3];
    evaluateCoef(Pm_syn, t, n, qbp->nf, qbp->cs_em.Pm, 3);

    //-------------------------------------------------------------------------------
    //Back to  INSEM coordinates
    //-------------------------------------------------------------------------------
    double Ps[6];
    EMvtoIN(t, Ps_syn, Ps, qbp);
    double Pe[6];
    EMvtoIN(t, Pe_syn, Pe, qbp);
    double Pm[6];
    EMvtoIN(t, Pm_syn, Pm, qbp);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);


    //-------------------------------------------------------------------------------
    // Update the vector field
    //-------------------------------------------------------------------------------
    if(me != 0) f[3] += - me / pow(qPe2, 3.0 / 2) * (y[0] - Pe[0]);
    if(mm != 0) f[3] += - mm / pow(qPm2, 3.0 / 2) * (y[0] - Pm[0]);
    if(ms != 0) f[3] += - ms / pow(qPs2, 3.0 / 2) * (y[0] - Ps[0]);

    if(me != 0) f[4] += - me / pow(qPe2, 3.0 / 2) * (y[1] - Pe[1]);
    if(mm != 0) f[4] += - mm / pow(qPm2, 3.0 / 2) * (y[1] - Pm[1]);
    if(ms != 0) f[4] += - ms / pow(qPs2, 3.0 / 2) * (y[1] - Ps[1]);

    if(me != 0) f[5] += - me / pow(qPe2, 3.0 / 2) * (y[2] - Pe[2]);
    if(mm != 0) f[5] += - mm / pow(qPm2, 3.0 / 2) * (y[2] - Pm[2]);
    if(ms != 0) f[5] += - ms / pow(qPs2, 3.0 / 2) * (y[2] - Ps[2]);

    //-------------------------------------------------------------------------------
    //For differential equations
    //-------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++) b[i] = 0.0;


    //For differential equations
    b[0] += Uij(y, Pe, qPe2, me, 1.0, 0, 0) + Uij(y, Pm, qPm2, mm, 1.0, 0, 0) + Uij(y, Ps, qPs2, ms, 1.0, 0, 0);
    b[1] += Uij(y, Pe, qPe2, me, 1.0, 0, 1) + Uij(y, Pm, qPm2, mm, 1.0, 0, 1) + Uij(y, Ps, qPs2, ms, 1.0, 0, 1);
    b[2] += Uij(y, Pe, qPe2, me, 1.0, 1, 1) + Uij(y, Pm, qPm2, mm, 1.0, 1, 1) + Uij(y, Ps, qPs2, ms, 1.0, 1, 1);
    b[3] += Uij(y, Pe, qPe2, me, 1.0, 0, 2) + Uij(y, Pm, qPm2, mm, 1.0, 0, 2) + Uij(y, Ps, qPs2, ms, 1.0, 0, 2);
    b[4] += Uij(y, Pe, qPe2, me, 1.0, 1, 2) + Uij(y, Pm, qPm2, mm, 1.0, 1, 2) + Uij(y, Ps, qPs2, ms, 1.0, 1, 2);
    b[5] += Uij(y, Pe, qPe2, me, 1.0, 2, 2) + Uij(y, Pm, qPm2, mm, 1.0, 2, 2) + Uij(y, Ps, qPs2, ms, 1.0, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);

    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5]);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);
    return 0;
}



//----------------------------------------------------------------
// NC coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
void set_vareq_matrix(gsl_matrix* Q, double b[], double alpha[])
{
    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set(Q, 0, 0,  alpha[1]);
    gsl_matrix_set(Q, 0, 1,  alpha[2]);
    gsl_matrix_set(Q, 0, 3,  alpha[0]);
    gsl_matrix_set(Q, 1, 0, -alpha[2]);
    gsl_matrix_set(Q, 1, 1,  alpha[1]);
    gsl_matrix_set(Q, 1, 4,  alpha[0]);
    gsl_matrix_set(Q, 2, 2,  alpha[1]);
    gsl_matrix_set(Q, 2, 5,  alpha[0]);
    gsl_matrix_set(Q, 3, 0,  b[0] + alpha[14]);
    gsl_matrix_set(Q, 3, 1,  b[1]);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1]);
    gsl_matrix_set(Q, 4, 1,  b[2] + alpha[14]);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5] + alpha[14]);
    gsl_matrix_set(Q, 3, 3, -alpha[1]);
    gsl_matrix_set(Q, 3, 4,  alpha[2]);
    gsl_matrix_set(Q, 4, 3, -alpha[2]);
    gsl_matrix_set(Q, 4, 4, -alpha[1]);
    gsl_matrix_set(Q, 5, 5, -alpha[1]);
}

/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm(const double y[], gsl_matrix* Q, double alpha[],
            double ps[], double pe[], double pm[],
            double qps2, double qpe2, double qpm2,
            double ms, double me, double mm,
            double gamma)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[5] / pow(gamma, 3.0);

    b[0] =     Uij(y, pe, qpe2, me, factor, 0, 0)
               + Uij(y, pm, qpm2, mm, factor, 0, 0)
               + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates.
 **/
int qbfbp_vfn_varnonlin(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp  = (QBCP_L*) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, qbp->nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, qbp->nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, qbp->nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0] - pe[0]) * (y[0] - pe[0]) + (y[1] - pe[1]) * (y[1] - pe[1]) + (y[2] - pe[2]) * (y[2] - pe[2]);
    double qps2 = (y[0] - ps[0]) * (y[0] - ps[0]) + (y[1] - ps[1]) * (y[1] - ps[1]) + (y[2] - ps[2]) * (y[2] - ps[2]);
    double qpm2 = (y[0] - pm[0]) * (y[0] - pm[0]) + (y[1] - pm[1]) * (y[1] - pm[1]) + (y[2] - pm[2]) * (y[2] - pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state(y, f, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);
    vfn_stm(y, Q, alpha, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------
// NC coordinates, (X, V)
//----------------------------------------------------------------
/**
 *  \brief Build the Variational Equation Matrix Q from the arrays b and alpha
 **/
void set_vareq_matrix_xv(gsl_matrix* Q, double b[], double alpha[], double alphad[])
{
    //-------------------------------------------------------------------------------
    //Intermediate coefficients
    //-------------------------------------------------------------------------------
    double alpha16 = alphad[1] - alphad[0] * alpha[1] / alpha[0] + alpha[1] * alpha[1] + alpha[2] * alpha[2]; //
    double alpha17 = alphad[2] - alphad[0] * alpha[2] / alpha[0]; //
    double alpha18 = alphad[0] / alpha[0]; //

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    gsl_matrix_set(Q, 0, 3,  1.0);
    gsl_matrix_set(Q, 1, 4,  1.0);
    gsl_matrix_set(Q, 2, 5,  1.0);

    gsl_matrix_set(Q, 3, 0,  b[0] + alpha16);
    gsl_matrix_set(Q, 3, 1,  b[1] + alpha17);
    gsl_matrix_set(Q, 3, 2,  b[3]);
    gsl_matrix_set(Q, 4, 0,  b[1] - alpha17);
    gsl_matrix_set(Q, 4, 1,  b[2] + alpha16);
    gsl_matrix_set(Q, 4, 2,  b[4]);
    gsl_matrix_set(Q, 5, 0,  b[3]);
    gsl_matrix_set(Q, 5, 1,  b[4]);
    gsl_matrix_set(Q, 5, 2,  b[5] + alpha16 - alpha[2]*alpha[2]);

    gsl_matrix_set(Q, 3, 3,  alpha18);
    gsl_matrix_set(Q, 3, 4,  2 * alpha[2]);
    gsl_matrix_set(Q, 4, 3, -2 * alpha[2]);
    gsl_matrix_set(Q, 4, 4,  alpha18);
    gsl_matrix_set(Q, 5, 5,  alpha18);
}



/**
 *  \brief Update the Normalized-Centered variational equation matrix
 **/
int vfn_stm_xv(const double y[], gsl_matrix* Q, double alpha[], double alphad[],
               double ps[], double pe[], double pm[],
               double qps2, double qpe2, double qpm2,
               double ms, double me, double mm,
               double gamma)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[0] * alpha[5] / pow(gamma, 3.0);

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
             + Uij(y, pm, qpm2, mm, factor, 0, 0)
             + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix_xv(Q, b, alpha, alphad);

    return GSL_SUCCESS;
}

/**
 *  \brief Vector field of the QBCP with EM units and Normalized-Centered coordinates. The state is here (x, xdot), and not (x, px)
 **/
int qbfbp_vfn_varnonlin_xv(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* Q    = gsl_matrix_calloc(6, 6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp  = (QBCP_L*) params_void;
    int noc      = qbp->numberOfCoefs;
    double ms    = qbp->us.ms;
    double me    = qbp->us.me;
    double mm    = qbp->us.mm;
    double n     = qbp->us.n;
    double gamma = qbp->cs.gamma;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas and their derivatives @ t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);
    double alphad[3];
    evaluateCoefDerivatives(alphad, t, n, qbp->nf, qbp->cs.coeffs, 3);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double ps[3];
    evaluateCoef(ps, t, n, qbp->nf, qbp->cs.ps, 3);
    double pe[3];
    evaluateCoef(pe, t, n, qbp->nf, qbp->cs.pe, 3);
    double pm[3];
    evaluateCoef(pm, t, n, qbp->nf, qbp->cs.pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qpe2 = (y[0] - pe[0]) * (y[0] - pe[0]) + (y[1] - pe[1]) * (y[1] - pe[1]) + (y[2] - pe[2]) * (y[2] - pe[2]);
    double qps2 = (y[0] - ps[0]) * (y[0] - ps[0]) + (y[1] - ps[1]) * (y[1] - ps[1]) + (y[2] - ps[2]) * (y[2] - ps[2]);
    double qpm2 = (y[0] - pm[0]) * (y[0] - pm[0]) + (y[1] - pm[1]) * (y[1] - pm[1]) + (y[2] - pm[2]) * (y[2] - pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vfn_state_xv(y, f, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(Q);
    vfn_stm_xv(y, Q, alpha, alphad, ps, pe, pm, qps2, qpe2, qpm2, ms, me, mm, gamma);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = Q * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Q, STM, 0.0, STMd);


    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(Q);

    return GSL_SUCCESS;
}


//----------------------------------------------------------------
// SYS coordinates, (X, P)
//----------------------------------------------------------------
/**
 *  \brief Update the variational equation matrix in system coordinates (non-normalized)
 **/
int vf_stm(const double y[], gsl_matrix* Q, double alpha[],
           double ps[], double pe[], double pm[],
           double qps2, double qpe2, double qpm2,
           double ms, double me, double mm)
{
    double b[6];

    //-------------------------------------------------------------------------------
    // Variational equations, nonlinear case:
    // The nonlinear terms of the equations of motion
    // are the derivatives of the potential of the primaries U: Ux, Uy, Uz.
    // The variational equations contain the second derivatives of this potential
    //
    // b1 = dUxx, b2 = dUxy, b4 = dUxz
    // b2 = dUyx, b3 = dUyy, b5 = dUyz
    // b4 = dUzx, b5 = dUzy, b6 = dUzz
    //
    //-------------------------------------------------------------------------------
    double factor = alpha[5];

    b[0] =   Uij(y, pe, qpe2, me, factor, 0, 0)
             + Uij(y, pm, qpm2, mm, factor, 0, 0)
             + Uij(y, ps, qps2, ms, factor, 0, 0);

    b[1] =   Uij(y, pe, qpe2, me, factor, 0, 1)
             + Uij(y, pm, qpm2, mm, factor, 0, 1)
             + Uij(y, ps, qps2, ms, factor, 0, 1);

    b[2] =   Uij(y, pe, qpe2, me, factor, 1, 1)
             + Uij(y, pm, qpm2, mm, factor, 1, 1)
             + Uij(y, ps, qps2, ms, factor, 1, 1);

    b[3] =   Uij(y, pe, qpe2, me, factor, 0, 2)
             + Uij(y, pm, qpm2, mm, factor, 0, 2)
             + Uij(y, ps, qps2, ms, factor, 0, 2);

    b[4] =   Uij(y, pe, qpe2, me, factor, 1, 2)
             + Uij(y, pm, qpm2, mm, factor, 1, 2)
             + Uij(y, ps, qps2, ms, factor, 1, 2);

    b[5] =   Uij(y, pe, qpe2, me, factor, 2, 2)
             + Uij(y, pm, qpm2, mm, factor, 2, 2)
             + Uij(y, ps, qps2, ms, factor, 2, 2);

    //-------------------------------------------------------------------------------
    // Build the variational equation matrix Q such that dot(STM) = Q * STM
    //-------------------------------------------------------------------------------
    set_vareq_matrix(Q, b, alpha);

    return GSL_SUCCESS;
}


/**
 *  \brief Vector field of the QBCP with EM units and EM reference frame. Variational equations included.
 **/
int qbfbp_vf_varnonlin(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Memory allocation
    //-------------------------------------------------------------------------------
    gsl_matrix* STM  = gsl_matrix_calloc(6, 6);
    gsl_matrix* STMd = gsl_matrix_calloc(6, 6);
    gsl_matrix* B    = gsl_matrix_calloc(6, 6);

    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;
    int noc   = qbp->numberOfCoefs;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);

    //-------------------------------------------------------------------------------
    //Phase space derivatives: x', y', z', px', py', pz'
    //-------------------------------------------------------------------------------
    vf_state(y, f, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //-------------------------------------------------------------------------------
    //Variationnal equations, nonlinear case
    //-------------------------------------------------------------------------------
    gsl_matrix_set_zero(B);
    vf_stm(y, B, alpha, Ps, Pe, Pm, qPs2, qPe2, qPm2, ms, me, mm);

    //-------------------------------------------------------------------------------
    //STM update & dot(STM) computation
    //-------------------------------------------------------------------------------
    //STM update
    gslc_vectorToMatrix(STM, y, 6, 6, 6);

    //dot(STM) = B * STM
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, B, STM, 0.0, STMd);

    //dot(STM) stored in f
    gslc_matrixToVector(f, STMd, 6, 6 , 6);

    //Memory release
    gsl_matrix_free(STMd);
    gsl_matrix_free(STM);
    gsl_matrix_free(B);

    return GSL_SUCCESS;
}


//========================================================================================
//
//          Reduced vector field
//
//========================================================================================
/**
 *  \brief Computes the reduce vector field dot(s) = f(s,t)
 *  \param t the current time
 *  \param y the complex state s (CCM), with real and imag parts separated (CCM8 form)
 *  \param f the complex reduced vector field f(s,t), with real and imag parts separated (CCM8 form)
 *  \param params_void a (void) pointer towards a RVF structure that contains:
 *           (1) the Fourier-Taylor expansions of the rvf
 *           (2) the order of the evaluation
 *           (3) the pulsation n of the system
 *           (4) an OFS temp variable.
 *
 *     Note that the use of void *params_void instead of the direct use of an RVF structure is to comply with GSL ODE integration tools,
 *     that require a vector field of the exact form : int vf(double t, const double y[], double f[], void *params_void)
 *
 **/
int qbfbp_fh(double t, const double y[], double f[], void* params_void)
{
    //-------------------------------
    //Initialization
    //-------------------------------
    RVF* rvf = (RVF*) params_void;

    //-------------------------------
    //Evaluation of the reduced vector field at order rvf->order
    //-------------------------------
    CCM8toRVF8(y, t, rvf->n, rvf->order, rvf->ofs_order, rvf->reduced_nv, *rvf->fh, *rvf->ofs, f);

    return GSL_SUCCESS;
}


//========================================================================================
//
// Hamiltonians
//
//========================================================================================
/**
 *  \brief Hamiltonian of the QBCP with SYS units and SYS coordinates. Note that alpha[14] (alpha15) is zero for the QBCP
 **/
double qbfbp_H(double t, const double y[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us.ms;
    double me = qbp->us.me;
    double mm = qbp->us.mm;
    double n  = qbp->us.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);

    //-------------------------------------------------------------------------------
    // Hamiltonian
    //-------------------------------------------------------------------------------
    double H = 0.5 * alpha[0] * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) + alpha[1] * (y[3] * y[0] + y[4] * y[1] + y[5] * y[2])
               + alpha[2] * (y[3] * y[1] - y[4] * y[0])
               + alpha[3] * y[0] + alpha[4] * y[1]
               - 0.5 * alpha[14] * (y[0] * y[0] + y[1] * y[1] + y[2] * y[2])
               - alpha[5] * ( me / pow(qPe2, 1.0 / 2)
                              + mm / pow(qPm2, 1.0 / 2)
                              + ms / pow(qPs2, 1.0 / 2) );

    return H;
}


/**
 *  \brief Hamiltonian of the QBCP with SEM units and SEM coordinates
 **/
double qbfbp_H_SEM(double t, const double y[], void* params_void)
{
    //-------------------------------------------------------------------------------
    // Misc parameters
    //-------------------------------------------------------------------------------
    //Retrieving the parameters
    QBCP_L* qbp = (QBCP_L*) params_void;
    int noc   = qbp->numberOfCoefs;
    double ms = qbp->us_sem.ms;
    double me = qbp->us_sem.me;
    double mm = qbp->us_sem.mm;
    double n  = qbp->us_sem.n;

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t
    //1, 3, 4, 6, 7 even
    //2, 5, 8 odd
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t, n, qbp->nf, qbp->cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the primaries positions @ t
    //-------------------------------------------------------------------------------
    double Ps[3];
    evaluateCoef(Ps, t, n, qbp->nf, qbp->cs_sem.Ps, 3);
    double Pe[3];
    evaluateCoef(Pe, t, n, qbp->nf, qbp->cs_sem.Pe, 3);
    double Pm[3];
    evaluateCoef(Pm, t, n, qbp->nf, qbp->cs_sem.Pm, 3);

    //-------------------------------------------------------------------------------
    // Distances to 2nd power
    //-------------------------------------------------------------------------------
    double qPe2 = (y[0] - Pe[0]) * (y[0] - Pe[0]) + (y[1] - Pe[1]) * (y[1] - Pe[1]) + (y[2] - Pe[2]) * (y[2] - Pe[2]);
    double qPs2 = (y[0] - Ps[0]) * (y[0] - Ps[0]) + (y[1] - Ps[1]) * (y[1] - Ps[1]) + (y[2] - Ps[2]) * (y[2] - Ps[2]);
    double qPm2 = (y[0] - Pm[0]) * (y[0] - Pm[0]) + (y[1] - Pm[1]) * (y[1] - Pm[1]) + (y[2] - Pm[2]) * (y[2] - Pm[2]);

    //-------------------------------------------------------------------------------
    // Hamiltonian
    //-------------------------------------------------------------------------------
    double H = 0.5 * alpha[0] * (y[3] * y[3] + y[4] * y[4] + y[5] * y[5]) + alpha[1] * (y[3] * y[0] + y[4] * y[1] + y[5] * y[2])
               + alpha[2] * (y[3] * y[1] - y[4] * y[0])
               + alpha[3] * y[0] + alpha[4] * y[1]
               - 0.5 * alpha[14] * (y[0] * y[0] + y[1] * y[1] + y[2] * y[2])
               - alpha[5] * ( me / pow(qPe2, 1.0 / 2)
                              + mm / pow(qPm2, 1.0 / 2)
                              + ms / pow(qPs2, 1.0 / 2) );

    return H;
}



//========================================================================================
//
// Testing the QBTBP
//
//========================================================================================
//-----------------------------------------------------------------------------
// Integrating the QBTBP
//-----------------------------------------------------------------------------
/**
 *  \brief Derivatives of the QBTBP. To plug into GSL integrator.
 *
 * Note: the use of gsl livrary forces us to use double variables
 * As a consequence, in the case of z and Z, we need to use real and imaginary parts
 * as separate variables
 */
int qbtbp_derivatives(double t, const double y[], double f[], void* params)
{
    //Parameters for the qbtbp
    double mu = * (double* ) params;
    double ms = * ((double* ) (params) + 1);

    //Reconstruction of z and Z
    cdouble z = y[0] + I * y[1];
    cdouble Z = y[2] + I * y[3];

    cdouble temp1 = Z - mu * z;
    temp1 = temp1 / pow(cabs(temp1), 3.0);

    cdouble temp2 = Z + (1 - mu) * z;
    temp2 = temp2 / pow(cabs(temp2), 3.0);

    cdouble zdd = +0.0 * I - z / pow(cabs(z), 3.0) + ms * (temp1 - temp2);
    cdouble Zdd = -(1 + ms) * (mu * temp2 + (1 - mu) * temp1);

    //-------------------------------------------------------------------------------
    //Phase space derivatives
    //-------------------------------------------------------------------------------
    f[0] = y[4];
    f[1] = y[5];
    f[2] = y[6];
    f[3] = y[7];
    f[4] = creal(zdd);
    f[5] = cimag(zdd);
    f[6] = creal(Zdd);
    f[7] = cimag(Zdd);

    return GSL_SUCCESS;
}

/**
 * \brief Computes the semi-analytical position of the Earth, the Moon, and the Sun from the coefficients alpha and delta of the EM and SEM vector field.
 **/
void primary_analytical_position(int prim, double tem, double yEM_to_IN[], double ySEM_to_IN[], QBCP_L& qbcp_l)
{

    //====================================================================
    // 1. Init
    //====================================================================
    int noc     = qbcp_l.numberOfCoefs;
    int nf      = qbcp_l.nf;
    double tsem = tem * qbcp_l.us_em.ns;   //new tem in SEM units

    //====================================================================
    // 2. Evaluate the Fourier coefficients
    //====================================================================
    //-------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //-------------------------------------------------------------------------------
    double delta[noc];
    evaluateCoef(delta, tsem, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);
    double deltad[noc];
    evaluateCoefDerivatives(deltad, tsem, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, tem, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);
    double alphad[noc];
    evaluateCoefDerivatives(alphad, tem, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);


    //====================================================================
    // 3. Evaluate the semi-analytical position of the Primary
    //====================================================================
    double ySEM[6], yEM[6];

    switch(prim)
    {
    case EARTH:
        //Earth position & momenta in SEM ref and SEM units
        ySEM[0] = delta[8];
        ySEM[1] = delta[9];
        ySEM[2] = 0.0;
        ySEM[3] = deltad[8];
        ySEM[4] = deltad[9];
        ySEM[5] = 0.0;
        //Earth position & momenta in EM ref and EM units
        yEM[0] = qbcp_l.us_em.mu_EM;
        yEM[1] = 0.0;
        yEM[2] = 0.0;
        yEM[3] = 0.0;
        yEM[4] = 0.0;
        yEM[5] = 0.0;
        break;

    case MOON:
        //Moon position & velocity in SEM ref and SEM units
        ySEM[0] = delta[10];
        ySEM[1] = delta[11];
        ySEM[2] = 0.0;
        ySEM[3] = deltad[10];
        ySEM[4] = deltad[11];
        ySEM[5] = 0.0;
        //Moon position & velocity in EM ref and EM units
        yEM[0] = qbcp_l.us_em.mu_EM - 1;
        yEM[1] = 0.0;
        yEM[2] = 0.0;
        yEM[3] = 0.0;
        yEM[4] = 0.0;
        yEM[5] = 0.0;
        break;

    case SUN:
        //Sun position & velocity in SEM ref and SEM units
        ySEM[0] = qbcp_l.us_sem.mu_SEM;
        ySEM[1] = 0.0;
        ySEM[2] = 0.0;
        ySEM[3] = 0.0;
        ySEM[4] = 0.0;
        ySEM[5] = 0.0;
        //Sun position & velocity in EM ref and EM units
        yEM[0] = alpha[6];
        yEM[1] = alpha[7];
        yEM[2] = 0.0;
        yEM[3] = alphad[6];
        yEM[4] = alphad[7];
        yEM[5] = 0.0;

        //            cout << "t*qbcp_l.us_em.n/(2*pi) = " << tem*qbcp_l.us_em.n/(2*M_PI) << endl;
        //            cout << "angle/(2*pi)  = " << atan2 (yEM[1],yEM[0])/(2*M_PI) << endl;
        //            cout << "Sun position = " << endl;
        //            vector_printf(yEM, 3);

        break;

    default:
        perror("Unknown primary");
    }



    //====================================================================
    // 4. Back in inertial coordinates,
    //    still with the use of the semi-analytical expressions
    //====================================================================
    double ytp[6];
    qbcp_coc(tsem, ySEM, ytp, VSEM, VEM);
    EMvtoIN(tem, ytp, ySEM_to_IN, &qbcp_l);
    EMvtoIN(tem, yEM, yEM_to_IN,  &qbcp_l);
}

/**
 * \brief Computes the "true" integrated position of the Earth, the Moon, and the Sun from the integrated bicircular movement contained in Z(t) and z(t).
 **/
void primary_integrated_position(int prim, double yIN[], cdouble z0, cdouble Z0, cdouble zdot0, cdouble Zdot0, QBCP_L& qbcp_l)
{
    //Initialization
    USYS us_em  = qbcp_l.us_em;

    //====================================================================
    // 1. Evaluate the integrated position of the Primary
    //====================================================================
    switch(prim)
    {
    case EARTH:
        yIN[0] = -us_em.ms / (1.0 + us_em.ms) * creal(Z0) + us_em.mu_EM * creal(z0);
        yIN[1] = -us_em.ms / (1.0 + us_em.ms) * cimag(Z0) + us_em.mu_EM * cimag(z0);
        yIN[2] = 0.0;
        yIN[3] = -us_em.ms / (1.0 + us_em.ms) * creal(Zdot0) + us_em.mu_EM * creal(zdot0);
        yIN[4] = -us_em.ms / (1.0 + us_em.ms) * cimag(Zdot0) + us_em.mu_EM * cimag(zdot0);
        yIN[5] = 0.0;
        break;

    case MOON:
        yIN[0] = -us_em.ms / (1.0 + us_em.ms) * creal(Z0) - (1 - us_em.mu_EM) * creal(z0);
        yIN[1] = -us_em.ms / (1.0 + us_em.ms) * cimag(Z0) - (1 - us_em.mu_EM) * cimag(z0);
        yIN[2] = 0.0;
        yIN[3] = -us_em.ms / (1.0 + us_em.ms) * creal(Zdot0) - (1 - us_em.mu_EM) * creal(zdot0);
        yIN[4] = -us_em.ms / (1.0 + us_em.ms) * cimag(Zdot0) - (1 - us_em.mu_EM) * cimag(zdot0);
        yIN[5] = 0.0;
        break;

    case SUN:
        yIN[0] = 1.0 / (1.0 + us_em.ms) * creal(Z0);
        yIN[1] = 1.0 / (1.0 + us_em.ms) * cimag(Z0);
        yIN[2] = 0.0;
        yIN[3] = 1.0 / (1.0 + us_em.ms) * creal(Zdot0);
        yIN[4] = 1.0 / (1.0 + us_em.ms) * cimag(Zdot0);
        yIN[5] = 0.0;
        break;

    default:
        perror("Unknown primary");
    }

}

/**
 *  \brief Test function to compare the analytical solution of the QBTBP to the numerical integration of the equations of motion.
 */
void qbtbp_test(double t1, QBCP_L& qbcp_l)
{
    cout << "---------------------------------------------------" << endl;
    cout << "                 Test of the                       " << endl;
    cout << "   Quasi-Bicircular Three-Body Problem (QBTBP)     " << endl;
    cout << "---------------------------------------------------" << endl;
    //Param
    double n   = qbcp_l.us_em.n;
    double ns  = qbcp_l.us_em.ns;
    double as  = qbcp_l.us_em.as;
    double ni  = qbcp_l.us_em.ni;
    double ai  = qbcp_l.us_em.ai;
    int noc    = qbcp_l.numberOfCoefs;
    int nf     = qbcp_l.nf;

    //====================================================================
    // 0. Gnuplot
    //====================================================================
    //----------------------------------------------------------
    //Gnuplot window
    //----------------------------------------------------------
    gnuplot_ctrl* h1;
    h1 = gnuplot_init();
    gnuplot_cmd(h1, "set title \"Variations of the error in the Primaries' state vector: EM vs IN\" ");
    gnuplot_cmd(h1, "set grid");
    gnuplot_set_xlabel(h1, (char*)"t [EM units]");
    gnuplot_set_ylabel(h1, (char*)"Error [x Tsem]");
    gnuplot_cmd(h1, "set logscale y");
    gnuplot_cmd(h1, "set format y \"1e\%%L\"");


    gnuplot_ctrl* h2;
    h2 = gnuplot_init();
    gnuplot_cmd(h2, "set title \"Variations of the error in the Primaries' state vector: SEM vs IN\" ");
    gnuplot_cmd(h2, "set grid");
    gnuplot_set_xlabel(h2, (char*)"t [EM units]");
    gnuplot_set_ylabel(h2, (char*)"Error [x Tsem]");
    gnuplot_cmd(h2, "set logscale y");
    gnuplot_cmd(h2, "set format y \"1e\%%L\"");


    //====================================================================
    // 1. Init the integration tools
    //====================================================================
    //Stepper
    const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
    OdeStruct ode_s;
    //Root-finding
    const gsl_root_fsolver_type* T_root = gsl_root_fsolver_brent; //Brent-Dekker root finding method
    //Ode solver parameters
    double param[2];
    param[0] = qbcp_l.us_em.mu_EM; //note that the qbtbp is computed in EM framework
    param[1] = qbcp_l.us_em.ms;    //note that the qbtbp is computed in EM framework
    //General structures
    init_ode_structure(&ode_s, T, T_root, PREC_ABS, PREC_REL, PREC_ROOT, PREC_DIFF,  8, PREC_HSTART,  qbtbp_derivatives, NULL, param);

    //====================================================================
    // 2. Initial conditions
    //====================================================================
    double t0  = 0.0;             //new t0 in EM units
    double t0c = t0 * qbcp_l.us_em.ns;   //new t0 in SEM units

    //z(0) and Z(0)
    cdouble z0    = evz(qbcp_l.cs_em.zt, t0, n, ni, ai);
    cdouble Z0    = evz(qbcp_l.cs_em.Zt, t0, n, ns, as);
    cdouble zdot0 = evzdot(qbcp_l.cs_em.zt, qbcp_l.cs_em.ztdot, t0, n, ni, ai);
    cdouble Zdot0 = evzdot(qbcp_l.cs_em.Zt, qbcp_l.cs_em.Ztdot, t0, n, ns, as);

    //Put the IC in real form
    double yv[8];
    yv[0] = creal(z0);
    yv[1] = cimag(z0);
    yv[2] = creal(Z0);
    yv[3] = cimag(Z0);
    yv[4] = creal(zdot0);
    yv[5] = cimag(zdot0);
    yv[6] = creal(Zdot0);
    yv[7] = cimag(Zdot0);

    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(15);
    cout << "Initial positions: " << endl;
    cout << "Internal motion: z(t=0.0) = " << creal(z0) << " + " << cimag(z0) << "i" <<  endl;
    cout << "External motion: Z(t=0.0) = " << creal(Z0) << " + " << cimag(Z0) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 2.2 Same for the Primaries
    //====================================================================
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);
    double yEM_to_IN[6], ySEM_to_IN[6], yIN[6];

    //-------------------------------------------------------------------------------
    //Initial Earth position
    //-------------------------------------------------------------------------------
    primary_analytical_position(EARTH, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(EARTH, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Earth position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    cout << "-------------------------------------------" << endl;

    //-------------------------------------------------------------------------------
    //Initial Moon position
    //-------------------------------------------------------------------------------
    primary_analytical_position(MOON, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(MOON, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Moon position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    cout << "-------------------------------------------" << endl;


    //-------------------------------------------------------------------------------
    //Initial Sun position
    //-------------------------------------------------------------------------------
    primary_analytical_position(SUN, t0, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(SUN, yIN, z0, Z0, zdot0, Zdot0, qbcp_l);
    cout << "Initial Sun position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 3. Initial conditions of the primaries
    //====================================================================
    //-------------------------------------------------------------------------------
    //Evaluate the deltas @t = t0c
    //-------------------------------------------------------------------------------
    double delta[noc];
    evaluateCoef(delta, t0c, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);
    double deltad[noc];
    evaluateCoefDerivatives(deltad, t0c, qbcp_l.us_sem.n, nf, qbcp_l.cs_sem.coeffs, noc);

    //-------------------------------------------------------------------------------
    //Evaluate the alphas @t = t0
    //-------------------------------------------------------------------------------
    double alpha[noc];
    evaluateCoef(alpha, t0, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);
    double alphad[noc];
    evaluateCoefDerivatives(alphad, t0, qbcp_l.us_em.n, nf, qbcp_l.cs_em.coeffs, noc);


    //====================================================================
    //Loop
    //====================================================================
    int Nplot = 200;
    double t = t0, ti;
    double tvec[Nplot + 1], yEM_vs_IN_E[Nplot + 1], yEM_vs_IN_M[Nplot + 1], yEM_vs_IN_S[Nplot + 1];
    double ySEM_vs_IN_E[Nplot + 1], ySEM_vs_IN_M[Nplot + 1], ySEM_vs_IN_S[Nplot + 1];

    for(int p = 0; p <= Nplot; p++)
    {
        //-------------------------------------------------------------------------------
        //Integration from t to ti
        //-------------------------------------------------------------------------------
        ti = t0 + (t1 - t0) * p / Nplot;
        gsl_odeiv2_driver_apply (ode_s.d, &t , ti, yv);
        tvec[p] = ti / qbcp_l.us_em.T;

        //-------------------------------------------------------------------------------
        //Current Earth position
        //-------------------------------------------------------------------------------
        primary_analytical_position(EARTH, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(EARTH, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);

        yEM_vs_IN_E[p]  = DENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_E[p] = DENorm(ySEM_to_IN, yIN, 6);

        //-------------------------------------------------------------------------------
        //Current Moon position
        //-------------------------------------------------------------------------------
        primary_analytical_position(MOON, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(MOON, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);

        yEM_vs_IN_M[p]  = DENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_M[p] = DENorm(ySEM_to_IN, yIN, 6);

        //-------------------------------------------------------------------------------
        //Current Sun position
        //-------------------------------------------------------------------------------
        primary_analytical_position(SUN, ti, yEM_to_IN, ySEM_to_IN, qbcp_l);
        primary_integrated_position(SUN, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);

        yEM_vs_IN_S[p]  = DENorm(yEM_to_IN, yIN, 6);
        ySEM_vs_IN_S[p] = DENorm(ySEM_to_IN, yIN, 6);
    }

    //-------------------------------------------------------------------------------
    //Plot: EM vs IN
    //-------------------------------------------------------------------------------
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_E, Nplot + 1, (char*)"Earth", "lines", "3", "2", 3);
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_M, Nplot + 1, (char*)"Moon", "lines", "3", "2", 4);
    gnuplot_plot_xy(h1, tvec, yEM_vs_IN_S, Nplot + 1, (char*)"Sun", "lines", "3", "2", 5);

    gnuplot_fplot_xy(tvec, yEM_vs_IN_E, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_EM_vs_IN_E.txt").c_str());
    gnuplot_fplot_xy(tvec, yEM_vs_IN_M, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_EM_vs_IN_M.txt").c_str());
    gnuplot_fplot_xy(tvec, yEM_vs_IN_S, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_EM_vs_IN_S.txt").c_str());


    //-------------------------------------------------------------------------------
    //Plot: SEM vs IN
    //-------------------------------------------------------------------------------
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_E, Nplot + 1, (char*)"Earth", "lines", "3", "2", 3);
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_M, Nplot + 1, (char*)"Moon", "lines", "3", "2", 4);
    gnuplot_plot_xy(h2, tvec, ySEM_vs_IN_S, Nplot + 1, (char*)"Sun", "lines", "3", "2", 5);

    gnuplot_fplot_xy(tvec, ySEM_vs_IN_E, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_SEM_vs_IN_E.txt").c_str());
    gnuplot_fplot_xy(tvec, ySEM_vs_IN_M, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_SEM_vs_IN_M.txt").c_str());
    gnuplot_fplot_xy(tvec, ySEM_vs_IN_S, Nplot + 1, (char*) (qbcp_l.cs.F_PLOT + "QBTBP_SEM_vs_IN_S.txt").c_str());

    //====================================================================
    //Final state
    //====================================================================
    cout << "-------------------------------------------" << endl;
    cout << "End of integration." << endl;
    cout << "Final t: " << t << endl;
    cout << "-------------------------------------------" << endl;
    cout << "Results from integration" << endl;
    cout << "Internal motion: z(t=t1) = " << yv[0] << " + " << yv[1] << "i" << endl;
    cout << "External motion: Z(t=t1) = " << yv[2] << " + " << yv[3] << "i" << endl;
    cout << "-------------------------------------------" << endl;

    //Analytical final state
    cdouble zfinal = evz(qbcp_l.cs_em.zt, t1, n, ni, ai);
    cdouble Zfinal = evz(qbcp_l.cs_em.Zt, t1, n, ns, as);

    cdouble zdotfinal = evzdot(qbcp_l.cs_em.zt, qbcp_l.cs_em.ztdot, t1, n, ni, ai);
    cdouble Zdotfinal = evzdot(qbcp_l.cs_em.Zt, qbcp_l.cs_em.Ztdot, t1, n, ns, as);


    cout << "Analytical results" << endl;
    cout << "Internal motion: z(t=t1) = " << creal(zfinal) << " + " << cimag(zfinal) << "i" <<  endl;
    cout << "External motion: Z(t=t1) = " << creal(Zfinal) << " + " << cimag(Zfinal) << "i" <<  endl;
    cout << "-------------------------------------------" << endl;


    cout << "Absolute delta between analytical and numerical results: " << endl;
    cout << "dz    = " << cabs(zfinal - yv[0] - I * yv[1])    << endl;
    cout << "dzdot = " << cabs(zdotfinal - yv[4] - I * yv[5]) << endl;
    cout << "dZ    = " << cabs(Zfinal - yv[2] - I * yv[3])    << endl;
    cout << "dZdot = " << cabs(Zdotfinal - yv[6] - I * yv[7]) << endl;
    cout << "-------------------------------------------" << endl;

    //====================================================================
    // 2.2 Same for the Primaries
    //====================================================================
    cout << std::showpos << setiosflags(ios::scientific)  << setprecision(5);

    //-------------------------------------------------------------------------------
    //Current Earth position
    //-------------------------------------------------------------------------------
    primary_analytical_position(EARTH, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(EARTH, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);
    cout << "Final Earth position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    cout << "-------------------------------------------" << endl;

    //-------------------------------------------------------------------------------
    //Current Moon position
    //-------------------------------------------------------------------------------
    primary_analytical_position(MOON, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(MOON, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);
    cout << "Final Moon position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    cout << "-------------------------------------------" << endl;


    //-------------------------------------------------------------------------------
    //Current Sun position
    //-------------------------------------------------------------------------------
    primary_analytical_position(SUN, t1, yEM_to_IN, ySEM_to_IN, qbcp_l);
    primary_integrated_position(SUN, yIN, yv[0] + I * yv[1], yv[2] + I * yv[3], yv[4] + I * yv[5], yv[6] + I * yv[7], qbcp_l);
    cout << "Final Sun position:" << endl;
    cout << "-------------------------------------------" << endl;
    cout << "EM to IN (1)     SEM to IN (2)     IN (3)         (1) - (3)       (2) - (3)" << endl;

    for(int i = 0; i < 6; i++) cout << yEM_to_IN[i] << "    "  << ySEM_to_IN[i] << "    "  << yIN[i] << "    "
                                        << yEM_to_IN[i] - yIN[i] << "    "  << ySEM_to_IN[i] - yIN[i] << endl;

    char ch;
    gnuplot_cmd(h1, "set logscale y");
    printf("Press ENTER to close the gnuplot window(s)\n");
    scanf("%c", &ch);
    gnuplot_close(h1);
    gnuplot_close(h2);
}


