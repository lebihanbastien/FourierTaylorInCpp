#include "Invman.h"

//========================================================================================
// LOG
//
// - Evaluation is validated for RTBP.
//
// @todo Update the CU/CS for QBCP EML1/EML2
// @todo Write things clearly @ the beginning of this file.
// @todo Add tests at the beginning of each routine that is using ofts_order and ofs_order as inputs (compare with inner values!)
//========================================================================================



//========================================================================================
// Constructor
//========================================================================================
/**
 *  \brief Constructor.
 *         If the graph style is used, the series are very sparse
 *         and some simplifications can be made.
 **/
Invman::Invman(int ofts_order_, int ofs_order_, CSYS& csys):
//----------------------------------------------------------------------------------------
// Initialization list
//----------------------------------------------------------------------------------------
    ofs(ofs_order_),
    cs(&csys),
    n(csys.us.n),
    ofts_order(ofts_order_),
    ofs_order(ofs_order_),
    pmType(csys.pmType),
    manType(csys.manType),
    reduced_nv(compRNV(csys)),
    fwrk(csys.fwrk),
    ncs(csys.fwrk==F_SEM?NCSEM:NCEM),
    scs(csys.fwrk==F_SEM?PSEM:PEM),
    Mcoc(6, 6, ofs_order_),
    MIcoc(6, 6, ofs_order_),
    Vcoc(6, Ofsc(ofs_order_)),
    DWh(),
    Wh(),
    W(),
    Hy(),
    mIn(6, reduced_nv, ofs_order),
    nc_B_SYS(6, Ofsc(ofs_order_)),
    nc_R_SYS(6, 6, ofs_order_),
    SYS_R_nc(6, 6, ofs_order_),
    VSYS_R_SYS(6, 6, ofs_order_),
    SYS_R_VSYS(6, 6, ofs_order_)
//----------------------------------------------------------------------------------------
// Body of the constructor
//----------------------------------------------------------------------------------------
{
    //====================================================================================
    // Parameterization in TFC form
    //====================================================================================
    Wh.reserve(6);
    //------------------------------------------------------------------------------------
    //Initialize with respect to the manifold type and parameterization style
    //------------------------------------------------------------------------------------
    switch(csys.pmType)
    {
        //--------------------------------------------------------------------------------
        //If we use the graph style, some simplification can be made
        //--------------------------------------------------------------------------------
    case PMS_GRAPH:
        switch(csys.manType)
        {
        case MAN_CENTER:
            //----------------------------------------------------------------------------
            //If it is a center manifold parameterized using the graph style,
            //the following assertions are true:
            //      - The dimension 0, 2, 3, and 5 are linear in the parameters (order 1)
            //      - The dimension 1 and 4 are full fourier-taylor series
            //----------------------------------------------------------------------------
            Wh.push_back(Oftsc(reduced_nv, 1         , OFS_NV, ofs_order)); //0
            Wh.push_back(Oftsc(reduced_nv, ofts_order, OFS_NV, ofs_order)); //1
            Wh.push_back(Oftsc(reduced_nv, 1         , OFS_NV, ofs_order)); //2
            Wh.push_back(Oftsc(reduced_nv, 1         , OFS_NV, ofs_order)); //3
            Wh.push_back(Oftsc(reduced_nv, ofts_order, OFS_NV, ofs_order)); //4
            Wh.push_back(Oftsc(reduced_nv, 1         , OFS_NV, ofs_order)); //5

            //----------------------------------------------------------------------------
            //Read from file
            //----------------------------------------------------------------------------
            readVOFTS_bin(Wh, csys.F_PMS+"W/Wh");
            break;

        case MAN_CENTER_S:
        case MAN_CENTER_U:
            //----------------------------------------------------------------------------
            //If it is a center-(un)stable manifold parameterized using the graph style,
            //the following assertions are true:
            //      - Wh[0] = Whc[0] + F1(s5)
            //      - Wh[2] = Whc[2]
            //      - Wh[3] = Whc[3] + F3(s5)
            //      - Wh[5] = Whc[5]
            // where Whc is the center-manifold parameterization, and F1 and F3 are two
            // one-dimensionnal full series in the variable s5.
            //
            // Moreover, regarding the two remaining dimensions:
            //      - Wh[4] (resp. Wh[1]) is a full 5-dim FT series for CU (resp. CS).
            //      - Wh[1] =  Whc[1] + F2(s5) for CU,
            //      - Wh[4] =  Whc[4] + F2(s5) for CS,
            // where F2 is a one-dimensionnal first order series in the variable s5.
            //
            //----------------------------------------------------------------------------
            //Hence, Wh is initialized as in the center case,
            //for the dimension 0, 2, 3, and 5 + one additionnal
            //that depends on the man type

            //----------------------------------------------------------------------------
            //First line
            //----------------------------------------------------------------------------
            Wh.push_back(Oftsc(4, 1         , OFS_NV, ofs_order));
            //Read from file. However, csys does not contain the right folder, so we need
            //to go and seek the "graph style center" folder (not PMS folder!)
            readOFTS_bin(Wh[0], (csys.F_GS+"W/Wh[0].bin"));

            //----------------------------------------------------------------------------
            //Second line
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_U)
            {
                Wh.push_back(Oftsc(4, ofts_order, OFS_NV, ofs_order));
                //Read from file in GS
                readOFTS_bin(Wh[1], (csys.F_GS+"W/Wh[1].bin"));
            }
            else
            {
                Wh.push_back(Oftsc(5, ofts_order, OFS_NV, ofs_order));
                //Read from file in PMS
                readOFTS_bin(Wh[1], (csys.F_PMS+"W/Wh[1].bin"));
            }

            //----------------------------------------------------------------------------
            //Third line
            //----------------------------------------------------------------------------
            Wh.push_back(Oftsc(4, 1         , OFS_NV, ofs_order));
            //Read from file. However, csys does not contain the right folder, so we need
            //to go and seek the "graph style center" folder (not PMS folder!)
            readOFTS_bin(Wh[2], (csys.F_GS+"W/Wh[2].bin"));

            //----------------------------------------------------------------------------
            //Fourth line
            //----------------------------------------------------------------------------
            Wh.push_back(Oftsc(4, 1         , OFS_NV, ofs_order));
            //Read from file. However, csys does not contain the right folder, so we need
            //to go and seek the "graph style center" folder (not PMS folder!)
            readOFTS_bin(Wh[3], (csys.F_GS+"W/Wh[3].bin"));


            //----------------------------------------------------------------------------
            //Fifth line
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_S)
            {
                Wh.push_back(Oftsc(4, ofts_order, OFS_NV, ofs_order));
                //Read from file in GS
                readOFTS_bin(Wh[4], (csys.F_GS+"W/Wh[4].bin"));
            }
            else
            {
                Wh.push_back(Oftsc(5, ofts_order, OFS_NV, ofs_order));
                //Read from file in PMS
                readOFTS_bin(Wh[4], (csys.F_PMS+"W/Wh[4].bin"));
            }

            //----------------------------------------------------------------------------
            //Sixth line
            //----------------------------------------------------------------------------
            Wh.push_back(Oftsc(4, 1         , OFS_NV, ofs_order));
            //Read from file. However, csys does not contain the right folder, so we need
            //to go and seek the "graph style center" folder (not PMS folder!)
            readOFTS_bin(Wh[5], (csys.F_GS+"W/Wh[5].bin"));


            //----------------------------------------------------------------------------
            //Hy is initialized as the following:
            //Hy[0] will contain either Wh[4] or Wh[1], i.e. the full 5-dim FT series.
            //Hy[1 to 3] will contain F1 to F3
            //----------------------------------------------------------------------------
            Hy.reserve(4);
            Hy.push_back(Oftsc(5, ofts_order, OFS_NV, ofs_order)); //0
            Hy.push_back(Oftsc(1, ofts_order, OFS_NV, ofs_order)); //1
            Hy.push_back(Oftsc(1,          1, OFS_NV, ofs_order)); //2
            Hy.push_back(Oftsc(1, ofts_order, OFS_NV, ofs_order)); //3

            //----------------------------------------------------------------------------
            //To get Hy[2], we read Wh[1] or Wh[4] at order one with Hy[0],
            //then update Hy[2]
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_U)
                readOFTS_bin(Hy[0], csys.F_PMS+"W/Wh[1].bin", 1);
            else
                readOFTS_bin(Hy[0], csys.F_PMS+"W/Wh[4].bin", 1);
            //Updating Hy[1]
            Hy[2].getCA(1,0)->ccopy(*Hy[0].getCA(1, 4));



            //----------------------------------------------------------------------------
            //Hy[0] is found at the corresponding position
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_U)
                readOFTS_bin(Hy[0], csys.F_PMS+"W/Wh[4].bin");
            else
                readOFTS_bin(Hy[0], csys.F_PMS+"W/Wh[1].bin");

            //----------------------------------------------------------------------------
            //Hy[1] and Hy[3] are found at the corresponding positions
            //----------------------------------------------------------------------------
            readOFTS_bin(Hy[1], csys.F_PMS+"W/F[0].bin");
            readOFTS_bin(Hy[3], csys.F_PMS+"W/F[1].bin");
            break;

        case MAN_CENTER_US:
            //----------------------------------------------------------------------------
            //Full fourier-taylor series, for now
            //----------------------------------------------------------------------------
            for(int i = 0; i < 6; i++)
                Wh.push_back(Oftsc(reduced_nv, ofts_order, OFS_NV, ofs_order));

            //----------------------------------------------------------------------------
            //Read from file
            //----------------------------------------------------------------------------
            readVOFTS_bin(Wh, csys.F_PMS+"W/Wh");
            break;
        }
        break;

    case PMS_MIXED:
    case PMS_NORMFORM:
        //----------------------------------------------------------------------------
        //Full fourier-taylor series, for now
        //----------------------------------------------------------------------------
        for(int i = 0; i < 6; i++)
            Wh.push_back(Oftsc(reduced_nv, ofts_order, OFS_NV, ofs_order));
        break;

        //--------------------------------------------------------------------------------
        //Read from file
        //--------------------------------------------------------------------------------
        readVOFTS_bin(Wh, csys.F_PMS+"W/Wh");
        break;
    }

    //====================================================================================
    // Parameterization in NC form. Uncomment if necessary
    //====================================================================================
    //Reserve
    W.reserve(6);
    //Initialize the FULL Fourier-Taylor series
    for(int i = 0; i < 6; i++) W.push_back(Oftsc(reduced_nv, ofts_order, OFS_NV, ofs_order));
    //Read from file
    readVOFTS_bin(W,  csys.F_PMS+"W/W");


    //====================================================================================
    // CCM_R_RCM
    //====================================================================================
    CCM_R_RCM = gsl_matrix_complex_calloc(reduced_nv, reduced_nv);
    switch(csys.manType)
    {
    case MAN_CENTER:
        rotmat_CC_R_RCM_CENTER(CCM_R_RCM);
        break;
    case MAN_CENTER_S:
    case MAN_CENTER_U:
        rotmat_CC_R_RCM_CENTER_HYP(CCM_R_RCM);
        break;
    case MAN_CENTER_US: //nothing is done here
        rotmat_CC_R_RCM_CENTER_6(CCM_R_RCM);
        break;
    }

    //====================================================================================
    // Jacobian
    //====================================================================================
    //Reserve the size
    DWh.reserve(6, reduced_nv);

    //------------------------------------------------------------------------------------
    //Initialize with respect to the manifold type and parameterization style
    //------------------------------------------------------------------------------------
    switch(csys.pmType)
    {
        //--------------------------------------------------------------------------------
        //If we use the graph style, some simplification can be made
        //--------------------------------------------------------------------------------
    case PMS_GRAPH:
        switch(csys.manType)
        {
        case MAN_CENTER:
            //----------------------------------------------------------------------------
            //If it is a center manifold parameterized using the graph style,
            //the following assertions are true:
            //      - The dimension 0, 2, 3, and 5 are linear in the parameters (order 1),
            //        therefore the jacobian is of order 0!
            //      - The dimension 1 and 4 are full fourier-taylor series
            // For the jacobian, the sequence of orders of the series is thus,
            // with n = ofts_order:
            //
            //                  | 0  0  0  0  0 |
            //                  | n  n  n  n  n |
            //                  | 0  0  0  0  0 |
            //  order(DWh) =    | 0  0  0  0  0 |
            //                  | n  n  n  n  n |
            //                  | 0  0  0  0  0 |
            //
            // all being 4-dimensional series.
            //----------------------------------------------------------------------------
            //First line is of order 0
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, 0, OFS_NV, ofs_order));

            //Second line is of full order
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, ofts_order-1, OFS_NV, ofs_order));

            //Third line is of order 0
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, 0, OFS_NV, ofs_order));

            //Fourth line is of order 0
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, 0, OFS_NV, ofs_order));

            //Fifth line is of full order
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, ofts_order-1, OFS_NV, ofs_order));

            //Sixth line is of order 0
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, 0, OFS_NV, ofs_order));

            //----------------------------------------------------------------------------
            //Read from file
            //----------------------------------------------------------------------------
            readMOFTS_bin(DWh, csys.F_PMS+"DWf/DWhc");
            break;

        case MAN_CENTER_S:
        case MAN_CENTER_U:
            //----------------------------------------------------------------------------
            //If it is a center-(un)stable manifold parameterized using the graph style,
            //the following assertions are true:
            //      - Wh[0] = Whc[0] + F1(s5)
            //      - Wh[2] = Whc[2]
            //      - Wh[3] = Whc[3] + F2(s5)
            //      - Wh[5] = Whc[5]
            // where Whc is the center-manifold parameterization, and F1 and F2 are two
            // one-dimensionnal series in the variable s5.
            //
            // Moreover, regarding the two remaining dimensions:
            //      - Wh[4] (resp. Wh[1]) is a full 5-dim FT series for CU (resp. CS).
            //      - Wh[1] = Whc[1] + c1*s5 for CU
            //      - Wh[4] = Whc[4] + c2*s5 for CS
            //
            // For the Jacobian, the sequence of orders of the series is thus,
            // with n = ofts_order:
            //
            //               | 0  0  0  0  0  n |          | 0  0  0  0  0  n |
            //               | n  n  n  n  n  0 |          | n  n  n  n  n  n |
            //               | 0  0  0  0  0  0 |          | 0  0  0  0  0  0 |
            //  order(DWh) = | 0  0  0  0  0  n | for CU,  | 0  0  0  0  0  n |, for CS.
            //               | n  n  n  n  n  n |          | n  n  n  n  n  0 |
            //               | 0  0  0  0  0  0 |          | 0  0  0  0  0  0 |
            //
            // The sequence of dimensions of the series is:
            //
            //               | 4  4  4  4  4  1 |          | 4  4  4  4  4  1 |
            //               | 4  4  4  4  4  4 |          | 5  5  5  5  5  5 |
            //               | 4  4  4  4  4  1 |          | 4  4  4  4  4  1 |
            //  order(DWh) = | 4  4  4  4  4  1 | for CU,  | 4  4  4  4  4  1 |, for CS.
            //               | 5  5  5  5  5  5 |          | 4  4  4  4  4  4 |
            //               | 4  4  4  4  4  1 |          | 4  4  4  4  4  1 |
            //
            //----------------------------------------------------------------------------

            //----------------------------------------------------------------------------
            //First line
            //----------------------------------------------------------------------------
            for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));
            DWh.mpush_back(Oftsc(1, ofts_order-1, OFS_NV, ofs_order));

            //Read from file
            readOFTS_bin(*DWh.getCA(0, 0), (csys.F_PMS+"DWf/DWhc[0][0].bin"));
            readOFTS_bin(*DWh.getCA(0, 4), (csys.F_PMS+"W/Fd[0].bin"));

            //----------------------------------------------------------------------------
            //Second line
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_U)
            {
                for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, ofts_order-1, OFS_NV, ofs_order));
                DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));

                //Read from file (F_GS)
                readOFTS_bin(*DWh.getCA(1, 0), (csys.F_GS+"DWf/DWhc[1][0].bin"));
                readOFTS_bin(*DWh.getCA(1, 1), (csys.F_GS+"DWf/DWhc[1][1].bin"));
                readOFTS_bin(*DWh.getCA(1, 2), (csys.F_GS+"DWf/DWhc[1][2].bin"));
                readOFTS_bin(*DWh.getCA(1, 3), (csys.F_GS+"DWf/DWhc[1][3].bin"));
                DWh.getCA(1, 4)->getCA(0,0)->ccopy(*Hy[2].getCA(1,0));
            }
            else
            {
                for(int j = 0; j < 5; j++)
                {
                    DWh.mpush_back(Oftsc(5, ofts_order-1, OFS_NV, ofs_order));
                    readOFTS_bin(*DWh.getCA(1, j), (csys.F_PMS+"DWf/DWhc[1]["+numTostring(j)+"].bin"));
                }
            }

            //----------------------------------------------------------------------------
            //Third line
            //----------------------------------------------------------------------------
            for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));
            DWh.mpush_back(Oftsc(1, 0, OFS_NV, ofs_order));

            //Read from file
            readOFTS_bin(*DWh.getCA(2, 1), (csys.F_PMS+"DWf/DWhc[2][1].bin"));

            //----------------------------------------------------------------------------
            //Fourth line
            //----------------------------------------------------------------------------
            for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));
            DWh.mpush_back(Oftsc(1, ofts_order-1, OFS_NV, ofs_order));
            //Read from file
            readOFTS_bin(*DWh.getCA(3, 2), (csys.F_PMS+"DWf/DWhc[3][2].bin"));
            readOFTS_bin(*DWh.getCA(3, 4), (csys.F_PMS+"W/Fd[1].bin"));

            //----------------------------------------------------------------------------
            //Fifth line
            //----------------------------------------------------------------------------
            if(csys.manType == MAN_CENTER_S)
            {
                for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, ofts_order-1, OFS_NV, ofs_order));
                DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));

                //Read from file (F_GS)
                readOFTS_bin(*DWh.getCA(4, 0), (csys.F_GS+"DWf/DWhc[4][0].bin"));
                readOFTS_bin(*DWh.getCA(4, 1), (csys.F_GS+"DWf/DWhc[4][1].bin"));
                readOFTS_bin(*DWh.getCA(4, 2), (csys.F_GS+"DWf/DWhc[4][2].bin"));
                readOFTS_bin(*DWh.getCA(4, 3), (csys.F_GS+"DWf/DWhc[4][3].bin"));
                DWh.getCA(4, 4)->getCA(0,0)->ccopy(*Hy[2].getCA(1,0));
            }
            else
            {
                for(int j = 0; j < 5; j++)
                {
                    DWh.mpush_back(Oftsc(5, ofts_order-1, OFS_NV, ofs_order));
                    readOFTS_bin(*DWh.getCA(4, j), (csys.F_PMS+"DWf/DWhc[4]["+numTostring(j)+"].bin"));
                }
            }

            //----------------------------------------------------------------------------
            //Sixth line
            //----------------------------------------------------------------------------
            for(int j = 0; j < 4; j++) DWh.mpush_back(Oftsc(4, 0, OFS_NV, ofs_order));
            DWh.mpush_back(Oftsc(1, 0, OFS_NV, ofs_order));
            //Read from file
            readOFTS_bin(*DWh.getCA(5, 3), (csys.F_PMS+"DWf/DWhc[5][3].bin"));


            break;

        case MAN_CENTER_US:
            //----------------------------------------------------------------------------
            //Full fourier-taylor series, for now
            //----------------------------------------------------------------------------
            for(int i = 0; i <6; i++)
                for(int j = 0; j <reduced_nv; j++)
                    DWh.mpush_back(Oftsc(reduced_nv, ofts_order-1, OFS_NV, ofs_order));


            //----------------------------------------------------------------------------
            //Read from file
            //----------------------------------------------------------------------------
            readMOFTS_bin(DWh, csys.F_PMS+"DWf/DWhc");
            break;
        }
        break;

    case PMS_MIXED:
    case PMS_NORMFORM:
        //--------------------------------------------------------------------------------
        //Full fourier-taylor series, for now
        //--------------------------------------------------------------------------------
        for(int i = 0; i <6; i++)
            for(int j = 0; j <reduced_nv; j++)
                DWh.mpush_back(Oftsc(reduced_nv, ofts_order-1, OFS_NV, ofs_order));

        //--------------------------------------------------------------------------------
        //Read from file
        //--------------------------------------------------------------------------------
        readMOFTS_bin(DWh, csys.F_PMS+"DWf/DWhc");
        break;
    }


    //====================================================================================
    // COC
    //====================================================================================
    initCOC(Mcoc, MIcoc, Vcoc, csys);

    //====================================================================================
    // For projection purposes
    //====================================================================================
    omega1 = cimag(Wh[0].getCA(1,0)->getCoef(0));
    omega3 = cimag(Wh[2].getCA(1,1)->getCoef(0));


    //====================================================================================
    // For COC: NCEM <-> NCSEM
    //====================================================================================
    Ofsc a1, a2, a3, ma3, ma2da1, a3da1, ma3da1, oda1, temp;

    //Read the vector field
    readOFS_txt(a1, cs->F_COEF+"alpha1_fft");
    readOFS_txt(a2, cs->F_COEF+"alpha2_fft");
    readOFS_txt(a3, cs->F_COEF+"alpha3_fft");

    //------------------------------------------------------------------------------------
    // nc_B_SYS = c1 * (1  0  0 -alpha_2/alpha_1  +alpha_3/alpha_1 0)^T
    //------------------------------------------------------------------------------------
    // nc_B_SYS[0] = c1
    nc_B_SYS[0].setCoef(cs->c1, 0);

    // nc_B_SYS[3] = -c1*alpha_2/alpha_1
    ma2da1.ofs_div(a2, a1, temp);
    ma2da1 *= -1.0;
    nc_B_SYS[3].ofs_smult(ma2da1, +cs->c1, ofs_order);

    // nc_B_SYS[4] = +c1*alpha_3/alpha_1
    a3da1.ofs_div(a3, a1, temp);

    nc_B_SYS[4].ofs_smult(a3da1, +cs->c1, ofs_order);

    //------------------------------------------------------------------------------------
    // nc_R_SYS = 1/gamma * diag(-1 -1 1 -1 -1 1)
    // SYS_R_nc = + gamma * diag(-1 -1 1 -1 -1 1)
    //------------------------------------------------------------------------------------
    nc_R_SYS.setCoef(-1.0/cs->gamma, 0, 0);
    nc_R_SYS.setCoef(-1.0/cs->gamma, 1, 1);
    nc_R_SYS.setCoef(+1.0/cs->gamma, 2, 2);

    nc_R_SYS.setCoef(-1.0/cs->gamma, 3, 3);
    nc_R_SYS.setCoef(-1.0/cs->gamma, 4, 4);
    nc_R_SYS.setCoef(+1.0/cs->gamma, 5, 5);


    SYS_R_nc.setCoef(-1.0*cs->gamma, 0, 0);
    SYS_R_nc.setCoef(-1.0*cs->gamma, 1, 1);
    SYS_R_nc.setCoef(+1.0*cs->gamma, 2, 2);

    SYS_R_nc.setCoef(-1.0*cs->gamma, 3, 3);
    SYS_R_nc.setCoef(-1.0*cs->gamma, 4, 4);
    SYS_R_nc.setCoef(+1.0*cs->gamma, 5, 5);


    //------------------------------------------------------------------------------------
    //                 |  1       0       0       0    0     0    |
    //                 |  0       1       0       0    0     0    |
    //                 |  0       0       1       0    0     0    |
    //  SYS_R_VSYS =   | -a2/a1  -a3/a1   0       1/a1   0   0    |
    //                 | +a3/a1  -a2/a1   0       0    1/a1  0    |
    //                 |  0       0      -a2/a1   0    0     1/a1 |
    //------------------------------------------------------------------------------------
    //oda1 = 1/a1;
    oda1.ofs_pows(a1, -1.0);
    //ma3da1 = -a3/a1
    ma3da1.ofs_smult(a3da1, -1.0, ofs_order);

    SYS_R_VSYS.setCoef(1.0, 0, 0);
    SYS_R_VSYS.setCoef(1.0, 1, 1);
    SYS_R_VSYS.setCoef(1.0, 2, 2);

    SYS_R_VSYS.setCoef(ma2da1, 3, 0);
    SYS_R_VSYS.setCoef(ma2da1, 4, 1);
    SYS_R_VSYS.setCoef(ma2da1, 5, 2);

    SYS_R_VSYS.setCoef(ma3da1, 3, 1);
    SYS_R_VSYS.setCoef(a3da1,  4, 0);

    SYS_R_VSYS.setCoef(oda1, 3, 3);
    SYS_R_VSYS.setCoef(oda1, 4, 4);
    SYS_R_VSYS.setCoef(oda1, 5, 5);

    //------------------------------------------------------------------------------------
    //                 |  1    0    0    0    0   0  |
    //                 |  0    1    0    0    0   0  |
    //                 |  0    0    1    0    0   0  |
    //  VSYS_R_SYS =   | +a2  +a3   0    a1   0   0  |
    //                 | -a3  +a2   0    0    a1  0  |
    //                 |  0    0   +a2   0    0   a1 |
    //------------------------------------------------------------------------------------
    //ma3 = -a3
    ma3.ofs_smult(a3, -1.0, ofs_order);

    VSYS_R_SYS.setCoef(1.0, 0, 0);
    VSYS_R_SYS.setCoef(1.0, 1, 1);
    VSYS_R_SYS.setCoef(1.0, 2, 2);

    VSYS_R_SYS.setCoef(a2, 3, 0);
    VSYS_R_SYS.setCoef(a2, 4, 1);
    VSYS_R_SYS.setCoef(a2, 5, 2);

    VSYS_R_SYS.setCoef(a3,  3, 1);
    VSYS_R_SYS.setCoef(ma3, 4, 0);

    VSYS_R_SYS.setCoef(a1, 3, 3);
    VSYS_R_SYS.setCoef(a1, 4, 4);
    VSYS_R_SYS.setCoef(a1, 5, 5);
}

//========================================================================================
// Destructor
//========================================================================================
/**
 *  \brief Destructor. Nothing to do here.
 **/
Invman::~Invman()
{
    //dtor
}

//========================================================================================
// Evaluate
//========================================================================================
/**
 *  \brief Evaluate the invariant manifold from RCM coordinates (double[]), to
 *         TFC coordinates (vector<Ofsc>, i.e. a vector of Fourier series).
 **/
void Invman::evalRCMtoTFC(double const st0[], vector<Ofsc>& zIn, const int ofts_order, const int ofs_order) const
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble s0[reduced_nv];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);

    //------------------------------------------
    // CCM to TFC
    //------------------------------------------
    this->evalCCMtoTFC(s0, zIn, ofts_order, ofs_order);
}


/**
 *  \brief Evaluate the invariant manifold from CCM coordinates (complex double[]), to
 *         TFC coordinates (vector<Ofsc>, i.e. a vector of Fourier series).
 **/
void Invman::evalCCMtoTFC(cdouble const s0[], vector<Ofsc>& zIn, const int ofts_order, const int ofs_order) const
{
    //------------------------------------------------------------------------------------
    //Basic checks
    //------------------------------------------------------------------------------------
    if(ofts_order > this->ofts_order)
    {
        cout << "Invman::evalCCMtoTFC. ofts_order is too big. return." << endl;
        return;
    }

    if(ofs_order > this->ofs_order)
    {
        cout << "Invman::evalCCMtoTFC. ofs_order is too big. return." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Basic evaluation
    //------------------------------------------------------------------------------------
    for(int p = 0; p < 6; p++)
    {
        Wh[p].evaluate(s0, zIn[p], ofts_order, ofs_order);
    }

    //------------------------------------------------------------------------------------
    //Specific case of PMS_GRAPH (graph style)
    //------------------------------------------------------------------------------------
    if(pmType == PMS_GRAPH)
    {
        switch(manType)
        {
        case MAN_CENTER:
            //nothing to do here, the manifold is well updated via the basic evaluation
            break;

        case MAN_CENTER_U:
            //----------------------------------------------------------------------------
            //Add the necessary values:
            // Wh[0] = Whc[0] + Hy[1](s5)
            // Wh[3] = Whc[3] + Hy[3](s5)
            // Wh[1] = Whc[1] + Hy[2](s5) for CU
            //----------------------------------------------------------------------------
            Hy[1].sevaluate(&s0[4], zIn[0], ofts_order, ofs_order);
            Hy[3].sevaluate(&s0[4], zIn[3], ofts_order, ofs_order);
            Hy[2].sevaluate(&s0[4], zIn[1], ofts_order, ofs_order);
            break;

        case MAN_CENTER_S:
            //----------------------------------------------------------------------------
            //Add the necessary values:
            // Wh[0] = Whc[0] + Hy[1](s5)
            // Wh[3] = Whc[3] + Hy[3](s5)
            // Wh[4] = Whc[4] + Hy[2](s5) for CS
            //----------------------------------------------------------------------------
            Hy[1].sevaluate(&s0[4], zIn[0], ofts_order, ofs_order);
            Hy[3].sevaluate(&s0[4], zIn[3], ofts_order, ofs_order);
            Hy[2].sevaluate(&s0[4], zIn[4], ofts_order, ofs_order);
            break;

        default:
            cout << "Invman::evalCCMtoTFC. Not defined for this type of manifold in graph mode." << endl;
            break;

        }
    }
}

/**
 *  \brief Evaluate the invariant manifold from CCM coordinates (complex double[]), to
 *         NC coordinates (double[]).
 **/
void Invman::evalCCMtoNC(cdouble const s0[], double const t, double z1[], const int ofts_order, const int ofs_order) const
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------
    // CCM to TFC
    //------------------------------------------
    this->evalCCMtoTFC(s0, zIn, ofts_order, ofs_order);

    //------------------------------------------
    // TFC to NC
    //------------------------------------------
    applyCOC(Mcoc, Vcoc, zIn, zOut);
    for(int p = 0; p < 6; p++) z1[p] = creal(zOut[p].evaluate(n*t, ofs_order));
}

/**
 *  \brief Evaluate the invariant manifold from RCM coordinates (double[]), to
 *         NC coordinates (double[]).
 **/
void Invman::evalRCMtoNC(double const st0[], double const t, double z1[], const int ofts_order, const int ofs_order) const
{
    //------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------
    cdouble s0[reduced_nv];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);

    //------------------------------------------
    // CCM to NC
    //------------------------------------------
    this->evalCCMtoNC(s0, t, z1, ofts_order, ofs_order);
}

//========================================================================================
// Evaluate the Jacobian
//========================================================================================
/**
 *  \brief Evaluate the jacobian of the invariant manifold from CCM coordinates
 *         (complex double[]), to TFC coordinates (matrix<Ofsc>, i.e. a matrix of
 *          Fourier series).
 **/
void Invman::evalDCCMtoTFC(cdouble const s0[], matrix<Ofsc>& mIn, const int ofts_order, const int ofs_order) const
{
    //------------------------------------------
    // 1. Check the size
    //------------------------------------------
    if(mIn.getSize(1) != DWh.getSize(1) || mIn.getSize(2) != DWh.getSize(2))
    {
        cout << "CCMtoTFC_JAC. Dimensions do not match. return." << endl;
        return;
    }

    //--------------------
    // 2. Zeroing
    //--------------------
    for(int i = 0; i < mIn.getSize(1); i++)
    {
        for(int j = 0; j < mIn.getSize(2); j++)
        {
            mIn.getCA(i,j)->zero();
        }
    }

    //Temp variables
    Ofsc ofs_temp(mIn.getCA(0,0)->getOrder());


    //--------------------
    // General computation
    //--------------------
    for(int i = 0; i < mIn.getSize(1); i++)
    {
        for(int j = 0; j < mIn.getSize(2); j++)
        {
            if(DWh.getCoef(i,j).getNV() == 1) //if there is only one variable, it is s0[reduced_nv-1]
                DWh.getCoef(i,j).evaluate(&s0[reduced_nv-1], ofs_temp, ofts_order, ofs_order);
            else
                DWh.getCoef(i,j).evaluate(s0, ofs_temp, ofts_order, ofs_order);

            mIn.getCA(i,j)->ccopy(ofs_temp);
        }
    }
}

/**
 *  \brief Evaluate the jacobian of the invariant manifold from RCM coordinates
 *         (double[]), to TFC coordinates (gsl_matrix_complex*, i.e. a matrix of
 *         complex number - gsl_complex format). note that this computation is only a
 *         partial computation, the real jacobian from RCM to TFC includes another matrix
 *         product. See evalDRCMtoNC for details.
 **/
void Invman::evalDRCMtoTFC_partial(double const st0[], double const t, gsl_matrix_complex* m1, const int ofts_order, const int ofs_order) const
{
    //------------------------------------------
    // Inner variables (CCM)
    //------------------------------------------
    cdouble s0[reduced_nv];

    //------------------------------------------
    // RCM to CCM
    //------------------------------------------
    RCMtoCCM(st0, s0, reduced_nv);

    //------------------------------------------
    // CCM to TFC
    //------------------------------------------
    this->evalDCCMtoTFC(s0, (matrix<Ofsc>&) mIn, ofts_order, ofs_order);

    //------------------------------------------
    // 3. Evaluate mIn in m2
    //------------------------------------------
    evaluate(n*t, mIn, m1);
}

/**
 *  \brief Evaluate the jacobian of the invariant manifold from RCM coordinates
 *         (double[]), to NC coordinates (gsl_matrix_complex*, i.e. a matrix of
 *          double).
 *
 *         The computation is as follows:
 *
 *         m1 = dzNC/dsRCM = Mcoc * dzTFC/dsCCM * dsCCM/dsRCM
 *
 *         where Mcoc is a private component of "this", dzTFC/dsCCM is the result of
 *         evalDRCMtoTFC_partial, and dsCCM/dsRCM = CCM_R_RCM is also a component of
 *         "this".
 *
 *         Although most of the computation is in complex form, the end result is real.
 *
 **/
void Invman::evalDRCMtoNC(double const st0[], double const t, gsl_matrix* m1, const int ofts_order, const int ofs_order) const
{
    //------------------------------------------------------------------------------------
    // Check sizes
    //------------------------------------------------------------------------------------
    if( (int) m1->size1 != 6 || (int) m1->size2 != reduced_nv)
    {
        cout << "Invman::evalDRCMtoNC. Dimension mismatch: " << endl;
        cout << "m1->size1 = " << m1->size1 << "instead of " << 6 << endl;
        cout << "m1->size2 = " << m1->size2 << "instead of " << reduced_nv << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Inner variables
    //------------------------------------------------------------------------------------
    gsl_matrix_complex* m1c = gsl_matrix_complex_calloc(6, reduced_nv);
    gsl_matrix_complex* PC  = gsl_matrix_complex_calloc(6, 6);
    gsl_matrix_complex* K1c = gsl_matrix_complex_calloc(6, reduced_nv);
    gsl_matrix_complex* K2c = gsl_matrix_complex_calloc(6, reduced_nv);

    //------------------------------------------------------------------------------------
    // RCM to TFC
    //------------------------------------------------------------------------------------
    evalDRCMtoTFC_partial(st0, t, m1c, ofts_order, ofs_order);

    //------------------------------------------------------------------------------------
    // Evaluate Mcoc
    //------------------------------------------------------------------------------------
    evaluate(t, n, Mcoc, PC);

    //------------------------------------------------------------------------------------
    //K2c = PC*m1c*CCM_R_RCM,
    // in C(6,6)*C(6,rnv)*C(rnv,rnv) = C(6,rnv)
    //------------------------------------------------------------------------------------
    //K1c = m1c*CCM_R_RCM
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), m1c, CCM_R_RCM, gslc_complex(0.0,0.0), K1c);
    //K2c = PC*m1c*CCM_R_RCM
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), PC, K1c, gslc_complex(0.0,0.0), K2c);

    //------------------------------------------------------------------------------------
    //m1 = real(K2c), in R(6,rnv)
    //------------------------------------------------------------------------------------
    gslc_matrix_complex_to_matrix(K2c, m1);

    //------------------------------------------------------------------------------------
    // Free
    //------------------------------------------------------------------------------------
    gsl_matrix_complex_free(m1c);
    gsl_matrix_complex_free(PC);
    gsl_matrix_complex_free(K1c);
    gsl_matrix_complex_free(K2c);
}

/**
 *  \brief Evaluate the jacobian of the invariant manifold from RCM coordinates
 *         (double[]), to coord_type coordinates (gsl_matrix_complex*, i.e. a matrix of
 *          double).
 *
 *         The computation uses evalDRCMtoNC to get the jacobian from RCM coordinates
 *         to NC coordinates. Then, it uses the routine rot_mat_coc to get the COC matrix
 *         to transform this matrix in coord_type coordinates. Then, it computes the
 *         product of this two matrices.
 *
 **/
void Invman::evalDRCMtoCOORD(double const st0[], double const t, gsl_matrix* m1, const int ofts_order, const int ofs_order, const int coord_type) const
{
    //------------------------------------------------------------------------------------
    // Check sizes
    //------------------------------------------------------------------------------------
    if( (int) m1->size1 != 6 || (int) m1->size2 != reduced_nv)
    {
        cout << "Invman::evalDRCMtoCOORD. Dimension mismatch: " << endl;
        cout << "m1->size1 = " << m1->size1 << "instead of "    << 6 << endl;
        cout << "m1->size2 = " << m1->size2 << "instead of "    << reduced_nv << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Computation: if the desired coord_type is different from ncs, which is the
    // normalized coordinate system (either NCSEM or NCEM), the the matrix product is
    // necessary. If it is equal, then we are back to the evalDRCMtoNC case.
    //------------------------------------------------------------------------------------
    if(ncs != coord_type)
    {
        //--------------------------------------------------------------------------------
        // Inner variables
        //--------------------------------------------------------------------------------
        gsl_matrix* NC_J_RCM = gsl_matrix_calloc(6, reduced_nv);
        gsl_matrix* COORD_R_NC = gsl_matrix_calloc(6,6);


        //--------------------------------------------------------------------------------
        // RCM to NC
        //--------------------------------------------------------------------------------
        evalDRCMtoNC(st0, t, NC_J_RCM, ofts_order, ofs_order);

        //--------------------------------------------------------------------------------
        // NC to coord_type, if necessary
        //--------------------------------------------------------------------------------
        //COORD_R_NC = COC matrix from ncs (either NCSEM or NCEM) to coord_type
        rot_mat_coc(t, COORD_R_NC, ncs, coord_type);

        //COORD_J_RCM = COORD_R_NC*NC_J_RCM, in R(6,6)*R(6,5) = R(6,5)
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, COORD_R_NC, NC_J_RCM, 0.0, m1);


        //--------------------------------------------------------------------------------
        // Free
        //--------------------------------------------------------------------------------
        gsl_matrix_free(NC_J_RCM);
        gsl_matrix_free(COORD_R_NC);
    }
    else
    {
        //--------------------------------------------------------------------------------
        //A direction computation is possible, since the result is desired in NC
        //coordinates
        //--------------------------------------------------------------------------------
        evalDRCMtoNC(st0, t, m1, ofts_order, ofs_order);
    }
}

/**
 *  \brief Evaluate the time derivative of the invariant manifold from RCM coordinates
 *         (double[]), to NC coordinates (double[]).
 **/
void Invman::evaldotRCMtoNC(double const st0[], double const t, double z1[], const int ofts_order, const int ofs_order) const
{
    //------------------------------------------------------------------------------------
    // Inner variables (CCM, TFC, NC)
    //------------------------------------------------------------------------------------
    vector<Ofsc> zIn(6), zOut(6);

    //------------------------------------------------------------------------------------
    // RCM to TFC
    //------------------------------------------------------------------------------------
    this->evalRCMtoTFC(st0, zIn, ofts_order, ofs_order);

    //------------------------------------------------------------------------------------
    // TFC to NC
    //------------------------------------------------------------------------------------
    applyCOC(Mcoc, Vcoc, zIn, zOut);

    //------------------------------------------------------------------------------------
    // zIn = dot(zOut)
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++)
    {
        zIn[i].zero();
        zIn[i].dot(zOut[i], n);
    }

    //------------------------------------------
    // Evaluate zIn in z1
    //------------------------------------------
    for(int p = 0; p < 6; p++) z1[p] = creal(zIn[p].evaluate(n*t, ofs_order));
}


/**
 *  \brief Evaluate the time derivative of an EM invariant manifold in NCSEM coordinates.
 *
 *         CAREFUL: this routine has not been properly checked, because it is not used for
 *         now.
 **/
void Invman::evaldotRCMEMtoNCSEM(double const st0[], double const t, gsl_vector *dzNCSEM,
                                 const int ofts_order, const int ofs_order,
                                 const Invman &invman_SEM) const
{
    //====================================================================================
    // Check that fwrk == F_EM, otherwise, no sense!
    //====================================================================================
    if(fwrk != F_EM)
    {
        cout << "Warning in Invman::evaldotRCMEMtoNCSEM.        " << endl;
        cout << "This routine should be used only on a invariant" << endl;
        cout << "manifold in the EM system (F_EM framework.     " << endl;
        cout << "The routine will be stopped.                   " << endl;
        cout << "Use evaldotRCMtoNC if you are in the F_SEM.    " << endl;
        return;
    }

    //====================================================================================
    // Inner variables (CCM, TFC, NC)
    //====================================================================================
    vector<Ofsc> zTFC(6), zEM(6), zVEM(6), dzVEM(6);

    //====================================================================================
    // First step: we compute d zINEM / dt, all in EM units
    //====================================================================================
    //------------------------------------------------------------------------------------
    // RCM to TFC
    //------------------------------------------------------------------------------------
    this->evalRCMtoTFC(st0, zTFC, ofts_order, ofs_order);

    //------------------------------------------------------------------------------------
    // TFC to NCEM
    //------------------------------------------------------------------------------------
    applyCOC(Mcoc, Vcoc, zTFC, zEM);

    //------------------------------------------------------------------------------------
    // NCEM to EM
    //------------------------------------------------------------------------------------
    //Step 1: translation
    for(int i = 0; i < 6; i++) zEM[i] -= nc_B_SYS[i];

    //Step 2: scaling
    zEM[0] *= -cs->gamma;
    zEM[1] *= -cs->gamma;
    zEM[2] *= +cs->gamma;
    zEM[3] *= -cs->gamma;
    zEM[4] *= -cs->gamma;
    zEM[5] *= +cs->gamma;

    //------------------------------------------------------------------------------------
    // PEM to VEM: zVEM = VSYS_R_SYS*zEM, in Fourier space
    //------------------------------------------------------------------------------------
    smvprod_ofs(VSYS_R_SYS, zEM, zVEM);

    //------------------------------------------------------------------------------------
    // dzVEM = dot(zVEM)
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 6; i++)
    {
        dzVEM[i].zero();
        dzVEM[i].dot(zVEM[i], n);
    }

    //------------------------------------------------------------------------------------
    // Then, evaluation
    //------------------------------------------------------------------------------------
    gsl_vector_complex *zVEMg  = gsl_vector_complex_alloc(6);
    gsl_vector_complex *dzVEMg = gsl_vector_complex_alloc(6);

    for(int p = 0; p < 6; p++)
    {
        gsl_vector_complex_set(zVEMg,  p, gslc_complex(zVEM[p].evaluate(n*t, ofs_order)));
        gsl_vector_complex_set(dzVEMg, p, gslc_complex(dzVEM[p].evaluate(n*t, ofs_order)));
    }


    //------------------------------------------------------------------------------------
    // Then, evaluation of the matrices INEM_R_VEM, dINEM_R_VEM and so on
    //------------------------------------------------------------------------------------
    gsl_matrix_complex* INEM_R_VEM  = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* dINEM_R_VEM = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* EM_R_INEM   = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* dEM_R_INEM  = gsl_matrix_complex_alloc(6,6);

    gsl_vector_complex *INEM_B_EM   = gsl_vector_complex_alloc(6);
    gsl_vector_complex *dINEM_B_EM  = gsl_vector_complex_alloc(6);

    eval_IN_B_SYS(t, INEM_B_EM, dINEM_B_EM, INEM_R_VEM, dINEM_R_VEM,
                     EM_R_INEM, dEM_R_INEM, ofs_order);

    //------------------------------------------------------------------------------------
    // Then, we evaluate zINEM and dzINEM/dt
    //------------------------------------------------------------------------------------
    gsl_vector_complex *zINEM   = gsl_vector_complex_alloc(6);
    gsl_vector_complex *dzINEM  = gsl_vector_complex_alloc(6);
    gsl_vector_complex *temp1   = gsl_vector_complex_alloc(6);
    gsl_vector_complex *temp2   = gsl_vector_complex_alloc(6);

    // zINEM = INEM_R_VEM * zVEMg - INEM_B_EM
    gsl_blas_zgemv(CblasNoTrans, gslc_complex(1.0,0.0), INEM_R_VEM, zVEMg, gslc_complex(0.0,0.0), zINEM);
    gsl_vector_complex_sub(zINEM, INEM_B_EM);

    //temp1 = dINEM_R_VEM * zVEMg
    gsl_blas_zgemv(CblasNoTrans, gslc_complex(1.0,0.0), dINEM_R_VEM, zVEMg, gslc_complex(0.0,0.0), temp1);

    //temp2 = INEM_R_VEM * dzVEMg
    gsl_blas_zgemv(CblasNoTrans, gslc_complex(1.0,0.0), INEM_R_VEM, dzVEMg, gslc_complex(0.0,0.0), temp2);

    //Then dzINEM = temp1 + temp2 - dINEM_B_EM
    gsl_vector_complex_memcpy(dzINEM, temp1);
    gsl_vector_complex_add(dzINEM, temp2);
    gsl_vector_complex_sub(dzINEM, dINEM_B_EM);

    //====================================================================================
    // Second step: we use the invariant manifold from SEM to finish the change of coordinates
    //====================================================================================

    //------------------------------------------------------------------------------------
    // First, we need to take into account the change of units
    //------------------------------------------------------------------------------------
    double tsem = t*SEML.us_em.ns; //now the time is in SEM units

    //------------------------------------------------------------------------------------
    // Then, INEM -> INSEM
    //------------------------------------------------------------------------------------
    gsl_vector_complex *zINSEM   = gsl_vector_complex_alloc(6);
    gsl_vector_complex *dzINSEM  = gsl_vector_complex_alloc(6);

    cdouble ctemp = 0.0;
    for(int i = 0; i < 6; i++)
    {
        //--------------------------------------------------------------------------------
        // zINEM
        //--------------------------------------------------------------------------------
        //Take the complex value
        ctemp = gslc_complex(gsl_vector_complex_get(zINEM, i));

        //Change of units
        if(i < 3) ctemp *= 1.0/SEML.us_em.as;
        else ctemp *= 1.0/(SEML.us_em.as*SEML.us_em.ns);

        //Update dzINSEM
        gsl_vector_complex_set(zINSEM, i, gslc_complex(ctemp));

        //--------------------------------------------------------------------------------
        // dzINEM
        //--------------------------------------------------------------------------------
        //Take the complex value
        ctemp = gslc_complex(gsl_vector_complex_get(dzINEM, i));

        //Change of units
        if(i < 3) ctemp *= 1.0/(SEML.us_em.as*SEML.us_em.ns);
        else ctemp *= 1.0/(SEML.us_em.as*SEML.us_em.ns*SEML.us_em.ns);

        //Update dzINSEM
        gsl_vector_complex_set(dzINSEM, i, gslc_complex(ctemp));
    }


    //------------------------------------------------------------------------------------
    // Then, we use the invariant manifold at SEML
    //------------------------------------------------------------------------------------
    gsl_matrix_complex* INSEM_R_VSEM  = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* dINSEM_R_VSEM = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* SEM_R_INSEM   = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* dSEM_R_INSEM  = gsl_matrix_complex_alloc(6,6);
    gsl_vector_complex *INSEM_B_SEM   = gsl_vector_complex_alloc(6);
    gsl_vector_complex *dINSEM_B_SEM  = gsl_vector_complex_alloc(6);

    invman_SEM.eval_IN_B_SYS(tsem, INSEM_B_SEM, dINSEM_B_SEM, INSEM_R_VSEM,
                                  dINSEM_R_VSEM, SEM_R_INSEM, dSEM_R_INSEM,
                                  ofs_order);


    //------------------------------------------------------------------------------------
    // Finally:
    //
    //  dzncsem = sem_R_SEM * ( dSEM_R_INSEM * zINSEM + SEM_R_INSEM * dzINSEM ) +  dSEM_B_sem
    //------------------------------------------------------------------------------------
    cdouble dzNCSEMc[6];

    //------------------------------------------------------------------------------------
    // SEM to NCSEM
    //------------------------------------------------------------------------------------
    gsl_vector_complex *dzSEMg  = gsl_vector_complex_alloc(6);

    // dzSEMg = dSEM_R_INSEM * zINSEM
    gsl_blas_zgemv(CblasNoTrans, gslc_complex(1.0,0.0), dSEM_R_INSEM, zINSEM, gslc_complex(0.0,0.0), dzSEMg);

    // dzSEMg += SEM_R_INSEM * dzINSEM
    gsl_blas_zgemv(CblasNoTrans, gslc_complex(1.0,0.0), SEM_R_INSEM, dzINSEM, gslc_complex(1.0,0.0), dzSEMg);

    //Copy into dzNCSEMc
    for(int i = 0; i < 6; i++) dzNCSEMc[i] = gslc_complex(gsl_vector_complex_get(dzSEMg, i));

    //------------------------------------------------------------------------------------
    // SEM to NCSEM: sem_R_SEM * dzSEM +  dSEM_B_sem
    //------------------------------------------------------------------------------------
    //Step 1: scaling
    dzNCSEMc[0] /= -invman_SEM.cs->gamma;
    dzNCSEMc[1] /= -invman_SEM.cs->gamma;
    dzNCSEMc[2] /= +invman_SEM.cs->gamma;
    dzNCSEMc[3] /= -invman_SEM.cs->gamma;
    dzNCSEMc[4] /= -invman_SEM.cs->gamma;
    dzNCSEMc[5] /= +invman_SEM.cs->gamma;

    //Step 2: translation (zTFC[0] is used as a temporary object)
    // dzNCSEM is finally computed.
    for(int i = 0; i < 6; i++)
    {
        zTFC[0].dot(invman_SEM.nc_B_SYS[i], invman_SEM.cs->us.n);
        dzNCSEMc[i] +=  zTFC[0].evaluate(n*t);  //note: n*t is adimensionalized
        gsl_vector_set(dzNCSEM, i, creal(dzNCSEMc[i]));
    }


    //------------------------------------------------------------------------------------
    // Free
    //------------------------------------------------------------------------------------
    gsl_vector_complex_free(zVEMg);
    gsl_vector_complex_free(dzVEMg);

    gsl_matrix_complex_free(INEM_R_VEM);
    gsl_matrix_complex_free(dINEM_R_VEM);
    gsl_matrix_complex_free(EM_R_INEM );
    gsl_matrix_complex_free(dEM_R_INEM);

    gsl_vector_complex_free(INEM_B_EM );
    gsl_vector_complex_free(dINEM_B_EM);

    gsl_vector_complex_free(zINEM );
    gsl_vector_complex_free(dzINEM);
    gsl_vector_complex_free(temp1 );
    gsl_vector_complex_free(temp2 );

    gsl_vector_complex_free(zINSEM );
    gsl_vector_complex_free(dzINSEM);

    gsl_matrix_complex_free(INSEM_R_VSEM);
    gsl_matrix_complex_free(dINSEM_R_VSEM);
    gsl_matrix_complex_free(SEM_R_INSEM );
    gsl_matrix_complex_free(dSEM_R_INSEM);
    gsl_vector_complex_free(INSEM_B_SEM );
    gsl_vector_complex_free(dINSEM_B_SEM);

    gsl_vector_complex_free(dzSEMg);
}

//========================================================================================
// Evaluate other matrices & vectors
//========================================================================================
/**
 *  \brief Evalute the vector & matrices IN_B_SYS that appears in the COC between the
 *         system coordinates and the inertial coordinates, in native units.
 *         IN_B_SYS is non null only in the F_EM (Earth-Moon framework).
 *         The time derivative of this matrix is also evaluated.
 **/
void Invman::eval_IN_B_SYS(double const t,
                           gsl_vector_complex* IN_B_SYS,
                           gsl_vector_complex* dIN_B_SYS,
                           gsl_matrix_complex* IN_R_VSYS,
                           gsl_matrix_complex* dIN_R_VSYS,
                           gsl_matrix_complex* SYS_R_IN,
                           gsl_matrix_complex* dSYS_R_IN,
                           const int ofs_order) const
{
    gsl_matrix_complex* VSYS_R_IN  = gsl_matrix_complex_alloc(6,6);
    gsl_matrix_complex* dVSYS_R_IN = gsl_matrix_complex_alloc(6,6);


    //------------------------------------------------------------------------------------
    // Initialisation and computing z & Z
    //------------------------------------------------------------------------------------
    double n  = cs->us.n;
    double ms = cs->us.ms;
    double mm = cs->us.mm;
    double me = cs->us.me;
    double ns = cs->us.ns;
    double as = cs->us.as;
    double ni = cs->us.ni;
    double ai = cs->us.ai;

    //z
    cdouble z    = evz(cs->zt, t, n, ni, ai);
    cdouble zd   = evzdot(cs->zt, cs->ztdot, t, n, ni, ai);
    cdouble zdd  = evzddot(cs->zt, cs->ztdot, cs->ztddot, t, n, ni, ai);

    //Z
    cdouble Z    = evz(cs->Zt, t, n, ns, as);
    cdouble Zd   = evzdot(cs->Zt, cs->Ztdot, t, n, ns, as);
    cdouble Zdd  = evzddot(cs->Zt, cs->Ztdot, cs->Ztddot, t, n, ns, as);

    //Inner variables
    double r1, r2, dr1, dr2, ddr1, ddr2;

    //------------------------------------------------------------------------------------
    // Switch between the framework
    //------------------------------------------------------------------------------------
    switch(fwrk)
    {
    case F_EM:
    {
        //--------------------------------------------------------------------------------
        //z is selected for the rest of the computation
        //--------------------------------------------------------------------------------
        r1 = creal(z);
        r2 = cimag(z);

        dr1 = creal(zd);
        dr2 = cimag(zd);

        ddr1 = creal(zdd);
        ddr2 = cimag(zdd);

        //--------------------------------------------------------------------------------
        //IN_B_SYS = ms/(mm + me + ms) * (R1 R2 0 dot(R1) dot(R2) 0)^T
        //--------------------------------------------------------------------------------
        double factor = ms/(ms + mm + me);
        gsl_vector_complex_set(IN_B_SYS, 0, gslc_complex(factor*creal(Z), 0.0));
        gsl_vector_complex_set(IN_B_SYS, 1, gslc_complex(factor*cimag(Z), 0.0));

        gsl_vector_complex_set(IN_B_SYS, 3, gslc_complex(factor*creal(Zd), 0.0));
        gsl_vector_complex_set(IN_B_SYS, 4, gslc_complex(factor*cimag(Zd), 0.0));

        //--------------------------------------------------------------------------------
        //dIN_B_SYS = ms/(mm + me + ms) * (dot(R1) dot(R2) 0 ddot(R1) ddot(R2) 0)^T
        //--------------------------------------------------------------------------------
        gsl_vector_complex_set(dIN_B_SYS, 0, gslc_complex(factor*creal(Zd), 0.0));
        gsl_vector_complex_set(dIN_B_SYS, 1, gslc_complex(factor*cimag(Zd), 0.0));

        gsl_vector_complex_set(dIN_B_SYS, 3, gslc_complex(factor*creal(Zdd), 0.0));
        gsl_vector_complex_set(dIN_B_SYS, 4, gslc_complex(factor*cimag(Zdd), 0.0));

        break;
    }

    case F_SEM:
    default:
    {
        //--------------------------------------------------------------------------------
        //Z is selected for the rest of the computation
        //--------------------------------------------------------------------------------
        r1 = creal(Z);
        r2 = cimag(Z);

        dr1 = creal(Zd);
        dr2 = cimag(Zd);

        ddr1 = creal(Zdd);
        ddr2 = cimag(Zdd);
        break;
    }


    }


    //------------------------------------------------------------------------------------
    //                 |  +b1   +b2   0    0    0   0  |
    //                 |  -b2   +b1   0    0    0   0  |
    //                 |   0     0    a6   0    0   0  |
    //  VSYS_R_IN  =   | +db1  +db2   0    b1  +b2  0  |
    //                 | -db2  +db1   0   -b2  +b1  0  |
    //                 |   0     0  +da6   0    0  +b1 |
    //
    // with
    //          b1 = r1/r^2
    //          b2 = r2/r^2
    //          a6 = 1/r  = alpha6 (or delta6)
    //
    //  and db1 = dot(b1)
    //------------------------------------------------------------------------------------
    double r  = sqrt(r1*r1+r2*r2);
    double dr = 1.0/r*(r1*dr1 + r2*dr2);

    double b1 = r1/(r*r);
    double b2 = r2/(r*r);

    double db1 = (dr1 * r - 2 * dr * r1)/(r*r*r);
    double db2 = (dr2 * r - 2 * dr * r2)/(r*r*r);

    double a6  = 1.0/r;
    double da6 = -dr/(r*r);

    gsl_matrix_complex_set(VSYS_R_IN, 0, 0, gslc_complex( b1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 0, 1, gslc_complex( b2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 1, 0, gslc_complex(-b2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 1, 1, gslc_complex( b1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 2, 2, gslc_complex( a6, 0.0));

    gsl_matrix_complex_set(VSYS_R_IN, 3, 0, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 3, 1, gslc_complex( db2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 0, gslc_complex(-db2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 1, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 5, 2, gslc_complex( da6, 0.0));

    gsl_matrix_complex_set(VSYS_R_IN, 3, 3, gslc_complex( b1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 3, 4, gslc_complex( b2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 3, gslc_complex(-b2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 4, gslc_complex( b1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 5, 5, gslc_complex( a6, 0.0));


    //------------------------------------------------------------------------------------
    //                 |  +db1   +db2   0     0    0    0  |
    //                 |  -db2   +db1   0     0    0    0  |
    //                 |   0      0    da6    0    0    0  |
    //  dVSYS_R_IN =   | +ddb1  +ddb2   0    db1  +db2  0  |
    //                 | -ddb2  +ddb1   0   -db2  +db1  0  |
    //                 |   0      0   +dda6   0    0  +db1 |
    //------------------------------------------------------------------------------------
    double ddr  = 1.0/(r*r)*( (dr1*dr1 + r1 *ddr1 + dr2*dr2 + r2*ddr2)*r - dr*(r1*dr1 + r2*dr2 ));
    double ddb1 = (r* r*ddr1 - 2*r*( 2*dr1*dr + r1*ddr ) + 6*r1* dr*dr)/(r* r* r*r);
    double ddb2 = (r* r*ddr2 - 2*r*( 2*dr2*dr + r2*ddr ) + 6*r2* dr*dr)/(r* r* r*r);
    double dda6 = (2*dr*dr - r*ddr)/(r* r*r);

    gsl_matrix_complex_set(dVSYS_R_IN, 0, 0, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 0, 1, gslc_complex( db2, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 1, 0, gslc_complex(-db2, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 1, 1, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 2, 2, gslc_complex( da6, 0.0));

    gsl_matrix_complex_set(dVSYS_R_IN, 3, 0, gslc_complex( ddb1, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 3, 1, gslc_complex( ddb2, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 4, 0, gslc_complex(-ddb2, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 4, 1, gslc_complex( ddb1, 0.0));
    gsl_matrix_complex_set(dVSYS_R_IN, 5, 2, gslc_complex( dda6, 0.0));

    gsl_matrix_complex_set(VSYS_R_IN, 3, 3, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 3, 4, gslc_complex( db2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 3, gslc_complex(-db2, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 4, 4, gslc_complex( db1, 0.0));
    gsl_matrix_complex_set(VSYS_R_IN, 5, 5, gslc_complex( da6, 0.0));


    //------------------------------------------------------------------------------------
    //                 |  +r1   -r2   0    0    0   0  |
    //                 |  +r2   +r1   0    0    0   0  |
    //                 |   0     0    r    0    0   0  |
    //  IN_R_VSYS  =   | +dr1  -dr2   0    r1  -r2  0  |
    //                 | +dr2  +dr1   0   +r2  +r1  0  |
    //                 |   0     0  +dr    0    0  +r  |
    //------------------------------------------------------------------------------------
    gsl_matrix_complex_set(IN_R_VSYS, 0, 0, gslc_complex( r1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 0, 1, gslc_complex(-r2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 1, 0, gslc_complex( r2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 1, 1, gslc_complex( r1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 2, 2, gslc_complex( r , 0.0));

    gsl_matrix_complex_set(IN_R_VSYS, 3, 0, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 3, 1, gslc_complex(-dr2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 4, 0, gslc_complex( dr2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 4, 1, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 5, 2, gslc_complex( dr , 0.0));

    gsl_matrix_complex_set(IN_R_VSYS, 3, 3, gslc_complex( r1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 3, 4, gslc_complex(-r2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 4, 3, gslc_complex( r2, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 4, 4, gslc_complex( r1, 0.0));
    gsl_matrix_complex_set(IN_R_VSYS, 5, 5, gslc_complex( r , 0.0));

    //------------------------------------------------------------------------------------
    //                 |  +dr1   -dr2   0     0     0    0  |
    //                 |  +dr2   +dr1   0     0     0    0  |
    //                 |   0      0     dr    0     0    0  |
    //  dIN_R_VSYS =   | +ddr1  -ddr2   0     dr1  -dr2  0  |
    //                 | +ddr2  +ddr1   0    +dr2  +dr1  0  |
    //                 |   0      0    +ddr   0     0   +dr |
    //------------------------------------------------------------------------------------
    gsl_matrix_complex_set(dIN_R_VSYS, 0, 0, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 0, 1, gslc_complex(-dr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 1, 0, gslc_complex( dr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 1, 1, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 2, 2, gslc_complex( dr , 0.0));

    gsl_matrix_complex_set(dIN_R_VSYS, 3, 0, gslc_complex( ddr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 3, 1, gslc_complex(-ddr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 4, 0, gslc_complex( ddr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 4, 1, gslc_complex( ddr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 5, 2, gslc_complex( ddr , 0.0));

    gsl_matrix_complex_set(dIN_R_VSYS, 3, 3, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 3, 4, gslc_complex(-dr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 4, 3, gslc_complex( dr2, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 4, 4, gslc_complex( dr1, 0.0));
    gsl_matrix_complex_set(dIN_R_VSYS, 5, 5, gslc_complex( dr , 0.0));


    //------------------------------------------------------------------------------------
    //  SYS_R_IN  = SYS_R_VSYS * VSYS_R_IN
    //------------------------------------------------------------------------------------
    // SYS_R_VSYSg = SYS_R_VSYS in GSL format
    gsl_matrix_complex* SYS_R_VSYSg = gsl_matrix_complex_alloc(6,6);
    evaluate(n*t, SYS_R_VSYS, SYS_R_VSYSg);

    //SYS_R_IN  = SYS_R_VSYS * VSYS_R_IN
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), SYS_R_VSYSg, VSYS_R_IN, gslc_complex(0.0,0.0), SYS_R_IN);

    //------------------------------------------------------------------------------------
    //  dSYS_R_IN  = dSYS_R_VSYS * VSYS_R_IN + SYS_R_VSYS * dVSYS_R_IN
    //------------------------------------------------------------------------------------
    Ofsc ofs_temp(ofs_order);
    gsl_matrix_complex* dSYS_R_VSYSg = gsl_matrix_complex_alloc(6,6);
    evaluatedot(t, n, SYS_R_VSYS, dSYS_R_VSYSg, ofs_temp);

    //dSYS_R_IN  = dSYS_R_VSYS * VSYS_R_IN
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), dSYS_R_VSYSg, VSYS_R_IN, gslc_complex(0.0,0.0), dSYS_R_IN);

    //dSYS_R_IN  += SYS_R_VSYS * dVSYS_R_IN
    gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, gslc_complex(1.0,0.0), SYS_R_VSYSg, dVSYS_R_IN, gslc_complex(1.0,0.0), dSYS_R_IN);

    gsl_matrix_complex_free( SYS_R_VSYSg  );
    gsl_matrix_complex_free( dSYS_R_VSYSg );
    gsl_matrix_complex_free( VSYS_R_IN    );
    gsl_matrix_complex_free( dVSYS_R_IN    );
}

//========================================================================================
//          Projection on (un)stable manifold
//========================================================================================
/**
 * \brief Projects in sti a given NC state on the corresponding center manifold given
 *        by this, seen as a center manifold in graph style. If the graph style
 *        is not detected, a warning message is displayed and nothing is done.
 **/
void Invman::NCprojCCMtoCM(double* yv, double tv, double sti[5])
{
    //------------------------------------------------------------------------------------
    //Check that the graph style is used in this->invman
    //------------------------------------------------------------------------------------
    if(this->pmType != PMS_GRAPH)
    {
        cout << "Invman::NCprojCCMtoCM. Impossible to use if the invariant manifold ";
        cout << "of the current orbit is not based on the graph style. return." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    //Get the closest point on the center manifold, in scp[4]
    //------------------------------------------------------------------------------------
    cdouble scp[4];
    NCprojCCM(yv, tv, this->n, OFS_ORDER, this->MIcoc,
              this->Vcoc, this->omega1, this->omega3, scp, 4);

    //------------------------------------------------------------------------------------
    //Get the correspondance in RCM coordinates
    //------------------------------------------------------------------------------------
    CCMtoRCM(scp, sti, 4);
}

//========================================================================================
// Getters
//========================================================================================
const vector<Oftsc>& Invman::getWh() const
{
    return Wh;
}

const vector<Oftsc>& Invman::getW() const
{
    return W;
}

const matrix<Oftsc>& Invman::getDWh() const
{
    return DWh;
}

const vector<Ofsc>& Invman::getVcoc() const
{
    return Vcoc;
}

const matrix<Ofsc>& Invman::getMcoc() const
{
    return Mcoc;
}

const matrix<Ofsc>& Invman::getMIcoc() const
{
    return MIcoc;
}

const vector<Oftsc>& Invman::getHy() const
{
    return Hy;
}

const CSYS* Invman::getCS() const
{
    return cs;
}

double Invman::getOmega1() const
{
    return omega1;
}

double Invman::getOmega3() const
{
    return omega3;
}

double Invman::getN() const
{
    return n;
}

int Invman::getPmType() const
{
    return pmType;
}

int Invman::getManType() const
{
    return manType;
}

int Invman::getRnv() const
{
    return reduced_nv;
}

int Invman::getFwrk() const
{
    return fwrk;
}

int Invman::getNCS() const
{
    return ncs;
}

int Invman::getSCS() const
{
    return scs;
}

int Invman::getOFTSORDER() const
{
    return ofts_order;
}


//========================================================================================
// Energy
//========================================================================================
/**
 *   \brief Update the array st0 with values from orbit.s0 and the value sr. Used in orbit_init_pmap and deltham.
 *
 *          The following patterns are followed (examples):
 *
 *             - if orbit.dim = 2: s2 = s4 = sr
 *                   st0[0] = sr;
 *                   st0[1] = refineH->st0[1];
 *                   st0[2] = sr;
 *                   st0[3] = refineH->st0[3];
 **/
void update_s0_direction(RefineH* refineH, double st0[], double sr)
{
    switch(refineH->vdim)
    {
    case 3:
    {
        st0[0] = refineH->st0[0];
        st0[1] = 0.0;
        st0[2] = -sr;
        st0[3] = 0.0;
        break;
    }

    case 1:
    {
        st0[0] = sr;
        st0[1] = refineH->st0[1];
        st0[2] = sr;
        st0[3] = refineH->st0[3];
        break;
    }
    case 2:
    case 4:
    {
        st0[0] = refineH->st0[0];
        st0[1] = sr;
        st0[2] = refineH->st0[2];
        st0[3] = 0.0;//sr;
        break;
    }

    default: //as case 2
    {
        st0[0] = refineH->st0[0];
        st0[1] = sr;
        st0[2] = refineH->st0[2];
        st0[3] = sr;
        break;
    }
    }

    //st0[4] = refineH->st0[4];
}

/**
 *   \brief Computes the hamiltonian at the position st0, in outputType coordinates and units.
 **/
double Invman::H_SYS(double st0[], double t0)
{
    //------------------------------------------------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------------------------------------------------
    double zNC[6];

    //------------------------------------------------------------------------------------
    // Evaluate the manifold
    //------------------------------------------------------------------------------------
    evalRCMtoNC(st0, t0, zNC, ofts_order, ofs_order);

    //------------------------------------------------------------------------------------
    // Energy
    //------------------------------------------------------------------------------------
    return qbcp_H_complete(t0, zNC, ncs, scs);
}

/**
 *   \brief Computes the hamiltonian at the position st0, in outputType coordinates and units, at order ofts_order_0
 **/
double Invman::H_SYS(double st0[], double t0, int ofts_order_0)
{
    //------------------------------------------------------------------------------------
    // Inner variables (NC, TFC)
    //------------------------------------------------------------------------------------
    double zNC[6];

    //------------------------------------------------------------------------------------
    // Evaluate the manifold
    //------------------------------------------------------------------------------------
    evalRCMtoNC(st0, t0, zNC, ofts_order_0, ofs_order);

    //------------------------------------------------------------------------------------
    // Energy
    //------------------------------------------------------------------------------------
    return qbcp_H_complete(t0, zNC, ncs, scs);
}

/**
 *   \brief  Computes the difference between an given Ham value and the state configuration defined in the routine update_s0 (see comments therein).
 *   \param  sr a double to complete the current tested configuration
 *   \param  params a pointer to the orbit with a given Ham value (why void? the idea is to a have a generic function, but might be useless at this point)
 *   \return the difference between the two hamiltonians
 **/
double deltham(double sr, void* params)
{
    RefineH* refineH = (RefineH*) params;
    double st0[5];
    double Hv;

    //Update st0
    update_s0_direction(refineH, st0, sr);

    //Return the hamiltonian value
    Hv = refineH->invman->H_SYS(st0, refineH->t0, refineH->order);

    return (Hv - refineH->Hv);
}

/**
 *   \brief Initialize an orbit wrt a Poincare map so that H(orbit.s0) = H(Pmap)
 **/
int init_s0_energy(RefineH* refineH, double st0[], double t0)
{

    //====================================================================================
    //Root bracketing
    //====================================================================================

    //------------------------------------------------------------------------------------
    // Initialization
    //------------------------------------------------------------------------------------
    gsl_function F;                  //called function in the root finding routine
    double s_low, s_high;            //variables for root bracketing
    double Hup, Hdown, r, fy;        //Energy values and root
    int status;
    int itermax = 100;

    //Initialization of F
    F.function = &deltham; //Energy delta with respect to target
    F.params   = refineH;  //refineH structure contains the parameters

    double snorm =  sqrt(st0[0]*st0[0] + st0[1]*st0[1] + st0[2]*st0[2] + st0[3]*st0[3]);

    //------------------------------------------------------------------------------------
    // We diminish the order of refineH during the bracketing,
    // because we do not need a precise value for the energy
    //------------------------------------------------------------------------------------
    int order0 = refineH->order;
    refineH->order = 5;

    //------------------------------------------------------------------------------------
    //Bracketing the root
    //------------------------------------------------------------------------------------
    s_low  = 0.0;
    s_high = 1e-6;  //arbitrary small value maybe multiplied by a given number?

    //Increasing s_high until a root is found (Hup*Hdown < 0)
    int iter = 0;
    do
    {
        //Evaluating the bracket
        Hup   = deltham (s_high , refineH);
        Hdown = deltham (s_low,   refineH);
        s_high *= 1.1;
    }
    while(Hup*Hdown > 0 && iter < itermax && s_high < 10*snorm);

    //Swap the two variables if needed
    if(s_low > s_high)
    {
        s_low  = s_low + s_high;
        s_high = s_low - s_high;
        s_low  = s_low - s_high;
    }

    // TURN OFF GSL ERROR ERROR HANDLER
    // (GSL JUST PRINT ERROR MESSAGE AND KILL THE PROGRAM IF FLAG IS ON)
    gsl_set_error_handler_off();


    //------------------------------------------------------------------------------------
    // We put back the original order of refineH
    //------------------------------------------------------------------------------------
    refineH->order = order0;

    //====================================================================================
    //Root finding
    //====================================================================================
    //If a root is found, refine root
    if(Hup*Hdown < 0)
    {
        //Setting the solver
        status = gsl_root_fsolver_set (refineH->s_root, &F, s_low, s_high);
        //Loop
        iter = 0;
        do
        {
            status = gsl_root_fsolver_iterate (refineH->s_root);         //updating the solver
            r      = gsl_root_fsolver_root (refineH->s_root);            //updating the root
            fy     = deltham(r, refineH);                                //Checking convergence
            status = gsl_root_test_residual (fy , refineH->eps_root);    //Checking convergence
        }
        while (status == GSL_CONTINUE && (++iter)< 50);

        if(status == GSL_SUCCESS)
        {
            //Update st0
            update_s0_direction(refineH, st0, r);
            return GSL_SUCCESS;
        }
        else
        {
            //cout <<  "init_s0_energy: No refined root was found (1). The IC are unchanged" << endl;
            return GSL_FAILURE;
        }
    }
    else
    {
        //cout <<  "init_s0_energy: No root was found (2). The IC are unchanged" << endl;
        return GSL_FAILURE;
    }
}


//========================================================================================
// Static
//========================================================================================
/**
 *  \brief Returns the number of variables associated to the coordinate system csys.
 *
 *         It is probably not the best way to compute this value that is fundamentally
 *         associated to the parameterization, and not to the coordinate system,
 *         but it allows to limit the number of inputs in the constructor.
 **/
int Invman::compRNV(CSYS& csys)
{
    switch(csys.manType)
    {
    case MAN_CENTER:
        return 4;
        break;
    case MAN_CENTER_S:
    case MAN_CENTER_U:
        return 5;
        break;
    case MAN_CENTER_US:
        return 6;
        break;
    default:
        cout << "Invman::compRNV. Unknown manifold type. 4 is returned." << endl;
        return 4;
    }
}


//========================================================================================
//
//          SUBROUTINES
//
//========================================================================================
/**
 *  \brief Set the matrix CCM_R_RCM_C  as the rotation matrix between CCM and RCM
 *         coordinates in a center manifold (4 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER(gsl_matrix_complex* CCM_R_RCM_C)
{
    //------------------------------------------------------------------------------------
    // Check that the dimensions are okay
    //------------------------------------------------------------------------------------
    if(CCM_R_RCM_C->size1 != 4 || CCM_R_RCM_C->size2 != 4)
    {
        cout << "rotmat_CC_R_RCM_CENTER. Dimension mismatch." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Update
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        gsl_matrix_complex_set(CCM_R_RCM_C, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    }

    gsl_matrix_complex_set(CCM_R_RCM_C, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_C, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_C, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_C, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));
}

/**
 *  \brief Set the matrix CCM_R_RCM_CH  as the rotation matrix between CCM and RCM
 *         coordinates in a center-hyperbolic manifold (5 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER_HYP(gsl_matrix_complex* CCM_R_RCM_CH)
{
    //------------------------------------------------------------------------------------
    // Check that the dimensions are okay
    //------------------------------------------------------------------------------------
    if(CCM_R_RCM_CH->size1 != 5 || CCM_R_RCM_CH->size2 != 5)
    {
        cout << "rotmat_CC_R_RCM_CENTER_HYP. Dimension mismatch." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Update
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        gsl_matrix_complex_set(CCM_R_RCM_CH, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    }
    gsl_matrix_complex_set(CCM_R_RCM_CH, 4, 4, gslc_complex(1.0, 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));
}

/**
 *  \brief Set the matrix CCM_R_RCM_CH  as the rotation matrix between CCM and RCM
 *         coordinates in the complete phase space manifold (6 dimensions).
 **/
void rotmat_CC_R_RCM_CENTER_6(gsl_matrix_complex* CCM_R_RCM_CH)
{
    //------------------------------------------------------------------------------------
    // Check that the dimensions are okay
    //------------------------------------------------------------------------------------
    if(CCM_R_RCM_CH->size1 != 5 || CCM_R_RCM_CH->size2 != 5)
    {
        cout << "rotmat_CC_R_RCM_CENTER_HYP. Dimension mismatch." << endl;
        return;
    }

    //------------------------------------------------------------------------------------
    // Update
    //------------------------------------------------------------------------------------
    for(int i = 0; i < 4; i++)
    {
        gsl_matrix_complex_set(CCM_R_RCM_CH, i, i, gslc_complex(1.0/sqrt(2), 0.0));
    }
    gsl_matrix_complex_set(CCM_R_RCM_CH, 4, 4, gslc_complex(1.0, 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 5, 5, gslc_complex(1.0, 0.0));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 0, 2, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 1, 3, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 2, 0, gslc_complex(0.0, -1.0/sqrt(2)));
    gsl_matrix_complex_set(CCM_R_RCM_CH, 3, 1, gslc_complex(0.0, -1.0/sqrt(2)));
}


//========================================================================================
// Tests
//========================================================================================
/**
 *  \brief Test of the routine Invman::evalCCMtoTFC. From CCM coordinates
 *         to TFC coordinates
 **/
void test_evalCCMtoTFC()
{
    //Reduced number of variables
    int reduced_nv = Invman::compRNV(*SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via Invman
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via initOFS
    //------------------------------------------------------------------------------------
    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc> CM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_TFC;     //center manifold in TFC coordinates
    //Init routine
    initOFTS(CM_NC, CM_TFC, reduced_nv, OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //------------------------------------------------------------------------------------
    //Initial conditions (random)
    //------------------------------------------------------------------------------------
    cdouble s0[reduced_nv];
    for(int i = 0; i < reduced_nv-1; i++) s0[i] = 0.01 * (rand() % 10);
    s0[reduced_nv-1] = 0.001 * (rand() % 10); //10% on the hyperbolic part

    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    vector<Ofsc> zIn(6, Ofsc(OFS_ORDER));
    cdouble zOld1[6], zOld0[6], zNew[6];

    cout << "---------------------------------------------------" << endl;
    cout << "   Comparison of old and new version of CCMtoTFC   " << endl;
    cout << "---------------------------------------------------" << endl;

    tic();
    CCMtoTFC(s0, OFTS_ORDER, OFS_ORDER, CM_TFC, zIn, 0);
    cout << "With old version, isGS = 0 (reference), in " << toc() << endl;
    for(int i = 0; i < 6; i++) zOld0[i] = zIn[i].evaluate(SEML.cs->us.n);
    //vector_complex_printf_prec(zOld0, 6);
    cout << endl;


    tic();
    CCMtoTFC(s0, OFTS_ORDER, OFS_ORDER, CM_TFC, zIn, 1);
    cout << "With old version, isGS = 1, in " << toc() << endl;
    for(int i = 0; i < 6; i++) zOld1[i] = zIn[i].evaluate(SEML.cs->us.n);
    //vector_complex_printf_prec(zOld1, 6);
    cout << "Delta with the reference = " << DENorm(zOld0, zOld1, 6) << endl;
    cout << endl;

    tic();
    invman.evalCCMtoTFC(s0, zIn, OFTS_ORDER, OFS_ORDER);
    cout << "With new version, in " << toc() << endl;
    for(int i = 0; i < 6; i++) zNew[i] = zIn[i].evaluate(SEML.cs->us.n);
    //vector_complex_printf_prec(zNew, 6);
    cout << "Delta with the reference = " << DENorm(zOld0, zNew, 6) << endl;
    cout << endl;
}

/**
 *  \brief Test of the routine Invman::evalDCCMtoTFC.
 **/
void test_evalDCCMtoTFC()
{
    //Reduced number of variables
    int reduced_nv = Invman::compRNV(*SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via Invman
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via initOFS
    //------------------------------------------------------------------------------------
    //--------------------------------------
    //DCM_TFC : Jacobian of CM_TFC
    //--------------------------------------
    matrix<Oftsc> DCM_TFC(6, reduced_nv, reduced_nv, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_TFC, SEML.cs->F_PMS+"DWf/DWhc");

    //------------------------------------------------------------------------------------
    //Initial conditions (random)
    //------------------------------------------------------------------------------------
    cdouble s0[reduced_nv];
    for(int i = 0; i < reduced_nv-1; i++) s0[i] = 0.01 * (rand() % 10);
    s0[reduced_nv-1] = 0.001 * (rand() % 10); //10% on the hyperbolic part


    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    matrix<Ofsc> mIn(6, reduced_nv, OFS_ORDER);
    cdouble** zOld1 = dcmatrix(0, 6, 0, reduced_nv-1);
    cdouble** zOld0 = dcmatrix(0, 6, 0, reduced_nv-1);
    cdouble** zNew  = dcmatrix(0, 6, 0, reduced_nv-1);

    cout << "---------------------------------------------------" << endl;
    cout << "   Comparison of old and new version of DCCMtoTFC  " << endl;
    cout << "---------------------------------------------------" << endl;

    tic();
    CCMtoTFC_JAC(s0, OFTS_ORDER, OFS_ORDER, DCM_TFC, mIn, 0);
    cout << "With old version, isGS = 0 (reference), in " << toc() << endl;
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < reduced_nv; j++)  zOld0[i][j] = mIn.getCoef(i,j).evaluate(SEML.cs->us.n);
    }
    //matrix_complex_printf((cdouble **)zOld0, 6, reduced_nv);
    cout << endl;


    tic();
    CCMtoTFC_JAC(s0, OFTS_ORDER, OFS_ORDER, DCM_TFC, mIn, 1);
    cout << "With old version, isGS = 1, in " << toc() << endl;
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < reduced_nv; j++)  zOld1[i][j] = mIn.getCoef(i,j).evaluate(SEML.cs->us.n);
    }
    //matrix_complex_printf((cdouble **)zOld1, 6, reduced_nv);
    cout << "Delta with the reference = " << DENorm((const cdouble**)zOld0, (const cdouble**)zOld1, 6, reduced_nv) << endl;
    cout << endl;

    tic();
    invman.evalDCCMtoTFC(s0, mIn, OFTS_ORDER, OFS_ORDER);
    cout << "With new version, in " << toc() << endl;
    for(int i = 0; i < 6; i++)
    {
        for(int j = 0; j < reduced_nv; j++)  zNew[i][j] = mIn.getCoef(i,j).evaluate(SEML.cs->us.n);
    }
    //matrix_complex_printf((cdouble **)zNew, 6, reduced_nv);
    cout << "Delta with the reference = " << DENorm((const cdouble**)zOld0, (const cdouble**)zNew, 6, reduced_nv) << endl;
    cout << endl;
}

/**
 *  \brief Test of the routine Invman::evalRCMtoNC. RCM to NC coordinates
 **/
void test_evalRCMtoNC()
{
    //Reduced number of variables
    int reduced_nv = Invman::compRNV(*SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via Invman
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via initOFS
    //------------------------------------------------------------------------------------
    //--------------------------------------
    // Center-((Un)stable) manifolds
    //--------------------------------------
    vector<Oftsc> CM_NC;      //center manifold in NC coordinates
    vector<Oftsc> CM_TFC;     //center manifold in TFC coordinates
    //Init routine
    initOFTS(CM_NC, CM_TFC, reduced_nv, OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //--------------------------------------
    // COC objects
    //--------------------------------------
    matrix<Ofsc>  Mcoc(6,6);    //COC matrix
    matrix<Ofsc>  Pcoc(6,6);    //COC matrix (Mcoc = Pcoc*Complex matrix)
    matrix<Ofsc>  MIcoc(6,6);   //COC matrix = inv(Mcoc)
    matrix<Ofsc>  PIcoc(6,6);   //COC matrix = inv(Pcoc)
    vector<Ofsc>  Vcoc(6);      //COC vector

    //Read from files
    initCOC(Pcoc, Mcoc, PIcoc, MIcoc, Vcoc, SEML);

    //--------------------------------------
    //DCM_TFC : Jacobian of CM_TFC
    //--------------------------------------
    matrix<Oftsc> DCM_TFC(6, 4, reduced_nv, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_TFC, SEML.cs->F_PMS+"DWf/DWhc");

    //------------------------------------------------------------------------------------
    //Initial conditions (random)
    //------------------------------------------------------------------------------------
    double s0[reduced_nv];
    for(int i = 0; i < reduced_nv-1; i++) s0[i] = 0.01 * (rand() % 10);
    s0[reduced_nv-1] = 0.001 * (rand() % 10); //10% on the hyperbolic part

    //Initial time
    double t = 0.01 * (rand() % 10);

    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    Ofsc ofs(OFS_ORDER);
    double zOld1[6], zOld0[6], zNew[6];

    cout << "---------------------------------------------------" << endl;
    cout << "   Comparison of old and new version of RCMtoNC    " << endl;
    cout << "---------------------------------------------------" << endl;

    tic();
    RCMtoNCbyTFC(s0, t, SEML.us->n, OFTS_ORDER, OFS_ORDER, reduced_nv, CM_TFC, Mcoc, Vcoc, zOld0, 0);
    cout << "With old version, isGS = 0 (reference), in " << toc() << endl;
    //vector_printf_prec(zOld0, 6);
    cout << endl;

    tic();
    RCMtoNCbyTFC(s0, t, SEML.us->n, OFTS_ORDER, OFS_ORDER, reduced_nv, CM_TFC, Mcoc, Vcoc, zOld1, 1);
    cout << "With old version, isGS = 1, in " << toc() << endl;
    //vector_printf_prec(zOld1, 6);
    cout << "Delta with the reference = " << DENorm(zOld0, zOld1, 6) << endl;
    cout << endl;

    tic();
    invman.evalRCMtoNC(s0, t, zNew, OFTS_ORDER, OFS_ORDER);
    cout << "With new version, in " << toc() << endl;
    //vector_printf_prec(zNew, 6);
    cout << "Delta with the reference = " << DENorm(zOld0, zNew, 6) << endl;
    cout << endl;
}

/**
 *  \brief Test of the routine Invman::evalDRCMtoTFC_partial.
 **/
void test_evalDRCMtoTFC()
{
    //Reduced number of variables
    int reduced_nv = Invman::compRNV(*SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via Invman
    //------------------------------------------------------------------------------------
    Invman invman(OFTS_ORDER, OFS_ORDER, *SEML.cs);

    //------------------------------------------------------------------------------------
    //Initialization of the manifold via initOFS
    //------------------------------------------------------------------------------------
    //--------------------------------------
    //DCM_TFC : Jacobian of CM_TFC
    //--------------------------------------
    matrix<Oftsc> DCM_TFC(6, reduced_nv, reduced_nv, OFTS_ORDER-1, OFS_NV, OFS_ORDER);
    readMOFTS_bin(DCM_TFC, SEML.cs->F_PMS+"DWf/DWhc");

    //------------------------------------------------------------------------------------
    //Initial conditions (random)
    //------------------------------------------------------------------------------------
    double s0[reduced_nv];
    for(int i = 0; i < reduced_nv-1; i++) s0[i] = 0.01 * (rand() % 10);
    s0[reduced_nv-1] = 0.001 * (rand() % 10); //10% on the hyperbolic part

    //Initial time
    double t = 0.01 * (rand() % 10);

    //------------------------------------------------------------------------------------
    //Comparison
    //------------------------------------------------------------------------------------
    matrix<Ofsc> ofs(6, reduced_nv, OFS_ORDER);
    gsl_matrix_complex* zOld0 = gsl_matrix_complex_alloc(6, reduced_nv);
    gsl_matrix_complex* zOld1 = gsl_matrix_complex_alloc(6, reduced_nv);
    gsl_matrix_complex* zNew  = gsl_matrix_complex_alloc(6, reduced_nv);

    cout << "---------------------------------------------------" << endl;
    cout << "   Comparison of old and new version of DRCMtoNC   " << endl;
    cout << "---------------------------------------------------" << endl;

    tic();
    RCMtoTFC_JAC(s0, t, SEML.us->n, OFTS_ORDER, OFS_ORDER, reduced_nv, DCM_TFC, ofs, zOld0, 0);
    cout << "With old version, isGS = 0 (reference), in " << toc() << endl;
    //gslc_matrix_complex_printf(zOld0);
    cout << endl;

    tic();
    RCMtoTFC_JAC(s0, t, SEML.us->n, OFTS_ORDER, OFS_ORDER, reduced_nv, DCM_TFC, ofs, zOld1, 1);
    cout << "With old version, isGS = 1, in " << toc() << endl;
    //gslc_matrix_complex_printf(zOld1);
    cout << "Delta with the reference = " << gslc_matrix_complex_diff_L2(zOld0, zOld1) << endl;
    cout << endl;

    tic();
    invman.evalDRCMtoTFC_partial(s0, t, zNew, OFTS_ORDER, OFS_ORDER);
    cout << "With new version, in " << toc() << endl;
    //gslc_matrix_complex_printf(zNew);
    cout << "Delta with the reference = " << gslc_matrix_complex_diff_L2(zOld0, zNew) << endl;
    cout << endl;
}
