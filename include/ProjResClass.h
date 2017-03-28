#ifndef PROJRESCLASS_H
#define PROJRESCLASS_H


//========================================================================================
//
//          Subroutines
//
//========================================================================================

//========================================================================================
//
//          Class
//
//========================================================================================
/**
 *  \struct ProjResClass
 *  \brief  Define a given structure to store the results from a projection
 **/
class ProjResClass
{
    //------------------------------------------------------------------------------------
    // Parameters
    //------------------------------------------------------------------------------------
private:

    //Data
    vector<double> t0_CMU_EM;
    vector<double> tf_CMU_EM;

    vector<double> s1_CMU_EM;
    vector<double> s2_CMU_EM;
    vector<double> s3_CMU_EM;
    vector<double> s4_CMU_EM;
    vector<double> s5_CMU_EM;

    vector<double> pmin_dist_SEM;
    vector<double> tf_man_SEM;

    vector<double> s1_CM_SEM;
    vector<double> s2_CM_SEM;
    vector<double> s3_CM_SEM;
    vector<double> s4_CM_SEM;

    vector<double> crossings_NCSEM;
    vector<double> collision_NCEM;

    //Sort vector
    vector<size_t> sortId;

    //Size
    int csize;

    //------------------------------------------------------------------------------------
    //Public routine. Rk: do we need a destructor???
    //------------------------------------------------------------------------------------
public:
    ProjResClass():csize(0) {}


    void push_back(ProjResClass& inSt, int kpor)
    {
        this->t0_CMU_EM.push_back(inSt.t0_CMU_EM[kpor]);
        this->tf_CMU_EM.push_back(inSt.tf_CMU_EM[kpor]);

        this->s1_CMU_EM.push_back(inSt.s1_CMU_EM[kpor]);
        this->s2_CMU_EM.push_back(inSt.s2_CMU_EM[kpor]);
        this->s3_CMU_EM.push_back(inSt.s3_CMU_EM[kpor]);
        this->s4_CMU_EM.push_back(inSt.s4_CMU_EM[kpor]);
        this->s5_CMU_EM.push_back(inSt.s5_CMU_EM[kpor]);

        this->pmin_dist_SEM.push_back(inSt.pmin_dist_SEM[kpor]);
        this->tf_man_SEM.push_back(inSt.tf_man_SEM[kpor]);


        this->s1_CM_SEM.push_back(inSt.s1_CM_SEM[kpor]);
        this->s2_CM_SEM.push_back(inSt.s2_CM_SEM[kpor]);
        this->s3_CM_SEM.push_back(inSt.s3_CM_SEM[kpor]);
        this->s4_CM_SEM.push_back(inSt.s4_CM_SEM[kpor]);


        this->crossings_NCSEM.push_back(inSt.crossings_NCSEM[kpor]);
        this->collision_NCEM.push_back(inSt.collision_NCEM[kpor]);

        this->csize =  this->t0_CMU_EM.size();
    }

    void pop_back()
    {
        t0_CMU_EM.pop_back();
        tf_CMU_EM.pop_back();
        s1_CMU_EM.pop_back();
        s2_CMU_EM.pop_back();
        s3_CMU_EM.pop_back();
        s4_CMU_EM.pop_back();
        s5_CMU_EM.pop_back();
        pmin_dist_SEM.pop_back();
        tf_man_SEM.pop_back();
        s1_CM_SEM.pop_back();
        s2_CM_SEM.pop_back();
        s3_CM_SEM.pop_back();
        s4_CM_SEM.pop_back();
        crossings_NCSEM.pop_back();
        collision_NCEM.pop_back();
        csize = t0_CMU_EM.size();
    }


    int size()
    {
        return csize;
    }

    vector<size_t>& getSortId()
    {
        return sortId;
    }

    void displayEntry(int k)
    {
        int ind = sortId[k];
        cout << "t0_EM    = "  << t0_CMU_EM[ind]  << endl;
        cout << "tf_EM    = "  << tf_CMU_EM[ind] << endl;
        cout << "pmin_SEM = "  << pmin_dist_SEM[ind] << endl;
        cout << "s_CM_EM  = (" << s1_CMU_EM[ind] << ", " << s2_CMU_EM[ind] << ", " << s3_CMU_EM[ind] << ", " << s4_CMU_EM[ind] << ")" << endl;
        cout << "s_CM_SEM = (" << s1_CM_SEM[ind] << ", " << s2_CM_SEM[ind] << ", " << s3_CM_SEM[ind] << ", " << s4_CM_SEM[ind] << ")" << endl;
    }

    void displayFirstEntry()
    {
        displayEntry(0);
    }
    void displayLasttEntry()
    {
        displayEntry(csize-1);
    }

    /**
     *  \brief Read in a data file the connections between EML1,2 and SEML1,2.
     **/
    int readIntProjCU_bin(string filename)
    {
        string fname = "readIntProjCU_bin";

        //====================================================================================
        //Check the existence of the filename
        //====================================================================================
        if(!fileExists(filename))
        {
            cout << fname << ". " << filename << ": " << strerror(errno) << endl;
            return FTC_ENOENT;
        }

        //====================================================================================
        //Get the number of columns
        //====================================================================================
        int ncol[5] = {61, 56, 55, 39, 37};
        int ncol0   = numberOfColumns(filename, ncol, 5);

        //cout << "Number of columns in " << filename << " is: " << ncol0 << endl;
        //pressEnter(true);

        //====================================================================================
        //Open and read datafile (the existence) has already been checked)
        //====================================================================================
        fstream filestream;
        filestream.open (filename.c_str(), ios::binary | ios::in);
        if (filestream.is_open())
        {
            double res;

            //--------------------------------------------------------------------------------
            // Read data
            //--------------------------------------------------------------------------------
            do
            {
                //0. Label if ncol >= 56
                if(ncol0 >= 56) filestream.read((char*) &res, sizeof(double));

                // 1. time grid in NCEM units
                filestream.read((char*) &res, sizeof(double));
                this->t0_CMU_EM.push_back(res);

                // 2-7. initial state in NCEM coordinates: NOT SAVED
                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));

                }

                // 8-13. initial state in SEM coordinates: NOT SAVED
                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                }

                // 14-18. initial state in RCM coordinates: SAVED
                filestream.read((char*) &res, sizeof(double));
                this->s1_CMU_EM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s2_CMU_EM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s3_CMU_EM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s4_CMU_EM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s5_CMU_EM.push_back(res);

                // 19. minimum distance of projection: SAVED
                filestream.read((char*) &res, sizeof(double));
                this->pmin_dist_SEM.push_back(res);

                // 20. associated dv: NOT SAVED
                filestream.read((char*) &res, sizeof(double));

                // 21. tf at SEM: SAVED
                filestream.read((char*) &res, sizeof(double));
                this->tf_man_SEM.push_back(res);

                //Moreover, we save the value in tf_CMU_EM
                this->tf_CMU_EM.push_back(res/SEML.us_em.ns);


                // 22-27. final_state_CMU_SEM state in SE coordinates: NOT SAVED
                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                }

                // 28-33. projected_state_CMU_SEM state in SE coordinates: NOT SAVED
                for (int k = 0; k < 6; k++)
                {
                    filestream.read((char*) &res, sizeof(double));
                }

                // 34-37. projected_state_CMU_RCM state in SE coordinates: SAVED
                filestream.read((char*) &res, sizeof(double));
                this->s1_CM_SEM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s2_CM_SEM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s3_CM_SEM.push_back(res);

                filestream.read((char*) &res, sizeof(double));
                this->s4_CM_SEM.push_back(res);

                //----------------------------------------------------------------------------
                // Last columns depend on the number ncol0
                //----------------------------------------------------------------------------
                switch(ncol0)
                {
                case 37:
                {
                    //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                    res = -1.0;
                    this->crossings_NCSEM.push_back(res);

                    //39. Collision flag, from NCEM flow
                    res = -1.0;
                    this->collision_NCEM.push_back(res);

                    break;
                }
                case 39:
                {
                    //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                    filestream.read((char*) &res, sizeof(double));
                    this->crossings_NCSEM.push_back(res);

                    //39. Collision flag, from NCEM flow
                    filestream.read((char*) &res, sizeof(double));
                    this->collision_NCEM.push_back(res);

                    break;
                }
                case 55:
                case 56:
                {
                    //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                    filestream.read((char*) &res, sizeof(double));
                    this->crossings_NCSEM.push_back(res);

                    //39. Collision flag, from NCEM flow
                    filestream.read((char*) &res, sizeof(double));
                    this->collision_NCEM.push_back(res);

                    //40-43: H0 at IC - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //44-47: H0 at emli - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //48-51: Hf - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //52-55: Hf at semli - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    break;
                }
                case 61:
                {
                    //38. Number of crossings of the x = -1 line (clock/counterclockwise)
                    filestream.read((char*) &res, sizeof(double));
                    this->crossings_NCSEM.push_back(res);

                    //39. Collision flag, from NCEM flow
                    filestream.read((char*) &res, sizeof(double));
                    this->collision_NCEM.push_back(res);

                    //40-43: H0 at IC - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //44-47: H0 at emli - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //48-51: Hf - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //52-55: Hf at semli - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));
                    filestream.read((char*) &res, sizeof(double));

                    //57. init_time_EM_seed - NOT SAVED FOR NOW
                    filestream.read((char*) &res, sizeof(double));

                    //58-61. initial seed in RCM coordinates - NOT SAVED FOR NOW
                    for (int k = 0; k < 4; k++) filestream.read((char*) &res, sizeof(double));

                    break;
                }
                }
            }
            while(!filestream.eof());

            filestream.close();
        }
        else
        {

            cout << "readClosestIntProjCU_bin. Unable to open the file. Check file name: " << endl;
            cout << filename << endl;
            return FTC_FAILURE;
        }


        //====================================================================================
        //Delete last element that is not a real value
        //====================================================================================
        this->pop_back();

        //====================================================================================
        //Update the size
        //====================================================================================
        this->csize = this->t0_CMU_EM.size();

        return FTC_SUCCESS;
    }


    /**
    *  \brief Read in a data file the connections between EML2 and SEML1,2.
    *         Interpolate in the data set to get the right desired t0 at EML2 departures.
    **/
    int readClosestIntProjCU_bin(string filename, double t0_des, int typeOfTimeSelection)
    {
        string fname = "readClosestIntProjCU_bin";

        //====================================================================================
        //Check the existence of the filename
        //====================================================================================
        if(!fileExists(filename))
        {
            cout << fname << ". " << filename << ": " << strerror(errno) << endl;
            return FTC_ENOENT;
        }


        //====================================================================================
        //Open and read datafile
        //====================================================================================
        ProjResClass readSt;
        readSt.readIntProjCU_bin(filename);

        //====================================================================================
        //Type of selection
        //====================================================================================
        coutmp();
        cout << fname << ". ";

        if(typeOfTimeSelection == TIME_SELECTION_ABSOLUTE)
            cout << " the initial time is searched in [0, Inf[ " << endl;
        else
            cout << " the initial time is searched in [0, T[ " << endl;

        int ti = 0;
        std::vector<size_t> indRes;
        switch(typeOfTimeSelection)
        {
        case TIME_SELECTION_ABSOLUTE:
        {
            //----------------------------------------------------------------------------
            //Get the unique elements in t0_CMU_EM
            //----------------------------------------------------------------------------
            //Copy t0_CMU_EM into t0_CMU_EM_UNIQUE
            vector<double> t0_CMU_EM_UNIQUE(readSt.t0_CMU_EM);
            //Get unique elements
            vector_getUnique(t0_CMU_EM_UNIQUE);

            //----------------------------------------------------------------------------
            // Print the range in t0_CMU_EM_UNIQUE
            //----------------------------------------------------------------------------
            this->sortId = sort_indexes(t0_CMU_EM_UNIQUE);
            cout << "--------------------------------------" << endl;
            cout << "There is " << t0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
            cout << "[" << t0_CMU_EM_UNIQUE[this->sortId[0]]/SEML.us_em.T << ", ";
            cout << t0_CMU_EM_UNIQUE[this->sortId[t0_CMU_EM_UNIQUE.size()-1]]/SEML.us_em.T;
            cout << "] x T." << endl;

            //----------------------------------------------------------------------------
            // Find the nearest t0 value, or let the user choose if t0_des < 0
            //----------------------------------------------------------------------------
            cout << "Desired time is  = " << t0_des/SEML.us_em.T;
            cout << " x T." <<  endl;
            if(t0_des >= 0)
            {
                double dmin = fabs(t0_CMU_EM_UNIQUE[0] - t0_des);
                for(int i = 1; i < (int) t0_CMU_EM_UNIQUE.size(); i++)
                {
                    if(fabs(t0_CMU_EM_UNIQUE[i] - t0_des) < dmin)
                    {
                        dmin = fabs(t0_CMU_EM_UNIQUE[i] - t0_des);
                        ti = i;
                    }
                }
                cout << "The index of the time closer to the desired one is " << ti;
                cout << ", which corresponds to:" << endl;
            }
            else
            {
                do
                {
                    cout << "Please enter a number between " << 0 << " and " << t0_CMU_EM_UNIQUE.size()-1;
                    cout << " to select a specific starting time" << endl;
                    scanf("%d", &ti);
                }
                while(ti < 0 || ti > (int) t0_CMU_EM_UNIQUE.size()-1);
                cout << "You have selected " << ti << ", which corresponds to:" << endl;
            }

            cout << "t0 = " << t0_CMU_EM_UNIQUE[ti]/SEML.us_em.T;
            cout << " x T." <<  endl;

            //----------------------------------------------------------------------------
            // Select the values that matches the desired time
            //----------------------------------------------------------------------------
            vector_getIndices_range(indRes, readSt.t0_CMU_EM, t0_CMU_EM_UNIQUE[ti]);

            break;
        }

        case TIME_SELECTION_RATIO:
        default:
        {
            //----------------------------------------------------------------------------
            //We get all the unique times, as %T, between [0 and T]
            //----------------------------------------------------------------------------
            vector<double> r0_CMU_EM(readSt.t0_CMU_EM);

            //Get the modulo[T], up to RATIO_DIGITS digits
            double factor = 1.0*pow(10, RATIO_DIGITS);
            for(int k = 0; k < (int) r0_CMU_EM.size(); k++)
            {
                r0_CMU_EM[k] = round(factor*fmod(r0_CMU_EM[k]/SEML.us_em.T, 1.0))/factor;
            }

            //Get unique elements
            vector<double> r0_CMU_EM_UNIQUE(r0_CMU_EM);
            vector_getUnique(r0_CMU_EM_UNIQUE);

            //----------------------------------------------------------------------------
            // Print the range in r0_CMU_EM_UNIQUE
            //----------------------------------------------------------------------------
            this->sortId = sort_indexes(r0_CMU_EM_UNIQUE);
            cout << "--------------------------------------" << endl;
            cout << "There is " << r0_CMU_EM_UNIQUE.size() << " different times in data, in the following range:" << endl;
            cout << "[" << r0_CMU_EM_UNIQUE[this->sortId[0]] << ", ";
            cout << r0_CMU_EM_UNIQUE[this->sortId[r0_CMU_EM_UNIQUE.size()-1]];
            cout << "] x T." << endl;

            //----------------------------------------------------------------------------
            // Find the nearest t0 value, or let the user choose if t0_des < 0
            //----------------------------------------------------------------------------
            double r0_des = round(factor*fmod(t0_des/SEML.us_em.T, 1.0))/factor;

            cout << "Desired time is  = " << r0_des;
            cout << " x T." <<  endl;

            if(t0_des >= 0)
            {
                double dmin = fabs(r0_CMU_EM_UNIQUE[0] - r0_des);

                for(int i = 1; i < (int) r0_CMU_EM_UNIQUE.size(); i++)
                {
                    if(fabs(r0_CMU_EM_UNIQUE[i] - r0_des) < dmin)
                    {
                        dmin = fabs(r0_CMU_EM_UNIQUE[i] - r0_des);
                        ti = i;
                    }
                }
                cout << "The index of the time closer to the desired one is " << ti;
                cout << ", which corresponds to:" << endl;
            }
            else
            {
                do
                {
                    cout << "Please enter a number between " << 0 << " and " << r0_CMU_EM_UNIQUE.size()-1;
                    cout << " to select a specific starting time" << endl;
                    scanf("%d", &ti);
                }
                while(ti < 0 || ti > (int) r0_CMU_EM_UNIQUE.size()-1);
                cout << "You have selected " << ti << ", which corresponds to:" << endl;
            }

            cout << "t0 = " << r0_CMU_EM_UNIQUE[ti];
            cout << " x T." <<  endl;

            //----------------------------------------------------------------------------
            // Select the values that matches the desired time
            //----------------------------------------------------------------------------
            vector_getIndices_range(indRes, r0_CMU_EM, r0_CMU_EM_UNIQUE[ti]);
            break;
        }


        }
        cout << "--------------------------------------" << endl;
        coutlp();


        //====================================================================================
        // Copy the selected results in the inputs
        //====================================================================================
        for(int ind = 0; ind < (int) indRes.size(); ind++)
        {
            //Times
            this->t0_CMU_EM.push_back(readSt.t0_CMU_EM[indRes[ind]]);
            this->tf_CMU_EM.push_back(readSt.tf_CMU_EM[indRes[ind]]);
            //CMU of EM
            this->s1_CMU_EM.push_back(readSt.s1_CMU_EM[indRes[ind]]);
            this->s2_CMU_EM.push_back(readSt.s2_CMU_EM[indRes[ind]]);
            this->s3_CMU_EM.push_back(readSt.s3_CMU_EM[indRes[ind]]);
            this->s4_CMU_EM.push_back(readSt.s4_CMU_EM[indRes[ind]]);
            this->s5_CMU_EM.push_back(readSt.s5_CMU_EM[indRes[ind]]);
            //Projection distance
            this->pmin_dist_SEM.push_back(readSt.pmin_dist_SEM[indRes[ind]]);
            //tf_man_SEM
            this->tf_man_SEM.push_back(readSt.tf_man_SEM[indRes[ind]]);
            //CM of SEM
            this->s1_CM_SEM.push_back(readSt.s1_CM_SEM[indRes[ind]]);
            this->s2_CM_SEM.push_back(readSt.s2_CM_SEM[indRes[ind]]);
            this->s3_CM_SEM.push_back(readSt.s3_CM_SEM[indRes[ind]]);
            this->s4_CM_SEM.push_back(readSt.s4_CM_SEM[indRes[ind]]);
            //Crossings
            this->crossings_NCSEM.push_back(readSt.crossings_NCSEM[indRes[ind]]);
            this->collision_NCEM.push_back(readSt.collision_NCEM[indRes[ind]]);
        }

        //====================================================================================
        //Update the size
        //====================================================================================
        this->csize = this->t0_CMU_EM.size();

        //====================================================================================
        // Sort data wrt to the tof (this is done here using the tf, since t0 is constant)
        //====================================================================================
        this->sortId = sort_indexes(this->tf_CMU_EM);
        int ind;

        cout << "--------------------------------------" << endl;
        cout << "The min and max time of flights are:" << endl;
        ind = this->sortId[0];
        cout << "min(tof_EM)  = " << this->tf_CMU_EM[ind]/SEML.us->T - this->t0_CMU_EM[ind]/SEML.us->T  << " x T" << endl;
        ind = this->sortId[this->csize -1];
        cout << "max(tof_EM)  = " << this->tf_CMU_EM[ind]/SEML.us->T - this->t0_CMU_EM[ind]/SEML.us->T  << " x T" << endl;

        //====================================================================================
        // Sort data wrt to the projection distance
        //====================================================================================
        this->sortId = sort_indexes(this->pmin_dist_SEM);
        coutlp();

        return FTC_SUCCESS;
    }

    bool push_back_conditional(ProjResClass& inSt, RefSt & refSt)
    {
        //====================================================================================
    // 3. Select given intervals for the inputs
    //====================================================================================
    double s1_CMU_EM_MIN, s1_CMU_EM_MAX;
    double s3_CMU_EM_MIN, s3_CMU_EM_MAX;
    double s2_CMU_EM_MIN, s2_CMU_EM_MAX;
    double s4_CMU_EM_MIN, s4_CMU_EM_MAX;
    double tof_MIN, tof_MAX;

    if(refSt.isLimUD)
    {
        cout << "-------------------------------------------------------" << endl;
        cout << "   selectemlisemli. User-defined interval of research  " << endl;
        cout << "-------------------------------------------------------" << endl;
        cout << "Enter a value for s1_CMU_EM_MIN: ";
        cin >> s1_CMU_EM_MIN;
        cout << "Enter a value for s1_CMU_EM_MAX: ";
        cin >> s1_CMU_EM_MAX;
        cout << "Enter a value for s3_CMU_EM_MIN: ";
        cin >> s3_CMU_EM_MIN;
        cout << "Enter a value for s3_CMU_EM_MAX: ";
        cin >> s3_CMU_EM_MAX;

        cout << "Enter a value for tof_MIN (% of T, -1 if no use): ";
        cin >> tof_MIN;
        tof_MIN *= SEML.us->T;
        cout << "Enter a value for tof_MAX (% of T, -1 if no use): ";
        cin >> tof_MAX;
        tof_MAX *= SEML.us->T;
    }
    else
    {
        s1_CMU_EM_MIN = refSt.s1_CMU_EM_MIN;
        s1_CMU_EM_MAX = refSt.s1_CMU_EM_MAX;
        s3_CMU_EM_MIN = refSt.s3_CMU_EM_MIN;
        s3_CMU_EM_MAX = refSt.s3_CMU_EM_MAX;
        tof_MIN = refSt.tof_MIN;
        tof_MAX = refSt.tof_MAX;
    }

    s2_CMU_EM_MIN = refSt.s2_CMU_EM_MIN;
    s2_CMU_EM_MAX = refSt.s2_CMU_EM_MAX;
    s4_CMU_EM_MIN = refSt.s4_CMU_EM_MIN;
    s4_CMU_EM_MAX = refSt.s4_CMU_EM_MAX;

        int kpor  = 0;
        bool cst  = 0;
        bool flag = 0;

        for(int kpos = 0; kpos < (int) inSt.size(); kpos++)
        {
            kpor = inSt.sortId[kpos];

            //--------------------------------------------------------------------------------
            // Limit in the reduced coordinates (s1, s2, s3, s4)
            //--------------------------------------------------------------------------------
            cst  = (inSt.s1_CMU_EM[kpor] >= s1_CMU_EM_MIN) & (inSt.s1_CMU_EM[kpor] <= s1_CMU_EM_MAX);
            cst  = cst & (inSt.s3_CMU_EM[kpor] >= s3_CMU_EM_MIN) & (inSt.s3_CMU_EM[kpor] <= s3_CMU_EM_MAX);

            if(refSt.is3D())
            {
                cst  = cst & (inSt.s2_CMU_EM[kpor] >= s2_CMU_EM_MIN) & (inSt.s2_CMU_EM[kpor] <= s2_CMU_EM_MAX);
                cst  = cst & (inSt.s4_CMU_EM[kpor] >= s4_CMU_EM_MIN) & (inSt.s4_CMU_EM[kpor] <= s4_CMU_EM_MAX);
            }

            //--------------------------------------------------------------------------------
            // Limits in the time of flight (TOF) if necessary
            //--------------------------------------------------------------------------------
            if(tof_MIN > 0) cst = cst & ( (inSt.tf_CMU_EM[kpor] - inSt.t0_CMU_EM[kpor]) > tof_MIN );
            if(tof_MAX > 0) cst = cst & ( (inSt.tf_CMU_EM[kpor] - inSt.t0_CMU_EM[kpor]) < tof_MAX );

            //--------------------------------------------------------------------------------
            // Limits in the crossings if necessary
            //--------------------------------------------------------------------------------
            if(refSt.crossings > 0)
            {
                cst = cst & (inSt.crossings_NCSEM[kpor] == refSt.crossings);
            }

            if(cst)
            {
                flag = 1;
                this->push_back(inSt, kpor);
            }
        }

        //--------------------------------------------------------------------------------
        // Update the size
        //--------------------------------------------------------------------------------
        this->csize = this->t0_CMU_EM.size();

        return flag;
    }


    void sort_pmin_dist_SEM()
    {
        this->sortId = sort_indexes(this->pmin_dist_SEM);
    }

    void update_ic(double st_EM[5], double st_SEM[5], double t_EM[2],
                    double* t0_SEM, double* pmin_dist_SEM_out, int k)
    {
        int kpos = this->sortId[k];

        //================================================================================
        // 6. Initialize local variables: EM
        //================================================================================
        //RCM coordinates
        st_EM[0] = this->s1_CMU_EM[kpos];
        st_EM[1] = this->s2_CMU_EM[kpos];
        st_EM[2] = this->s3_CMU_EM[kpos];
        st_EM[3] = this->s4_CMU_EM[kpos];
        st_EM[4] = PROJ_EPSILON;
        //Time
        t_EM[0]  = this->t0_CMU_EM[kpos];
        t_EM[1]  = this->tf_CMU_EM[kpos];


        //================================================================================
        // 7. Initialize local variables: SEM
        //================================================================================
        //RCM coordinates
        st_SEM[0] = this->s1_CM_SEM[kpos];
        st_SEM[1] = this->s2_CM_SEM[kpos];
        st_SEM[2] = this->s3_CM_SEM[kpos];
        st_SEM[3] = this->s4_CM_SEM[kpos];
        st_SEM[4] = 0.0;
        //Initial time in SEM units
        *t0_SEM = this->tf_CMU_EM[kpos]*SEML.us_em.ns;

        //Minimum projection distance
        *pmin_dist_SEM_out = this->pmin_dist_SEM[kpos];


    }

};



#endif // PROJRESCLASS_H
