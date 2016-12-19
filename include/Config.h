#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
using namespace std;

class Config
{
    private:

        //--------------------------------------------------------------------------------
        //Constants for ODE structure (see ode.h & cpp)
        //--------------------------------------------------------------------------------
        double PREC_ABS;
        double PREC_REL;
        double PREC_ROOT;
        double PREC_DIFF;
        double PREC_HSTART;

        //--------------------------------------------------------------------------------
        //Constants for differential correction procedures
        //--------------------------------------------------------------------------------
        int DC_ITERMAX;

    public:
        static Config& configManager()
        {
            static Config instance;
            return instance;
        }

        //--------------------------------------------------------------------------------
        //Getters
        //--------------------------------------------------------------------------------
        double G_PREC_ABS(){return PREC_ABS;}
        double G_PREC_REL(){return PREC_REL;}
        double G_PREC_ROOT(){return PREC_ROOT;}
        double G_PREC_DIFF(){return PREC_DIFF;}
        double G_PREC_HSTART(){return PREC_HSTART;}

        int    G_DC_ITERMAX(){return DC_ITERMAX;}

        //--------------------------------------------------------------------------------
        //Setters
        //--------------------------------------------------------------------------------
        /**
         *  \brief Hard precision for integration purposes
         **/
        void C_PREC_HARD();

        /**
         *  \brief Soft precision for integration purposes, to fasten correction loops
         *         It is the default case.
         **/
        void C_PREC_SOFT();

    private:
        Config();
        Config(Config const&);
        ~Config();
        void operator=(Config const&);
};

#endif // CONFIG_H
