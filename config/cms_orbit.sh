# Configuration file for FourierTaylorInCpp

#-----------------------------------------------------
# Source the parameters
#-----------------------------------------------------
source config/constants.sh

#=====================================================
#  ---- General parameters ----
#=====================================================

#-----------------------------------------------------
# TYPE OF COMPUTATION (COMP_TYPE)
#-----------------------------------------------------
COMP_TYPE=$COMP_CM_EML2_TO_CMS_SEML

#-----------------------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------
MODEL=$M_QBCP 

#-----------------------------------------------------
# DEFAULT LIBRATION POINT FOR EM & SEM SYSTEM
#-----------------------------------------------------
LI_EM=2
LI_SEM=2

#-----------------------------------------------------
# Orders for the semi-numerical expansions
#-----------------------------------------------------
# Order of the Fourier series
OFS_ORDER=30

# Order of the Taylor series, 
# depending on where the computation takes place
if [ $SERVER == 1 ]; then
	OFTS_ORDER=20
else
	OFTS_ORDER=20
fi


#-----------------------------------------------------
# Parameters for parallel computation
#-----------------------------------------------------
# Number of parallel threads
if [ $SERVER == 1 ]; then
	NUM_THREADS=50
else
	NUM_THREADS=4
fi

# Parallel computation is on by default
ISPAR=0		

#-----------------------------------------------------
# NOHUP condition
#-----------------------------------------------------
ISNOHUP=0

#-----------------------------------------------------
# Position of the poincaré section in NCEM coordinates
# if -1, 3BSOI: about 159198km of radius around the Moon, 
# hence 2.46761593 in NCEM coordinates about EML2
#-----------------------------------------------------
RPS=-1 

#-----------------------------------------------------
# I/O Handling
#-----------------------------------------------------
IO_HANDLING=$IO_BASH

#=====================================================
#  ---- Refinement parameters ----
#=====================================================

#-----------------------------------------------------
# Parameters that change often
#-----------------------------------------------------
REFST_TYPE=$REF_ORBIT       # Type of refinement - either REF_ORBIT or REF_CONT_ORBIT
REFST_DIM=$REF_PLANAR       # Type of dimensions planar or 3d?
REFST_T0_DES=0.0            # Initial time - given as %T, with T the SEM period   

# Domain of search (min s1, max s1, min s2, max s2, min s3, max s3, min s4, max s4) for the first guess
REFST_SI_CMU_EM_LIM=(-1e5 1e5 -1e5 +1e5 -1e5 +1e5 -1e5 +1e5)
# Or, if we want the user to define such domain:
REFST_ISLIMUD=0

# Domain of search (min s1, max s1, min s3, max s3) for the seed of first guess
REFST_SI_SEED_EM_LIM=(10 10 -40 +40 -40 +40 -40 +40)

# Maximum projection distance allowed during subselection
REFST_PMAX_DIST_SEM=6e-4

# Limits for the time of flight during transfers - not used if -1
REFST_TOF_LIM=(-1 -1)

# Values for crossings
REFST_CROSSINGS=-1

# Number of steps in the continuation procedure
REFST_CONT_STEP_MAX=+450;       # with fixed times
REFST_CONT_STEP_MAX_VT=+50;     # with variable times

# Initial step in the continuation procedure
if [ $LI_EM == 1 ]; then
	REFST_FIXED_TIME_DS0=3e-2 # for fixed times
else
	REFST_FIXED_TIME_DS0=5e-1 # for fixed times
fi
REFST_VAR_TIME_DS0=8e-2	          # for variable times


# Desired number of iterations in Newton's method in the continuation procedure
REFST_FIXED_TIME_NU0=2           # for fixed times
REFST_VAR_TIME_NU0=4 	         # for variable times

# Direction of the continuation procedure
REFST_ISDIRUD=0			 # is it user defined?
REFST_DIR=-1    		 # if not, +1 or -1
 
# User parameters
REFST_ISFLAGON=0   	         # do we have steps in the procedure - asking the user to press enter to go on?
REFST_ISPLOTTED=0   		 # do we plot the results during the computation?
REFST_ISSAVED=1     		 # do we save the results in data files?
REFST_ISFROMSERVER=0		 # does the raw data comes from server files?

# Maximum angle around SEMLi if REF_COND_T is used (in degrees)
REFST_THETAMAX=180           # should be a multiple of 90°

# Filenames (used only if IO_HANDLING==$IO_BASH)
FILE_PCU="projcu_order_20_dest_L2_Orbit_10_eps_1e-5.bin" #"Serv/projcu_order_20_dest_L2_orbit.bin" #"Serv/projcu_order_20_dest_L2_QHalo.bin" 

FILE_CONT="cont_atf_order_20_dest_L2_TEST.txt"
FILE_CONT_RES="cont_atf_traj_order_20_dest_L2_TEST.bin"

FILE_TRAJ_FROM_W="traj_from_w_order_20_dest_L2_TEST.bin"
FILE_TRAJ_FROM_C="traj_from_c_order_20_dest_L2_TEST.bin" 

FILE_JPL="/cont_jpl_order_20_dest_L2_TEST_local.bin"

#-----------------------------------------------------
# Parameters that are stable
#-----------------------------------------------------
REFST_ISDEBUG=0			         # if yes, additionnal tests are performed
REFST_GRIDSIZE=20        	     # number of points on the refinement grid
REFST_MPLOT=200        	         # number of points per plot between to pach points (e.g. total plot points is REFST_MPLOT*REFST_GRIDSIZE)

REFST_TIME=$REF_VAR_TN		 # type of constraints on the times in REF_CONT
REFST_GRID=$REF_FIXED_GRID	 # type of grid
REFST_TERMINATION=$REF_COND_T   # Termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
REFST_COORD_TYPE=$NCSEM		 # coordinates system in the refinement procedure

REFST_XPS=0.65			  # position of the poincaré section in NCSEM coordinates
REFST_XPS_NCEM=-1         # position of the poincaré section in NCEM coordinates 3BSOI: about 159198km of radius around the Moon, hence 2.46761593 in NCEM coordinates about EML2
REFST_ISJPL=1		      # is the JPL refinement performed when possible?
REFST_DJPLCOORD=$NJ2000   # coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with $NJ2000
REFST_SIDIM=0		      # 0 or 2 - component of s0 that stays constant when t0 is free

# Sampling frequencies in REF_COMP (complete trajectory) in days
REFST_SF_EML2=2			 # orbit at EML2	
REFST_SF_MAN=5			 # transfer leg
REFST_SF_SEML2=10		 # orbit at SEML2

# Integration window for each orbit
REFST_TSPAN_EM=10    	 # given as %T, where T is the SEM period, in EM units
REFST_TSPAN_SEM=10		 # given as %T, where T is the SEM period, in SEM units

# Type of computation for each orbit
REFST_COMP_ORB_EM=$INT_TRY_BOTH  #may use INT_RED_COORD directly for bigger orbits (QHalo?)
REFST_COMP_ORB_SEM=$INT_TRY_BOTH
	
# Storing the orbits at each step?
REFST_ISSAVED_EM=0               # 0: don't save, 1: save using projection method
REFST_ISSAVED_SEM=0              # 0: don't save, 1: save using projection method, 2: save using integration in reduced coordinates

# Type of time selection
REFST_TYPE_OF_T_SEL=$TIME_SELECTION_RATIO
