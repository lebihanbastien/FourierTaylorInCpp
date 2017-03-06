# Configuration file for FourierTaylorInCpp
#
# BLB 2017

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
COMP_TYPE=$COMP_CM_EML2_TO_CM_SEML

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
	OFTS_ORDER=16
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
ISPAR=1		

#-----------------------------------------------------
# NOHUP condition
#-----------------------------------------------------
ISNOHUP=0

#=====================================================
#  ---- Projection parameters ----
#=====================================================

#-----------------------------------------------------
# Parameters that change often
#-----------------------------------------------------
# Time grid: min, max and number of points on the grid
TMIN=0.995     # (given as %T, with T the SEM period)
TMAX=0.25      # (given as %T, with T the SEM period)
TSIZE=0	  

# Configuration (s1, s2, s3, s4) grid
GLIM_S1=(-35 +35)
GLIM_S2=(-12 +12)
GLIM_S3=(-35 +35)
GLIM_S4=(-12 +12)
GSIZE_SI=(+100 +6 +100 +6)

#-----------------------------------------------------
# Parameters that are stable
#-----------------------------------------------------
TM=12 	      # Maximum integration time (given as %T, with T the SEM period)
MSIZE=500     # Number of points on each trajectory
NSMIN=20      # Number of sorted solutions
YNMAX=0.6     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
SNMAX=0.6     # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
NOD=6         # Number of dimensions on which we compute the norm of the projection error

#====================================================================================
# Computation parameters: Projection CMU EML2 to CM SEMLi
#
#  Used for producing the complete sets:
#
#    #TIME
#    double TM    = 12*SEML.us.T;
#    double TMIN  = 0.75*SEML.us.T;
#    double TMAX  = 1*SEML.us.T;
#    int    TSIZE = 50;
#
#    #UNSTABLE MANIFOLD AT EML2
#    double GMIN_S1  = -35;
#    double GMAX_S1  = +35;
#    int    GSIZE_S1 = 200;
#
#    double GMIN_S3  = -35;
#    double GMAX_S3  = +35;
#    int    GSIZE_S3 = 200;
#
#    #TRAJECTORTY REFINMENT
#    int MSIZE = 500;    #number of points on each trajectory
#    int NSMIN = 20;      #number of sorted solutions
#
#====================================================================================
