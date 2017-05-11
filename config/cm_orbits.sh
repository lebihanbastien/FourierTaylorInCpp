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
COMP_TYPE=$COMP_ORBIT_EML2_TO_CM_SEML

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

#-----------------------------------------------------
# I/O Handling
#-----------------------------------------------------
IO_HANDLING=$IO_BASH

#=====================================================
#  ---- Projection parameters ----
#=====================================================

#-----------------------------------------------------
# Parameters that change often
#-----------------------------------------------------
# Time grid: min, max and number of points on the grid
TMIN=0.0    # (given as %T, with T the SEM period)
TMAX=1.0    # (given as %T, with T the SEM period)
TSIZE=200	  

# Configuration (s1, s2, s3, s4) grid
GLIM_S1=(+10 +40)
GLIM_S2=(+0 +0)
GLIM_S3=(-10 +35)
GLIM_S4=(+0 +0)
GSIZE_SI=(+0 +0 +10 +0)

# Primary family - keep in mind that the first minimum rule is used now for the primary!
PRIMARY=1

# Fixed delta of energy (-1 if not used)
DHD=-1

# Time frequency, in %T
DT=0.001
	
# Hyperbolic components
HYP_EPSILON_EML2=1e-5	# Hyperbolic component at eml2
HYP_EPSILON_SEML2=1e-6	# Hyperbolic component at seml2

# Filenames (used only if IO_HANDLING==$IO_BASH)
FILE_CU="cu_orbit.bin"
FILE_PCU="projcu_order_20_dest_L2_Orbit_10_eps_1e-5.bin"


#-----------------------------------------------------
# Parameters that are stable
#-----------------------------------------------------
TM=12 	      # Maximum integration time (given as %T, with T the SEM period)
MSIZE=500     # Number of points on each trajectory
NSMIN=20      # Number of sorted solutions
YNMAX=0.6     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
SNMAX=0.6     # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
NOD=6         # Number of dimensions on which we compute the norm of the projection error
