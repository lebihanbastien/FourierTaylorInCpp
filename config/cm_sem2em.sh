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
COMP_TYPE=$COMP_CMU_SEMLi_TO_CM_EMLj

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
ISPAR=0

#-----------------------------------------------------
# NOHUP condition
#-----------------------------------------------------
ISNOHUP=0

#-----------------------------------------------------
# I/O Handling
#-----------------------------------------------------
IO_HANDLING=$IO_DIALOG

#=====================================================
#  ---- Projection parameters ----
#=====================================================

#-----------------------------------------------------
# Parameters that change often
#-----------------------------------------------------
# Time grid: min, max and number of points on the grid
TMIN=0.2       # (given as %T, with T the SEM period)
TMAX=0.5       # (given as %T, with T the SEM period)
TSIZE=0	  

# Configuration (s1, s2, s3, s4) grid
GLIM_S1=(-0.2 0.2)
GLIM_S2=(0 0)
GLIM_S3=(-0.2 0.2)
GLIM_S4=(0 0)
GSIZE_SI=(+40 +0 +40 +0)

# Primary family - keep in mind that the first minimum rule is used now for the primary!
PRIMARY=0

# Fixed delta of energy (-1 if not used)
DHD=-1

# Hyperbolic components
HYP_EPSILON_EML2=1e-6	# Hyperbolic component at eml1,2
HYP_EPSILON_SEML2=1e-6	# Hyperbolic component at seml2

# Time frequency, in %T, just here to avoid a warning in main.sh, but not used in the sequel
DT=0.005

# Filenames (used only if IO_HANDLING==$IO_BASH)
FILE_CU="cu_orbit.bin"
FILE_PCU="projcu_orbit.bin"

#-----------------------------------------------------
# Parameters that are stable
#-----------------------------------------------------
TM=15 	      # Maximum integration time (given as %T, with T the SEM period)
MSIZE=500     # Number of points on each trajectory
NSMIN=20      # Number of sorted solutions
YNMAX=0.6     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of  EMLi
SNMAX=56.0    # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of  EMLi
NOD=6         # Number of dimensions on which we compute the norm of the projection error
