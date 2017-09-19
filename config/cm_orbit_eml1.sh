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
COMP_TYPE=$COMP_SINGLE_ORBIT_EML2_TO_CM_SEML

#-----------------------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------
MODEL=$M_QBCP 

#-----------------------------------------------------
# DEFAULT LIBRATION POINT FOR EM & SEM SYSTEM
#-----------------------------------------------------
LI_EM=1
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
ISPAR=1

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
#  ---- Projection parameters ----
#=====================================================

#-----------------------------------------------------
# Parameters that change often
#-----------------------------------------------------
# Time grid: min, max and number of points on the grid
TMIN=0.0    # (given as %T, with T the SEM period)
TMAX=1.0    # (given as %T, with T the SEM period)
TSIZE=100  

# Configuration (s1, s2, s3, s4) grid
GLIM_S1=(+1.5 +10)
GLIM_S2=(+0.6 +0)
GLIM_S3=(-1.5 -10)
GLIM_S4=(+0.6 +0)
GSIZE_SI=(+0 +0 +10 +0)

# Values for QHalo orbit (small one)
# GLIM_S1=(+28 +30)
# GLIM_S2=(+1.74347452709299 +0)
# GLIM_S3=(36 +35)
# GLIM_S4=(1.74347452709299 +0)
# GSIZE_SI=(+0 +0 +10 +0)

# Values for QHalo orbit (medium one)
# GLIM_S1=(+28 +30)
# GLIM_S2=(+2.81422620372967 +0)
# GLIM_S3=(36 +35)
# GLIM_S4=(2.81422620372967 +0)
# GSIZE_SI=(+0 +0 +10 +0)

# Values for QHalo orbit (big one)
# GLIM_S1=(+41.7272727272727 +30)
# GLIM_S2=(1.55435533128842 +0)
# GLIM_S3=(13.9090909090909 +35)
# GLIM_S4=(1.55435533128842 +0)
# GSIZE_SI=(+0 +0 +10 +0)

# Primary family
PRIMARY=0

# Fixed delta of energy (-1 if not used)
DHD=-1

# Time frequency, in %T
DT=0.001
	
# Hyperbolic components
HYP_EPSILON_EML2=-1e-5	# Hyperbolic component at eml2
HYP_EPSILON_SEML2=1e-6	# Hyperbolic component at seml2

# Filenames (used only if IO_HANDLING==$IO_BASH)
# Filenames (used only if IO_HANDLING==$IO_BASH)
if [ $SERVER == 1 ]; then
	FILE_CU="cu_orbit.bin"
	FILE_PCU="projcu_order_20_dest_L2_Orbit_Liss_s1_30_s2_5.bin"
else
	FILE_CU="cu_orbit.bin"
	FILE_PCU="Serv/projcu_order_20_dest_L2_Orbit_s1_15_s2_6_All_TM_8.bin"
fi

#-----------------------------------------------------
# Parameters that are stable
#-----------------------------------------------------
TM=8 	      # Maximum integration time (given as %T, with T the SEM period)
MSIZE=500     # Number of points on each trajectory
NSMIN=20      # Number of sorted solutions
YNMAX=0.6     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
SNMAX=0.6     # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
NOD=6         # Number of dimensions on which we compute the norm of the projection error
