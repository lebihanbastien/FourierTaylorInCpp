# Constants for FourierTaylorInCpp
#
# The values stored in this file 
# should match the constants defined
# in the header file via #define
#=====================================================

#=====================================================
# Find the number of cores to detect 
# if we are on the server
#=====================================================
NCORES=`nproc --all`

SERVER=0
if [ $NCORES -gt 10 ]; then
	echo "Use of the server is detected"
	SERVER=1
fi

#=====================================================
#  ---- General parameters ----
#
#	COMP_TYPE
#	MODEL_TYPE
#
#=====================================================

#-----------------------------------------------------
# TYPE OF COMPUTATION (COMP_TYPE)
#-----------------------------------------------------
export COMP_CM_EML2_TO_CM_SEML=0    	     #EML2 Center Manifold to SEMLi Center Manifold
export COMP_CM_EML2_TO_CMS_SEML=1   	     #EML2 Center Manifold to SEMLi Center-Stable Manifold
export COMP_SINGLE_ORBIT=2    		         #Just some example of solutions
export COMP_CM_EML2_TO_CMS_SEML_READ=3       #Read
export COMP_VF_TEST=4    		             #Test of the QBCP vector field. Should probably be put in OOFTDA
export COMP_CM_EML2_TO_CM_SEML_REFINE=5      #EML2 Center Manifold to SEMLi Center Manifold
export COMP_EPHEMERIDES_TEST=6    	         #Ephemerides test
export COMP_CM_EML2_TO_CM_SEML_3D=7 	     #EML2 Center Manifold to SEMLi Center Manifold in 3D
export COMP_VOFTS_TO_VOFTS=8     	         #Store CS/CU into one-dimensionnal series to gain memory
export COMP_test_INVMAN=9    		         #Test of the new invariant manifold implementation
export COMP_REF_JPL=10  		             #Refine to JPL ephemerides
export COMP_CM_EML2_TO_CM_SEML_H=11          #Planar EMLi Center Manifold to SEMLi Center Manifold, at a given energy
export COMP_ORBIT_EML2_TO_CM_SEML=12         #EMLi Center Manifold to SEMLi Center Manifold, on a given set of orbits
export COMP_SINGLE_ORBIT_EML2_TO_CM_SEML=13  #EMLi Center Manifold to SEMLi Center Manifold, on a given orbit
export COMP_CMU_SEMLi_TO_CM_EMLj=14          #SEMLi Center-Unstable Manifold to EMLj Center Manifold

#-----------------------------------------------------
# MODEL
# M_RTBP = 0; M_QBCP = M_1; M_BCP = 2; M_ERTBP = 3
#-----------------------------------------------------
export M_RTBP=0
export M_QBCP=1
export M_BCP=2  
export M_ERTBP=3  

#-----------------------------------------------------
# Hyperbolic component (default values)
#-----------------------------------------------------
export HYP_EPSILON_EML2_DEFAULT=1e-5
export HYP_EPSILON_SEML2_DEFAULT=1e-6


#=====================================================
#  ---- Refinement parameters ----
#
#	REFST_TYPE
#	REFST_DIM
#	REFST_TIME
#	REFST_GRID
#	REFST_TERMINATION
#	REFST_COORD_TYPE
#	REFST_DJPL_COORD
#
#=====================================================

#-----------------------------------------------------
# REFST_TYPE
#-----------------------------------------------------
# Planar or 3D solutions
export REF_PLANAR=0
export REF_3D=1
export REF_MIXED=101

# Only one connection refined
export REF_SINGLE=2

# Only some orbit connections refined
export REF_ORBIT=21

# Continuation procedures
export REF_CONT=30		#simple continuation, defined with REF_FIXED_TIME/REF_VAR_TN/REF_VAR_TIME and REF_FIXED_GRID/REF_VAR_GRID/REF_GIVEN_GRID
export REF_CONT_D=31		#double continuation: (1) first a continuation with variable final time, (2) second continuation with fixed time 
export REF_CONT_D_HARD_CASE=32  #same as REF_CONT_D but with additionnal refinement procedures between the two continuations
export REF_CONT_ORBIT=33		# only some orbit connections refined, but a continuation procedure is performed on the final time

# Refinement of a complete trajectory: orbit at EML2 + transfer + orbit at SEMLi
export REF_COMP=4

# Type of constraints on the times in REF_CONT
export REF_FIXED_TIME=5		#all times are fixed
export REF_VAR_TN=6		#the final time (at SEMLi) is free
export REF_VAR_TIME=7		#all times except the first one are free

# Type of grid in REF_CONT
export REF_FIXED_GRID=8		#fixed time grid
export REF_VAR_GRID=9		#variable time grid, given by the numerical integration
export REF_GIVEN_GRID=10	#the time grid has been computed in a previous continuation

# Termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
export REF_COND_S5=11    	#the continuation is stopped when s5 = 0 at SEMLi
export REF_COND_T=12	    #the continuation is stopped when "enough turns are performed at SEMLi

#-----------------------------------------------------
# COORDINATES SYSTEMS
#-----------------------------------------------------
# Synodical QBCP
export NCSEM=0
export NCEM=1
export VNCSEM=2
export VNCEM=3
export PSEM=4
export PEM=5
export VSEM=6
export VEM=7
# Inertial QBCP
export INEM=8
export INSEM=9
export ECISEM=10
# Ephemerides
export VECLI=11
export J2000=12
export NJ2000=13
export VSYNEM=14
export VSYNSEM=15

#-----------------------------------------------------
# For selection of the times
#-----------------------------------------------------
export TIME_SELECTION_ABSOLUTE=1
export TIME_SELECTION_RATIO=0

#-----------------------------------------------------
# I/O handling
#-----------------------------------------------------
export IO_DEFAULT=0
export IO_BASH=1
export IO_DIALOG=2

#-----------------------------------------------------
# Type of computation for each orbit
#-----------------------------------------------------
export INT_NO_COMP=0
export INT_PROJ_CHECK=1
export INT_RED_COORD=2
export INT_PROJ_FREE=3
export INT_TRY_BOTH=4
