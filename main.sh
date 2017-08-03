# main script file for FourierTaylorInCpp
# RECALL: make it executable with "$ chmod +x main.sh"

#--------------------------------------------------------------
# Source the configuration file passed as argument
#--------------------------------------------------------------
source $1

#--------------------------------------------------------------
# Get the current directory
#--------------------------------------------------------------
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DIR=$DIR"/"


#--------------------------------------------------------------
# Function to update unset parameters
#--------------------------------------------------------------
function set_param
{
	# argAry2 contains the arguments after $1, 
	# which may correspond to either an array or a scalar
	declare -a argAry2=("${!2}")
	local size=${#argAry2[@]}
    
	
	# If the size of argAry2 is equal to 1, we the second argument is a scalar.
	# if not, it is an array and we use a specific call
	echo "#------------------------------------------#"
	if [ "$size" == "1" ]; then
	
		#Evaluate
		eval "$1=$2"
		
		# Warn the user
		echo 'WARNING: the variable' $1 'is not set.'
		echo $2 'is used by default.'
	
	
	else
	
		#Evaluate
		eval "$1=( ${argAry2[@]} )"  
		
		# Warn the user
		echo 'WARNING: the variable' $1 'is not set.'
		echo '(' ${argAry2[@]} ') is used by default.'
	fi
	echo "#------------------------------------------#"
}

#--------------------------------------------------------------
# Title
#--------------------------------------------------------------
echo
echo "#------------------------------------------#"
echo "#     FourierTaylorInCpp configuration     #"
echo "#------------------------------------------#"
echo "Current directory: " "$DIR"
echo "Source file      : " "$1"
echo ''

#=====================================================
#  ---- General parameters ----
#=====================================================

#-----------------------------------------------------
# TYPE OF COMPUTATION (COMP_TYPE)
#-----------------------------------------------------
echo "#------------------------------------------#"
#Check that the variable $COMP_TYPE exist
if [ -z ${COMP_TYPE+x} ]; then
	echo 'WARNING: the variable COMP_TYPE is not set.'
	echo 'No computation.'
	return
else
	echo 'Type of computation:'
	case $COMP_TYPE in
		$COMP_CM_EML2_TO_CM_SEML_H)      echo 'COMP_TYPE    = COMP_CM_EML2_TO_CM_SEML_H'
		;;
		$COMP_CM_EML2_TO_CM_SEML)        echo 'COMP_TYPE    = COMP_CM_EML2_TO_CM_SEML'
		;;
		$COMP_CM_EML2_TO_CMS_SEML) 	     echo 'COMP_TYPE    = COMP_CM_EML2_TO_CMS_SEML'
		;;
		$COMP_SINGLE_ORBIT)              echo 'COMP_TYPE    = COMP_SINGLE_ORBIT'
		;;
		$COMP_CM_EML2_TO_CMS_SEML_READ)  echo 'COMP_TYPE    = COMP_CM_EML2_TO_CMS_SEML_READ'
		;;
		$COMP_VF_TEST)                   echo 'COMP_TYPE    = COMP_VF_TEST'
		;;
		$COMP_CM_EML2_TO_CM_SEML_REFINE) echo 'COMP_TYPE    = COMP_CM_EML2_TO_CM_SEML_REFINE'
		;;
		$COMP_EPHEMERIDES_TEST)          echo 'COMP_TYPE    = COMP_EPHEMERIDES_TEST'
		;;
		$COMP_CM_EML2_TO_CM_SEML_3D)     echo 'COMP_TYPE    = COMP_CM_EML2_TO_CM_SEML_3D'
		;;
		$COMP_VOFTS_TO_VOFTS)            echo 'COMP_TYPE    = COMP_VOFTS_TO_VOFTS'
		;;
		$COMP_test_INVMAN) 		         echo 'COMP_TYPE    = COMP_test_INVMAN'
		;;
		$COMP_REF_JPL) 			         echo 'COMP_TYPE    = COMP_REF_JPL'
		;;
		$COMP_ORBIT_EML2_TO_CM_SEML) 	 echo 'COMP_TYPE    = COMP_ORBIT_EML2_TO_CM_SEML'
		;;
		$COMP_SINGLE_ORBIT_EML2_TO_CM_SEML) echo 'COMP_TYPE    = COMP_SINGLE_ORBIT_EML2_TO_CM_SEML'
		;;
		$COMP_CMU_SEMLi_TO_CM_EMLj) 	 echo 'COMP_TYPE    = COMP_CMU_SEMLi_TO_CM_EMLj'
		;;
		*)     				             echo "COMP_TYPE    = "$COMP_TYPE". Unknown type."
	esac
fi



#=====================================================
#  ---- General parameters ----
#=====================================================
echo "#------------------------------------------#"
echo ''
echo "#------------------------------------------#"
echo 'Current set of general parameters:'
echo ''

#-----------------------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------
#Check that the variable $MODEL exist
if [ -z ${MODEL+x} ]; then
	set_param "MODEL" $M_QBCP
else
	if [ "$MODEL" == "$RTBP" ]; then
		OFS_ORDER=0
	else  
		OFS_ORDER=30
	fi
fi

# display MODEL
case $MODEL in
	$M_RTBP)  echo 'MODEL                  = M_RTBP'
	;;
	$M_QBCP)  echo 'MODEL                  = M_QBCP'
	;;
	$M_BCP)   echo 'MODEL                  = M_BCP'
	;;
	$M_ERTBP) echo 'MODEL                  = M_ERTBP'
	;;
	*)        echo "MODEL                  = "$MODEL". Unknown type."
esac

#-----------------------------------------------------
# DEFAULT LIBRATION POINT FOR EM & SEM SYSTEM
#-----------------------------------------------------
echo "LI_EM  		       =" $LI_EM
echo "LI_SEM 		       =" $LI_SEM

#-----------------------------------------------------
# Orders for the semi-numerical expansions
#-----------------------------------------------------
echo "OFS_ORDER  	       =" $OFS_ORDER
echo "OFTS_ORDER 	       =" $OFTS_ORDER

#-----------------------------------------------------
# Parameters for parallel computation
#-----------------------------------------------------
echo "ISPAR                  =" $ISPAR
echo "NUM_THREADS            =" $NUM_THREADS


#-----------------------------------------------------
# Hyperbolic component
#-----------------------------------------------------
# At EML2
if [ -z ${HYP_EPSILON_EML2+x} ]; then
		set_param "HYP_EPSILON_EML2" $HYP_EPSILON_EML2_DEFAULT
fi
echo "HYP_EPSILON_EML2       =" $HYP_EPSILON_EML2

# At SEML2
if [ -z ${HYP_EPSILON_SEML2+x} ]; then
		set_param "HYP_EPSILON_SEML2" $HYP_EPSILON_SEML2_DEFAULT
fi
echo "HYP_EPSILON_SEML2      =" $HYP_EPSILON_SEML2
echo ''


#-----------------------------------------------------
# Position of the poincaré section in NCEM coordinates
# if -1, 3BSOI: about 159198km of radius around the Moon, 
# hence 2.46761593 in NCEM coordinates about EML2
#-----------------------------------------------------
# At EML2
if [ -z ${RPS+x} ]; then
		set_param "RPS" -1
fi
echo "RPS       =" $RPS

#-----------------------------------------------------
# I/O Handling
#-----------------------------------------------------
# At EML2
if [ -z ${IO_HANDLING+x} ]; then
		set_param "IO_HANDLING" $IO_DEFAULT
fi

case $IO_HANDLING in
		$IO_DEFAULT)  echo 'IO_HANDLING  = IO_DEFAULT. Use of the default filenames.'
		;;
		$IO_BASH)     echo 'IO_HANDLING  = IO_BASH. Use of the filames contained in the config files'
		;;
		$IO_DIALOG)   echo 'IO_HANDLING  = IO_DIALOG. A file dialog window pops up each time a filename is needed.'
		;;
		*)     		  echo "IO_HANDLING  = "$IO_HANDLING". Unknown type."
esac
echo ''

		
#=====================================================
#  ---- Projection parameters ----
#=====================================================
# We check that TMIN, the first projection parameter
# is set or not. If not, we set default values for 
# each projection parameters

echo "#------------------------------------------#"
if [ -z ${TMIN+x} ]; then
	
	#-----------------------------------------------------
	# Warn the user
	#-----------------------------------------------------
	echo 'The projection parameters are not set.'
	echo 'Default values are used.'

	#-----------------------------------------------------
	# Parameters that change often
	#-----------------------------------------------------
	# Time grid: min, max and number of points on the grid
	TMIN=0.00     # (given as %T, with T the SEM period)
	TMAX=0.25     # (given as %T, with T the SEM period)
	TSIZE=0	  

	# Configuration (s1, s2, s3, s4) grid
	GLIM_S1=(-35 +35)
	GLIM_S2=(+0  +10)
	GLIM_S3=(-35 +35)
	GLIM_S4=(+0  +10)
	GSIZE_SI=(+50 +5 +50 +5)
	
	# Primary family
	PRIMARY=0
	
	# Fixed delta of energy (-1 if not used)
	DHD=-1
	
	# Time frequency, in %T
	DT=0.005
	
	#-----------------------------------------------------
	# Parameters that are stable
	#-----------------------------------------------------
	TM=12 	      # Maximum integration time (given as %T, with T the SEM period)
	MSIZE=500     # Number of points on each trajectory
	NSMIN=20      # Number of sorted solutions
	YNMAX=0.6     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
	SNMAX=0.6     # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
	NOD=6         # Number of dimensions on which we compute the norm of the projection error
	
	#-----------------------------------------------------
	# Filenames (used only if IO_HANDLING==$IO_BASH
	#-----------------------------------------------------
	if [ -z ${FILE_CU+x} ]; then
		set_param "FILE_CU" "cu.bin"
	fi
	
	if [ -z ${FILE_PCU+x} ]; then
		set_param "FILE_PCU" "projcu.bin"
	fi
else

	#-----------------------------------------------------
	# Primary Family
	#-----------------------------------------------------
	if [ -z ${PRIMARY+x} ]; then
		set_param "PRIMARY" 0
	fi
	
	#-----------------------------------------------------
	# Fixed delta of energy (-1 if not used)
	#-----------------------------------------------------
	if [ -z ${DHD+x} ]; then
		set_param "DHD" -1
	fi
	
	#-----------------------------------------------------
	# Fixed delta of energy (-1 if not used)
	#-----------------------------------------------------
	if [ -z ${DT+x} ]; then
		set_param "DT" 0.005
	fi
	
	#-----------------------------------------------------
	# Display current set of parameters
	#-----------------------------------------------------
	echo 'Current set of projection parameters:'
	echo ''
	echo "TMIN     =" $TMIN
	echo "TMAX     =" $TMAX
	echo "TSIZE    =" $TSIZE
	echo ''
	echo "GLIM_S1  = ["${GLIM_S1[*]}"]"
	echo "GLIM_S2  = ["${GLIM_S2[*]}"]"
	echo "GLIM_S3  = ["${GLIM_S3[*]}"]"
	echo "GLIM_S4  = ["${GLIM_S4[*]}"]"
	echo "GSIZE_SI = ["${GSIZE_SI[*]}"]"
	echo ''
	echo "TM       =" $TM
	echo "MSIZE    =" $MSIZE
	echo "NSMIN    =" $NSMIN
	echo "YNMAX    =" $YNMAX
	echo "SNMAX    =" $SNMAX
	echo "NOD      =" $NOD
	echo ''
	echo "PRIMARY  =" $PRIMARY
	echo ''
	echo "DHD      =" $DHD
	echo "DT       =" $DT
	echo ''
	
	#-----------------------------------------------------
	# Filenames (used only if IO_HANDLING==$IO_BASH)
	#-----------------------------------------------------
	if [ -z ${FILE_CU+x} ]; then
		set_param "FILE_CU" "cu.bin"
	fi
	
	if [ -z ${FILE_PCU+x} ]; then
		set_param "FILE_PCU" "projcu.bin"
	fi
	
	echo "FILE_CU  =" $FILE_CU
	echo "FILE_PCU =" $FILE_PCU
	echo ''
fi

#=====================================================
#  ---- Refinement parameters ----
#=====================================================
# We check that REFST_TYPE, the first refinement parameter
# is set or not. If not, we set default values for 
# each refinement parameters
echo "#------------------------------------------#"
echo ''
echo "#------------------------------------------#"
if [ -z ${REFST_TYPE+x} ]; then
	
	#-----------------------------------------------------
	# Warn the user
	#-----------------------------------------------------
	echo 'The refinement parameters are not set.'
	echo 'Default values are used.'

	#-----------------------------------------------------
	# Parameters that change often
	#-----------------------------------------------------
	REFST_TYPE=$REF_CONT_D           # Type of refinement - rk: set REF_CONT_D_HARD_CASE for difficult cases with REF_CONT_D (ex: EML2-SEMLi via SEML1...)
	REFST_DIM=$REF_PLANAR            # Type of dimensions planar or 3d?
	REFST_T0_DES=0.00                # Initial time - given as %T, with T the SEM period   

	# Domain of search (min s1, max s1, min s3, max s3) for the first guess
	REFST_SI_CMU_EM_LIM=(-35 +35 0 0 -35 +35 0 0)
	# Or, if we want the user to define such domain:
	REFST_ISLIMUD=0
	
	# Domain of search (min s1, max s1, min s3, max s3) for the seed of first guess
	REFST_SI_SEED_EM_LIM=(-40 +40 0 0 -40 +40 0 0)

	# Limits for the time of flight during transfers - not used if -1
	REFST_TOF_LIM=(-1 -1)
	
	# Values for crossings
	REFST_CROSSINGS=-1
	
	# Maximum projection distance allowed during subselection
    REFST_PMAX_DIST_SEM=5e-4

	# Number of steps in the continuation procedure
	REFST_CONT_STEP_MAX=+450;        # with fixed times
	REFST_CONT_STEP_MAX_VT=+50;      # with variable times

	# Initial step in the continuation procedure
	REFST_VAR_TIME_DS0=8e-2	         # for variable times
	if [ $LI_EM == 1 ]; then
		REFST_FIXED_TIME_DS0=3e-2   # for fixed times
	else
		REFST_FIXED_TIME_DS0=5e-1   # for fixed times
	fi

	# Desired number of iterations in Newton's method in the continuation procedure
	REFST_VAR_TIME_NU0=4 	# for variable times
	REFST_FIXED_TIME_NU0=2  # for fixed times

	# Direction of the continuation procedure
	REFST_ISDIRUD=0			 # is it user defined?
	REFST_DIR=-1    		 # if not, +1 or -1

	# User parameters
	REFST_ISFLAGON=1   	         # do we have steps in the procedure - asking the user to press enter to go on?
	REFST_ISPLOTTED=1   		 # do we plot the results during the computation?
	REFST_ISSAVED=1     		 # do we save the results in data files?
	REFST_ISFROMSERVER=1		 # does the raw data comes from server files?
	
	# Maximum angle around SEMLi if REF_COND_T is used (in degrees)
	REFST_THETAMAX=240               # should be a multiple of 90°

	#-----------------------------------------------------
	# Parameters that are stable
	#-----------------------------------------------------
	REFST_ISDEBUG=0			 # if yes, additionnal tests are performed
	REFST_GRIDSIZE=20        	 # number of points on the refinement grid
	REFST_MPLOT=200        	         # number of points per plot between to pach points (e.g. total plot points is REFST_MPLOT*REFST_GRIDSIZE)
	REFST_FPLOT=5                    # frequency of plotting in hours
	
	REFST_TIME=$REF_VAR_TN		 # type of constraints on the times in REF_CONT
	REFST_GRID=$REF_FIXED_GRID	 # type of grid
	REFST_TERMINATION=$REF_COND_T    # Termination condition in the continuation with variable final time (either REF_VAR_TN/REF_VAR_TIME)
	REFST_COORD_TYPE=$NCSEM		 # coordinates system in the refinement procedure

	REFST_XPS=0.7			     # position of the poincaré section in NCSEM coordinates
	REFST_XPS_NCEM=0.8			 # position of the poincaré section in NCEM coordinates
	REFST_ISJPL=1		         # is the JPL refinement performed when possible?
	REFST_DJPLCOORD=-1		     # coordinate system used during the JPL refinement (if -1, it is user defined) Best results obtained with $NJ2000
	REFST_SIDIM=0		         # 0 or 2 - component of s0 that stays constant when t0 is free

	# Sampling frequencies in REF_COMP (complete trajectory) in days
	REFST_SF_EML2=2			 # orbit at EML2	
	REFST_SF_MAN=5			 # transfer leg
	REFST_SF_SEML2=10		 # orbit at SEML2

	# Integration window for each orbit
	REFST_TSPAN_EM=10  		 # given as %T, where T is the SEM period, in EM units
	REFST_TSPAN_SEM=10 		 # given as %T, where T is the SEM period, in SEM units
	
	# Storing the orbits at each step?
	REFST_ISSAVED_EM=0   # 0: don't save, 1: save using projection method
	REFST_ISSAVED_SEM=0  # 0: don't save, 1: save using projection method, 2: save using integration in reduced coordinates
	
	# Type of computation for each orbit
	REFST_COMP_ORB_EM=$INT_TRY_BOTH
	REFST_COMP_ORB_SEM=$INT_TRY_BOTH
	
	# Type of time selection
	REFST_TYPE_OF_T_SEL=$TIME_SELECTION_RATIO
	
	# Desired number of solutions
	REFST_NREF=-1
	
	#-----------------------------------------------------
	# Filenames (used only if IO_HANDLING==$IO_BASH)
	#-----------------------------------------------------
	FILE_CONT="cont_atf.txt"
	FILE_CONT_RES="cont_atf_traj.bin"
	FILE_TRAJ_FROM_W="traj_from_w.bin"
	FILE_TRAJ_FROM_C="traj_from_c.bin"
	FILE_JPL="cont_jpl.bin"
	
else
	#-----------------------------------------------------
	# Display current set of parameters
	#-----------------------------------------------------
	echo 'Current set of refinement parameters:'
	echo ''
	
	#-----------------------------------------------------
	# Crossings
	#-----------------------------------------------------
	if [ -z ${REFST_CROSSINGS+x} ]; then
		set_param "REFST_CROSSINGS" -1
	fi
	
	#-----------------------------------------------------
	# Crossings
	#-----------------------------------------------------
	if [ -z ${REFST_PMAX_DIST_SEM+x} ]; then
		set_param "REFST_PMAX_DIST_SEM" 5e-4
	fi
	
	#-----------------------------------------------------
	# Type of time selection
	#-----------------------------------------------------
	if [ -z ${REFST_TYPE_OF_T_SEL+x} ]; then
		set_param "REFST_TYPE_OF_T_SEL" $TIME_SELECTION_RATIO
	fi
	
	#-----------------------------------------------------
	# Domain of search (min s1, max s1, min s3, max s3) for the seed of first guess
	#-----------------------------------------------------
	if [ -z ${REFST_SI_SEED_EM_LIM+x} ]; then
		TEMP=(-40 +40 0 0 -40 +40 0 0)
		set_param "REFST_SI_SEED_EM_LIM" TEMP[@]
	fi
	
	#-----------------------------------------------------
	# position of the poincaré section in NCSEM & NCEM coordinates
	#-----------------------------------------------------
	if [ -z ${REFST_XPS+x} ]; then
		set_param "REFST_XPS" 0.7
	fi
	
	if [ -z ${REFST_XPS_NCEM+x} ]; then
		set_param "REFST_XPS_NCEM" 0.8
	fi
	
	#-----------------------------------------------------
	# Type of computation for each orbit
	#-----------------------------------------------------
	if [ -z ${REFST_COMP_ORB_EM+x} ]; then
		set_param "REFST_COMP_ORB_EM" $INT_TRY_BOTH
	fi
	
	if [ -z ${REFST_COMP_ORB_SEM+x} ]; then
		set_param "REFST_COMP_ORB_SEM" $INT_TRY_BOTH
	fi
	

	#-----------------------------------------------------
	# Parameters that change often
	#-----------------------------------------------------
	#REFST_TYPE
	case $REFST_TYPE in
		$REF_SINGLE)            echo 'REFST_TYPE             = REF_SINGLE'
		;;
		$REF_ORBIT)             echo 'REF_ORBIT              = REF_ORBIT'
		;;
		$REF_CONT) 	            echo 'REFST_TYPE             = REF_CONT'
		;;
		$REF_CONT_D)            echo 'REFST_TYPE             = REF_CONT_D'
		;;
		$REF_CONT_D_HARD_CASE)  echo 'REFST_TYPE             = REF_CONT_D_HARD_CASE'
		;;
		$REF_COMP)  		    echo 'REFST_TYPE             = REF_COMP'
		;;
		*)     			        echo "REFST_TYPE             = "$REFST_TYPE". Unknown type."
	esac

	#REFST_DIM
	case $REFST_DIM in
		$REF_PLANAR) echo 'REFST_DIM              = REF_PLANAR'
		;;
		$REF_3D)     echo 'REFST_DIM              = REF_3D'
		;;
		$REF_MIXED)  echo 'REFST_DIM              = REF_MIXED'
		;;
		*)     	     echo "REFST_DIM              = "$REFST_DIM". Unknown type."
	esac

	echo "REFST_T0_DES           =" $REFST_T0_DES
	echo ''
	echo "REFST_SI_CMU_EM_LIM    = ["${REFST_SI_CMU_EM_LIM[*]}"]"
	echo "REFST_ISLIMUD          =" $REFST_ISLIMUD
	echo "REFST_SI_SEED_EM_LIM   = ["${REFST_SI_SEED_EM_LIM[*]}"]"
	echo "REFST_TOF_LIM          = ["${REFST_TOF_LIM[*]}"]"
	echo "REFST_CROSSINGS        =" $REFST_CROSSINGS
	echo "REFST_PMAX_DIST_SEM    =" $REFST_PMAX_DIST_SEM
	echo ''
	echo "REFST_CONT_STEP_MAX    =" $REFST_CONT_STEP_MAX
	echo "REFST_CONT_STEP_MAX_VT =" $REFST_CONT_STEP_MAX_VT
	echo ''
	echo "REFST_VAR_TIME_DS0     =" $REFST_VAR_TIME_DS0
	echo "REFST_FIXED_TIME_DS0   =" $REFST_FIXED_TIME_DS0
	echo ''
	echo "REFST_VAR_TIME_NU0     =" $REFST_VAR_TIME_NU0
	echo "REFST_FIXED_TIME_NU0   =" $REFST_FIXED_TIME_NU0
	echo ''
	echo "REFST_ISDIRUD          =" $REFST_ISDIRUD
	echo "REFST_DIR              =" $REFST_DIR
	echo ''
	echo "REFST_ISFLAGON         =" $REFST_ISFLAGON
	echo "REFST_ISPLOTTED        =" $REFST_ISPLOTTED
	echo "REFST_ISSAVED          =" $REFST_ISSAVED
	echo "REFST_ISFROMSERVER     =" $REFST_ISFROMSERVER
	echo ''
	echo "REFST_THETAMAX         =" $REFST_THETAMAX
	echo ''

	#-----------------------------------------------------
	# Parameters that are stable
	#-----------------------------------------------------
	#REFST_TIME
	case $REFST_TIME in
		$REF_FIXED_TIME) echo 'REFST_TIME            = REF_FIXED_TIME'
		;;
		$REF_VAR_TN) 	 echo 'REFST_TIME            = REF_VAR_TN'
		;;
		$REF_VAR_TIME) 	 echo 'REFST_TIME            = REF_VAR_TIME'
		;;
		*)     		 echo "REFST_TIME            = "$REFST_TIME". Unknown type."
	esac
	
	#REFST_GRID
	case $REFST_GRID in
		$REF_FIXED_GRID) echo 'REFST_GRID            = REF_FIXED_GRID'
		;;
		$REF_VAR_GRID) 	 echo 'REFST_GRID            = REF_VAR_GRID'
		;;
		$REF_GIVEN_GRID) echo 'REFST_GRID            = REF_GIVEN_GRID'
		;;
		*)     		     echo "REFST_GRID            = "$REFST_GRID". Unknown type."
	esac

	#REFST_TERMINATION
	case $REFST_TERMINATION in
		$REF_COND_S5) echo 'REFST_TERMINATION     = REF_COND_S5'
		;;
		$REF_COND_T)  echo 'REFST_TERMINATION     = REF_COND_T'
		;;
		*)     	      echo "REFST_TERMINATION     = "$REFST_TERMINATION". Unknown type."
	esac
	
	#REFST_COORD_TYPE
	case $REFST_COORD_TYPE in
		$NCSEM)  echo 'REFST_COORD_TYPE      = NCSEM'
		;;
		$NCEM)   echo 'REFST_COORD_TYPE      = NCEM'
		;;
		$VNCSEM) echo 'REFST_COORD_TYPE      = VNCSEM'
		;;
		$VNCEM)  echo 'REFST_COORD_TYPE      = VNCEM'
		;;
		$PSEM)   echo 'REFST_COORD_TYPE      = PSEM'
		;;
		$PEM)    echo 'REFST_COORD_TYPE      = PEM'
		;;
		$VSEM)   echo 'REFST_COORD_TYPE      = VSEM'
		;;
		$VEM)    echo 'REFST_COORD_TYPE      = VEM'
		;;
		*)       echo "REFST_COORD_TYPE      = "$REFST_COORD_TYPE". Unknown type."
	esac
	echo ''
	
	#-----------------------------------------------------
	# REFST_FPLOT
	#-----------------------------------------------------
	if [ -z ${REFST_FPLOT+x} ]; then
		set_param "REFST_FPLOT" 5
	fi

	echo "REFST_ISDEBUG         =" $REFST_ISDEBUG
	echo "REFST_GRIDSIZE        =" $REFST_GRIDSIZE
	echo "REFST_MPLOT           =" $REFST_MPLOT
	echo "REFST_FPLOT           =" $REFST_FPLOT
	echo ''
	echo "REFST_XPS             =" $REFST_XPS
	echo "REFST_XPS_NCEM        =" $REFST_XPS_NCEM
	echo "REFST_ISJPL           =" $REFST_ISJPL
	echo "REFST_DJPLCOORD       =" $REFST_DJPLCOORD
	echo "REFST_SIDIM           =" $REFST_SIDIM
	echo ''
	echo "REFST_SF_EML2         =" $REFST_SF_EML2
	echo "REFST_SF_MAN          =" $REFST_SF_MAN
	echo "REFST_SF_SEML2        =" $REFST_SF_SEML2
	echo ''
	echo "REFST_TSPAN_EM        =" $REFST_TSPAN_EM
	echo "REFST_TSPAN_SEM       =" $REFST_TSPAN_SEM
	echo ''
	echo "REFST_ISSAVED_EM      =" $REFST_ISSAVED_EM
	echo "REFST_ISSAVED_SEM     =" $REFST_ISSAVED_SEM
	echo ''
	echo "REFST_COMP_ORB_EM     =" $REFST_COMP_ORB_EM
	echo "REFST_COMP_ORB_SEM    =" $REFST_COMP_ORB_SEM
	echo ''
	
	
	
	#-----------------------------------------------------
	# Type of time selection
	#-----------------------------------------------------
	#REFST_TYPE_OF_T_SEL
	case $REFST_TYPE_OF_T_SEL in
		$TIME_SELECTION_ABSOLUTE) echo 'REFST_TYPE_OF_T_SEL   = TIME_SELECTION_ABSOLUTE'
		;;
		$TIME_SELECTION_RATIO)  echo 'REFST_TYPE_OF_T_SEL   = TIME_SELECTION_RATIO'
		;;
		*) echo "REFST_TYPE_OF_T_SEL   = "$REFST_TYPE_OF_T_SEL". Unknown type."
	esac
	
	#-----------------------------------------------------
	# Desired number of solutions
	#-----------------------------------------------------
	if [ -z ${REFST_NREF+x} ]; then
		set_param "REFST_NREF" -1
	fi
	echo "REFST_NREF            =" $REFST_NREF
	echo ''
	
	#-----------------------------------------------------
	# Filenames (used only if IO_HANDLING==$IO_BASH)
	#-----------------------------------------------------
	if [ -z ${FILE_PCU+x} ]; then
		set_param "FILE_PCU" "projcu.bin"
	fi
	
	if [ -z ${FILE_CONT+x} ]; then
		set_param "FILE_CONT" "cont_atf.txt"
	fi
	
	if [ -z ${FILE_CONT_RES+x} ]; then
		set_param "FILE_CONT_RES" "cont_atf_traj.bin"
	fi
	
	if [ -z ${FILE_TRAJ_FROM_W+x} ]; then
		set_param "FILE_TRAJ_FROM_W" "traj_from_w.bin"
	fi
	
	if [ -z ${FILE_TRAJ_FROM_C+x} ]; then
		set_param "FILE_TRAJ_FROM_C" "traj_from_c.bin"
	fi
	
	if [ -z ${FILE_JPL+x} ]; then
		set_param "FILE_JPL" "cont_jpl.bin"
	fi
	
	echo "FILE_PCU              =" $FILE_PCU
	echo "FILE_CONT             =" $FILE_CONT
	echo "FILE_CONT_RES         =" $FILE_CONT_RES
	echo "FILE_TRAJ_FROM_W      =" $FILE_TRAJ_FROM_W
	echo "FILE_TRAJ_FROM_C      =" $FILE_TRAJ_FROM_C
	echo "FILE_JPL              =" $FILE_JPL
	echo ''
fi


#-----------------------------------------------------------------------------------------
# Go on with the implementation?
#-----------------------------------------------------------------------------------------
echo "#------------------------------------------#"
echo ''
echo -e "Do you want to go on with the computation (y/n)? \c "
read  ans
# bash check the answer
if [ "$ans" == "y" ]; then

	#-----------------------------------------------------------------------------
	#Build the array of coefficients
	#-----------------------------------------------------------------------------

	# The general parameters are common to all types of computation
	COEFFS=($OFTS_ORDER $OFS_ORDER $COMP_TYPE $MODEL $LI_EM $LI_SEM $ISPAR $NUM_THREADS $HYP_EPSILON_EML2 $HYP_EPSILON_SEML2 $RPS $IO_HANDLING)

	# Then, for each type, we add some parameters
	case $COMP_TYPE in
		$COMP_CM_EML2_TO_CM_SEML | $COMP_CM_EML2_TO_CM_SEML_3D | $COMP_CM_EML2_TO_CM_SEML_H | $COMP_ORBIT_EML2_TO_CM_SEML | $COMP_SINGLE_ORBIT_EML2_TO_CM_SEML | $COMP_CMU_SEMLi_TO_CM_EMLj)        
			COEFFS=(${COEFFS[*]}  $TMIN $TMAX $TM $TSIZE)
			COEFFS=(${COEFFS[*]}  ${GLIM_S1[*]} ${GLIM_S2[*]} ${GLIM_S3[*]} ${GLIM_S4[*]} ${GSIZE_SI[*]})
			COEFFS=(${COEFFS[*]}  $PRIMARY $MSIZE $NSMIN $YNMAX $SNMAX $NOD $DHD $DT $FILE_CU $FILE_PCU)
		;;
		$COMP_CM_EML2_TO_CMS_SEML | $COMP_REF_JPL)
			COEFFS=(${COEFFS[*]}  $REFST_TYPE $REFST_DIM $REFST_T0_DES)
		    COEFFS=(${COEFFS[*]}  $REFST_ISDIRUD $REFST_DIR)
		    COEFFS=(${COEFFS[*]}  ${REFST_SI_CMU_EM_LIM[*]} $REFST_ISLIMUD)
		    COEFFS=(${COEFFS[*]}  ${REFST_SI_SEED_EM_LIM[*]})
			COEFFS=(${COEFFS[*]}  ${REFST_TOF_LIM[*]} $REFST_CROSSINGS $REFST_PMAX_DIST_SEM)
		  	COEFFS=(${COEFFS[*]}  $REFST_CONT_STEP_MAX $REFST_CONT_STEP_MAX_VT)
			COEFFS=(${COEFFS[*]}  $REFST_FIXED_TIME_DS0 $REFST_VAR_TIME_DS0)
			COEFFS=(${COEFFS[*]}  $REFST_FIXED_TIME_NU0 $REFST_VAR_TIME_NU0)
			COEFFS=(${COEFFS[*]}  $REFST_ISFLAGON $REFST_ISPLOTTED $REFST_ISSAVED $REFST_ISFROMSERVER)
			COEFFS=(${COEFFS[*]}  $REFST_THETAMAX)

			COEFFS=(${COEFFS[*]}  $REFST_ISDEBUG $REFST_GRIDSIZE $REFST_MPLOT $REFST_FPLOT)
			COEFFS=(${COEFFS[*]}  $REFST_TIME $REFST_GRID $REFST_TERMINATION $REFST_COORD_TYPE)
			COEFFS=(${COEFFS[*]}  $REFST_XPS $REFST_ISJPL $REFST_DJPLCOORD $REFST_SIDIM)
			COEFFS=(${COEFFS[*]}  $REFST_SF_EML2 $REFST_SF_MAN $REFST_SF_SEML2)
			COEFFS=(${COEFFS[*]}  $REFST_TSPAN_EM $REFST_TSPAN_SEM)
			COEFFS=(${COEFFS[*]}  $REFST_ISSAVED_EM $REFST_ISSAVED_SEM)
			COEFFS=(${COEFFS[*]}  $REFST_COMP_ORB_EM $REFST_COMP_ORB_SEM)
			COEFFS=(${COEFFS[*]}  $REFST_TYPE_OF_T_SEL $REFST_NREF)
			
			COEFFS=(${COEFFS[*]}   $FILE_PCU $FILE_CONT $FILE_CONT_RES $FILE_TRAJ_FROM_W $FILE_TRAJ_FROM_C $FILE_JPL)
		;;
		$COMP_SINGLE_ORBIT)              

		;;
		$COMP_CM_EML2_TO_CMS_SEML_READ)  

		;;
		$COMP_VF_TEST)                   

		;;
		$COMP_CM_EML2_TO_CM_SEML_REFINE) 

		;;
		$COMP_EPHEMERIDES_TEST)          

		;;
		$COMP_VOFTS_TO_VOFTS)            

		;;
		$COMP_test_INVMAN) 		 

		;;
		*)  echo "COMP_TYPE    = "$COMP_TYPE". Unknown type. Return."
			exit 0
	esac
	 

	#-------------------------------
	#Call software
	#-------------------------------
	#Check the NOHUP condition
	if [ "$ISNOHUP" == "1" ]; then
			#If true, check that an output file has been set
			if [ -z ${OUT+x} ]; then
				echo 'WARNING: OUT variable has not been set.'
				echo 'OUT = default.out by default.'
				OUT='default.out'
			fi
			#Call software
			nohup bin/Release/FourierTaylorInCpp ${COEFFS[*]} > $OUT &
	else
			bin/Release/FourierTaylorInCpp ${COEFFS[*]}
	fi

	notify-send "Computation is done in FourierTaylorInCpp"
	
else  
	echo "Stop. No computation."
fi 

