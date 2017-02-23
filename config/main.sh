# main script file for FourierTaylorInCpp
# RECALL: make it executable with "$ chmod +x main.sh"

#--------------------------------------------------------------
# Source the configuration file passed as argument
#--------------------------------------------------------------
source $1

#--------------------------------------------------------------
# Title
#--------------------------------------------------------------
echo
echo "#------------------------------------------#"
echo "#          OOFTDA configuration            #"
echo "#------------------------------------------#"
echo "Source file: " "$1"

#=====================================================
#  ---- General parameters ----
#=====================================================


#-----------------------------------------------------
# TYPE OF COMPUTATION (COMP_TYPE)
#-----------------------------------------------------
#Check that the variable $COMP_TYPE exist
if [ -z ${COMP_TYPE+x} ]; then
	echo 'WARNING: the variable COMP_TYPE is not set.'
	echo 'No computation.'
	return
fi


#-----------------------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#-----------------------------------------------------
#Check that the variable $MODEL exist
if [ -z ${MODEL+x} ]; then
	#If not, default values to QBCP
	echo 'WARNING: the variable MODEL is not set.'
	echo 'QBCP is used by default.'
	MODEL=$QBCP
else
	if [ "$MODEL" == "$RTBP" ]; then
		OFS_ORDER=0
	else  
		OFS_ORDER=30
	fi
fi


#=====================================================
#  ---- Projection parameters ----
#=====================================================
# We check that TMIN, the first projection parameter
# is set or not. If not, we set default values for 
# each projection parameters

if [ -z ${TMIN+x} ]; then
	
	# Warn the user
	echo 'The projection parameters are not set.'
	echo 'Default values are used.'

	# Time grid: min, max and number of points on the grid
	TMIN=0.00     # (given as %T, with T the SEM period)
	TMAX=0.25     # (given as %T, with T the SEM period)
	TSIZE=0	  

	# Configuration (s1, s2, s3, s4) grid
	GLIM_S1=(-35 +35)
	GLIM_S2=(+0  +10)
	GLIM_S3=(-35 +35)
	GLIM_S4=(+0  +10)
	
	# Parameters that are stable
	TM=12 	      # Maximum integration time (given as %T, with T the SEM period)
	MSIZE=500     # Number of points on each trajectory
	NSMIN=20      # Number of sorted solutions
	YNMAX=0.5     # The maximum norm (in SEM normalized units) for a projection to occur on the CM_NC of SEMLi
	SNMAX=0.6     # The maximum norm (in RCM normalized units) for a projection to occur on the CM_NC of SEMLi
	NOD=6         # Number of dimensions on which we compute the norm of the projection error
fi



