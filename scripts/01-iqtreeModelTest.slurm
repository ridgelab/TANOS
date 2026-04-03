#! /bin/bash

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=${SLURM_CPUS_ON_NODE:-1}

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#	Some handy variables
#${SLURM_MEM_PER_CPU}
#${SLURM_MEM_PER_NODE}
#${SLURM_JOB_NAME}
#${SLURM_NTASKS}
#${SLURM_JOB_NUM_NODES}

if [ -n "$SLURM_JOB_ID" ] # basically, if this is managed by slurm vs being run locally
then
	if [ -n "$SLURM_JOB_NUM_NODES" ] && [ $SLURM_JOB_NUM_NODES -ne 1 ]
	then
		printf "%s\n" "This job is meant to be run with a single node" 1>&2
		exit 1
	elif [ -n "$SLURM_MEM_PER_CPU" ]
	then
		MEM_TASK_IN_MB=${SLURM_MEM_PER_CPU}
		MEM_JOB_IN_MB=$((${MEM_TASK_IN_MB}*${SLURM_NTASKS}))
		MEM_JOB_IN_GB=$((${MEM_JOB_IN_MB}/1024))
	elif [ -n "$SLURM_MEM_PER_NODE" ]
	then
		MEM_JOB_IN_MB=$((${SLURM_MEM_PER_NODE}*${SLURM_JOB_NUM_NODES}))
		MEM_JOB_IN_GB=$((${MEM_JOB_IN_MB}/1024))
		MEM_TASK_IN_MB=$(bc <<< "${MEM_JOB_IN_MB}/${SLURM_NTASKS}")
	else
		printf "%s\n" '$SLURM_MEM_PER_NODE and $SLURM_MEM_PER_CPU not specificed.' 1>&2
		exit 1
	fi
fi

#	move into the correct place
if [ -n "${SLURM_SUBMIT_DIR}" ]
then
	cd "$SLURM_SUBMIT_DIR"
else
	SLURM_SUBMIT_DIR=.
fi

#	manage job cleanup
cleanup()
{
	# cleanup tmp dir
	if [ -n $SLURM_JOB_ID ] && [ -e /tmp/${SLURM_JOB_ID} ]
	then
		rm -rf /tmp/${SLURM_JOB_ID} &> /dev/null
	elif [ -e /tmp/${$} ]
	then
		rm -rf /tmp/${$} &> /dev/null
	fi
}

control_c()
{
	cleanup
	exit 1
}

trap control_c SIGHUP SIGINT SIGTERM SIGQUIT

# 	load modules
module purge
# module load iqtree # we used version 1.6.12

#	setup variables for the job
#	# check that enough input was given
if [ $# -ne 2 ]
then
	printf "%s\n" "ERROR: Expected 2 arguments to this slurm file." 1>&2
	cleanup
	exit 1
fi

INPUT_ALN="${1}"
OUTPUT_PFX="${2}"
OUTPUT_FILES=("${OUTPUT_PFX}".{log,{ckp,model}.gz,treefile,iqtree})

# 	check for existence of input file(s)
#		We assume iqtree is capable of recognizing whether the file(s)
#		it requires exist(s).

# 	check for existence of expected output file(s)
declare -a OUTPUT_DIRS
declare -a ALREADY_EXISTING_OUTPUT_FILES
for OUTPUT_FILE in "${OUTPUT_FILES[@]}"
do
	OUTPUT_DIRS+=($(readlink -m `dirname "${OUTPUT_FILE}"`))

	if [ -e "${OUTPUT_FILE}" ]
	then
		ALREADY_EXISTING_OUTPUT_FILES+=("${OUTPUT_FILE}")
	fi
done
OUTPUT_DIRS=($(printf "%s\n" "${OUTPUT_DIRS[@]}" | sort | uniq | tr '\n' ' '))

if [ ${#ALREADY_EXISTING_OUTPUT_FILES[@]} -eq ${#OUTPUT_FILES[@]} ]
then
	printf "%s\n" "INFO: All ${#OUTPUT_FILES[@]} expected output files already exist! We assume that means this process can quit without running the intended command. Bye!" 1>&2
	cleanup
	exit 0

elif [ ${#ALREADY_EXISTING_OUTPUT_FILES[@]} -gt 0 ]
then
	printf "%s\n" "ERROR: ${#ALREADY_EXISTING_OUTPUT_FILES[@]} of ${#OUTPUT_FILES[@]} expected output files already exist. This is strange; please investigate. You might consider simply running the following command(s):" 1>&2
	printf "\trm -f \"%s\"\n" "${ALREADY_EXISTING_OUTPUT_FILES[@]}" 1>&2
	cleanup
	exit 1
fi
unset ALREADY_EXISTING_OUTPUT_FILES

#	create output directories, if needed
mkdir -p "${OUTPUT_DIRS[@]}" &> /dev/null
unset OUTPUT_DIRS

# sets defaults for threads and mem in case the file is not being run on the cluster
THREADS=${SLURM_NTASKS:-4}
MEM=${MEM_JOB_IN_GB:-8}

#	run the program of interest
time iqtree2 \
	-nt ${THREADS} \
	-mem ${MEM}G \
	-s "${INPUT_ALN}" \
	-pre "${OUTPUT_PFX}" \
	-m TESTONLY

EXIT_CODE=$?

#	cleanup and exit
if [ ${EXIT_CODE} -eq 0 ]
then
	chmod 444 "${OUTPUT_FILES[@]}" &> /dev/null
else
	rm -f "${OUTPUT_FILES[@]}" &> /dev/null
fi

cleanup
exit ${EXIT_CODE}

