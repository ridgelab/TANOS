#! /bin/bash

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

#	Some handy variables
#${SLURM_MEM_PER_CPU}
#${SLURM_MEM_PER_NODE}
#${SLURM_JOB_NAME}
#${SLURM_NTASKS}
#${SLURM_JOB_NUM_NODES}
#${SLURM_JOB_ID}
#${SLURM_ARRAY_JOB_ID}
#${SLURM_ARRAY_TASK_ID}
#${SLURM_ARRAY_TASK_COUNT}
#${SLURM_ARRAY_TASK_MIN}
#${SLURM_ARRAY_TASK_MAX}

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

if [ -z "${SLURM_ARRAY_TASK_ID}" ]
then
	printf "%s\n" "ERROR: SLURM_ARRAY_TASK_ID not defined." 1>&2
	control_c
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

	rm -rf /tmp/${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID} &> /dev/null
}

control_c()
{
	kill -SIGINT `jobs -p`
	cleanup
	exit 1
}

trap control_c SIGHUP SIGINT SIGTERM SIGQUIT

moveTempFilesBackToNetworkStorage()
{
	if [ -n $WORK_DIR ] && [ -e $WORK_DIR ] && [ -d $WORK_DIR ]
	then
		cp -v "${WORK_DIR}"/* `dirname "${OUTPUT_PFX}"` 1>&2
	fi
}

outOfTime()
{
	printf "%s\n" "This job ran out of time! SLURM sent signal USR1 and now we're trying to quite gracefully. (fingers crossed!)" 1>&2
	kill -SIGINT `jobs -p`

	printf "%s\n" "Attempting to move files from /tmp to network storage..." 1>&2
	moveTempFilesBackToNetworkStorage

	printf "%s\n" "Now using 'cleanup' function with status 'success'. Be advised: this process ran out of time- you will need to run this again with more time (unless it checkpoints, in which case it will still need to run again, just not necessarily with more time)." 1>&2
	cleanup

	exit 10 # SIGUSR1 == 10
}

trap outOfTime USR1

# 	load modules
module purge
# module load iqtree # we used version 1.6.12

#	check that enough input was given
if [ $# -ne 3 ]
then
	printf "%s\n" "ERROR: Expected 3 arguments to this slurm file." 1>&2
	cleanup
	exit 1
fi

#	setup variables for the job
SPLIT_NUM="${SLURM_ARRAY_TASK_ID}"

INPUT_ALN=$(readlink -m "${1}")
MODEL="${2}"
OUTPUT_PFX="${3}-${SPLIT_NUM}"

OUTPUT_FILES=("${OUTPUT_PFX}".{bionj,log,ckp.gz,treefile,iqtree,mldist})

# 	check for existence of input file(s)
#		We assume iqtree is capable of recognizing whether the file(s)
#		it requires exist(s).

# 	check for existence of expected output file(s)
# 		We will not do this because iqtree can resume where it
#		left off and, if already completed, will not obliterate
#		output from a previous command unless told to do so.

# 	Since we are working in /tmp space, we'll need to copy output
#	files from a previous run, if they exist to the /tmp space.
#	We will also need to create the output directories in the
#	compute space, if they don't already exist, and in /tmp space.

#		Identify output directories and files that already exist
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

#		create output directories, if needed
mkdir -p "${OUTPUT_DIRS[@]}" &> /dev/null
unset OUTPUT_DIRS

#		create /tmp output directory and copy existing files
WORK_DIR="/tmp/${SLURM_ARRAY_JOB_ID}-${SPLIT_NUM}"
mkdir -p "${WORK_DIR}" &> /dev/null

if [ ${#ALREADY_EXISTING_OUTPUT_FILES[@]} -gt 0 ]
then
	cp -v "${ALREADY_EXISTING_OUTPUT_FILES[@]}" "${WORK_DIR}" 1>&2
fi
unset ALREADY_EXISTING_OUTPUT_FILES

# write tmp info for later cleanup, if needed
#	get the lock
dotlockfile -l "${SLURM_SUBMIT_DIR}/.cleanup.tsv.lock"
printf "%s\t%s\n" "${SLURM_JOB_NODELIST}" "${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}" >> "${SLURM_SUBMIT_DIR}/cleanup.tsv"
dotlockfile -u "${SLURM_SUBMIT_DIR}/.cleanup.tsv.lock"

#	run the program of interest
printf "%s\n" "Running the following now:" 1>&2
printf "\t%s\n" "time iqtree -nt ${SLURM_NTASKS} -mem ${MEM_JOB_IN_GB}G -s ${INPUT_ALN} -t RANDOM -pre ${WORK_DIR}/`basename ${OUTPUT_PFX}` -m ${MODEL}" 1>&2
time iqtree2 \
	-nt ${SLURM_NTASKS} \
	-mem ${MEM_JOB_IN_GB}G \
	-s "${INPUT_ALN}" \
	-t "RANDOM" \
	-pre "${WORK_DIR}/`basename ${OUTPUT_PFX}`" \
	-m "${MODEL}" &

wait `jobs -p`
EXIT_CODE=$?

moveTempFilesBackToNetworkStorage

cleanup

exit ${EXIT_CODE}

