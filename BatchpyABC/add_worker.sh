#!/bin/sh
# Master script

PORT=${1}
N_NODES=${2}
PARTITION=${3}
TIME=${4}
CPUSPERTASK=${5}

if [ -z "${N_NODES}" ] # -z == string is null, that is, has zero length
then
        N_NODES=1
fi

# Get Redis JOB ID
JOB_ID=$(squeue --Format="jobid,name" | grep redis | grep ${PORT} | awk '{print $1}')

# Retrieve compute node of that job ID
COMP_NODE=$(squeue | grep $TEST | awk '{ print $8 }')

# Retrieve IP of compute node
REDIS_IP=$(host ${COMP_NODE} | awk '{ print $4 }')

# Find already working workers on that PORT
MAX_WORKER_ID=$(squeue --Format="jobid,name" | grep ${PORT} | grep worker | awk '{print $2}' | cut -f3 -d'_' | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print max}')

if [ -z "${MAX_WORKER_ID}" ]
then
        MAX_WORKER_ID=0
fi

NEW_START=$((${MAX_WORKER_ID}+1))
NEW_END=$((${N_NODES}+${NEW_START}-1))

for i in $(seq ${NEW_START} ${NEW_END})
do 
	NAME=worker_${PORT}_${i}
	sbatch --time=${TIME} --job-name=${NAME} --cpus-per-task ${CPUSPERTASK} --partition=${PARTITION} --account=fitmulticell submit_worker.sh ${REDIS_IP} ${PORT} ${TIME} ${CPUSPERTASK}
done
