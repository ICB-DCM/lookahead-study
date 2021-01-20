#!/bin/sh

# Kills all jobs on BQ_CLuster
PORT=${1}

# Kill redis-first
REDIS_JOB_ID=$(squeue --Format="jobid,name" | grep redis_${PORT} | awk '{print $1}')
scancel ${REDIS_JOB_ID}

# Kill python script
PYTHON_JOB_ID=$(squeue --Format="jobid,name" | grep python_${PORT} | awk '{print $1}')
scancel ${PYTHON_JOB_ID}

# Kill all workers
WORKER_JOB_IDS=$(squeue --Format="jobid,name" | grep worker_${PORT} | awk '{ print $1 }')

for i in ${WORKER_JOB_IDS}
do 
	scancel ${i}
done

echo "Redis, python and Jobs killed"
