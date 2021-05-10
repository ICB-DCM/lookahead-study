#!/bin/sh
# Master script


PORT=${1}
N_NODES=${2}
PARTITION=${3}
TIME=${4}
CPUSPERTASK=${5}
PYTHONFILE=${6}
if [ -z "${N_NODES}" ]
then
        N_NODES=1
fi

# Submit redis 
JOB_NAME=$(sbatch --time=${TIME} --job-name=redis_${PORT} --partition=${PARTITION} submit_redis.sh ${PWD} ${PORT})
JOB_ID=$(echo $JOB_NAME | cut -f4 -d' ') # Output sth like 8670323

echo "Job_name: $JOB_NAME "
echo "JOB_ID: $JOB_ID "
# Wait for redis to start
# On PBS/Torque I had a while loop here that checked when the status of the job switched to RUNNING
while [ $(squeue --Format="jobid,name,state" | grep redis_${PORT} | awk '{print $3}') != "RUNNING" ]
do
        echo 'Waiting for Redis to start'
        sleep 10
done
# Retrieve compute node of that job ID
COMP_NODE=$(squeue | grep ${JOB_ID} | awk '{ print $8 }')

# Retrieve IP of compute node
REDIS_IP=$(host ${COMP_NODE}i | awk '{ print $4 }')
echo 'Redis IP:' ${REDIS_IP}

# Start worker with IP and port
for i in $(seq 1 ${N_NODES})
do 
	NAME=worker_${PORT}_${i}
	sbatch --time=${TIME} --job-name=${NAME}  --dependency=after:$JOB_ID --cpus-per-task=${CPUSPERTASK} --partition=${PARTITION}  submit_worker.sh ${REDIS_IP} ${PORT} ${TIME} ${CPUSPERTASK} ${PWD}
done

echo 'last worker id=worker_'${PORT}_${N_NODES}
while [ $(squeue --Format="jobid,name,state" | grep worker_${PORT}_${N_NODES} | awk '{print $3}') != "RUNNING" ]
do
        echo 'Waiting for worker to start'
        sleep 10
done

# Start python script pass IP and port
sbatch --time=${TIME} --job-name=python_${PORT} --dependency=after:$JOB_ID --partition=${PARTITION} submit_python.sh ${PWD} ${REDIS_IP} ${PORT} ${PYTHONFILE}


