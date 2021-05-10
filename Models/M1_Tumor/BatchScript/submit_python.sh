#!/bin/sh

cur_PWD=${1}
IP=${2}
PORT=${3}
N_NODES=${4}
PYHTONFILE=${5}
cd ${cur_PWD}
echo '#### This is the python script job ####'

# Soruce Modules
source ./load_module.sh

# Start python script
python ${PYHTONFILE} --port ${PORT} --ip ${IP} --nodes ${N_NODES}> out.txt 2> errpy_${PORT}.txt
