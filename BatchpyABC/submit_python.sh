#!/bin/sh

cur_PWD=${1}
IP=${2}
PORT=${3}
PYHTONFILE=${4}
cd ${cur_PWD}
echo '#### This is the python script job ####'

# Soruce Modules
source ./load_module.sh

# Start python script
python ${PYHTONFILE} --port ${PORT} --ip ${IP} > out.txt 2> errpy.txt
