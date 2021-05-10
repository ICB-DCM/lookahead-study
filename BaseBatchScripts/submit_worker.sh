#!/bin/sh





HOST_IP=${1}
PORT=${2}
TIME=${3}
CPUSPERTASK=${4}
cur_PWD=${5}

cd ${cur_PWD}

echo '#### This is the worker script job ####'

# Source module
source ./load_module.sh

# Define associative array with port and password
declare -A PW_DICT

# Start redis-worker
/home/ealamoodi/.local/bin/abc-redis-worker --host=${HOST_IP} --port ${PORT} --runtime ${TIME:0:2}h --processes ${CPUSPERTASK}
