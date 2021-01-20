#!/bin/sh

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

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
/p/home/jusers/reck1/juwels/keksenv/bin/abc-redis-worker --host=${HOST_IP} --port ${PORT} --runtime ${TIME:0:2}h --processes ${CPUSPERTASK} --daemon=False
