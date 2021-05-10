#!/bin/sh

# Store passed port
cur_PWD=${1}
PORT=${2}

cd ${cur_PWD}
# Source Modules
source ./load_module.sh
 
# On torque sometimes we had problems that the job wouldn't know on which node it is currently
sleep 5 

# Find IP of host
HOSTNAME=$(hostname)
HOST_IP=$(host ${HOSTNAME} | awk '{ print $4 }')

# Start Redis server with IP and Port 
/p/home/jusers/reck1/juwels/miniconda3/bin/redis-server --port ${PORT} --protected-mode no

