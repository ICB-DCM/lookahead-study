# Batch_pyABC

Batch scripts to allow dynamically allocate master/worker nodes for the redis implementation in pyABC

## Start a parallel pyABC run:
run `submit_job.sh` with the follwing parameters:
* port number
* Number of nodes
* queue name
* job time 
* CPUSPERTASK
* python script
## To kill all running jobs:
run `kill_all.sh`
