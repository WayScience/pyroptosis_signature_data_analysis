#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --output=sample-%j.out



# This script is used to call all of the training jobs
# this is a cheeky implementation of a job queueing system
# the max number of jobs allowed for the HPC system that I am using is 1000
# so I will submit the first 750 jobs and then wait for the queue to lower to 248 before submitting the next 750 jobs
# 248 due to the fact that slurm indexes the jobs from 1-1000 thus 2 jobs are always running of one is pending

# get the array of channels
filename="../../0.split_data/cytokine_list/channel_splits.txt"
# read all lines of the file to an array
readarray -t channels < $filename

# set the maximum number of pending jobs allowed prior to submitting new jobs
max_pending_jobs=248

# loop through the data sets
for data_set in "${channels[@]}"; do
    echo "data_set: $data_set"

    while true; do
        # check the number of pending jobs in the queue
        num_pending_jobs=$(squeue -u $USER -t PENDING | wc -l)
        num_running_jobs=$(squeue -u $USER -t RUNNING | wc -l)
        num_active_jobs=$((num_pending_jobs + num_running_jobs))

        # submit a new job only if the number of pending jobs is less than the maximum allowed
        if [ $num_active_jobs -lt $max_pending_jobs ]; then
            sbatch train_regression_call_w_channel_splits.sh "$data_set"
            echo "Submitted new job to SLURM"
            break
        else
            sleep 120
        fi
    done
done
