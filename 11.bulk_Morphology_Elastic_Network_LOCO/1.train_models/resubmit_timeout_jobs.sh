#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --output=sample_resub-%j.out

jids_file="job_ids.txt"

timeout_jobs_file="list_timeout_jobs.txt"

# read all lines of the file to an array
while IFS= read -r line; do
    # Assuming the line contains job_id, cell_type, shuffle, feature_combination, cytokine
    job_id=$(echo "$line" | awk '{print $1}')
    cell_type=$(echo "$line" | awk '{print $2}')
    shuffle=$(echo "$line" | awk '{print $3}')
    feature_combination=$(echo "$line" | awk '{print $4}')
    cytokine=$(echo "$line" | awk '{print $5}')
    # get the number of jobs for the user
    number_of_jobs=$(squeue -u $USER | wc -l)
    while [ $number_of_jobs -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u $USER | wc -l)
    done
    # resubmit the job
    new_jid=$(sbatch --parsable --time=2:00:00 train_regression_call_cheeky_child.sh $cell_type $shuffle $feature_combination $cytokine)
    echo "'$new_jid' '$cell_type' '$shuffle' '$feature_combination' '$cytokine'" >> $timeout_jobs_file
done < $jids_file

# run the check_job_status.sh script
sbatch check_job_status.sh

echo "Jobs resubmitted"
