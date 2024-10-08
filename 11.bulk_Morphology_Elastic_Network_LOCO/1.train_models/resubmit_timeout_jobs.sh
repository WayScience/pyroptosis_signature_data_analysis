#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --output=sample_resub-%j.out

jids_file="job_ids.txt"
touch $jids_file
timeout_jobs_file="list_timeout_jobs.txt"

# read all lines of the file to an array
while IFS= read -r line; do
    # Assuming the line contains job_id, cell_type, shuffle, feature_combination, cytokine
    echo $line
    job_id=$(echo "$line" | awk -F"'" '{print $2}')
    cell_type=$(echo "$line" | awk -F"'" '{print $4}')
    shuffle=$(echo "$line" | awk -F"'" '{print $6}')
    feature_combination=$(echo "$line" | awk -F"'" '{print $8}')
    cytokine=$(echo "$line" | awk -F"'" '{print $10}')
    echo "'$job_id' '$cell_type' '$shuffle' '$feature_combination' '$cytokine'"
    # get the number of jobs for the user
    number_of_jobs=$(squeue -u $USER | wc -l)
    while [ $number_of_jobs -gt 990 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u $USER | wc -l)
    done
    # resubmit the job
    new_jid=$(sbatch --parsable --time=2:00:00 train_regression_call_cheeky_child.sh "$cell_type" "$shuffle" "$feature_combination" "$cytokine")

    echo "'$new_jid' '$cell_type' '$shuffle' '$feature_combination' '$cytokine'" >> $jids_file
done < $timeout_jobs_file

echo "Jobs resubmitted"
