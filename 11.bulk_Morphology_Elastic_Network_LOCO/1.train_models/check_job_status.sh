#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=5:00
#SBATCH --output=job_status_check-%j.out

# get the file with the job ids
jids_file="job_ids.txt"

# define each job status file to be created
completed_jobs_file="list_completed_jobs.txt"
failed_jobs_file="list_failed_jobs.txt"
timeout_jobs_file="lsit_timeout_jobs.txt"
other_jobs_file="list_other_jobs.txt"


# check if the files exist, delete them if they do
# create if they do not exist
if [ -f "$completed_jobs_file" ]; then
    rm $completed_jobs_file
fi
if [ -f "$failed_jobs_file" ]; then
    rm $failed_jobs_file
fi
if [ -f "$timeout_jobs_file" ]; then
    rm $timeout_jobs_file
fi
if [ -f "$other_jobs_file" ]; then
    rm $other_jobs_file
fi

[[ ! -f "$completed_jobs_file" ]] && touch $completed_jobs_file
[[ ! -f "$failed_jobs_file" ]] && touch $failed_jobs_file
[[ ! -f "$timeout_jobs_file" ]] && touch $timeout_jobs_file
[[ ! -f "$other_jobs_file" ]] && touch $other_jobs_file



# counters
total_counter=0
completed_counter=0
failed_counter=0
timeout_counter=0
other_counter=0

# read all lines of the file to an array
readarray -t job_ids < $jids_file

while IFS= read -r line; do
    # Assuming the line contains job_id, cell_type, shuffle, feature_combination, cytokine
    job_id=$(echo "$line" | awk '{print $1}')
    cell_type=$(echo "$line" | awk '{print $2}')
    shuffle=$(echo "$line" | awk '{print $3}')
    feature_combination=$(echo "$line" | awk '{print $4}')
    cytokine=$(echo "$line" | awk '{print $5}')

    total_counter=$((total_counter+1))
    status=$(sacct -j "$job_id" --format=State --noheader | awk 'NR==1{print $1}')
    echo "Job ID: $job_id"
    echo "Status: $status"
    # Display the status of the job
    # Check if the job has completed successfully or failed
    if [[ "$status" == "COMPLETED" ]]; then
        completed_counter=$((completed_counter+1))
        echo "$job_id" >> "$completed_jobs_file"
    elif [[ "$status" == "FAILED" ]]; then
        failed_counter=$((failed_counter+1))
        echo "$job_id" >> "$failed_jobs_file"
    elif [[ "$status" == "TIMEOUT" ]]; then
        timeout_counter=$((timeout_counter+1))
        echo "$job_id" >> "$timeout_jobs_file"
    else
        other_counter=$((other_counter+1))
        echo "$job_id" >> "$other_jobs_file"
    fi
done < "$jids_file"

echo "Total jobs: $total_counter"
echo "Completed jobs: $completed_counter"
echo "Failed jobs: $failed_counter"
echo "Timeout jobs: $timeout_counter"
echo "Other jobs: $other_counter"

echo "Completed job status check"
