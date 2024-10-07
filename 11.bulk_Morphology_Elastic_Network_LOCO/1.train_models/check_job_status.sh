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
completed_jobs_file="completed_jobs.txt"
failed_jobs_file="failed_jobs.txt"
timeout_jobs_file="timeout_jobs.txt"
other_jobs_file="other_jobs.txt"
touch $completed_jobs_file $failed_jobs_file $timeout_jobs_file $other_jobs_file


# read all lines of the file to an array
readarray -t job_ids < $jids_file

# counters
total_counter=0
completed_counter=0
failed_counter=0
timeout_counter=0
other_counter=0

for job_id in "${job_ids[@]}"; do
    total_counter=$((total_counter+1))
    status=$(sacct -j "$job_id" --format=State --noheader | awk '{print $1}')

    # Display the status of the job
    echo "$status"
    # Check if the job has completed successfully or failed
    if [[ "$status" == "COMPLETED" ]]; then
    completed_counter=$((completed_counter+1))
    echo "$job_id" >> $completed_jobs_file
    elif [[ "$status" == "FAILED" ]]; then
    failed_counter=$((failed_counter+1))
    echo "$job_id" >> $failed_jobs_file
    elif [[ "$status" == "TIMEOUT" ]]; then
    timeout_counter=$((timeout_counter+1))
    echo "$job_id" >> $timeout_jobs_file
    else
    other_counter=$((other_counter+1))
    echo "$job_id" >> $other_jobs_file
    fi
done

echo "Total jobs: $total_counter"
echo "Completed jobs: $completed_counter"
echo "Failed jobs: $failed_counter"
echo "Timeout jobs: $timeout_counter"
echo "Other jobs: $other_counter"

echo "Completed job status check"
