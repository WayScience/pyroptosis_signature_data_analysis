#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --time=72:00:00
#SBATCH --output=sample_resub-%j.out


module load anaconda
conda init bash
conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb

cd models/
ls regression/PBMC/aggregated_with_nomic/ > PBMC.txt
ls regression/SHSY5Y/aggregated_with_nomic/ > SHSY5Y.txt


cd ../scripts
python check_which_models_need_to_run.py
cd ../

manual_jobs_file="manual_jobs.txt"
jids_file="job_ids.txt"

# read all lines of the file to an array
while IFS= read -r line; do
    # Assuming the line contains job_id, cell_type, shuffle, feature_combination, cytokine
    echo $line
    cell_type=$(echo "$line" | awk -F"'" '{print $2}')
    shuffle=$(echo "$line" | awk -F"'" '{print $4}')
    feature_combination=$(echo "$line" | awk -F"'" '{print $6}')
    cytokine=$(echo "$line" | awk -F"'" '{print $8}')
    echo " '$cell_type' '$shuffle' '$feature_combination' '$cytokine'"
    # get the number of jobs for the user
    number_of_jobs=$(squeue -u $USER | wc -l)
    while [ $number_of_jobs -gt 3 ]; do
        sleep 1s
        number_of_jobs=$(squeue -u $USER | wc -l)
    done
    resubmit the job
    new_jid=$(sbatch train_regression_call_cheeky_child.sh "$cell_type" "$shuffle" "$feature_combination" "$cytokine")
    new_jid=$(echo $new_jid | awk '{print $4}')
    echo "'$new_jid' '$cell_type' '$shuffle' '$feature_combination' '$cytokine'" >> $jids_file
done < $manual_jobs_file


echo "Jobs resubmitted"
