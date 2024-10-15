#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --output=sample_parent-%j.out

# 32 channel combination * 2 cell types * 2 shuffles * 187 cytokines = 23936

module load anaconda
conda init bash
conda activate Interstellar_python

# get the array of cytokiens
feature_combination_key_file="../0.split_data/results/feature_combinations_keys.txt"
# get the array of cytokines
filename="../0.split_data/cytokine_list/cytokine_list.txt"

# read all lines of the file to an array
readarray -t cytokine_array < $filename
# get all of the feature combinations
# read all lines of the file to an array
readarray -t feature_combination_keys < $feature_combination_key_file


jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb

shuffles=( True False )
cell_types=( SHSY5Y PBMC )

jobs_submitted_counter=0

# make a file to store the job ids
touch job_ids.txt

for cell_type in "${cell_types[@]}"
do
    for shuffle in "${shuffles[@]}"
    do
        for feature_combination in "${feature_combination_keys[@]}"
        do
            for cytokine in "${cytokine_array[@]}"
            do

                # get the number of jobs for the user
                number_of_jobs=$(squeue -u $USER | wc -l)
                while [ $number_of_jobs -gt 990 ]; do
                    sleep 1s
                    number_of_jobs=$(squeue -u $USER | wc -l)
                done

                job_id=$(sbatch train_regression_call_cheeky_child.sh "$cell_type" "$shuffle" "$feature_combination" "$cytokine")
		        echo "$cell_type $shuffle $feature_combination '${cytokine}'"
                # append the job id to the file
                job_id=$(echo $job_id | awk '{print $4}')
                echo " '$job_id' '$cell_type' '$shuffle' '$feature_combination' '$cytokine'" >> job_ids.txt
	        let jobs_submitted_counter++
	    done
        done
    done
done

echo "$jobs_submitted_counter"

echo "Array complete"

# end this job once reaching this point
exit 0
