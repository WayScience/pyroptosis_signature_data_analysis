#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=10:00
#SBATCH --output=sample_parent-%j.out
#SBATCH --array=1-2%2

# 32 channel combination * 2 cell types * 2 shuffles = 128
# run the array at 16 tasks per node at one node for 72 hours

module load anaconda

conda init bash

conda activate Interstellar_python

# get the array of cytokiens
feature_combination_file="../../0.split_data/results/feature_combinations.toml"

# get all of the feature combinations
# read all lines of the file to an array
readarray -t feature_combination_keys < $feature_combination_key_file


jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb

shuffles=( True False )
cell_types=( SHSY5Y PBMC )
# calculate the number of jobs
# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
feature_combination_keys_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#feature_combination_keys[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
feature_combination=${feature_combination_keys[$feature_combination_keys_idx]}

sbatch train_regression_call.sh $cell_type $shuffle $feature_combination
