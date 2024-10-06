#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --mem=250M
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=1:30:00
#SBATCH --output=sample_train-%j.out
#SBATCH --array=1-187%25

# 1 cell type * 1 shuffle type * 187 cytokines = 187 per child script
# run the array at 50 tasks per node at one node for 15 minutes

module load anaconda
conda init bash
conda activate Interstellar_python

# get the args
cell_type=$1
shuffle=$2
feature_combinations_key=$3
feature_combinations_file=$4

# get the array of cytokines
filename="../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename

jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb

# calculate the number of jobs
# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cytokine_idx=$(((job_id % ${#cytokine_array[@]})))

cytokine=${cytokine_array[$cytokine_idx]}

cd scripts/

command="python 1.train_regression_multi_output.py"

echo "$cell_type $shuffle $cytokine"

$command --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine" --feature_combinations_key "$feature_combinations_key"

cd ../

echo "Complete"
