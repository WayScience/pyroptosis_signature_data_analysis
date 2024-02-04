#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --partition=amem
#SBATCH --mem= 500G
#SBATCH --qos=mem
#SBATCH --output=sample-%j.out
#SBATCH --time=48:00:00

module load anaconda

conda activate Interstellar

cell_types=( PBMC )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))

cell_type=${cell_types[$cell_type_idx]}


echo "Cell type: $cell_type"

command="python 2.combine_regression_tests.py"

$command --cell_type "$cell_type"

echo "Done"
