#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --output=sample-%j.out
#SBATCH --array=1-2%2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=200G

module load anaconda

conda activate Interstellar

cell_types=( SHSY5Y PBMC )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))

cell_type=${cell_types[$cell_type_idx]}


echo "Cell type: $cell_type"

command="python 2.combine_regression_tests.py"

$command --cell_type "$cell_type"

echo "Done"
