#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1

module load anaconda

conda activate Interstellar

# define the search parameters
cell_types=( "PBMC" )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))

cell_type=${cell_types[$cell_type_idx]}

# use papermill to run the notebooks with injected parameters
command="python 0.split_data_regression.py"

echo $cell_type

$command --cell_type "$cell_type"
