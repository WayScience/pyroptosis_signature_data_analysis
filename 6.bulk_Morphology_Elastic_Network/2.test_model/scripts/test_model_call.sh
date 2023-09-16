#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-4%4

module load anaconda

conda activate Interstellar

shuffles=( True False )
cell_types=( SHSY5Y PBMC )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))
shuffle_idx=$((job_id / ${#cell_types[@]}))

cell_type=${cell_types[$cell_type_idx]}
shuffle=${shuffles[$shuffle_idx]}

echo "Cell type: $cell_type" "Shuffle: $shuffle"

command="python 1.test_regression_multi_output.py"

$command --cell_type "$cell_type" --shuffle "$shuffle"

echo "Done"
