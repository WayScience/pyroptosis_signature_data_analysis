#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-8%4

module load anaconda

conda activate Interstellar

shuffles=(True False)
cell_types=( SHSY5Y PBMC )
aggregation=( True False )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
aggregation_idx=$(((job_id / (${#shuffles[@]} * ${#cell_types[@]})) % ${#aggregation[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
aggregation=${aggregation[$aggregation_idx]}

command="papermill \
    1.train_regression_multi_ouput.ipynb \
    1.train_regression_multi_ouput.ipynb"

echo $cell_type $aggregation $shuffle

$command \
    -p cell_type $cell_type \
    -p aggregation $aggregation \
    -p shuffle $shuffle \
    -p aggregation True \
    -p nomic True \
    -p flag True \
