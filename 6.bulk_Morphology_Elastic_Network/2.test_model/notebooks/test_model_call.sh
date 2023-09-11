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

cell_types=( SHSY5Y PBMC )
aggregations=( True False )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
aggregation_idx=$(((job_id / (${#shuffles[@]} * ${#cell_types[@]})) % ${#aggregations[@]}))

cell_type=${cell_types[$cell_type_idx]}
aggregation=${aggregations[$aggregation_idx]}

command1="papermill \
    2.visualize_regression_multi_ouput.ipynb \
    2.visualize_regression_multi_ouput.ipynb"

command2="Rscript ../scripts/2.visualize_regression_multi_ouput.r"

$command1 \
-p cell_type $cell_type \
        -p aggregation $aggregation \
        -p nomic True \
        -p flag True \
        -p shuffle $shuffle

$command2 \
    --aggregation $aggregation \
    --cell_type $cell_type
