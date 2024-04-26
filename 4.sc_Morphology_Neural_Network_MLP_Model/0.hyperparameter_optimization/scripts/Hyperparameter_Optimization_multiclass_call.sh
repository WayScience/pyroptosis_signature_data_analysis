#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --account=amc-general
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-2%2

module load anaconda

conda activate Interstellar_python

cell_types=( PBMC )
model_names=( MultiClass_MLP )

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

job_id=$((SLURM_ARRAY_TASK_ID - 1))
model_name_idx=$((job_id % ${#model_names[@]}))
cell_type_idx=$(((job_id / ${#model_names[@]}) % ${#cell_types[@]}))

cell_type=${cell_types[$cell_type_idx]}
model_name=${model_names[$model_name_idx]}

command="python Hyperparameter_Optimization_model_multiclass.py"

$command --cell_type "$cell_type" --model_name "$model_name"

echo "Done"
