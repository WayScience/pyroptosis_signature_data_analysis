#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --account=amc-general
#SBATCH --time=36:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-4%4

module load anaconda

conda activate Interstellar_python

cell_types=(  PBMC SHSY5Y )
model_names=( MultiClass_MLP )
shuffles=( True False )
jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
model_name_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#model_names[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
model_name=${model_names[$model_name_idx]}

command="python train_multiclass_model.py"

$command --CELL_TYPE "$cell_type" --MODEL_NAME "$model_name" --SHUFFLE "$shuffle"

echo "Done"
