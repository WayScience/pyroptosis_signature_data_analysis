#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-8%8

module load anaconda

conda activate Interstellar

cell_types=( SHSY5Y PBMC )
model_names=( MultiClass_MLP_h202_remove MultiClass_MLP )
shuffles=( True False )
jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb

shuffles=( True False )
cell_types=( SHSY5Y PBMC )
model_names=( MultiClass_MLP_h202_remove MultiClass_MLP )
# calculate the number of jobs
# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
model_name_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#model_names[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
model_name=${model_names[$model_name_idx]}

command="python train_multiclass_model.py"

$command --cell_type "$cell_type" --model_name "$model_name" --shuffle "$shuffle"


echo "Done"
