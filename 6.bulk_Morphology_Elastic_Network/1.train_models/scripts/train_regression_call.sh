#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=48:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-750%20

module load anaconda

conda activate Interstellar

# get the array of cytokiens
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

shuffles=( True False )
cell_types=( SHSY5Y PBMC )
# calculate the number of jobs
# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
cytokine_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#cytokine_array[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
cytokine=${cytokine_array[$cytokine_idx]}

command="python 1.train_regression_multi_output.py"

$command --cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle"
