#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --output=sample-%j.out
#SBATCH --array=1-750%150
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=200G

module load anaconda

conda activate Interstellar

# get the array of cytokiens
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename


shuffles=( "True" "False" )
cell_types=( PBMC SHSY5Y )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))
shuffle_idx=$(((job_id / ${#cell_types[@]}) % ${#shuffles[@]}))
cytokine_idx=$(((job_id / ${#cell_types[@]} / ${#shuffles[@]}) % ${#cytokine_array[@]}))

cell_type=${cell_types[$cell_type_idx]}
shuffle=${shuffles[$shuffle_idx]}
cytokine=${cytokine_array[$cytokine_idx]}

echo "Cell type: $cell_type" "Shuffle: $shuffle" "Cytokine: $cytokine"

command="python 1.test_regression_multi_output.py"

$command --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine"

echo "Done"
