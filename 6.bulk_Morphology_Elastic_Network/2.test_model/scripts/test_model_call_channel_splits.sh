#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --output=sample-%j.out
#SBATCH --array=1-750%100
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=200G

module load anaconda

conda activate Interstellar

# get the array of cytokiens
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename
channel_filename="../../0.split_data/cytokine_list/channel_splits.txt"
readarray -t channel_array < $channel_filename

shuffles=( "True" "False" )
cell_types=( SHSY5Y PBMC )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
cell_type_idx=$((job_id % ${#cell_types[@]}))
shuffle_idx=$(((job_id / ${#cell_types[@]}) % ${#shuffles[@]}))
cytokine_idx=$(((job_id / ${#cell_types[@]} / ${#shuffles[@]}) % ${#cytokine_array[@]}))
channel_idx=$(((job_id / ${#cell_types[@]} / ${#shuffles[@]} / ${#cytokine_array[@]}) % ${#channel_array[@]}))

cell_type=${cell_types[$cell_type_idx]}
shuffle=${shuffles[$shuffle_idx]}
cytokine=${cytokine_array[$cytokine_idx]}
channel=${channel_array[$channel_idx]}

echo "Cell type: $cell_type" "Shuffle: $shuffle" "Cytokine: $cytokine" "Channel: $channel"

command="python 1.test_regression_multi_output_channel_split.py"

$command --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine" --data "$channel"

echo "Done"
