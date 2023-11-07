#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=10:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-750%20

module load anaconda

conda activate Interstellar

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

# get the array of cytokines
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename
channel_filename="../../0.split_data/cytokine_list/channel_splits.txt"
readarray -t channel_array < $channel_filename
shuffles=( True False )
cell_types=( SHSY5Y PBMC )

run_parallel() {
    local job_id=$((SLURM_ARRAY_TASK_ID - 1))
    local shuffle_idx=$((job_id % ${#shuffles[@]}))
    local cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
    local cytokine_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#cytokine_array[@]}))
    local channel_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]} / ${#cytokine_array[@]}) % ${#channel_array[@]}))

    local shuffle="${shuffles[$shuffle_idx]}"
    local cell_type="${cell_types[$cell_type_idx]}"
    local cytokine="${cytokine_array[$cytokine_idx]}"

    local command="python 1.train_regression_multi_output.py"

    # use parallel to run the command for each channel
    parallel -j 4 $command --cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle" --data {} ::: "${channel_array[@]}"
}
# export the function so it can be used by parallel
export -f run_parallel

# run the parallel function for the current array index
run_parallel "$SLURM_ARRAY_TASK_ID"

# shuffles=( True False )
# cell_types=( SHSY5Y PBMC )
# # calculate the number of jobs
# # calculate the number of jobs
# job_id=$((SLURM_ARRAY_TASK_ID - 1))
# shuffle_idx=$((job_id % ${#shuffles[@]}))
# cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
# cytokine_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#cytokine_array[@]}))
# channel_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]} / ${#cytokine_array[@]}) % ${#channel_array[@]}))

# shuffle=${shuffles[$shuffle_idx]}
# cell_type=${cell_types[$cell_type_idx]}
# cytokine=${cytokine_array[$cytokine_idx]}
# channel=${channel_array[$channel_idx]}

# command="python 1.train_regression_multi_output.py"

# $command --cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle" --data "$channel"




