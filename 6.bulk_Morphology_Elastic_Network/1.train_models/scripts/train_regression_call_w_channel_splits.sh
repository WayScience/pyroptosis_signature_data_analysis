#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --qos=normal
#SBATCH --time=10:00:00
#SBATCH --output=sample-%j.out
#SBATCH --array=1-750%50

module load anaconda

conda activate Interstellar

module load gnu_parallel

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

# pass through the data set from the call script
data="$1"

# get the array of cytokines
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename
shuffles=( True False )
cell_types=( SHSY5Y PBMC )

# calculate the number of jobs
job_id=$((SLURM_ARRAY_TASK_ID - 1))
shuffle_idx=$((job_id % ${#shuffles[@]}))
cell_type_idx=$(((job_id / ${#shuffles[@]}) % ${#cell_types[@]}))
cytokine_idx=$(((job_id / ${#shuffles[@]} / ${#cell_types[@]}) % ${#cytokine_array[@]}))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
cytokine=${cytokine_array[$cytokine_idx]}

echo "cell_type: $cell_type cytokine: $cytokine shuffle: $shuffle data: $data"

# command="python 2.train_regression_multi_output_channel_selection.py"
# $command --cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle" --data "$data"





