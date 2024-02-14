#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --qos=long
#SBATCH --time=168:00:00
#SBATCH --output=main_parent-%j.out

module load anaconda

conda activate Interstellar_python

# get the array of channels
filename="../../0.split_data/cytokine_list/channel_splits.txt"
# read all lines of the file to an array
readarray -t channels < $filename

# get the array of cytokines
filename="../../0.split_data/cytokine_list/cytokine_list1.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename
shuffles=( True False )
cell_types=( PBMC )

channels=( Gasdermin_ER_Mito_PM.parquet )

# set the maximum number of jobs to run at once
MAX_JOBS=500
progress_counter=1
# get the length of the arrays
cytokine_array_length=${#cytokine_array[@]}
channels_length=${#channels[@]}
shuffles_length=${#shuffles[@]}
cell_types_length=${#cell_types[@]}
# calculate the total number of jobs
total_jobs=$((cytokine_array_length*channels_length*shuffles_length*cell_types_length))
echo "total_jobs: $total_jobs"



for cell_type in "${cell_types[@]}"; do
	for channel in "${channels[@]}"; do
		for shuffle in "${shuffles[@]}"; do
			for cytokine in "${cytokine_array[@]}"; do
			while true; do
				NUM_SLURMS=$(squeue -u "$USER" | wc -l)
				if [ "$NUM_SLURMS" -lt $MAX_JOBS ]; then
					echo "cell_type: $cell_type shuffle: $shuffle cytokine: $cytokine data: $channel"
					sbatch child_enet_LOCO_test.sh "$cell_type" "$shuffle" "$cytokine" "$channel"
					progress_counter=$((progress_counter+1))
					# calculate the progress
					progress=$((progress_counter*100/total_jobs))
					echo "progress: $progress%"
					break
				else
					sleep 120
				fi
				done
			done
		done
	done
done
