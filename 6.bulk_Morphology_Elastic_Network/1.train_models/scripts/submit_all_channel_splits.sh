#!/bin/bash
# This script is used to call all of the training jobs

# get the array of channels
filename="../../0.split_data/cytokine_list/channel_splits.txt"
# read all lines of the file to an array
readarray -t channels < $filename

for data_set in "${channels[@]}"; do
    echo "data_set: $data_set"
    sbatch train_regression_call_w_channel_splits.sh "$data_set"
done
