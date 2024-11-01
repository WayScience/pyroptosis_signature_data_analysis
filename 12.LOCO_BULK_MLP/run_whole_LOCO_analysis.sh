#!/bin/bash

# This script is used to call the training and testing scripts

conda activate Interstellar_python

convert=True

if [ "$convert" = True ]; then

cd 0.data_splits || exit
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
cd ../1.hyperparameter_optimization || exit
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
cd ../2.train_models || exit
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
cd .. || exit

fi

cell_types=( "SHSY5Y" "PBMC" )
shuffles=( "True" "False" )
model_name="MLP"
# read the feature combination key file
feature_combination_key_file="0.data_splits/results/feature_combinations_keys.txt"
readarray -t feature_combinations_keys < "$feature_combination_key_file"

for cell_type in "${cell_types[@]}"; do
    cd 0.data_splits || exit
    source call_data_split.sh "$cell_type" "$model_name"
    cd .. || exit
    for channel_combination_key in "${feature_combinations_keys[@]}"; do
        cd 1.hyperparameter_optimization || exit
        source call_hyperparameter_optimization.sh "$cell_type" "$model_name" "$channel_combination_key"
        cd .. || exit
        for shuffle in "${shuffles[@]}"; do
            cd 2.train_models || exit
            source call_training.sh "$cell_type" "$model_name" "$channel_combination_key" "$shuffle"
            cd .. || exit
        done
    done
done

echo "Training completed"

