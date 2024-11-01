#!/bin/bash

# This script is used to call the hyperparameter optimization

conda activate Interstellar_python

cell_type=$1
model_name=$2
channel_combination_key=$3
convert=False

if [ "$convert" = True ]; then
    jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
fi

cd scripts/ || exit

python 0.Hyperparameter_Optimization_model_multiclass.py --cell_type "$cell_type" --model_name "$model_name" --channel_combination_key "$channel_combination_key"

cd .. || exit

echo "Hyperparameter optimization completed"
