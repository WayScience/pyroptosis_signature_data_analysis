#!/bin/bash

# This script is used to call the training and testing scripts

conda activate Interstellar_python

cell_type=$1
model_name=$2
channel_combination_key=$3
shuffle=$4
convert=False

if [ "$convert" = True ]; then
    jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
fi

cd scripts/ || exit

python 0.train_multiclass_model.py --cell_type "$cell_type" --model_name "$model_name" --channel_combination_key "$channel_combination_key" --shuffle "$shuffle"
cd .. || exit
echo "Training completed"
