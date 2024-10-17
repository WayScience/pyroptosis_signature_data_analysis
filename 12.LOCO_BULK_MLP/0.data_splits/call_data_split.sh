#!/bin/bash

# This script is used to split the data for machine learning

conda activate Interstellar_python

cell_type=$1
model_name=$2
convert=False

if [ "$convert" = True ]; then
    jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb
fi

cd scripts/ || exit

python 0.LOCO_data_split.py --cell_type "$cell_type"

python 1.data_splits_model_multiclass.py --cell_type "$cell_type" --model_name "$model_name"

cd .. || exit

echo
