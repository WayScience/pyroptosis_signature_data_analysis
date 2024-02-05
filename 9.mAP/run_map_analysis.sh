#!/bin/bash

# init conda
conda init bash

# activate conda environment
conda activate map

# convert the .ipynb to .py
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts notebooks/*.ipynb

cd scripts || exit # change to the scripts directory but exit if it fails

# run the .py files (the map analysis)
python 0.generate_map_scores_class_level.py
python 1.aggregate_map_scores_class_level.py

# move back to the main directory
cd ../ || exit
