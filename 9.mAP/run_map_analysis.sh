#!/bin/bash

# init conda
conda init bash

# activate conda environment
conda activate map

# convert the .ipynb to .py
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts notebooks/*.ipynb

cd scripts || exit # change to the scripts directory but exit if it fails

# run the .py files (the map analysis)

python 0.generate_map_scores_morphology.py
python 1.aggregate_map_scores_morphology.py
python 2.generate_map_scores_secretome.py
python 3.aggregate_map_scores_secretome.py
python 4.generate_map_scores_morphology_treatment.py
python 5.aggregate_map_scores_morphology_treatment.py
python 6.generate_map_scores_secretome_treatment.py
python 7.aggregate_map_scores_secretome_treatment.py

# move back to the main directory
cd ../ || exit
