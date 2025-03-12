#!/bin/bash

# init conda
conda init bash

# activate conda environment
conda activate map

# convert the .ipynb to .py
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts notebooks/*.ipynb

cd scripts || exit # change to the scripts directory but exit if it fails

# run the .py files (the map analysis)

python 0.generate_map_scores_morphology.py --shuffle
python 0.generate_map_scores_morphology.py
python 1.generate_map_scores_secretome.py --shuffle
python 1.generate_map_scores_secretome.py

conda deactivate

conda activate Interstellar_R

Rscript 2.visualize_map_scores.r

conda deactivate

# move back to the main directory
cd ../ || exit
