#!/bin/bash

conda activate map

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

echo "All notebooks have been converted to scripts."

echo "Running all mAPs..."

python 0.generate_map_scores_morphology.py
python 1.aggregate_map_scores_morphology.py
python 2.generate_map_scores_secretome.py
python 3.aggregate_map_scores_secretome.py
python 4.generate_map_scores_morphology_treatment.py
python 5.aggregate_map_scores_morphology_treatment.py
python 6.generate_map_scores_secretome_treatment.py
python 7.aggregate_map_scores_secretome_treatment.py

echo "All mAPs have been generated and aggregated."
