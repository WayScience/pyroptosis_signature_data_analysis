#!/bin/bash

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

CELL_TYPE="PBMC"

filename="features/${CELL_TYPE}_feature_index.txt"
# read all lines of the file to an array
readarray -t features < $filename

echo "Cell type: $CELL_TYPE"

cd scripts/

# loop through the array
for feature in "${features[@]}"; do
    echo "Feature: $feature"
    time python 8.1_anova_all_groupings.py --feature $feature --cell_type $CELL_TYPE
done

cd ../
echo "All features submitted for ANOVA analysis"

