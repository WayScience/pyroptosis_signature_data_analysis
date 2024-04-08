#!/bin/bash
# This script is used to train the regression models for the elastic network

conda activate Interstellar

cell_types=( PBMC )

for cell_type in "${cell_types[@]}"; do
    echo "Cell type: $cell_type"
    time python 2.combine_regression_tests.py --cell_type "$cell_type"
done

echo "Done"
