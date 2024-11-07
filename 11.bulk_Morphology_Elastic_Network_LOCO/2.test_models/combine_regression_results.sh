#!/bin/bash
# This script is used to combine the results of the regression models

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb
cd scripts || exit

cell_types=( SHSY5Y PBMC )

for cell_type in "${cell_types[@]}"
do
    echo "$cell_type"
    python 2.combine_regression_tests.py --cell_type "$cell_type"
done

cd .. || exit

echo "All model results have been combined"


