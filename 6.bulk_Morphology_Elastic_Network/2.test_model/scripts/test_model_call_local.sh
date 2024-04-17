#!/bin/bash
# This script is used to train the regression models for the elastic network

 jupyter nbconvert --to python --output-dir=. ../notebooks/*.ipynb

conda activate Interstellar_python

# get the array of cytokiens
filename="../../0.split_data/cytokine_list/cytokine_list.txt"
# read all lines of the file to an array
readarray -t cytokine_array < $filename


shuffles=( "True" "False" )
cell_types=( PBMC )

# calculate total iterations
total_iterations=$((${#cytokine_array[@]} * ${#shuffles[@]} * ${#cell_types[@]}))


# loop through the array of cytokines
for cytokine in "${cytokine_array[@]}"; do
    for shuffle in "${shuffles[@]}"; do
        for cell_type in "${cell_types[@]}"; do
            echo "Cell type: $cell_type" "Shuffle: $shuffle" "Cytokine: $cytokine"
            python 1.test_regression_multi_output.py --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine"
        done
    done
done


echo "Done"
