#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

shuffles=(True False)
cell_types=( SHSY5Y PBMC )

for shuffle in "${shuffles[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        echo $cell_type $shuffle
        python 1.test_regression_multi_output.py--cell_type "$cell_type" --shuffle "$shuffle"
    done
done
echo "Done"
