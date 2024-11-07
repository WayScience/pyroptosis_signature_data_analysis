#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=10:00
#SBATCH --output=sample_test-%j.out

# 32 channel combination * 2 cell types * 2 shuffles * 187 cytokines = 23936
module load anaconda
conda activate Interstellar_python

cell_type=$1
shuffle=$2
feature_combinations_key=$3
cytokine=$4
data_split=$5

echo "$cell_type $shuffle $feature_combinations_key $cytokine"

cd scripts/ || exit

python 1.test_regression_multi_output.py --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine" --feature_combinations_key "$feature_combinations_key" --data_split "$data_split"

cd ../

echo "Complete"

