#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --mem=250M
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=1:30:00
#SBATCH --output=sample_train-%j.out

# 1 cell type * 1 shuffle type * 187 cytokines = 187 per child script
# run the array at 50 tasks per node at one node for 15 minutes

module load anaconda
conda init bash
conda activate Interstellar_python

# get the args
cell_type=$1
shuffle=$2
feature_combinations_key=$3
cytokine=$4

echo "$cell_type $shuffle $feature_combination_key $cytokine"

jupyter nbconvert --to=script --FilesWriter.build_directory=./scripts/ ./notebooks/*.ipynb

cd scripts/

#python 1.train_regression_multi_output.py --cell_type "$cell_type" --shuffle "$shuffle" --cytokine "$cytokine" --feature_combinations_key "$feature_combinations_key"

cd ../

echo "Complete"
