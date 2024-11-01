#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=1:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_python

cell_type=$1

cd scripts/ || exit

echo $cell_type
python 0.split_data_regression.py --cell_type "$cell_type"
python 1.get_cytokine_list.py
python 2.LOCO_data_split.py --cell_type "$cell_type"

cd ../

echo "Complete"
