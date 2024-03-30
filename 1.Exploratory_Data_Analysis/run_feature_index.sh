#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amem
#SBATCH --mem=600G
#SBATCH --qos=mem
#SBATCH --output=sample-%j.out
#SBATCH --time=24:00:00

module load anaconda

conda activate Interstellar_python

echo "Running feature index for all cell types"

cd scripts/

python 8.0_create_feature_index.py --cell_type "PBMC"

python 8.0_create_feature_index.py --cell_type "SHSY5Y"

cd ../

echo "Complete"
