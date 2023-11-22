#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amem
#SBATCH --mem=500G
#SBATCH --qos=mem
#SBATCH --output=sample-%j.out
#SBATCH --time=24:00:00

module load anaconda

conda activate Interstellar_python

feature=$1
CELL_TYPE=$2

echo "Feature: $feature" "Cell type: $CELL_TYPE"

python 8.1_anova_all_groupings.py --feature $feature --cell_type $CELL_TYPE

echo "Complete"
