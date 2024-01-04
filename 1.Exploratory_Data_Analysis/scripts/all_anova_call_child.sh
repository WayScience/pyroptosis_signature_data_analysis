#!/bin/bash
"""
This is a script to run the anova analysis on all the cell types and features.
This is the child process that gets spun up by the parent process all_anova_call_parent.sh.
I needed to do this to get around the limit of 1000 jobs that can be submitted to the cluster.

One would call this a "cheeky" solution.
"""
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
