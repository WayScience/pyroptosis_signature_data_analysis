#!/bin/bash
"""
This script combines the intermediate result files from
each of the child processes into one file.
The result is one file.
"""
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out


module load anaconda

conda activate Interstellar_python

python 8.2_combine_files.py

echo "Complete"
