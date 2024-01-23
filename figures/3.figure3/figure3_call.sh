#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_R

jupyter nbconvert --to=script --FilesWriter.build_directory=. *.ipynb

Rscript figure3.r
