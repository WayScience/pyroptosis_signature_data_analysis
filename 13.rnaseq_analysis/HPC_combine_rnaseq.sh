#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=400G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb

module load anaconda

conda activate rnaseq_r_env
cd scripts || exit

Rscript search_whole_tissue.r

cd ../ || exit

conda deactivate
