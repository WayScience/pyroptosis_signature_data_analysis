#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=256G
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=4:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "SHSY5Y"
papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "PBMC"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
