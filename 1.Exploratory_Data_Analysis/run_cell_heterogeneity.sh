#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

papermill 6.heterogeneity_analnysis.ipynb 6.heterogeneity_analnysis.ipynb -p celltype "SHSY5Y"
papermill 6.heterogeneity_analnysis.ipynb 6.heterogeneity_analnysis.ipynb -p celltype "PBMC"


jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
