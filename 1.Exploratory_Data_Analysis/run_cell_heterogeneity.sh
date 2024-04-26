#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cd notebooks/

papermill 6.heterogeneity_analysis.ipynb 6.heterogeneity_analysis.ipynb -p cell_type "SHSY5Y"
papermill 6.heterogeneity_analysis.ipynb 6.heterogeneity_analysis.ipynb -p cell_type "PBMC"


jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
