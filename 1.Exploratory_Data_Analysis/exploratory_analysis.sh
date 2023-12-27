#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_python

papermill notebooks/4.cell_count_analysis.ipynb noteboooks/4.cell_count_analysis.ipynb -p celltype "SHSY5Y"
papermill notebooks/1.umap_analysis_plate2.ipynb notebooks/1.umap_analysis_plate2.ipynb -p celltype "SHSY5Y"

#papermill notebooks/4.cell_count_analysis.ipynb notebooks/4.cell_count_analysis.ipynb -p celltype "PBMC"
#papermill notebooks/1.umap_analysis_plate2.ipynb notebooks/1.umap_analysis_plate2.ipynb -p celltype "PBMC"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts notebooks/*.ipynb

