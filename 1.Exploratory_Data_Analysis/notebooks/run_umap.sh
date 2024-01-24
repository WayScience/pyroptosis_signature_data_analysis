#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=48:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cell_types=( PBMC )
samples=( True False )
for cell_type in ${cell_types[@]}; do

    echo "cell_type: $cell_type"
    papermill 9.subset_umap.ipynb 9.subset_umap.ipynb -p cell_type "$cell_type"
    for sample in ${samples[@]}; do
        echo "sample: $sample"
        papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p cell_type "$cell_type" -p sample "$sample"
    done
done


jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
