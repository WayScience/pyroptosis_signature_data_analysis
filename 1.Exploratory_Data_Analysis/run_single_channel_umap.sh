#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cell_types=( SHSY5Y PBMC )

for cell_type in ${cell_types[@]}; do

    echo "cell_type: $cell_type"
    papermill 2.single_channel_umap_analysis.ipynb 2.single_channel_umap_analysis.ipynb -p cell_type "$cell_type"

done

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
