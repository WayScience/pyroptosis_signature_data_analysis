#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cell_types=( SHSY5Y PBMC )
model_names=( MultiClass_MLP_h202_remove MultiClass_MLP )
jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb*


for cell_type in "${cell_types[@]}"; do
    for model_name in "${model_names[@]}"; do
        Rscript \
        multi_classification_visualization.r \
        --cell_type $cell_type \
        --model_name $model_name
    done
done

echo "Done"
