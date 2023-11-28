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
model_names=(
    DMSO_0.025_vs_LPS_100
    DMSO_0.025_vs_Thapsigargin_10
    DMSO_0.025_vs_LPS_10
    DMSO_0.025_vs_LPS_1
    DMSO_0.025_vs_LPS_0.1
    DMSO_0.025_vs_LPS_0.01
    DMSO_0.025_vs_Thapsigargin_1
    )
selected_treatment_comparisons=(
    "DMSO_0.100_%_DMSO_0.025_% vs LPS_100.000_ug_per_ml_DMSO_0.025_%,DMSO_0.100_%_DMSO_0.025_% vs Thapsigargin_1.000_uM_DMSO_0.025_%,DMSO_0.100_%_DMSO_0.025_% vs Thapsigargin_10.000_uM_DMSO_0.025_%")

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb*


for cell_type in "${cell_types[@]}"; do
    for model_name in "${model_names[@]}"; do
        for selected_treatment_comparison in "${selected_treatment_comparisons[@]}"; do
            Rscript \
            binary_classification_testing_visualization.r \
            --celltype $cell_type \
            --model_name $model_name \
            --selected_treatment_comparisons "${selected_treatment_comparison}"
        done
    done
done
