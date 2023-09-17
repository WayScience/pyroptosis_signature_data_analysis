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
control_names=(
    DMSO_0.100_DMSO_0.025
    Flagellin_0.100_DMSO_0.025
    Flagellin_1.000_DMSO_0.025
    LPS_Nigericin_1.000_1.0_DMSO_0.025
    LPS_Nigericin_1.000_3.0_DMSO_0.025
    LPS_Nigericin_1.000_10.0_DMSO_0.025
    LPS_100.000_DMSO_0.025
    LPS_10.000_DMSO_0.025
    LPS_1.000_DMSO_0.025
    LPS_0.100_DMSO_0.025
    LPS_0.010_DMSO_0.025
    Thapsigargin_10.000_DMSO_0.025
    Thapsigargin_1.000_DMSO_0.025
    H2O2_100.000_DMSO_0.025 )
treatment_names=(
    DMSO_0.100_DMSO_0.025
    Flagellin_0.100_DMSO_0.025
    Flagellin_1.000_DMSO_0.025
    LPS_Nigericin_1.000_1.0_DMSO_0.025
    LPS_Nigericin_1.000_3.0_DMSO_0.025
    LPS_Nigericin_1.000_10.0_DMSO_0.025
    LPS_100.000_DMSO_0.025
    LPS_10.000_DMSO_0.025
    LPS_1.000_DMSO_0.025
    LPS_0.100_DMSO_0.025
    LPS_0.010_DMSO_0.025
    Thapsigargin_10.000_DMSO_0.025
    Thapsigargin_1.000_DMSO_0.025
    H2O2_100.000_DMSO_0.025 )

# generate model names from control_names and treatment_names
model_names=()
for control_name in "${!control_names[@]}"; do
    for treatment_name in "${!treatment_names[@]}"; do
        # if treatment_name is not equal to control_name
        if [[ $treatment_name != "$control_name" ]]; then
            model_names+=("${control_names[$control_name]}_vs_${treatment_names[$treatment_name]}")
        fi
    done
done

# save notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb*

# loop through all cell types, model names, and selected_treatment_comparisons
for cell_type in "${cell_types[@]}"; do
    for model_name in "${model_names[@]}"; do
        Rscript \
        binary_classification_training_visualization.r \
        --celltype "$cell_type" \
        --model_name "$model_name"
    done
done
