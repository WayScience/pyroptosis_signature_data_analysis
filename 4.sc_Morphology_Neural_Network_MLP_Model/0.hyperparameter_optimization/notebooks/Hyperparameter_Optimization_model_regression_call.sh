#!/bin/bash

# this file is used to run the Hyperparameter_Optimization_model_regression.ipynb file

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda
conda activate Interstellar

cell_types=( "SHSY5Y" "PBMC" )
control_names=( "DMSO_0.100_DMSO_0.025" )
treatment_names=(
    "LPS_100.000_DMSO_0.025"
    "Thapsigargin_10.000_DMSO_0.025"
    "LPS_10.000_DMSO_0.025"
    "LPS_1.000_DMSO_0.025"
    "LPS_0.100_DMSO_0.025"
    "LPS_0.010_DMSO_0.025"
    "Thapsigargin_1.000_DMSO_0.025" )
model_names=(
    "DMSO_0.025_vs_LPS_100"
    "DMSO_0.025_vs_Thapsigargin_10"
    "DMSO_0.025_vs_LPS_10"
    "DMSO_0.025_vs_LPS_1"
    "DMSO_0.025_vs_LPS_0.1"
    "DMSO_0.025_vs_LPS_0.01"
    "DMSO_0.025_vs_Thapsigargin_1" )

# for loop to run all the combinations of cell types, control_names, treatment_names, and model_names
for cell_type in "${cell_types[@]}"; do
    for control_name in "${control_names[@]}"; do
        for treatment_name in "${treatment_names[@]}"; do
            for model_name in "${model_names[@]}"; do
                echo "$cell_type $control_name" "${treatment_name}" "${model_name}"
                papermill \
                Hyperparameter_Optimization_model_regression.ipynb \
                Hyperparameter_Optimization_model_regression.ipynb \
                -p CELL_TYPE "$cell_type" \
                -p CONTROL_NAME "$control_name" \
                -p TREATMENT_NAME "$treatment_name" \
                -p MODEL_NAME "$model_name"
            done
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts ./*.ipynb


