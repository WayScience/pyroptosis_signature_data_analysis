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

shuffled_datas=( True False )
cell_types=( SHSY5Y PBMC )
control_names=( "DMSO_0.100_%_DMSO_0.025_%" )
treatment_names=(
    "LPS_100.000_ug_per_ml_DMSO_0.025_%"
    "Thapsigargin_10.000_uM_DMSO_0.025_%"
    "LPS_10.000_ug_per_ml_DMSO_0.025_%"
    "LPS_1.000_ug_per_ml_DMSO_0.025_%"
    "LPS_0.100_ug_per_ml_DMSO_0.025_%"
    "LPS_0.010_ug_per_ml_DMSO_0.025_%"
    "Thapsigargin_1.000_uM_DMSO_0.025_%"
    )
model_names=(
    DMSO_0.025_vs_LPS_100
    DMSO_0.025_vs_Thapsigargin_10
    DMSO_0.025_vs_LPS_10
    DMSO_0.025_vs_LPS_1
    DMSO_0.025_vs_LPS_0.1
    DMSO_0.025_vs_LPS_0.01
    DMSO_0.025_vs_Thapsigargin_1
    )
shuffles=( True False )

for shuffled_data in "${shuffled_datas[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        for control_name in "${control_names[@]}"; do
             for treatment_name in "${!treatment_names[@]}"; do
                for shuffle in "${shuffles[@]}"; do
                    papermill \
                    binary_classification_testing.ipynb \
                    binary_classification_testing.ipynb \
                    -p SHUFFLE $shuffled_data \
                    -p CELL_TYPE $cell_type \
                    -p CONTROL_NAME $control_name \
                    -p TREATMENT_NAME ${treatment_names[$treatment_name]} \
                    -p MODEL_NAME ${model_names[$treatment_name]} \
                    -p SHUFFLE_DATA $shuffle
                done
            done
        done
    done
done




