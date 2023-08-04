#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=aa100
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

shuffled_datas=( True False )
cell_types=( SHSY5Y )
control_names=( DMSO_0.100_DMSO_0.025 )
treatment_names=( LPS_100.000_DMSO_0.025 Thapsigargin_10.000_DMSO_0.025 )
model_names=( DMSO_0.025_vs_LPS_100 DMSO_0.025_vs_Thapsigargin_10 )

for shuffled_data in "${shuffled_datas[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        for control_name in "${control_names[@]}"; do
            for treatment_name in "${treatment_names[@]}"; do
                for model_name in "${model_names[@]}"; do
                    papermill \
                    Hyperparameter_Optimization_model_binary.ipynb \
                    Hyperparameter_Optimization_model_binary.ipynb \
                    -p SHUFFLE_DATA $shuffled_data \
                    -p CELL_TYPE $cell_type \
                    -p CONTROL_NAME $control_name \
                    -p TREATMENT_NAME $treatment_name \
                    -p MODEL_NAME $model_name
                done
            done
        done
    done
done

