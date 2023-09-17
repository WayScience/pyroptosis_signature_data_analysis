#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --output=sample-%j.out

# module load anaconda
#
# conda activate Interstellar

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

# get the number of combination of control_names treatment_names and cell_types
# this is the number of jobs
num_jobs=$(( ${#control_names[@]} * ${#treatment_names[@]} * ${#cell_types[@]} ))
echo "num_jobs: $num_jobs"
job=1

for cell_type in "${cell_types[@]}"; do
    for control_name in "${!control_names[@]}"; do
        for treatment_name in "${!treatment_names[@]}"; do
            echo "cell_type: $cell_type" control_name: "${control_names[$control_name]}" treatment_name: "${treatment_names[$treatment_name]}"
            job=$(( $job + 1 ))
            progress=$(( $job * 100 / $num_jobs ))
            echo "progress: $progress"
            # if treatment_name is not equal to control_name
            if [[ $treatment_name != "$control_name" ]]; then

                papermill \
                Hyperparameter_Optimization_model_binary.ipynb \
                Hyperparameter_Optimization_model_binary.ipynb \
                -p CELL_TYPE "$cell_type" \
                -p CONTROL_NAME "${control_names[$control_name]}" \
                -p TREATMENT_NAME "${treatment_names[$treatment_name]}"
            fi
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb
