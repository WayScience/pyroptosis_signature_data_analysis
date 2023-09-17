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
shuffles=( True False )
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

num_jobs=$(( ${#control_names[@]} * ${#treatment_names[@]} * ${#cell_types[@]} * ${#shuffles[@]} ))
echo "num_jobs: $num_jobs"
job=1


for cell_type in "${cell_types[@]}"; do
    for control_name in "${!control_names[@]}"; do
        for treatment_name in "${!treatment_names[@]}"; do
            for shuffle in "${!shuffles[@]}"; do
                echo "cell_type: $cell_type" control_name: "${control_names[$control_name]}" treatment_name: "${treatment_names[$treatment_name]}" shuffle: "${shuffles[$shuffle]}"
                job=$(( $job + 1 ))
                progress=$(( $job * 100 / $num_jobs ))
                echo "progress: $progress"
                # if treatment_name is not equal to control_name
                if [[ $treatment_name != "$control_name" ]]; then
                    papermill \
                    train_binary_model.ipynb \
                    train_binary_model.ipynb \
                    -p CELL_TYPE "$cell_type" \
                    -p CONTROL_NAME "${control_names[$control_name]}" \
                    -p TREATMENT_NAME "${treatment_names[$treatment_name]}"\
                    -p SHUFFLE "${shuffles[$shuffle]}"
                fi
            done
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb
