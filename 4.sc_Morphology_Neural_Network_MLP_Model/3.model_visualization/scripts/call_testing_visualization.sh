#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

# module load anaconda

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

# generate selected_treatment_comparisons from model_names (2 by 2)
selected_treatment_comparisons=()
for model_name in "${model_names[@]}"; do
    for model_name_2 in "${model_names[@]}"; do
	echo $model_name $model_name_2
	 if [[ $model_name != "$model_name_2" ]]; then
            selected_treatment_comparisons+=("${model_name},${model_name_2}")
        fi
    done
done


# save notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb*

# get the number of combination of control_names treatment_names and cell_types
# this is the number of jobs
num_jobs=$(( ${#cell_types[@]} * ${#model_names[@]} * ${#selected_treatment_comparisons[@]} ))
echo "num_jobs: $num_jobs"

job_id=$(( SLURM_ARRAY_TASK_ID - 1 ))
cell_type_idx=$(( job_id % "${#cell_types[@]}" ))
model_name_idx=$(( job_id / "${#cell_types[@]}" % "${model_names[@]}" ))
treatment_comparisons_idx=$(( job_id / "${cell_types[@]}" / "${model_names[@]}" % "${treatment_comparisons[@]}" ))

cell_type=${cell_types[$cell_type_idx]}
model_name=${model_names[$model_name_idx]}
treatment_comparison=${treatment_comparisons[$treatment_comparisons_idx]}

echo $cell_type $model_name $treatment_comparison

#command="Rscript binary_classification_testing_visualization.r"
#$command --cell_type "$cell_type" --model_name "$model_name" --selected_treatment_comparisons "{$selected_treatment_comparison}"


