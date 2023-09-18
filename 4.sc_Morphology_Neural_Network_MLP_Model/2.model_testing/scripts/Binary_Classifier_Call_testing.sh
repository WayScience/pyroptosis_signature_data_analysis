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

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

job_id=$(( SLURM_ARRAY_TASK_ID - 1 ))
shuffle_idx=$(( job_id % "${#shuffles[@]}" ))
cell_type_idx=$(( job_id / "${#shuffles[@]}" % "${#cell_types[@]}" ))
control_name_idx=$(( job_id / "${#shuffles[@]}" / "${#cell_types[@]}" % "${#control_names[@]}" ))
treatment_name_idx=$(( job_id / "${#shuffles[@]}" / "${#cell_types[@]}" / "${#control_names[@]}" % "${#treatment_names[@]}" ))

shuffle=${shuffles[$shuffle_idx]}
cell_type=${cell_types[$cell_type_idx]}
control_name=${control_names[$control_name_idx]}
treatment_name=${treatment_names[$treatment_name_idx]}

echo cell_type: $cell_type control_name: $control_name treatment_name: $treatment_name shuffle: $shuffle

command="python binaary_classification_testing.py"
$command --cell_type "$cell_type" --control_name "$control_name" --treatment_name "$treatment_name" --shuffle "$shuffle"

echo "completed"

