#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=700G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%

module load anaconda

conda activate Interstellar

treatment1s=( DMSO_0.100_DMSO_0.025 )
treatment2s=( LPS_100.000_DMSO_0.025 LPS_10.000_DMSO_0.025 LPS_1.000_DMSO_0.025 LPS_0.100_DMSO_0.025 LPS_0.010_DMSO_0.025)
treatment3s=( Thapsigargin_10.000_DMSO_0.025 Thapsigargin_1.000_DMSO_0.025 H202_100_Z-VAD-FMK_100 H202_100_DMSO_0.025)
cell_types=( SHSY5Y PBMC )

for treatment1 in ${treatment1s[@]}; do
    for treatment2 in ${treatment2s[@]}; do
        for treatment3 in ${treatment3s[@]}; do
            for cell_type in ${cell_types[@]}; do
                echo "treatment1: $treatment1 treatment2: $treatment2 treatment3: $treatment3 cell_type: $cell_type"
                papermill 5.anova.ipynb \
                5.anova.ipynb \
                -p cell_type $cell_type \
                -p treatment1 $treatment1 \
                -p treatment2 $treatment2 \
                -p treatment3 $treatment3
            done
        done
    done
done

# convert the notebook to script
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
