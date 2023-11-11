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
model_names=( MultiClass_MLP_h202_remove MultiClass_MLP )
shuffles=( True False )

for cell_type in "${cell_types[@]}"; do
    for shuffle in "${shuffles[@]}"; do
        for model_name in "${model_names[@]}"; do
            papermill \
            Hyperparameter_Optimization_model_binary.ipynb \
            Hyperparameter_Optimization_model_binary.ipynb \
            -p CELL_TYPE $cell_type \
            -p MODEL_NAME $model_name \
            -p SHUFFLE $shuffle
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb

echo "Done"
