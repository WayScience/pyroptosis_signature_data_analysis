#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=400G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=48:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cell_types=( SHSY5Y PBMC )
model_names=( MultiClass_MLP )

for cell_type in "${cell_types[@]}"; do
    for model_name in "${model_names[@]}"; do
        papermill \
        Hyperparameter_Optimization_model_multiclass.ipynb \
        Hyperparameter_Optimization_model_multiclass.ipynb \
        -p CELL_TYPE $cell_type \
        -p MODEL_NAME $model_name
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb

echo "Done"
