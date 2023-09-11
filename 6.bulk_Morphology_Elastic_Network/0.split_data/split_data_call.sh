#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

shuffles=(True False)
cell_types=( SHSY5Y PBMC )
# loop through the define cell types and if the data should be shuffled or not
for shuffle in "${shuffles[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        # use papermill to run the notebooks with injected parameters
        papermill \
            # input notebook
            0.split_data_regression.ipynb \
            # output notebook
            0.split_data_regression.ipynb \
            # parameters by name followed by the value
            -p cell_type $cell_type \
            -p aggregation True \
            -p nomic True \
            -p flag True \
            -p control "DMSO_0.100_DMSO_0.025" \
            -p treatment "LPS_100.000_DMSO_0.025" \
            -p shuffle $shuffle
    done
done
