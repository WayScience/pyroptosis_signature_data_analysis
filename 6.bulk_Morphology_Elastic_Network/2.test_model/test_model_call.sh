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
aggregations=( True False )


test_run() {
    shuffle=$1
    cell_type=$2
    aggregation=$3

    echo "shuffle: $shuffle" "cell_type: $cell_type" aggregation: $aggregation
    papermill \
        # input notebook
        1.test_regression_multi_ouput.ipynb \
        # output notebook
        1.test_regression_multi_ouput.ipynb \
        -p cell_type $cell_type \
        -p aggregation $aggregation \
        -p nomic True \
        -p flag True \
        -p shuffle $shuffle

    papermill \
        2.visualize_regression_multi_ouput.ipynb \
        2.visualize_regression_multi_ouput.ipynb \
        -p cell_type $cell_type \
        -p aggregation $aggregation \
        -p nomic True \
        -p flag True \
        -p shuffle $shuffle
}

export -f test_run

parallel --jobs 4 test_run ::: "${shuffles[@]}" ::: "${cell_types[@]}" ::: "${aggregations[@]}"

# loop through the define cell types and if the data should be shuffled or not
for shuffle in "${shuffles[@]}"; do
    for cell_type in "${cell_types[@]}"; do
        # use papermill to run the notebooks with injected parameters
        papermill \
            # input notebook
            1.test_regression_multi_ouput.ipynb \
            # output notebook
            1.test_regression_multi_ouput.ipynb \
            # parameters by name followed by the value
            -p cell_type $cell_type \
            -p aggregation True \
            -p nomic True \
            -p flag True \
            -p control "DMSO_0.100_DMSO_0.025" \
            -p treatment "LPS_100.000_DMSO_0.025" \
            -p shuffle $shuffle

        # use papermill to run the notebooks with injected parameters
        papermill \
        2.visualize_regression_multi_ouput.ipynb \
        2.visualize_regression_multi_ouput.ipynb \
        -p cell_type $cell_type \
        -p aggregation True \
        -p nomic True \
        -p flag True \
        -p control "DMSO_0.100_DMSO_0.025" \
        -p treatment "LPS_100.000_DMSO_0.025" \
        -p shuffle $shuffle

    done
done
