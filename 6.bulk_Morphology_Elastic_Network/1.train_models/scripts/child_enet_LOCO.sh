#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --partition=amem
#SBATCH --account=amc-general
#SBATCH --qos=mem
#SBATCH --time=72:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

# pass through the cytokine and channels from the call script
cell_type="$1"
cytokine="$2"
shuffle="$3"
data="$4"

echo "cell_type: $cell_type cytokine: $cytokine shuffle: $shuffle data: $data"

# call the python script
python 2.train_regression_multi_output_channel_selection.py \
--cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle" --data "$data"

