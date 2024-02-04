#!/bin/bash
# This script is used to train the regression models for the elastic network

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=sample-%j.out

module load anaconda
conda activate Interstellar_python
jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb

# pass through the cytokine and channels from the call script
cell_type="$1"
cytokine="$2"
shuffle="$3"
channel="$4"

echo "cell_type: $cell_type cytokine: $cytokine shuffle: $shuffle data: $channel"

# call the python script
python 4.test_regression_multi_output_channel_split.py --cell_type "$cell_type" --cytokine "$cytokine" --shuffle "$shuffle" --data "$channel"

echo "job finished"
