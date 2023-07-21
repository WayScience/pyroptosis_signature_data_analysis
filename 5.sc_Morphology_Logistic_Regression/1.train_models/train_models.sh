#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=48:00:00
#SBATCH --output=sample-%j.out

# shell loop to run notebook with multiple args
# set args to loop through
cell_types=(SHSY5Y PBMC)
aggregations=(True False)
nomics=(True False)
flag=(True)
controls=(DMSO_0.100_DMSO_0.025)
treatments=(LPS_100.000_DMSO_0.025 LPS_10.000_DMSO_0.025 LPS_1.000_DMSO_0.025 Thapsigargin_10.000_DMSO_0.025 Thapsigargin_1.000_DMSO_0.025)

for cell_type in ${cell_types[@]}; do
    for aggregation in ${aggregations[@]}; do
        for nomic in ${nomics[@]}; do
            for control in ${controls[@]}; do
                for treatment in ${treatments[@]}; do
                    echo "cell_type: $cell_type, aggregation: $aggregation, nomic: $nomic, flag: $flag, control: $control, treatment: $treatment"
                    papermill train_single_class.ipynb train_single_class.ipynb -p cell_type "$cell_type" -p aggregation "$aggregation" -p nomic "$nomic" -p flag "$flag" -p control "$control" -p treatment "$treatment"
                done
            done
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
