#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=noraml
#SBATCH --output=anova_parent-%j.out
#SBATCH --time=1:00:00


module load anaconda

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/ notebooks/*.ipynb

# change the directory to the scripts directory
cd scripts/ || exit

CELL_TYPE="PBMC"

filename="../features/${CELL_TYPE}_feature_index.txt"
# read all lines of the file to an array
readarray -t features < $filename

echo "Feature: $feature" "Cell type: $CELL_TYPE"
shuffles=( "True" "False" )

NUM_SLURMS=$(squeue -u $USER | wc -l)
# loop through the array
for feature in "${features[@]}"; do
    for shuffle in "${shuffles[@]}"; do
        echo "Feature: $feature" "Cell type: $CELL_TYPE" "Shuffle: $shuffle"
        # check if the number of slurms is less than 30
        while [ "$NUM_SLURMS" -gt 995 ]; do
           sleep 1s
            NUM_SLURMS=$(squeue -u $USER | wc -l)
        done
        sbatch all_anova_call_child.sh $feature $CELL_TYPE $shuffle
    done
done

cd ../ || exit
echo "All features submitted for ANOVA analysis"
