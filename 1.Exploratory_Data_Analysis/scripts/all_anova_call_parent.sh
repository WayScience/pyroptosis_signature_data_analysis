#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --output=anova_parent-%j.out
#SBATCH --time=168:00:00


# module load anaconda

# conda activate Interstellar_python

# jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb


CELL_TYPE="PBMC"

filename="../features/feature_index.txt"
# read all lines of the file to an array
readarray -t features < $filename

echo "Feature: $feature" "Cell type: $CELL_TYPE"

# loop through the array
for feature in "${features[@]}"; do
    # echo "Feature: $feature" "Cell type: $CELL_TYPE"
    # check if the number of slurms is less than 30
    while true; do
        NUM_SLURMS=$(squeue -u $USER | wc -l)
        if [ "$NUM_SLURMS" -lt 30 ]; then
            echo "Feature: $feature" "Cell type: $CELL_TYPE"
            sbatch anova_call_test.sh $feature $CELL_TYPE
            break
        else
            sleep 10
        fi
    done
done
