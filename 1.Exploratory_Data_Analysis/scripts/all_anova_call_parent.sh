#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amilan
#SBATCH --qos=long
#SBATCH --output=anova_parent-%j.out
#SBATCH --time=168:00:00


module load anaconda

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=. ../notebooks/*.ipynb


CELL_TYPE="PBMC"

filename="../features/feature_index.txt"
# read all lines of the file to an array
readarray -t features < $filename

echo "Feature: $feature" "Cell type: $CELL_TYPE"

# Loop until the number of jobs is less than 25
while true; do
    # Count the number of running and pending jobs
    total_jobs=$(squeue -h -u $USER | wc -l)

    # If the total number of jobs is less than 1000, submit the next job
    if [ $total_jobs -lt 30 ]; then
        # Get the next item from the list
        feature=${features[$SLURM_ARRAY_TASK_ID-1]}

        # Check if there are still items to process
        if [ -n "$feature" ]; then
            echo "Submitting job for item: $feature"
            sbatch all_anova_call_child.sh $feature $CELL_TYPE
        else
            echo "All items processed. Exiting."
            break
        fi
    else
        echo "Waiting for job slots to become available..."
        sleep 30  # Adjust the sleep duration as needed
    fi
done

