#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --output=sample-%j.out
#SBATCH --array=1-999
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=250G


module load anaconda

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

CELL_TYPE="SHSY5Y"

echo "creating features list"
papermill notebooks/8.0_create_feature_index.ipynb notebooks/8.0_create_feature_index_out.ipynb -p cell_type "$CELL_TYPE"

filename="./features/${CELL_TYPE}_features_index.txt"
# read all lines of the file to an array
readarray -t features < $filename

echo "Feature: $feature" "Cell type: $CELL_TYPE"

# Loop until the number of jobs is less than 1000
while true; do
    # Count the number of running and pending jobs
    total_jobs=$(squeue -h -u $USER | wc -l)

    # If the total number of jobs is less than 1000, submit the next job
    if [ $total_jobs -lt 990 ]; then
        # Get the next item from the list
        feature=${features[$SLURM_ARRAY_TASK_ID-1]}

        # Check if there are still items to process
        if [ -n "$feature" ]; then
            echo "Submitting job for item: $feature"
            srun python scripts.8.1_anova_all_groupings.py --feature $feature --cell_type $CELL_TYPE
        else
            echo "All items processed. Exiting."
            break
        fi
    else
        echo "Waiting for job slots to become available..."
        sleep 60  # Adjust the sleep duration as needed
    fi
done

