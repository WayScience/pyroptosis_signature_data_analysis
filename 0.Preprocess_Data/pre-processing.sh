#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

cell_types=( PBMC SHSY5Y )
aggragates=( True False )
nomics=( True False )

for cell_type in "${cell_types[@]}"; do
    papermill 1.preprocessing_morphology_data.ipynb 1.preprocessing_morphology_data.ipynb -p cell_type $cell_type
    papermill 2.preprocessing_nELISA_data.ipynb 2.preprocessing_nELISA_data.ipynb -p cell_type $cell_type
    for aggragate in "${aggragates[@]}"; do
        for nomic in "${nomics[@]}"; do
            if [ $aggragate == "False" ] && [ $nomic == "False" ]; then
                continue
            fi
            echo $cell_type $aggraggitate $nomic
            papermill 3.data_aggregation.ipynb 3.data_aggregation.ipynb -p cell_type $cell_type -p aggregation $aggragate -p nomic $nomic
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

