#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=700G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_python

cell_types=( SHSY5Y PBMC )
aggragates=( True )
nomics=( True )

for cell_type in ${cell_types[@]}; do
    echo $cell_type
    papermill 0.preprocess_profiles.ipynb 0.preprocess_profiles.ipynb -p cell_type $cell_type
    for aggragate in ${aggragates[@]}; do
        for nomic in ${nomics[@]}; do
            echo $cell_type $aggragate $nomic
                papermill 1.data_curation.ipynb 1.data_curation.ipynb -p cell_type $cell_type -p aggregation $aggragate -p nomic $nomic
        done
    done
done

# create the processed file
papermill 2.create_processed.ipynb 2.create_processed.ipynb
# create the library file
papermill 3.create_library.ipynb 3.create_library.ipynb

# convert the notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb
