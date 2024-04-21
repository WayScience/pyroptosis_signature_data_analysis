#!/bin/bash
# this is the script that should be run if using a local machine
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
