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

cell_types=( PBMC SHSY5Y)
aggragates=( True False )
nomics=( True False )

for cell_type in $cell_types; do
    for aggragate in $aggragates; do
        for nomic in $nomics; do
            echo $cell_type $aggragate $nomic
            papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p cell_type $cell_type -p aggregation $aggragate -p nomic $nomic
            papermill 2.data_aggregation.ipynb 2.data_aggregation.ipynb -p cell_type $cell_type -p aggregation $aggragate -p nomic $nomic
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
