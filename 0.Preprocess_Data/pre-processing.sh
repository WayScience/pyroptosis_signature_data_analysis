#!/bin/bash

conda activate Interstellar

papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "SHSY5Y"
papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "PBMC"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
