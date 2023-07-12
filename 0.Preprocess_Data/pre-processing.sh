#!/bin/bash

conda activate Interstellar

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "SHSY5Y"
# papermill 1.preprocessing.ipynb 1.preprocessing.ipynb -p celltype "PBMC"
