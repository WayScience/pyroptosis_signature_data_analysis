#!/bin/bash

conda activate Interstellar

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

papermill umap_analysis_plate2.ipynb umap_analysis_plate2.ipynb -p celltype "SHSY5Y"
