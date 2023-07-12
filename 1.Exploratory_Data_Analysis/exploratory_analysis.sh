#!/bin/bash

conda activate Interstellar

papermill 0.cell_count_analysis.ipynb 0.cell_count_analysis.ipynb -p celltype "SHSY5Y"
papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p celltype "SHSY5Y"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
