#!/bin/bash

# Run redundancy analysis for SHSY5Y and PBMC cell types with and without shuffling

papermill 0.redundancy_analysis.ipynb 0.redundancy_analysis.ipynb -p cell_type "SHSY5Y" -p Shuffle True
papermill 0.redundancy_analysis.ipynb 0.redundancy_analysis.ipynb -p cell_type "SHSY5Y" -p Shuffle False

papermill 0.redundancy_analysis.ipynb 0.redundancy_analysis.ipynb -p cell_type "PBMC" -p Shuffle True
papermill 0.redundancy_analysis.ipynb 0.redundancy_analysis.ipynb -p cell_type "PBMC" -p Shuffle False

# convert notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=../scripts *.ipynb

# run scripts to generate plots for redundancy analysis
Rscript ../scripts/1.redundancy_visualization.r --cell_type "SHSY5Y"
Rscript ../scripts/1.redundancy_visualization.r --cell_type "PBMC"
