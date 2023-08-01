#!/bin/bash


# Use papermill to run the notebooks with different cell types

# preprocess the data
papermill 0.nomic_cleanup.ipynb 0.nomic_cleanup.ipynb -p cell_type "SHSY5Y"
papermill 0.nomic_cleanup.ipynb 0.nomic_cleanup.ipynb -p cell_type "PBMC"

# Run the correlation of the nomic features
papermill 2a.correlation_nomic_features.ipynb 2a.correlation_nomic_features.ipynb -p cell_type "SHSY5Y"
papermill 2a.correlation_nomic_features.ipynb 2a.correlation_nomic_features.ipynb -p cell_type "PBMC"

# Run the correlation of the aggregated morphology features
papermill 2b.correlation_agg_morphology_features.ipynb 2b.correlation_agg_morphology_features.ipynb -p cell_type "SHSY5Y"
papermill 2b.correlation_agg_morphology_features.ipynb 2b.correlation_agg_morphology_features.ipynb -p cell_type "PBMC"

# Run the correlation of the aggregated morphology and nomic features
papermill 2c.correlation_nomic_and_agg_morphology.ipynb 2c.correlation_nomic_and_agg_morphology.ipynb -p cell_type "SHSY5Y"
papermill 2c.correlation_nomic_and_agg_morphology.ipynb 2c.correlation_nomic_and_agg_morphology.ipynb -p cell_type "PBMC"

# convert the notebooks to scripts
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

# define the cell types, levels and groups to loop through for the correlation visualization
celltypes=( SHSY5Y PBMC )
levels=( nomic aggregated_morphology aggregated_morphology_and_nomic )
groups=( wells treatments selected_treatments )

for celltype in "${celltypes[@]}"; do
    for level in "${levels[@]}"; do
        for group in ${groups[@]}; do
            # run the correlation visualization script
            Rscript scripts/3.correlation_visualization.r --cell_type="$celltype" --level="$level" --group="$group"
        done
    done
done
