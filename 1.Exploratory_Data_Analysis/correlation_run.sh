#!/bin/bash


papermill 2a.correlation_nomic_features.ipynb 2a.correlation_nomic_features.ipynb -p cell_type "SHSY5Y"
papermill 2a.correlation_nomic_features.ipynb 2a.correlation_nomic_features.ipynb -p cell_type "PBMC"

papermill 2b.correlation_agg_morphology_features.ipynb 2b.correlation_agg_morphology_features.ipynb -p cell_type "SHSY5Y"
papermill 2b.correlation_agg_morphology_features.ipynb 2b.correlation_agg_morphology_features.ipynb -p cell_type "PBMC"

papermill 2c.correlation_nomic_and_agg_morphology.ipynb 2c.correlation_nomic_and_agg_morphology.ipynb -p cell_type "SHSY5Y"
papermill 2c.correlation_nomic_and_agg_morphology.ipynb 2c.correlation_nomic_and_agg_morphology.ipynb -p cell_type "PBMC"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

celltypes=( SHSY5Y PBMC )
levels=( nomic aggregated_morphology aggregated_morphology_and_nomic )
groups=( wells_corr treatments_corr selected_treatments_corr )

for celltype in "${celltypes[@]}"; do
    for level in "${levels[@]}"; do
        for group in ${groups[@]}; do
            Rscript scripts/3.correlation_visualization.r --cell_type="$celltype" --level="$level" --group="$group"
        done
    done
done
