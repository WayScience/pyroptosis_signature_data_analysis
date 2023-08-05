#!/bin/bash


treatment1s=(DMSO_0.100_DMSO_0.025)
treatment2s=(LPS_100.000_DMSO_0.025)
treatment3s=(Thapsigargin_10.000_DMSO_0.025)
cell_types=(SHSY5Y)

for treatment1 in ${treatment1s[@]}; do
    for treatment2 in ${treatment2s[@]}; do
        for treatment3 in ${treatment3s[@]}; do
            for cell_type in ${cell_types[@]}; do
            echo "treatment1: $treatment1 treatment2: $treatment2 treatment3: $treatment3 cell_type: $cell_type"
                papermill 5.anova.ipynb \
                5.anova.ipynb \
                -p cell_type SHSY5Y \
                -p treatment1 DMSO_0.100_DMSO_0.025 \
                -p treatment2 LPS_100.000_DMSO_0.025 \
                -p treatment3 Thapsigargin_10.000_DMSO_0.025
            done
        done
    done
done

# conver the notebook to script
jupyter nbconvert --to=script --FilesWriter.build_directory=/scripts *.ipynb
