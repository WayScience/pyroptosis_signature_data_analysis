#!/bin/bash

conda activate Interstellar

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb


papermill 1a.fit_linear_model_1beta.ipynb 1a.fit_linear_model_1beta.ipynb -p celltype 'SHSY5Y'
Rscript scripts/2a.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'

papermill 1b.fit_linear_model_2beta.ipynb 1b.fit_linear_model_2beta.ipynb -p celltype 'SHSY5Y'
Rscript scripts/2b.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'

papermill 1c.fit_linear_model_3beta.ipynb 1c.fit_linear_model_3beta.ipynb -p celltype 'SHSY5Y'
Rscript scripts/2c.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'

papermill 1d.fit_linear_model_4beta.ipynb 1d.fit_linear_model_4beta.ipynb -p celltype 'SHSY5Y'
Rscript scripts/2d.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'


jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
