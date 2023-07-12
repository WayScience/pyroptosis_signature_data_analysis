#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=48:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

papermill 1a.fit_linear_model_1beta.ipynb 1a.fit_linear_model_1beta.ipynb -p celltype 'SHSY5Y'
papermill 1b.fit_linear_model_2beta.ipynb 1b.fit_linear_model_2beta.ipynb -p celltype 'SHSY5Y'
papermill 1c.fit_linear_model_3beta.ipynb 1c.fit_linear_model_3beta.ipynb -p celltype 'SHSY5Y'
papermill 1d.fit_linear_model_4beta.ipynb 1d.fit_linear_model_4beta.ipynb -p celltype 'SHSY5Y'

papermill 1a.fit_linear_model_1beta.ipynb 1a.fit_linear_model_1beta.ipynb -p celltype 'PBMC'
papermill 1b.fit_linear_model_2beta.ipynb 1b.fit_linear_model_2beta.ipynb -p celltype 'PBMC'
papermill 1c.fit_linear_model_3beta.ipynb 1c.fit_linear_model_3beta.ipynb -p celltype 'PBMC'
papermill 1d.fit_linear_model_4beta.ipynb 1d.fit_linear_model_4beta.ipynb -p celltype 'PBMC'

# copnvert to scripts after running papermill and conver r notebooks to scripts prior to running
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

Rscript scripts/2a.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'
Rscript scripts/2b.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'
Rscript scripts/2c.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'
Rscript scripts/2d.Linear_Modeling_Visualization.r --celltype 'SHSY5Y'

Rscript scripts/2a.Linear_Modeling_Visualization.r --celltype 'PBMC'
Rscript scripts/2b.Linear_Modeling_Visualization.r --celltype 'PBMC'
Rscript scripts/2c.Linear_Modeling_Visualization.r --celltype 'PBMC'
Rscript scripts/2d.Linear_Modeling_Visualization.r --celltype 'PBMC'
