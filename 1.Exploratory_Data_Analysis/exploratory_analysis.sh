#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=256G
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --output=sample-%j.out

module load anaoconda

conda activate Interstellar

papermill 0.cell_count_analysis.ipynb 0.cell_count_analysis.ipynb -p celltype "SHSY5Y"
papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p celltype "SHSY5Y"

papermill 0.cell_count_analysis.ipynb 0.cell_count_analysis.ipynb -p celltype "PBMC"
papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p celltype "PBMC"

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
