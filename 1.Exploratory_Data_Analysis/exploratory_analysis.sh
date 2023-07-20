#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar

papermill 0.cell_count_analysis.ipynb 0.cell_count_analysis.ipynb -p celltype "SHSY5Y"
papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p celltype "SHSY5Y"
papermill 2.Corelation_between_nomic_and_sc.ipynb 2.Corelation_between_nomic_and_sc.ipynb -p celltype "SHSY5Y"

papermill 0.cell_count_analysis.ipynb 0.cell_count_analysis.ipynb -p celltype "PBMC"
papermill 1.umap_analysis_plate2.ipynb 1.umap_analysis_plate2.ipynb -p celltype "PBMC"
papermill 2.Corelation_between_nomic_and_sc.ipynb 2.Corelation_between_nomic_and_sc.ipynb -p celltype "PBMC"


jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
