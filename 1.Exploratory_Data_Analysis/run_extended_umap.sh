#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

# module load anaconda

conda activate Interstellar

cell_types=(SHSY5Y PBMC)
controls=(DMSO_0.100_DMSO_0.025)
treatments=(LPS_100.000_DMSO_0.025 LPS_10.000_DMSO_0.025 LPS_1.000_DMSO_0.025 Thapsigargin_10.000_DMSO_0.025 Thapsigargin_1.000_DMSO_0.025)

for cell_type in ${cell_types[@]}; do
    for control in ${controls[@]}; do
        for treatment in ${treatments[@]}; do
            echo "cell_type: $cell_type, control: $control, treatment: $treatment"
            papermill 2.extended_umap_analysis.ipynb 2.extended_umap_analysis.ipynb -p cell_type "$cell_type" -p control "$control" -p treatment "$treatment"
        done
    done
done

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb
