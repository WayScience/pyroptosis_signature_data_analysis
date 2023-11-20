#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=700G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out

module load anaconda

conda activate Interstellar_python

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts *.ipynb

papermill 8.anova_all_groupings.ipynb 8.anova_all_groupings.ipynb

echo "Complete"
