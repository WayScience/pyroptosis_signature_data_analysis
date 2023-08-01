# This folder hosts jupyter notebooks for exploratory data analysis (EDA) for multiple 'waves' of data. See this repo's main page for a link to the Interstellar preprocessing page.

Where 'waves' of data are data from the same plate but using different channels for different analysis such as wavelength bleed through or antibody validation.

```Exploratory_Data_Analysis.ipynb```
is for clustering of cells across all features.
The goal is to aim to see if treatments cluster or not.
```cell_count_analysis.ipynb``` is for qc analysis to assess if there are cell count drop offs for a given treatment.

```umap_analysis.ipynb``` runs umap analysis on all image features across all treatments to see if replicates seperate or not.

```2.extended_umap_analysis.ipynb``` runs umap analysis on all image features across some treatments of interest (LPS, Thapsigargin, DMSO, etc) but creates one umap per image channel to decipher which channels are the most variable across single cells.

All figures are saved under `Figures` folder.


Small note about the `bash` files: these are meant to be run on a SLURM cluster.
Our PBMC data is large and requires a lot of memory for most analysis and one of our `128GiB` RAM workstations cannot handel these analysis.

This Repository will not go into the details of how to run these analysis on a SLURM cluster. Please see the [SLURM documentation](https://slurm.schedmd.com/documentation.html) for more information.
