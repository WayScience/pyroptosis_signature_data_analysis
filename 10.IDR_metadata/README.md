# Image Data Resource (IDR) Data Curation

This module contains the necessary scripts to adequately curate aggregated image-based profiles, nELISA secreted profiles, and metadata.
These data will be uploaded to IDR with the raw image data.

## Usage
There are four main notebooks in this module:
1. `0.preprocess_profiles.ipynb`: This notebooks curates and preprocesses the data for IDR.
Specifically, it gets the normalized, non-feature selected data ready for IDR upload.
2. `1.data_curation.ipynb`: This noteboook aggregates the data from the previous notebooks and creates the final aggregated dataset, which includes concatenated image-based profiles and nELISA data per well.
Here we curate the data for IDR.
3. `2.create_processed.ipynb`: This notebook creates the [processed](IDR_metadata/screenA/idr0000-screenA-processed.txt) metadata file for IDR.
This file will include the aggregated Image-based profiles, the nELISA profiles, and the metadata for the samples.
3. `3.create_library.ipynb`: This notebook creats the [library](IDR_metadata/screenA/idr0000-screenA-library.txt) metadata file needed fop IDR upload.

Each of these notebooks can be run via the [run_preprocessing_and_aggregation_local.sh](notebooks/run_preprocessing_and_aggregation_local.sh) script.
Note that there are two versions of the script.
One for a local machine and another for a slurm based cluster.
This can be accomplished by running the following command in the terminal:
For a local machine:
```bash
cd notebooks
bash run_preprocessing_and_aggregation_local.sh
```
For a slurm based cluster:
```bash
cd notebooks
bash run_preprocessing_and_aggregation_cluster.sh
```

Note for the slurm based cluster, the user might need to update the slurm parameters in the script depending on the cluster configuration.
In this script we use these parameters:
```bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=700G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out
```
Where `--nodes` is the number of nodes, `--ntasks` is the number of tasks, `--mem` is the memory, `--partition` is the partition, `--qos` is the quality of service, `--account` is the account, `--time` is the max wall time, and `--output` is the standard out output file.
