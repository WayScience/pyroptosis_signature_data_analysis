# Training Modules

Using `Binary_Classifier_Call.sh` to train the models on a slurm HPC system.
The PBMC dataset is too large to train on a local machine.

I submit the job via:
```bash
sbatch Binary_Classifier_Call.sh
```

With the following job parameters:
```bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=300G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --time=24:00:00
#SBATCH --output=sample-%j.out
```
which translates to:
- 1 node
- 16 tasks
- 300GB of memory
- partition: amem (the high memory partition)
- qos: mem (the high memory quality of service)
- 24 hours of runtime
- output file: sample-%j.out (where %j is the job id)
