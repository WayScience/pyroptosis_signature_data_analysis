#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=amem
#SBATCH --mem=400G
#SBATCH --qos=mem
#SBATCH --output=sample-%j.out
#SBATCH --time=24:00:00


feature=$1
CELL_TYPE=$2

echo "Feature: $feature" "Cell type: $CELL_TYPE"
