#!/bin/bash

echo "Running feature index for all cell types"

cd scripts/

python 8.0_create_feature_index.py --cell_type "PBMC"

python 8.0_create_feature_index.py --cell_type "SHSY5Y"

cd ../

echo "Complete"
