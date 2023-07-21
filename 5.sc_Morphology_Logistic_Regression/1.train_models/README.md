# Training models

Each logistic regression model will be run using:
- Model 1: Logistic regression - Binary Classifier
- Model 2: Logistic regression - Multiclass Classifier

Each config has the follwing parameters:
- `MODEL_TYPE`: This is the type of model to run. Either `Binary` or `Multiclass`
- `control`: The control treatment group in the column specified in the code `DMSO_0.100_DMSO_0.025`
- `treatments`: The treatment group(s) in the column specified in the code `LPS_100.000_DMSO_0.025`
- `aggregation`: A boolean that dictates if the single cell data should be aggregated by well or not `True` or `False`
- `nomic`: A boolean that dictates if the Nomic nELISA data will be used in the model or not `True` or `False`
- `cell_type`: A string for the cell type that is to be analyzed in this analaysis the options are `SHSY5Y` or `PBMC`

## Outputs
Each model will output:
* A model object
* A Confusion matrix
* A `.toml` file with the parameters used to run the model for reproducibility
* Two `.joblib` files with the model and scaler objects for reproducibility
    * One `.joblib` file for the trained model
    * One `.joblib` file for the shuffled model

