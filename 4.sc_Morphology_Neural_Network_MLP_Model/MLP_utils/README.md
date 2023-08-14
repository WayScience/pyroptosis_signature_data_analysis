# MLP_utils
## `Exceptions.py`
This file contains the exceptions rules for the MLP model.

## `parameters.py`
This file contains the parameters dataclass used for the mlp_params in:
* Hyperparameter Optimization
* Training
* Testing

## `utils.py`
This file contains the utility/helper functions for the MLP model.

## trained_models
This directory contains the trained models from the hyperparameter optimization and training.
This directory contains two subdirectories:
* architectures
    * Saves the model hyperparameters and architecture of the neural network
* model_save_states
    * saves the model weights and biases from the trained model

## in the `config.toml` file
### Binary Classification Configurations in `binary_config.toml`

Each part of the configuration file is explained below:

MODEL_TYPE: Binary Classification or Multi-Class Classification
DATA_SUBSET_OPTION: True or False for using a subset of the data
DATA_SUBSET_NUMBER: The size of the subset of the data
BATCH_SIZE: The batch size for the training
TRAIN_PROPORTION_SPLIT: The proportion of the data used for training (fraction of 1 e.g. 0.8)
VALIDATION_PROPORTION_SPLIT: The proportion of the data used for validation (fraction of 1 e.g. 0.1)
TEST_PROPORTION_SPLIT: The proportion of the data used for testing (fraction of 1 e.g. 0.1)
OPTIM_EPOCHS: The number of epochs for the hyperparameter optimization
N_TRIALS: The number of trials for the hyperparameter optimization
TRAIN_EPOCHS: The number of epochs for the training of model
MIN_Layers: The minimum number of layers for the hyperparameter optimization
MAX_Layers: The maximum number of layers for the hyperparameter optimization
LAYER_UNITS_MIN: The minimum number of units per layer for the hyperparameter optimization
LAYER_UNITS_MAX: The maximum number of units per layer for the hyperparameter optimization
DROP_OUT_MIN: The minimum dropout rate for the hyperparameter optimization
DROP_OUT_MAX: The maximum dropout rate for the hyperparameter optimization
DROP_OUT_STEP: The step size for the dropout rate for the hyperparameter optimization
LEARNING_RATE_MIN: The minimum learning rate for the hyperparameter optimization
LEARNING_RATE_MAX: The maximum learning rate for the hyperparameter optimization
OPTIMIZER_LIST: The list of optimizers for the hyperparameter optimization
METRIC: The metric for the hyperparameter optimization
DIRECTION: The direction for the hyperparameter optimization to use
MODEL_NAME: The name of the model
CONTROL_NAME: The name of the control group
TREATMENT_NAME: The name of the treatment group
CELL_TYPE: The cell type of the data

