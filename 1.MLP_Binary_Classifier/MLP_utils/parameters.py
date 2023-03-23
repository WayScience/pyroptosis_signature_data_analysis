"""Parametes for Machine Learning Model Training and Testing"""
from dataclasses import dataclass

import torch


@dataclass
class Parameters:
    """Class for keeping constants"""

    # Constants
    # Subset option yes or no? if no SUBSET_NUMBER won't be used
    DATA_SUBSET_OPTION = True

    # number of rows to subset main df for
    # in this case each row is 1 cell
    DATA_SUBSET_NUMBER = 1500

    # Batch of data to load into data loader (1 is equivalent to 1 row or 1 cell in this case)
    BATCH_SIZE = 100

    # number of epochs to use for model optimization
    OPTIM_EPOCHS = 100
    # number of trials to use for model optimization
    N_TRIALS = 50

    # number of epochs to use for optimized model
    TRAIN_EPOCHS = 1000

    # device use
    # defined as global for use in the optimizer function and training function
    # global DEVICE
    DEVICE = torch.device("cuda" if (torch.cuda.is_available()) else "cpu")

    MIN_LAYERS = 1
    MAX_LAYERS = 10

    LAYER_UNITS_MIN = 2
    LAYER_UNITS_MAX = 50

    DROPOUT_MIN = 0.1
    DROPOUT_MAX = 0.5
    DROP_OUT_STEP = 0.05

    LEARNING_RATE_MIN = 1e-5
    LEARNING_RATE_MAX = 1

    OPTIMIZER_LIST = ["Adam", "RMSprop", "SGD"]

    METRIC = "accuracy"
    DIRECTION = "maximize"
