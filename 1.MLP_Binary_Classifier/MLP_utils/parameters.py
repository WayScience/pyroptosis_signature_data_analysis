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
