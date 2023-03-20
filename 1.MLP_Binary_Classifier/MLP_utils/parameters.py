from dataclasses import dataclass

import torch


@dataclass
class Parameters:
    """Class for keeping constants"""

    def __init__(self):
        # Constants
        # Subset option yes or no? if no SUBSET_NUMBER won't be used
        self.SUBSET_OPTION = False

        # number of rows to subset main df for
        # in this casse each row is 1 cell
        self.SUBSET_NUMBER = 15000

        # Batch of data to load into data loader (1 is equivalent to 1 row or 1 cell in this case)
        self.BATCH_SIZE = 100

        # number of epochs to use for model optimization
        self.OPTIM_EPOCHS = 100
        # number of trials to use for model optimization
        self.N_TRIALS = 500

        # number of epochs to use for optimized model
        self.TRAIN_EPOCHS = 1000

        # device use
        # defined as global for use in the optimizer function and training function
        # global DEVICE
        self.DEVICE = torch.device("cuda" if (torch.cuda.is_available()) else "cpu")
