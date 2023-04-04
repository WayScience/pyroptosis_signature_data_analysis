"""Parameters for Machine Learning Model Training and Testing"""
from dataclasses import dataclass, field

import torch


@dataclass
class Parameters:
    """Class for keeping constants"""

    # Constants
    # Subset option yes or no? if no SUBSET_NUMBER won't be used
    DATA_SUBSET_OPTION: bool = True

    # number of rows to subset main df for
    # in this case each row is 1 cell
    DATA_SUBSET_NUMBER: int = 1500

    # Batch of data to load into data loader (1 is equivalent to 1 row or 1 cell in this case)
    BATCH_SIZE: int = 100

    # Data proportionality splits for training, validation, and testing
    TRAIN_PROPORTION_SPLIT: float = 0.8
    VALIDATION_PROPORTION_SPLIT: float = 0.1
    TEST_PROPORTION_SPLIT: float = 0.1

    # number of epochs to use for model optimization
    OPTIM_EPOCHS: int = 100
    # number of trials to use for model optimization
    N_TRIALS: int = 50

    # number of epochs to use for optimized model
    TRAIN_EPOCHS: int = 1000

    # device use
    # defined as global for use in the optimizer function and training function
    # global DEVICE
    DEVICE = torch.device("cuda" if (torch.cuda.is_available()) else "cpu")

    MIN_LAYERS: int = 1
    MAX_LAYERS: int = 10

    LAYER_UNITS_MIN: int = 2
    LAYER_UNITS_MAX: int = 50

    DROPOUT_MIN: float = 0.1
    DROPOUT_MAX: float = 0.5
    DROP_OUT_STEP: float = 0.05

    LEARNING_RATE_MIN: float = 1e-5
    LEARNING_RATE_MAX: float = 1

    OPTIMIZER_LIST: list[str] = field(
        default_factory=lambda: ["Adam", "RMSprop", "SGD"]
    )

    METRIC: str = "accuracy"
    DIRECTION: str = "maximize"
