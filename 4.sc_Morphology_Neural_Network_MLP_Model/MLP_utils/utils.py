"""
Functions for Machine Learning Model optimization, training and testing.
These are helper functions meant to be called in a separate notebook or script
"""

import json
import pathlib
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import seaborn as sns
import toml
import torch
import torch.nn as nn
import torch.optim as optim
from MLP_utils.exceptions import (
    ModelNameError,
    ModelTypeError,
    OptimizationMetricError,
    TrainingValidationTestingSplitError,
    YDataTypeError,
)
from MLP_utils.parameters import Parameters
from sklearn.metrics import (
    auc,
    classification_report,
    confusion_matrix,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize


def parameter_set(params: Parameters, config: toml) -> object:
    """reset parameter class defaults by updating from config

    Parameters
    ----------
    params : Parameters
        param class holding parameter information
    config: toml
        config class

    Returns
    -------
    object
        param class holding updated parameter information
    """
    params.MODEL_TYPE = config["MACHINE_LEARNING_PARAMS"]["MODEL_TYPE"]
    params.MODEL_NAME = config["MACHINE_LEARNING_PARAMS"]["MODEL_NAME"]
    params.DATA_SUBSET_OPTION = config["MACHINE_LEARNING_PARAMS"]["DATA_SUBSET_OPTION"]
    params.DATA_SUBSET_NUMBER = int(
        config["MACHINE_LEARNING_PARAMS"]["DATA_SUBSET_NUMBER"]
    )
    params.BATCH_SIZE = int(config["MACHINE_LEARNING_PARAMS"]["BATCH_SIZE"])
    params.TRAIN_PROPORTION_SPLIT = float(
        config["MACHINE_LEARNING_PARAMS"]["TRAIN_PROPORTION_SPLIT"]
    )
    params.VALIDATION_PROPORTION_SPLIT = float(
        config["MACHINE_LEARNING_PARAMS"]["VALIDATION_PROPORTION_SPLIT"]
    )
    params.TEST_PROPORTION_SPLIT = float(
        config["MACHINE_LEARNING_PARAMS"]["TEST_PROPORTION_SPLIT"]
    )
    params.OPTIM_EPOCHS = int(config["MACHINE_LEARNING_PARAMS"]["OPTIM_EPOCHS"])
    params.N_TRIALS = int(config["MACHINE_LEARNING_PARAMS"]["N_TRIALS"])
    params.TRAIN_EPOCHS = int(config["MACHINE_LEARNING_PARAMS"]["TRAIN_EPOCHS"])
    params.MIN_LAYERS = int(config["MACHINE_LEARNING_PARAMS"]["MIN_LAYERS"])
    params.MAX_LAYERS = int(config["MACHINE_LEARNING_PARAMS"]["MAX_LAYERS"])
    params.LAYER_UNITS_MIN = int(config["MACHINE_LEARNING_PARAMS"]["LAYER_UNITS_MIN"])
    params.LAYER_UNITS_MAX = int(config["MACHINE_LEARNING_PARAMS"]["LAYER_UNITS_MAX"])
    params.DROPOUT_MIN = float(config["MACHINE_LEARNING_PARAMS"]["DROPOUT_MIN"])
    params.DROPOUT_MAX = float(config["MACHINE_LEARNING_PARAMS"]["DROPOUT_MAX"])
    params.DROP_OUT_STEP = float(config["MACHINE_LEARNING_PARAMS"]["DROP_OUT_STEP"])
    params.LEARNING_RATE_MIN = float(
        config["MACHINE_LEARNING_PARAMS"]["LEARNING_RATE_MIN"]
    )
    params.LEARNING_RATE_MAX = float(
        config["MACHINE_LEARNING_PARAMS"]["LEARNING_RATE_MAX"]
    )
    params.OPTIMIZER_LIST = config["MACHINE_LEARNING_PARAMS"]["OPTIMIZER_LIST"]
    params.METRIC = config["MACHINE_LEARNING_PARAMS"]["METRIC"]
    params.DIRECTION = config["MACHINE_LEARNING_PARAMS"]["DIRECTION"]
    params.CONTROL_NAME = config["MACHINE_LEARNING_PARAMS"]["CONTROL_NAME"]
    params.TREATMENT_NAME = config["MACHINE_LEARNING_PARAMS"]["TREATMENT_NAME"]
    params.CELL_TYPE = config["MACHINE_LEARNING_PARAMS"]["CELL_TYPE"]
    return params


def data_split(
    X_vals: pd.DataFrame,
    y_vals: pd.Series,
    train_proportion: float = 0.8,
    val_proportion: float = 0.1,
    test_proportion: float = 0.1,
    seed: int = 1,
    params: Parameters = True,
) -> Tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    """Split the Data into train, validation, and test data in both X and Y

    Parameters
    ----------
    X_vals : _type_
        pandas DataFrame of the X values
        the X values are the data being used to predict (an) outcome(s)
    y_vals : _type_
        pandas DataFrame of the Y values
        the Y values are the values being used to validate the model
        the Y values are what we are aiming to predict
    train_proportion : float, optional
        the split of data for the training data, by default 0.8
    val_proportion : float, optional
        the split of data for the validation data, by default 0.1
    test_proportion : float, optional
        the split of data for the testing data, by default 0.1
    seed : int, optional
        the random state set for the random splits
        this is to ensure reproducibility during the development phase, by default 1
    params: Parameters, required
        Parameters class object

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]
        X_train, X_test, X_val, Y_train, Y_test, Y_val
        each training, validation, and testing dataset in both X and Y dimensions properly split

    Raises
    ------
    Exception
        If Training, validation, and testing split ratios do not add up to 1 raise error
        Each ratio is a percentage of the data and thus all should add up to 1
        Ex. training is 0.8, validation is 0.1, and testing is 0.1
    """

    if (train_proportion + test_proportion + val_proportion) != 1:
        raise TrainingValidationTestingSplitError

    if params.MODEL_TYPE == "Regression":
        # split data into train-test
        # TODO add a binning function to be able to properly stratify the continuos predictor feature
        X_train, X_val, Y_train, Y_val = train_test_split(
            X_vals,
            y_vals,
            test_size=val_proportion,
            random_state=seed,
        )
    else:
        X_train, X_val, Y_train, Y_val = train_test_split(
            X_vals, y_vals, test_size=val_proportion, random_state=seed, stratify=y_vals
        )

    # splitting from the train dataset a second time without replacement:
    # we need to adjust the ratio see docstring for more explanation
    test_proportion = test_proportion / (1 - val_proportion)
    if params.MODEL_TYPE == "Regression":
        # split data into train-test
        # TODO add a binning function to be able to properly stratify the continuos predictor feature
        X_train, X_test, Y_train, Y_test = train_test_split(
            X_train,
            Y_train,
            test_size=test_proportion,
            random_state=seed,
        )
    else:
        # split train data into train-validate
        X_train, X_test, Y_train, Y_test = train_test_split(
            X_train,
            Y_train,
            test_size=test_proportion,
            random_state=seed,
            stratify=Y_train,
        )

    # reset the index to avoid downstream errors
    X_train = X_train.reset_index(drop=True)
    Y_train = Y_train.reset_index(drop=True)
    X_val = X_val.reset_index(drop=True)
    Y_val = Y_val.reset_index(drop=True)
    X_test = X_test.reset_index(drop=True)
    Y_test = Y_test.reset_index(drop=True)

    return X_train, X_test, X_val, Y_train, Y_test, Y_val


# Class for x and y data
class Dataset_formatter:
    """
    A class for formatting data for a data loader

    Attributes:
    ----------
    X : Pandas DataFrame
        the X dimension of data (features)
    Y : Pandas DataFrame
        the Y dimension of data (predictor)

    Methods
    -------
    __len__:
        returns the length of the X dimension
    -------
    __getitem__:
        returns a row of the X and Y dimension given an index
    """

    def __init__(
        self,
        X,
        Y,
    ):
        self.X = X
        self.Y = Y

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.Y[idx]


# based on https://www.kaggle.com/code/ludovicocuoghi/pytorch-pytorch-lightning-w-optuna-opt
def build_model_custom(
    trial: optuna.study,
    params: Parameters,
) -> torch.nn.Sequential:
    """Generate a flexible pytorch Neural Network Model that allows for
    optuna hyperparameter optimization

    Parameters
    ----------
    trial : optuna.study
        an iteration of the optuna study object
    params : Parameters
        parameter dataclass object to import constants and params

    Returns
    -------
    torch.nn.Sequential
        this returns in a dict the architecture of the model with optimized parameters
    """

    # number of hidden layers
    # suggest.int takes into account the defined search space
    n_layers = trial.suggest_int("n_layers", params.MIN_LAYERS, params.MAX_LAYERS)

    #  layers will be added to this list and called upon later
    layers = []
    in_features = params.IN_FEATURES

    for i in range(n_layers):
        # the number of units within a hidden layer
        out_features = trial.suggest_int(
            "n_units_l{}".format(i), params.LAYER_UNITS_MIN, params.LAYER_UNITS_MAX
        )

        layers.append(nn.Linear(in_features, out_features))
        # activation function
        layers.append(nn.ReLU())
        # dropout rate
        p = trial.suggest_float(
            (f"dropout_{i}"), params.DROPOUT_MIN, params.DROPOUT_MAX
        )

        layers.append(nn.Dropout(p))
        in_features = out_features

    # final layer append
    layers.append(nn.Linear(in_features, params.OUT_FEATURES))

    # add layers to the model
    return nn.Sequential(*layers)


def train_n_validate(
    model: build_model_custom,
    optimizer: str,
    criterion: nn,
    train_acc: list,
    train_loss: list,
    valid_acc: list,
    valid_loss: list,
    total_step: int,
    total_step_val: int,
    params: Parameters,
    train_loader: torch.utils.data.DataLoader,
    valid_loader: torch.utils.data.DataLoader,
) -> Tuple[build_model_custom, str, object, list, list, list, list, int, int, int, int]:
    """This function trains and validates a machine learning neural network
    the output is used as feedback for optuna hyper parameter optimization


    Parameters
    ----------
    model : build_model_custom
        initialized class containing model
    optimizer : str
        optimizer type
    criterion : nn
        criterion function to be used to calculate loss
    train_acc : list
        list for adding training accuracy values
    train_loss : list
        list for adding training loss values
    valid_acc : list
        list for adding validation accuracy values
    valid_loss : list
        list for adding validation loss values
    total_step : int
        the length of the number of data points in training dataset
    total_step_val : int
        the length of the number of data points in validation dataset
    params : Parameters
        Dataclass containing constants and parameter spaces
    train_loader : torch.utils.data.DataLoader
        DataLoader for train data integration to pytorch
    valid_loader : torch.utils.data.DataLoader
        DataLoader for validation data integration to pytorch


    Returns
    -------
    Tuple[build_model_custom, str, object, list, list, list, list, int, int, int, int]
        model : build_model_custom
            initialized class containing model
        optimizer : str
            optimizer type
        criterion : nn
            criterion function to be used to calculate loss
        train_acc : list
            list for adding training accuracy values
        train_loss : list
            list for adding training loss values
        valid_acc : list
            list for adding validation accuracy values
        valid_loss : list
            list for adding validation loss values
        correct : int
            the number of correctly trained data points in training data
        total : int
            the total number of data points in the train data
        correct_v : int
            the number of correctly validated data points in validation data
        total_v : int
            the total number of correctly validated data points in validation data
    """

    running_loss = 0
    correct = 0
    total = 0

    # TRAINING

    model.train()

    # loop steps through data and optimizes the training weights via
    # finding the loss and accuracy
    # _ isn't used but is needed to "absorb" the data index
    for _, (X_train_batch, y_train_batch) in enumerate(train_loader):
        if params.MODEL_TYPE == "Multi_Class":
            y_train_batch = y_train_batch.type(torch.LongTensor)
        elif params.MODEL_TYPE == "Binary_Classification":
            pass
        elif params.MODEL_TYPE == "Regression":
            pass
        else:
            raise YDataTypeError

        X_train_batch, y_train_batch = X_train_batch.to(
            params.DEVICE
        ), y_train_batch.to(params.DEVICE)
        optimizer.zero_grad()
        output = model(X_train_batch)
        # print(output)
        if params.MODEL_TYPE == "Multi_Class":
            y_pred = torch.log_softmax(output, dim=1)
            _, y_pred = torch.max(y_pred, dim=1)
            loss = criterion(output, y_train_batch)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()  # sum loss for every batch
            correct += (y_pred == y_train_batch).sum().item()
            total += y_train_batch.size(0)
        elif params.MODEL_TYPE == "Binary_Classification":
            y_pred = torch.round(torch.sigmoid(output))
            loss = criterion(output, y_train_batch.unsqueeze(1))
            loss.backward()
            optimizer.step()
            running_loss += loss.item()  # sum loss for every batch
            correct += torch.sum(y_pred == y_train_batch.unsqueeze(1)).item()
            total += y_train_batch.size(0)
        elif params.MODEL_TYPE == "Regression":
            loss = criterion(output, y_train_batch.unsqueeze(1))
            loss.backward()
            optimizer.step()
            running_loss += loss.item()
        else:
            raise ModelTypeError
    if params.MODEL_TYPE == "Regression":
        pass
    else:
        train_acc.append(
            100 * correct / total
        )  # calculate accuracy among all entries in the batches
    train_loss.append(
        running_loss / total_step
    )  # get average loss among all batches dividing total loss by the number of batches

    # VALIDATION
    correct_v = 0
    total_v = 0
    batch_loss = 0
    with torch.no_grad():
        model.eval()
        for _, (X_valid_batch, y_valid_batch) in enumerate(valid_loader):
            if params.MODEL_TYPE == "Multi_Class":
                y_valid_batch = y_valid_batch.type(torch.LongTensor)
            elif params.MODEL_TYPE == "Binary_Classification":
                pass
            elif params.MODEL_TYPE == "Regression":
                pass
            else:
                raise YDataTypeError

            X_valid_batch, y_valid_batch = X_valid_batch.to(
                params.DEVICE
            ), y_valid_batch.to(params.DEVICE)
            # PREDICTION
            output = model(X_valid_batch)

            if params.MODEL_TYPE == "Multi_Class":
                y_pred = torch.log_softmax(output, dim=1)
                _, y_pred = torch.max(y_pred, dim=1)
                loss_v = criterion(output, y_valid_batch)
                batch_loss += loss_v.item()
                correct_v += (y_pred == y_valid_batch).sum().item()
                total_v += y_valid_batch.size(0)
            elif params.MODEL_TYPE == "Binary_Classification":
                y_pred = torch.round(torch.sigmoid(output))
                # LOSS
                loss_v = criterion(output, y_valid_batch.unsqueeze(1))
                batch_loss += loss_v.item()
                # ACCURACY
                correct_v += torch.sum(y_pred == y_valid_batch.unsqueeze(1)).item()
                total_v += y_valid_batch.size(0)
            elif params.MODEL_TYPE == "Regression":
                loss = criterion(output, y_valid_batch.unsqueeze(1))
                batch_loss += loss.item()
            else:
                raise ModelTypeError
        if params.MODEL_TYPE == "Regression":
            pass
        else:
            valid_acc.append(100 * correct_v / total_v)
        valid_loss.append(batch_loss / total_step_val)
    return (
        model,
        optimizer,
        criterion,
        train_acc,
        train_loss,
        valid_acc,
        valid_loss,
        correct,
        total,
        correct_v,
        total_v,
    )


# function for training and tracking model
def objective_model_optimizer(
    train_loader: torch.utils.data.DataLoader,
    valid_loader: torch.utils.data.DataLoader,
    trial: object = optuna.create_study,
    params: Parameters = False,
    metric: str = "loss",
    return_info: bool = False,
) -> str | int:
    """Optimizes the hyperparameter based on search space defined
    returns metrics of how well a given model is doing this is a helper
    function for the optuna optimizer

        Parameters
    ----------
    train_loader : torch.utils.data.DataLoader
        DataLoader for train data integration to pytorch
    valid_loader : torch.utils.data.DataLoader
        DataLoader for validation data integration to pytorch
    trial : trial from optuna
        hyperparameter optimization trial from optuna
    params : optuna.create_study
        Dataclass containing constants and parameter spaces
    metric : str
        metric to be tracked for model optimization
        valid options: 'accuracy' or 'loss'
    return_info : bool, optional
        the option to return more than one metric
        this is best to be False inside of the 'study.optimize' function
        as this function requires only one output metric, by default False


    Returns
    -------
    str | int
        str: returns printed statements of metrics
        when return_info=True
        int: returns metric specified in metric arg
        when return_info=False

    Raises
    ------
    optuna.exceptions.TrialPruned
        exception is raised if/when a trial is pruned due to poor intermediate values
    Exception
        raised when metric value is nor 'accuracy' or 'loss'
    """

    # calling model function
    model = build_model_custom(trial, params)

    # param dictionary for optimization
    optimization_params = {
        "learning_rate": trial.suggest_float(
            "learning_rate", params.LEARNING_RATE_MIN, params.LEARNING_RATE_MAX
        ),
        "optimizer": trial.suggest_categorical("optimizer", params.OPTIMIZER_LIST),
    }

    # param optimizer pick
    optimizer = getattr(optim, optimization_params["optimizer"])(
        model.parameters(), lr=optimization_params["learning_rate"]
    )
    # loss function

    # for binary model use different for multi-class
    if params.MODEL_TYPE == "Multi_Class":
        criterion = nn.CrossEntropyLoss()
    elif params.MODEL_TYPE == "Binary_Classification":
        criterion = nn.BCEWithLogitsLoss()
    elif params.MODEL_TYPE == "Regression":
        criterion = nn.MSELoss()
    else:
        raise ModelTypeError

    # send model to device(cuda)
    model = model.to(params.DEVICE)
    criterion = criterion.to(params.DEVICE)

    # train set accuracy and loss
    train_acc = []
    train_loss = []

    # validation set accuracy and loss
    valid_acc = []
    valid_loss = []

    # total number of data to pass through
    total_step = len(train_loader)
    total_step_val = len(valid_loader)

    for epoch in range(params.OPTIM_EPOCHS):
        (
            model,
            optimizer,
            criterion,
            train_acc,
            train_loss,
            valid_acc,
            valid_loss,
            correct,
            total,
            correct_v,
            total_v,
        ) = train_n_validate(
            model,
            optimizer,
            criterion,
            train_acc,
            train_loss,
            valid_acc,
            valid_loss,
            total_step,
            total_step_val,
            params,
            train_loader,
            valid_loader,
        )
        if metric == "accuracy":
            trial.report(np.mean(valid_acc), epoch)
        elif metric == "loss":
            trial.report(np.mean(valid_loss), epoch)
        else:
            raise OptimizationMetricError

        # Handle pruning based on the intermediate value
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    # I want information returned but only 1 metric required for the optimize function called by study.optimize
    # with out this conditional statement the optimization will fail
    if return_info == True:
        if params.MODEL_TYPE == "Regression":
            print(f"Validation Loss: {np.mean(valid_loss)}")
            print(f"Training Loss: {np.mean(train_loss)}")
            return (
                np.mean(valid_loss),
                np.mean(train_loss),
            )
        else:
            print(f"Validation Accuracy: {np.mean(valid_acc)}")
            print(f"Validation Loss: {np.mean(valid_loss)}")
            print(f"Training Accuracy: {np.mean(train_acc)}")
            print(f"Training Loss: {np.mean(train_loss)}")
            return (
                np.mean(valid_acc),
                np.mean(valid_loss),
                np.mean(train_acc),
                np.mean(train_loss),
            )
    elif metric == "accuracy":
        return np.mean(valid_acc)
    elif metric == "loss":
        return np.mean(valid_loss)
    else:
        raise OptimizationMetricError


def extract_best_trial_params(
    best_params: optuna.study, MLP_params: Parameters, model_name: str
) -> dict:
    """This function extracts the best parameters from the best trial.
    These extracted parameters will be used to train a new model.

    Parameters
    ----------
    best_params : optuna.study.best_params
        returns the best paramaters from the best study from optuna
    MLP_params : Parameters
        dataclass containing constants and parameter spaces
    model_name : str
        name of the model to be created

    Returns
    -------
    dict
        dictionary of all of the params for the best trial model
    """

    units = []
    dropout = []
    n_layers = best_params["n_layers"]
    optimizer = best_params["optimizer"]
    lr = best_params["learning_rate"]
    for i in range(best_params["n_layers"]):
        units.append(best_params[f"n_units_l{i}"])
        dropout.append(best_params[f"dropout_{i}"])
        param_dict = {
            "units": units,
            "dropout": dropout,
            "n_layers": n_layers,
            "optimizer": optimizer,
            "learning_rate": lr,
        }

    # write model architecture to file

    if MLP_params.MODEL_TYPE == "Multi_Class":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Multi_Class/{MLP_params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(architecture_path).mkdir(parents=True, exist_ok=True)
        with open(
            f"{architecture_path}/Multi_Class_{model_name}.json",
            "w",
        ) as f:
            json.dump(param_dict, f, indent=4)
        f.close()

    elif MLP_params.MODEL_TYPE == "Binary_Classification":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Binary_Classification/{MLP_params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(architecture_path).mkdir(parents=True, exist_ok=True)
        with open(
            f"{architecture_path}/Binary_Classification_{model_name}.json",
            "w",
        ) as f:
            json.dump(param_dict, f, indent=4)
        f.close()

    elif MLP_params.MODEL_TYPE == "Regression":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Regression/{MLP_params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(architecture_path).mkdir(parents=True, exist_ok=True)
        with open(
            f"{architecture_path}/Regression_{model_name}.json",
            "w",
        ) as f:
            json.dump(param_dict, f, indent=4)
        f.close()

    else:
        raise ModelTypeError

    return param_dict


# function for new optimized model
def optimized_model_create(
    # parameter_dict: dict,
    params: Parameters,
    model_name: str,
) -> torch.nn.Sequential:
    """creates the pytorch model architecture from the best trial
    from optuna hyperparameter optimization

    Parameters
    ----------
    # parameter_dict : dict
    #     dictionary of optimized model hyperparameters
    params : Parameters
        dataclass to store data of hyperparameters and parameters
    model_name : str
        name of the model to be created

    Returns
    -------
    torch.nn.Sequential
        this returns in a dict the architecture of the model with optimized parameters
    """
    # load in model architecture from saved model architecture
    if params.MODEL_TYPE == "Multi_Class":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Multi_Class/{params.CELL_TYPE}"
        ).resolve(strict=True)
        with open(
            f"{architecture_path}/Multi_Class_{model_name}.json",
            "r",
        ) as f:
            parameter_dict = json.load(f)
        f.close()
    elif params.MODEL_TYPE == "Binary_Classification":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Binary_Classification/{params.CELL_TYPE}"
        ).resolve(strict=True)
        with open(
            f"{architecture_path}/Binary_Classification_{model_name}.json",
            "r",
        ) as f:
            parameter_dict = json.load(f)
        f.close()
    elif params.MODEL_TYPE == "Regression":
        architecture_path = pathlib.Path(
            f"../../trained_models/architectures/Regression/{params.CELL_TYPE}"
        ).resolve(strict=True)
        with open(
            f"{architecture_path}/Regression_{model_name}.json",
            "r",
        ) as f:
            parameter_dict = json.load(f)
        f.close()
    else:
        raise ModelTypeError

    n_layers = parameter_dict["n_layers"]

    layers = []
    in_features = params.IN_FEATURES
    # loop through each layer
    for i in range(n_layers):
        # for each layer access the correct hyper paramter
        out_features = parameter_dict["units"][i]

        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = parameter_dict["dropout"][i]
        layers.append(nn.Dropout(p))
        in_features = out_features
    layers.append(nn.Linear(in_features, params.OUT_FEATURES))
    # output new model to train and test
    return nn.Sequential(*layers)


# Model Training
def train_optimized_model(
    EPOCHS: int,
    train_loader: torch.utils.data.DataLoader,
    valid_loader: torch.utils.data.DataLoader,
    parameter_dict: dict,
    params: Parameters,
    model_name: str,
) -> tuple[float, float, float, float, int, torch.nn.Sequential]:
    """This function trains the optimized model on the training dataset


    Parameters
    ----------
    EPOCHS : int
        The number of epochs to train the model for
    train_loader : torch.utils.data.DataLoader
        DataLoader for train data integration to pytorch
    valid_loader : torch.utils.data.DataLoader
        DataLoader for train data integration to pytorch
    parameter_dict : dict
        dictionary of optimized model hyperparameters
    params : Parameters
        Dataclass containing constants and parameter spaces
    model_name : str
        name of the model to be added to the save name

    Returns
    -------
    tuple[float, float, float, float, int, object]
        train_loss: float
        train_acc: float
        valid_loss: float
        valid_acc: float
        epochs_ran: int
        model: torch.nn.Sequential

    """
    model = optimized_model_create(params, model_name)
    model = model.to(params.DEVICE)
    # criterion is the method in which we measure our loss
    # isn't defined as loss as it doesn't represent the loss value but the method
    if params.MODEL_TYPE == "Multi_Class":
        criterion = nn.CrossEntropyLoss()
    elif params.MODEL_TYPE == "Binary_Classification":
        criterion = nn.BCEWithLogitsLoss()
    elif params.MODEL_TYPE == "Regression":
        criterion = nn.MSELoss()
    else:
        raise ModelTypeError

    optim_method = parameter_dict["optimizer"].strip("'")
    print(optim_method)

    optimizer = f'optim.{optim_method}(model.parameters(), lr={parameter_dict["learning_rate"]})'

    optimizer = eval(optimizer)

    early_stopping_patience = 15
    early_stopping_counter = 0

    train_acc = []
    train_loss = []

    valid_acc = []
    valid_loss = []

    total_step = len(train_loader)
    total_step_val = len(valid_loader)

    valid_loss_min = np.inf

    epochs_ran = []

    for epoch in range(EPOCHS):
        epochs_ran.append(epoch + 1)

        (
            model,
            optimizer,
            criterion,
            train_acc,
            train_loss,
            valid_acc,
            valid_loss,
            correct,
            total,
            correct_v,
            total_v,
        ) = train_n_validate(
            model,
            optimizer,
            criterion,
            train_acc,
            train_loss,
            valid_acc,
            valid_loss,
            total_step,
            total_step_val,
            params,
            train_loader,
            valid_loader,
        )

        if np.mean(valid_loss) <= valid_loss_min:
            if params.MODEL_TYPE == "Multi_Class":
                save_state_path = pathlib.Path(
                    f"../../trained_models/model_save_states/Multi_Class/{params.CELL_TYPE}"
                ).resolve(strict=True)
                pathlib.Path(save_state_path).mkdir(parents=True, exist_ok=True)
                torch.save(
                    model.state_dict(),
                    f"{save_state_path}/Multi_Class_{model_name}.pt",
                )
            elif params.MODEL_TYPE == "Binary_Classification":
                save_state_path = pathlib.Path(
                    f"../../trained_models/model_save_states/Binary_Classification/{params.CELL_TYPE}"
                ).resolve(strict=True)
                pathlib.Path(save_state_path).mkdir(parents=True, exist_ok=True)
                torch.save(
                    model.state_dict(),
                    f"{save_state_path}/Binary_Classification_{model_name}.pt",
                )
            elif params.MODEL_TYPE == "Regression":
                save_state_path = pathlib.Path(
                    f"../../trained_models/model_save_states/Regression/{params.CELL_TYPE}"
                ).resolve(strict=True)
                pathlib.Path(save_state_path).mkdir(parents=True, exist_ok=True)
                torch.save(
                    model.state_dict(),
                    f"{save_state_path}/Regression_{model_name}.pt",
                )
            else:
                raise ModelTypeError

            print(
                f"Epoch {epoch + 0:01}: Validation loss decreased ({valid_loss_min:.6f} --> {np.mean(valid_loss):.6f}).  Saving model ..."
            )
            valid_loss_min = np.mean(valid_loss)
            early_stopping_counter = 0  # reset counter if validation loss decreases
        else:
            print(f"Epoch {epoch + 0:01}: Validation loss did not decrease")
            early_stopping_counter += 1

        if early_stopping_counter > early_stopping_patience:
            print("Early stopped at epoch :", epoch)
            break
        if params.MODEL_TYPE == "Regression":
            pass
        else:
            print(
                f"\t Train_Loss: {np.mean(train_loss):.4f} Train_Acc: {(100 * correct / total):.3f} Val_Loss: {np.mean(valid_loss):.4f}  BEST VAL Loss: {valid_loss_min:.4f}  Val_Acc: {(100 * correct_v / total_v):.3f}\n"
            )
    return train_loss, train_acc, valid_loss, valid_acc, epochs_ran, model


# Plot training data over epochs
def plot_metric_vs_epoch(
    df: pd.DataFrame,
    x: str,
    y1: str,
    y2: str,
    title: str,
    x_axis_label: str,
    y_axis_label: str,
    params: Parameters,
    model_name: str,
    shuffle: bool = False,
) -> None:
    """Plot x vs y1 and x vs y2 using seaborn.

    Parameters
    ----------
    df : pd.DataFrame
        pandas data frame of choice
    x : str
        x variable in this case epoch number
    y1 : str
        y1 variable in this case train metric
    y2 : str
        y1 variable in this case validation metric
    title : str
        string of the title assigned for the given graph
    x_axis_label: str
        label for the x axis
    y_axis_label: str
        label for the y axis
    params : Parameters
        Dataclass containing constants and parameter spaces
    model_name : str
        name of the model to be added to the save name
    shuffle : bool, optional
        whether or not the data was shuffled, by default False
    """
    # sns.lineplot(x=x, y=y1, data=df)
    # sns.lineplot(x=x, y=y2, data=df)
    sns.lineplot(x=df[x], y=df[y1], palette="blue", label="Train")
    sns.lineplot(x=df[x], y=df[y2], palette="orange", label="Validation")
    plt.title(title)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    plt.legend()
    # create graph directory for this model
    graph_path = pathlib.Path(
        f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
    ).resolve(strict=True)
    pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)

    if shuffle:
        graph_path = pathlib.Path(
            f"{graph_path}/{y_axis_label}_graph_shuffled_data.png"
        ).resolve(strict=True)
    elif not shuffle:
        graph_path = pathlib.Path(f"{graph_path}/{y_axis_label}_graph.png").resolve(
            strict=True
        )
    else:
        raise ModelNameError

    plt.tight_layout()
    plt.savefig(graph_path)


def test_optimized_model(
    model: torch.nn.Sequential,
    test_loader: torch.utils.data.DataLoader,
    params: Parameters,
    model_name: str,
) -> Tuple[list, list]:
    """test the trained model on test data

    Parameters
    ----------
    model : torch.nn.Sequential
        pytorch model to us
    test_loader : torch.utils.data.DataLoader
        DataLoader for test data integration to pytorch
    params : Parameters
        Dataclass containing constants and parameter spaces
    model_name : str
        name of the model to be used for loading


    Returns
    -------
    Tuple[list, list]
        y_pred_list: list of predicted values for Y data

        y_pred_prob_list: list of probabilities of
        those predicted values
    """
    model = model.to(params.DEVICE)
    if params.MODEL_TYPE == "Multi_Class":
        save_state_path = pathlib.Path(
            f"../../trained_models/model_save_states/Multi_Class/{params.CELL_TYPE}"
        ).resolve(strict=True)
        model.load_state_dict(
            torch.load(f"{save_state_path}/Multi_Class_{model_name}.pt")
        )
    elif params.MODEL_TYPE == "Binary_Classification":
        save_state_path = pathlib.Path(
            f"../../trained_models/model_save_states/Binary_Classification/{params.CELL_TYPE}"
        ).resolve(strict=True)
        model.load_state_dict(
            torch.load(f"{save_state_path}/Binary_Classification_{model_name}.pt")
        )
    elif params.MODEL_TYPE == "Regression":
        save_state_path = pathlib.Path(
            f"../../trained_models/model_save_states/Regression/{params.CELL_TYPE}"
        ).resolve(strict=True)
        model.load_state_dict(
            torch.load(f"{save_state_path}/Regression_{model_name}.pt")
        )
    else:
        raise ModelTypeError
    y_pred_prob_list = []
    y_pred_list = []

    with torch.no_grad():
        model.eval()
        for _, (X_test_batch, _) in enumerate(test_loader):
            X_test_batch = X_test_batch.to(params.DEVICE)
            # PREDICTION
            output = model(X_test_batch)

            if params.MODEL_TYPE == "Multi_Class":
                _, y_pred = torch.max(output, dim=1)
                y_pred_list.append(y_pred.cpu().numpy())
            elif params.MODEL_TYPE == "Binary_Classification":
                y_pred_prob = torch.sigmoid(output)
                y_pred_prob_list.append(y_pred_prob.cpu().numpy())
                y_pred = torch.round(y_pred_prob)
                y_pred_list.append(y_pred.cpu().numpy())
            elif params.MODEL_TYPE == "Regression":
                y_pred_list.append(output.cpu().numpy())
            else:
                raise ModelTypeError
    y_pred_list = [a.squeeze().tolist() for a in y_pred_list]
    if params.MODEL_TYPE == "Multi_Class":
        return y_pred_list
    elif params.MODEL_TYPE == "Binary_Classification":
        y_pred_prob_list = [a.squeeze().tolist() for a in y_pred_prob_list]
        return y_pred_list, y_pred_prob_list
    elif params.MODEL_TYPE == "Regression":
        return y_pred_list
    else:
        raise ModelTypeError


# If output list is nested
def un_nest(lst) -> list:
    """returns an un-nested list from a nested list

    Parameters
    ----------
    lst : _type_
        a list of lists of any data type

    Returns
    -------
    list: list
        flattened list
    """

    new_lst = []
    for i in lst:
        for j in i:
            new_lst.append(j)
    return new_lst


def results_output(
    prediction_list: list,
    test_data: list,
    params: Parameters,
    prediction_probability_list: list = False,
    test_name: str = "test",
    model_name: str = "model",
    title: str = "Test Results",
    shuffle: bool = False,
) -> None:
    """Function outputs visualization of testing the model

    Parameters
    ----------
    prediction_list : list
        lis of predicted values
    test_data : list
        actual values (true values)
    params: Parameters
        dataclass of parameters
    prediction_probability_list : list, optional
        list of probabilities of 0 or 1, by default False
    test_name : str, optional
        name of the test, by default "test"
    model_name : str, optional
        name of the model, by default "model"
    title : str, optional
        title of the graph, by default "Test Results"
    shuffle : bool, optional
        whether the data was shuffled or not, by default False


    Raises
    ------
    Exception
        raised if proper model type is not specified
    """

    # Classification report

    if params.MODEL_TYPE == "Multi_Class":
        print(classification_report(test_data, prediction_list))
        confusion_matrix_df = pd.DataFrame(confusion_matrix(test_data, prediction_list))
        ax = sns.heatmap(confusion_matrix_df, annot=True, fmt="d")
        ax.invert_xaxis()
        ax.invert_yaxis()
        plt.title(f"Confusion Matrix for Binary Classifier \n {title}", fontsize=20)
        plt.xlabel("Actual Values", size=15)
        plt.ylabel("Predicted Values", size=15)

        # create graph directory for this model
        graph_path = pathlib.Path(
            f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
        if shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/confusion_matrix_graph_{test_name}_shuffled_data.png"
            ).resolve(strict=True)
        elif not shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/confusion_matrix_graph_{test_name}.png"
            ).resolve(strict=True)
        else:
            raise ModelNameError

        plt.tight_layout()
        plt.savefig(graph_path)
        for i in range(params.OUT_FEATURES):
            class_label = (
                i  # the class for which you want to calculate precision and recall
            )
            precision = precision_score(
                test_data, prediction_list, average="weighted", labels=[class_label]
            )
            recall = recall_score(
                test_data, prediction_list, average="weighted", labels=[class_label]
            )

            print(f"Precision for class {class_label}: {precision}")
            print(f"Recall for class {class_label}: {recall}")

        n_classes = (
            params.OUT_FEATURES
        )  # the number of classes in your multi-class problem

        print(n_classes)

        # Binarize the labels for the multiclass ROC calculation
        y_true_binarized = label_binarize(test_data, classes=np.arange(n_classes))
        y_score_binarized = label_binarize(
            prediction_list, classes=np.arange(n_classes)
        )

        # Compute ROC curve and ROC area for each class
        fpr = {}
        tpr = {}
        roc_auc = {}
        for i in range(n_classes):
            fpr[i], tpr[i], _ = roc_curve(
                y_true_binarized[:, i], y_score_binarized[:, i]
            )
            roc_auc[i] = auc(fpr[i], tpr[i])

        # Compute micro-average ROC curve and ROC area
        fpr["micro"], tpr["micro"], _ = roc_curve(
            y_true_binarized.ravel(), y_score_binarized.ravel()
        )
        roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

        # Plot ROC curves for each class
        plt.figure()
        lw = 2
        colors = ["blue", "red", "green", "purple"]
        for i, color in zip(range(n_classes), colors):
            plt.plot(
                fpr[i],
                tpr[i],
                color=color,
                lw=lw,
                label="ROC curve (AUC = %0.2f) for class %d" % (roc_auc[i], i),
            )

        # Plot micro-average ROC curve
        plt.plot(
            fpr["micro"],
            tpr["micro"],
            label="Micro-average ROC curve (AUC = {0:0.2f})".format(roc_auc["micro"]),
            color="deeppink",
            linestyle=":",
            linewidth=4,
        )

        plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"Receiver Operating Characteristic (ROC) Curve \n {title}")
        plt.legend(loc="lower right")
        # create graph directory for this model
        graph_path = pathlib.Path(
            f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
        if shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}_shuffled_data.png"
            ).resolve(strict=True)
        elif not shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}.png"
            ).resolve(strict=True)
        else:
            raise ModelNameError

        plt.tight_layout()
        plt.savefig(graph_path)
        plt.show()
        return confusion_matrix_df
    elif params.MODEL_TYPE == "Binary_Classification":
        print(classification_report(test_data, prediction_list))
        confusion_matrix_df = pd.DataFrame(confusion_matrix(test_data, prediction_list))
        confusion_matrix_df.columns = ["Negative", "Positive"]
        confusion_matrix_df.index = ["Negative", "Positive"]
        ax = sns.heatmap(confusion_matrix_df, annot=True, fmt="d")
        ax.invert_xaxis()
        ax.invert_yaxis()
        plt.title(f"Confusion Matrix for Binary Classifier \n {title}", fontsize=20)
        plt.xlabel("Actual Values", size=15)
        plt.ylabel("Predicted Values", size=15)
        # create graph directory for this model
        graph_path = pathlib.Path(
            f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
        if shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/confusion_matrix_graph_{test_name}_shuffled_data.png"
            ).resolve(strict=True)
        elif not shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/confusion_matrix_graph_{test_name}.png"
            ).resolve(strict=True)
        else:
            raise ModelNameError

        plt.tight_layout()
        plt.savefig(graph_path)
        # AUC graph of accuracy and false positive rates
        plt.figure(figsize=(5.5, 4))
        fpr, tpr, _ = roc_curve(test_data, prediction_probability_list)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, "b", label="AUC = %0.2f" % roc_auc)
        plt.plot([0, 1], [0, 1], "r--")
        plt.title(f"ROC curve \n {title}", fontsize=25)
        plt.ylabel("True Positive Rate", fontsize=18)
        plt.xlabel("False Positive Rate", fontsize=18)
        plt.legend(
            loc="lower right",
            fontsize=24,
            fancybox=True,
            shadow=True,
            frameon=True,
            handlelength=0,
        )
        # create graph directory for this model
        graph_path = pathlib.Path(
            f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
        if shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}_shuffled_data.png"
            ).resolve(strict=True)
        elif not shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}.png"
            ).resolve(strict=True)
        else:
            raise ModelNameError

        plt.tight_layout()
        plt.savefig(graph_path)
        plt.show()

    elif params.MODEL_TYPE == "Regression":
        mse = mean_squared_error(test_data, prediction_list)
        r_square = r2_score(test_data, prediction_list)
        print("Mean Squared Error :", mse)
        print("R^2 :", r_square)
        a, b = np.polyfit(test_data, prediction_list, 1)
        # plt.figure(figsize=(5.5, 4))
        plt.scatter(test_data, prediction_list, s=25)
        plt.plot(
            test_data,
            a * test_data + b,
            color="red",
            label="R2={0:0.2f}".format(r_square),
        )
        plt.title(
            f"Regression Nerual Network Prediction vs. True \n {title}", fontsize=25
        )
        plt.ylabel("Predicted", fontsize=18)
        plt.xlabel("Target", fontsize=18)

        plt.legend(
            loc="upper left",
            fontsize=24,
            fancybox=True,
            shadow=True,
            frameon=True,
            handlelength=0,
        )
        # create graph directory for this model
        graph_path = pathlib.Path(
            f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/{params.CELL_TYPE}"
        ).resolve(strict=True)
        pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
        if shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}_shuffled_data.png"
            ).resolve(strict=True)
        elif not shuffle:
            graph_path = pathlib.Path(
                f"{graph_path}/ROC_graph_{test_name}.png"
            ).resolve(strict=True)
        else:
            raise ModelNameError

        plt.tight_layout()
        plt.savefig(graph_path)
        plt.show()
    else:
        raise ModelTypeError
