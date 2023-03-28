"""
Functions for Machine Learning Model optimization, training and testing.
These are helper functions meant to be called in a separate notebook or script
"""


import ast
from configparser import ConfigParser
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import seaborn as sns
import torch
import torch.nn as nn
import torch.optim as optim
from MLP_utils.parameters import Parameters
from sklearn.metrics import (
    accuracy_score,
    auc,
    classification_report,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split

config = ConfigParser()
config.optionxform = str
config.read("MLP_utils/config.ini")

params = Parameters()


def parameter_set(params: Parameters, config: object) -> object:
    """reset parameter class defaults by updating from config

    Parameters
    ----------
    params : object
        param class holding parameter information
    config: object
        config class

    Returns
    -------
    object
        param class holding updated parameter information
    """
    params.DATA_SUBSET_OPTION = config["DEFAULT"]["DATA_SUBSET_OPTION"]
    params.DATA_SUBSET_NUMBER = int(config["DEFAULT"]["DATA_SUBSET_NUMBER"])
    params.BATCH_SIZE = int(config["DEFAULT"]["BATCH_SIZE"])
    params.OPTIM_EPOCHS = int(config["DEFAULT"]["OPTIM_EPOCHS"])
    params.N_TRIALS = int(config["DEFAULT"]["N_TRIALS"])
    params.TRAIN_EPOCHS = int(config["DEFAULT"]["TRAIN_EPOCHS"])
    params.MIN_LAYERS = int(config["DEFAULT"]["MIN_LAYERS"])
    params.MAX_LAYERS = int(config["DEFAULT"]["MAX_LAYERS"])
    params.LAYER_UNITS_MIN = int(config["DEFAULT"]["LAYER_UNITS_MIN"])
    params.LAYER_UNITS_MAX = int(config["DEFAULT"]["LAYER_UNITS_MAX"])
    params.DROPOUT_MIN = float(config["DEFAULT"]["DROPOUT_MIN"])
    params.DROPOUT_MAX = float(config["DEFAULT"]["DROPOUT_MAX"])
    params.DROP_OUT_STEP = float(config["DEFAULT"]["DROP_OUT_STEP"])
    params.LEARNING_RATE_MIN = float(config["DEFAULT"]["LEARNING_RATE_MIN"])
    params.LEARNING_RATE_MAX = float(config["DEFAULT"]["LEARNING_RATE_MAX"])
    params.OPTIMIZER_LIST = config["DEFAULT"]["OPTIMIZER_LIST"]
    params.METRIC = config["DEFAULT"]["METRIC"]
    params.DIRECTION = config["DEFAULT"]["DIRECTION"]
    return params


params = parameter_set(params, config)


def data_split(
    X_vals: pd.DataFrame,
    y_vals: pd.Series,
    train_proportion: float = 0.8,
    val_proportion: float = 0.1,
    test_proportion: float = 0.1,
    seed: int = 1,
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
        raise Exception("Train, val, and test sets must add to 1")

    # split data into train-test
    X_train, X_val, Y_train, Y_val = train_test_split(
        X_vals, y_vals, test_size=val_proportion, random_state=seed, stratify=y_vals
    )

    # splitting from the train dataset a second time without replacement:
    # we need to adjust the ratio see docstring for more explanation
    test_proportion = test_proportion / (1 - val_proportion)
    # split train data into train-validate
    X_train, X_test, Y_train, Y_test = train_test_split(
        X_train, Y_train, test_size=test_proportion, random_state=seed, stratify=Y_train
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
    trial: optuna.study, in_features: int, final_layer_out_features: int, params: object
) -> torch.nn.Sequential:
    """Generate a flexible pytorch Neural Network Model that allows for
    optuna hyperparameter optimization

    Parameters
    ----------
    trial : optuna.study
        an iteration of the optuna study object
    in_features : int
        the number of in features for the initial network layer
    final_layer_out_features : int
        the number of out features for the network final layer
    params : object
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
            ("dropout_{}".format(i)), params.DROPOUT_MIN, params.DROPOUT_MAX
        )

        layers.append(nn.Dropout(p))
        in_features = out_features

    # final layer append
    layers.append(nn.Linear(in_features, final_layer_out_features))

    # add layers to the model
    return nn.Sequential(*layers)


def train_n_validate(
    model: build_model_custom,
    optimizer: str,
    criterion: object,
    train_acc: list,
    train_loss: list,
    valid_acc: list,
    valid_loss: list,
    total_step: int,
    total_step_val: int,
    params: object,
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
    criterion : object
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
    params : object
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
        criterion : object
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
        X_train_batch, y_train_batch = X_train_batch.to(
            params.DEVICE
        ), y_train_batch.to(params.DEVICE)
        optimizer.zero_grad()
        output = model(X_train_batch)
        y_pred = torch.round(torch.sigmoid(output))
        # LOSS
        loss = criterion(output, y_train_batch.unsqueeze(1))
        loss.backward()
        optimizer.step()
        running_loss += loss.item()  # sum loss for every batch
        # ACCURACY
        correct += torch.sum(y_pred == y_train_batch.unsqueeze(1)).item()
        total += y_train_batch.size(0)
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
            X_valid_batch, y_valid_batch = X_valid_batch.to(
                params.DEVICE
            ), y_valid_batch.to(params.DEVICE)
            # PREDICTION
            output = model(X_valid_batch)
            y_pred = torch.round(torch.sigmoid(output))
            # LOSS
            loss_v = criterion(output, y_valid_batch.unsqueeze(1))
            batch_loss += loss_v.item()
            # ACCURACY
            correct_v += torch.sum(y_pred == y_valid_batch.unsqueeze(1)).item()
            total_v += y_valid_batch.size(0)
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
    trial: object = object,
    in_features: int = 1,
    out_features: int = 1,
    params: object = params,
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
    in_features : int
        the number of in features for the initial network layer
    out_features : int
        the number of out features for the final network layer
    params : object
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
    model = build_model_custom(trial, in_features, out_features, params)

    # param dictionary for optimization
    optimization_params = {
        "learning_rate": trial.suggest_float(
            "learning_rate", params.LEARNING_RATE_MIN, params.LEARNING_RATE_MAX
        ),
        "optimizer": trial.suggest_categorical(
            "optimizer", ast.literal_eval(params.OPTIMIZER_LIST)
        ),
    }

    # param optimizer pick
    optimizer = getattr(optim, optimization_params["optimizer"])(
        model.parameters(), lr=optimization_params["learning_rate"]
    )
    # loss function

    # for binary model use different for multi-class
    criterion = nn.BCEWithLogitsLoss()

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

        trial.report(np.mean(valid_loss), epoch)

        # Handle pruning based on the intermediate value
        if trial.should_prune():
            raise optuna.exceptions.TrialPruned()

    # I want information returned but only 1 metric required for the optimize function called by study.optimize
    # with out this conditional statement the optimization will fail
    if return_info == True:
        print(f"Validation Accuracy: {np.mean(valid_acc)}")
        print(f"Validation Loss: {np.mean(valid_loss)}")
        print(f"Training Accuracy: {np.mean(train_acc)}")
        print(f"Training Loss: {np.mean(train_loss)}")
    elif metric == "accuracy":
        return np.mean(valid_acc)
    elif metric == "loss":
        return np.mean(valid_loss)
    else:
        raise Exception("metric not defined as accuracy or loss")


def extract_best_trial_params(best_params: optuna.study) -> dict:
    """This function extracts the best parameters from the best trial.
    These extracted parameters will be used to train a new model.

    Parameters
    ----------
    best_params : optuna.study.best_params
        returns the best parmaters from the best study from optuna

    Returns
    -------
    dict
        dictionary of all of the params for the best trial model
    """

    best_params
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
    return param_dict


# function for new optimized model
def optimized_model_create(
    in_features: int, final_layer_out_features: int, parameter_dict: dict
) -> torch.nn.Sequential:
    """creates the pytorch model architecture from the best trial
    from optuna hyperparameter optimization

    Parameters
    ----------
    in_features : int
        the number of in features for the initial network layer
    final_layer_out_features : int
        the number of out features for the network final layer
    parameter_dict : dict
        dictionary of optimized model hyperparameters

    Returns
    -------
    torch.nn.Sequential
        this returns in a dict the architecture of the model with optimized parameters
    """

    n_layers = parameter_dict["n_layers"]

    layers = []

    # loop through each layer
    for i in range(n_layers):

        # for each layer access the correct hyper paramter
        out_features = parameter_dict["units"][i]

        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = parameter_dict["dropout"][i]
        layers.append(nn.Dropout(p))
        in_features = out_features
    layers.append(nn.Linear(out_features, final_layer_out_features))
    # output new model to train and test
    return nn.Sequential(*layers)


# Model Training
def train_optimized_model(
    EPOCHS: int,
    train_loader: torch.utils.data.DataLoader,
    valid_loader: torch.utils.data.DataLoader,
    in_features: int,
    out_features: int,
    parameter_dict: dict,
    params: object,
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
    in_features : int
        the number of in features for the initial network layer
    out_features : int
        the number of out features for the final network layer
    parameter_dict : dict
        dictionary of optimized model hyperparameters
    params : object
        Dataclass containing constants and parameter spaces

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

    model = optimized_model_create(in_features, out_features, parameter_dict)
    model = model.to(params.DEVICE)
    # criterion is the method in which we measure our loss
    # isn't defined as loss as it doesn't represent the loss value but the method
    criterion = nn.BCEWithLogitsLoss()
    optim_method = parameter_dict["optimizer"].strip("'")
    print(optim_method)

    optimizer = f'optim.{optim_method}(model.parameters(), lr={parameter_dict["learning_rate"]})'

    # optimizer = (
    # f'optim.{optim_method}(model.parameters(), lr=parameter_dict["learning_rate"])'
    # )
    # print(optimizer)
    optimizer = eval(optimizer)
    # optimizer = ast.literal_eval(optimizer)
    print(optimizer)

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
            torch.save(model.state_dict(), f"./state_dict.pt")
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

        print(
            f"\t Train_Loss: {np.mean(train_loss):.4f} Train_Acc: {(100 * correct / total):.3f} Val_Loss: {np.mean(valid_loss):.4f}  BEST VAL Loss: {valid_loss_min:.4f}  Val_Acc: {(100 * correct_v / total_v):.3f}\n"
        )
    return train_loss, train_acc, valid_loss, valid_acc, epochs_ran, model


# Plot training data over epochs
def plot_metric_vs_epoch(df: pd.DataFrame, x: str, y1: str, y2: str) -> None:
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
    """
    sns.lineplot(x=x, y=y1, data=df)
    sns.lineplot(x=x, y=y2, data=df)


def test_optimized_model(
    model: torch.nn.Sequential,
    test_loader: torch.utils.data.DataLoader,
    in_features: int,
    out_features: int,
    parameter_dict: dict,
    params: object,
) -> Tuple[list, list]:
    """test the trained model on test data

    Parameters
    ----------
    model : torch.nn.Sequential
        pytorch model to us
    test_loader : torch.utils.data.DataLoader
        DataLoader for test data integration to pytorch
    in_features : int
        the number of in features for the initial network layer
    out_features : int
        the number of out features for the final network layer
    parameter_dict : dict
        dictionary of optimized model hyperparameters
    params : object
        Dataclass containing constants and parameter spaces


    Returns
    -------
    Tuple[list, list]
        y_pred_list: list of predicted values for Y data

        y_pred_prob_list: list of probabilities of
        those predicted values
    """

    model = optimized_model_create(in_features, out_features, parameter_dict)
    model = model.to(params.DEVICE)
    criterion = nn.BCEWithLogitsLoss()
    optim_method = parameter_dict["optimizer"].strip("'")
    optimizer = (
        f'optim.{optim_method}(model.parameters(), lr=parameter_dict["learning_rate"])'
    )
    optimizer = eval(optimizer)

    y_pred_prob_list = []
    y_pred_list = []

    with torch.no_grad():
        model.eval()
        for _, (X_test_batch, y_test_batch) in enumerate(test_loader):
            X_test_batch = X_test_batch.to(params.DEVICE)
            # PREDICTION
            output = model(X_test_batch)
            y_pred_prob = torch.sigmoid(output)
            y_pred_prob_list.append(y_pred_prob.cpu().numpy())
            y_pred = torch.round(y_pred_prob)
            y_pred_list.append(y_pred.cpu().numpy())
    y_pred_prob_list = [a.squeeze().tolist() for a in y_pred_prob_list]
    y_pred_list = [a.squeeze().tolist() for a in y_pred_list]

    return y_pred_list, y_pred_prob_list


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
    prediction_list: list, prediction_probability_list: list, test_data: list
) -> None:
    """Function outputs visualization of testing the model


    Parameters
    ----------
    prediction_list : list
        list of predicted values from model
    prediction_probability_list : list
        list of probabilities of predicted values from model
    test_data : list
        input data to model actual labels

    Return:
        classification report
        confusion matrix
        AUC graph of accuracy and false positive rates
    """

    # Classification report
    print(classification_report(test_data, prediction_list))

    # confusion matrix
    confusion_matrix(test_data, prediction_list)

    # AUC graph of accuracy and false positive rates
    plt.figure(figsize=(5.5, 4))
    fpr, tpr, _ = roc_curve(test_data, prediction_probability_list)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, "b", label="AUC = %0.2f" % roc_auc)
    plt.plot([0, 1], [0, 1], "r--")
    plt.title("ROC curve", fontsize=25)
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
    plt.show()
