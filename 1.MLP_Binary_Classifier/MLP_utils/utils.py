"""
Functions for Machine Learning Model optimization, training and testing.
These are helper functions meant to be called in a separate notebook or script
"""


from configparser import ConfigParser
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import optuna
import seaborn as sns
import torch
import torch.nn as nn
import torch.optim as optim
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
config.read("Interstellar_Analysis/1.MLP_Binary_Classifier/MLP_utils/config.ini")
config["DEFAULT"]


# Set parameters for model optimization and training/testing
params = Parameters(
    DATA_SUBSET_OPTION=config["DEFAULT"]["DATA_SUBSET_OPTION"],
    DATA_SUBSET_NUMBER=config["DEFAULT"]["DATA_SUBSET_NUMBER"],
    BATCH_SIZE=config["DEFAULT"]["BATCH_SIZE"],
    OPTIM_EPOCHS=config["DEFAULT"]["OPTIM_EPOCHS"],
    N_TRIALS=config["DEFAULT"]["N_TRIALS"],
    TRAIN_EPOCHS=config["DEFAULT"]["TRAIN_EPOCHS"],
    MIN_LAYERS=config["DEFAULT"]["MIN_LAYERS"],
    MAX_LAYERS=config["DEFAULT"]["MAX_LAYERS"],
    LAYER_UNITS_MIN=config["DEFAULT"]["LAYER_UNITS_MIN"],
    LAYER_UNITS_MAX=config["DEFAULT"]["LAYER_UNITS_MAX"],
    DROPOUT_MIN=config["DEFAULT"]["DROPOUT_MIN"],
    DROPOUT_MAX=config["DEFAULT"]["DROPOUT_MAX"],
    DROP_OUT_STEP=config["DEFAULT"]["DROP_OUT_STEP"],
    LEARNING_RATE_MIN=config["DEFAULT"]["LEARNING_RATE_MIN"],
    LEARNING_RATE_MAX=config["DEFAULT"]["LEARNING_RATE_MAX"],
    OPTIMIZER_LIST=config["DEFAULT"]["OPTIMIZER_LIST"],
)


def pandas_reset_index(df):
    df = df.reset_index(drop=True)
    return df


def data_split(X_vals, y_vals, train, val, test, seed=1):
    """
    This function splits data into
    training, validation and testing data
    Note that ratios of splits must be adjusted to achieve true percentage splits

    Example: 80:10:10 : Train:Validation:Testing split
    after taking 10% of the data for validation remaining is 90% of original data

    Now taking 10% of the remaining 90% does not yield 10% of the original data:
    We must adjust the ratios by diving the 10% wanted split from the original data
    by the adjusted 1 - 0.1 = 0.9 ratio of data remaining

    thus finally we have 0.1 / (1 - 0.1) = 0.11111111111111112

    So ~11% of the remaining 90% of data yields 10% of the original data. Cheers Mates.

    Parameters:
        x_vals
        y_vals
        train
        val
        test
        params
    """

    if train + test + val != 1:
        return

    # split data into train-test
    X_train, X_val, Y_train, Y_val = train_test_split(
        X_vals, y_vals, test_size=val, random_state=seed, stratify=y_vals
    )

    # splitting from the train dataset a second time without replacement:
    # we need to adjust the ratio see docstring for more explanation
    test = test / (1 - val)
    # split train data into train-validate
    X_train, X_test, Y_train, Y_test = train_test_split(
        X_train, Y_train, test_size=test, random_state=seed, stratify=Y_train
    )

    # reset the index to avoid downstream errors
    for i in [X_train, Y_train, X_val, Y_val, X_test, Y_test]:
        i = pandas_reset_index(i)
    X_train = X_train.reset_index(drop=True)
    Y_train = Y_train.reset_index(drop=True)
    X_test = X_test.reset_index(drop=True)
    Y_test = Y_test.reset_index(drop=True)
    X_val = X_val.reset_index(drop=True)
    Y_val = Y_val.reset_index(drop=True)

    return X_train, X_test, X_val, Y_train, Y_test, Y_val


# Data class for x and y data
class Dataset:
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
def build_model_custom(trial, in_features, final_layer_out_features, params):
    """
    This function lays out the general architecture of a Neural Network.
    There are variables throughout to optimize the hyperparameter of this model.
    This function is meant to be used with optuna to optimize functions.

    Parameters:
        trial : optuna object
            an optuna object foe which optimization trial to input for
            what parameters to use in the defined search space
        in_features : int
            the number of input features to define the shape of the model

    Return:
        nn.Sequential(*layers) : dict
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
            "n_units_l{}".format(i), params.LAYER_UNITS_MIN, params.AYER_UNITS_MAX
        )

        layers.append(nn.Linear(in_features, out_features))
        # activation function
        layers.append(nn.ReLU())

        # dropout rate
        p = trial.suggest_float(
            "dropout_{}".format(i),
            params.DROPOUT_MIN,
            params.DROPOUT_MAX,
            params.DROP_OUT_STEP,
        )
        layers.append(nn.Dropout(p))
        in_features = out_features

    # final layer append
    layers.append(nn.Linear(in_features, final_layer_out_features))

    # add layers to the model
    return nn.Sequential(*layers)


def train_n_validate(
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
):
    """_summary_

    Parameters
    ----------
    model : _type_
        _description_
    optimizer : _type_
        _description_
    criterion : _type_
        _description_
    train_acc : _type_
        _description_
    train_loss : _type_
        _description_
    valid_acc : _type_
        _description_
    valid_loss : _type_
        _description_
    total_step : _type_
        _description_
    total_step_val : _type_
        _description_
    params : _type_
        _description_
    train_loader : _type_
        _description_
    valid_loader : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    running_loss = 0
    correct = 0
    total = 0

    # TRAINING

    model.train()

    # loop steps through data and optimizes the training weights via
    # finding the loss and accuracy
    for (X_train_batch, y_train_batch) in enumerate(train_loader):
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
        for (X_valid_batch, y_valid_batch) in enumerate(valid_loader):
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
def objective(
    trial, in_features, train_loader, valid_loader, params, return_info=False
):
    """
    This function trains the model and tests it on validation data.
    The accuracy and loss output is how the success of the model is tracked.

    Parameters:
        trial : optuna object

        in_features : int
            the number of in features for a given dataset

        train_loader : Object
            object from torch an iterable dataclass

        return_info : bool
            If set to False only one metric will be returned to optimize the model
            If set to True multiple metrics will be returned via printing
            This is required as the optimization of the model is tracked by one output metric
            and returning more than one metric will cause the optimization to fail but after optimization
            it is nice to know about what the other output metrics are

    Return:
        metric(s) : str or float
            if return_info == True:
                Mean Validation Accuracy
                Mean Validation Loss
                Mean Training Accuracy
                Mean Training Loss
            if return_info == False:
                return the mean validation accuracy

    """

    # calling model function
    model = build_model_custom(trial, in_features)

    # param dictionary for optimization
    optimization_params = {
        "learning_rate": trial.suggest_float("learning_rate", 1e-5, 1e-1),
        "optimizer": trial.suggest_categorical("optimizer", ["Adam", "RMSprop", "SGD"]),
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
    else:
        return np.mean(valid_loss)


def extract_best_trial_params(best_params):
    """
    This function extracts the best parameters from the best trial.

    These extracted parameters will be used to train a new model.

    Parameters:
        best_params : obj
            study.best_params
            the best_params function of the study object
    Return:
        param_dict : dict
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
def optimized_model(in_features, parameter_dict):
    """
    This function uses the extracted optimized functions to create a new model

    Parameters:
        in_features : int
            this is the number of in features to used for the model
        parameter_dict : dict
            this is a dictionary returned from the extract_best_trial_params function

    Return:
        nn.Sequential(*layers) : dict
            this returns in a dict the architecture of the model with optimized parameters
    """
    n_layers = parameter_dict["n_layers"]

    layers = []

    for i in range(n_layers):

        out_features = parameter_dict["units"][i]

        layers.append(nn.Linear(in_features, out_features))
        layers.append(nn.ReLU())
        p = parameter_dict["dropout"][i]
        layers.append(nn.Dropout(p))
        in_features = out_features
    layers.append(nn.Linear(out_features, 1))

    return nn.Sequential(*layers)


# Model Training
def train_optimized_model(
    EPOCHS,
    train_loader,
    valid_loader,
    in_features,
    parameter_dict,
    params,
    return_info=False,
):
    """
    This function trains the optimized model on the training dataset

    Parameters:
        EPOCHS : int
            the number of epochs to train the model for
    Return:
        training metrics : str

    """

    model = optimized_model(in_features, parameter_dict)
    model = model.to(params.DEVICE)
    # criterion is the method in which we measure our loss
    # isn't defined as loss as it doesn't represent the loss value but the method
    criterion = nn.BCEWithLogitsLoss()
    optim_method = parameter_dict["optimizer"].strip("'")
    optimizer = (
        f'optim.{optim_method}(model.parameters(), lr=parameter_dict["learning_rate"])'
    )
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
            torch.save(model.state_dict(), f"./{date.today()}_state_dict.pt")
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
def plot_metric_vs_epoch(df, x, y1, y2):

    sns.lineplot(x=x, y=y1, data=df)
    sns.lineplot(x=x, y=y2, data=df)


def test_optimized_model(model, test_loader, in_features, parameter_dict, params):
    """
    This function tests the trained optimized model on test dataset

    Parameters:
        None
    Return:
        y_pred_list : list
            list of predicted values
        y_pred_prob_list : list
            list of probabilities for predicted values
    """

    model = optimized_model(in_features, parameter_dict)
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
        for batch_idx, (X_test_batch, y_test_batch) in enumerate(test_loader):
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
def un_nest(lst):
    """
    returns an un-nested list from a nested list

    Parameters:
        lst : list
            a list of lists
    """
    new_lst = []
    for i in lst:
        for j in i:
            new_lst.append(j)
    return new_lst


def results_output(prediction_list, prediction_probability_list, test_data):
    """
    Function outputs visualization of testing the model

    Parameters:
        prediction_list : list of predicted values from model
        prediction_probability_list : list of probabilities of predicted values from model
        test_data : input data to model

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
