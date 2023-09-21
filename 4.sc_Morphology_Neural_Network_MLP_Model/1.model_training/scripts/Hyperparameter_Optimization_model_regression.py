#!/usr/bin/env python
# coding: utf-8
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Hyperparameter tuning via Optuna

# %% [markdown]
# %% [markdown]
# %% [markdown]
# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.
# ### Here I will be predicting if a cell received a treatment or not

# %%
import pathlib
import sys

import numpy as np
import optuna
import pandas as pd
import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

sys.path.append("../..")

from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    data_split,
    extract_best_trial_params,
    objective_model_optimizer,
    parameter_set,
    plot_metric_vs_epoch,
    results_output,
    test_optimized_model,
    train_optimized_model,
    un_nest,
)

sys.path.append("../../..")
from utils.utils import df_stats

# %% [markdown]
# #### Set up Data to be compatible with model

# %% [markdown]
# ##### Regression Model Data Wrangling and Set Up
# comment out if not using regression

# %%
# Parameters
CELL_TYPE = "SHSY5Y"
CONTROL_NAME = "DMSO_0.100_DMSO_0.025"
TREATMENT_NAME = "LPS_100.000_DMSO_0.025"
MODEL_NAME = "DMSO_0.025_vs_LPS_100"
SHUFFLE = False

# %%
ml_configs_file = pathlib.Path("../../MLP_utils/regression_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)

# overwrite params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
mlp_params.CONTROL_NAME = CONTROL_NAME
mlp_params.TREATMENT_NAME = TREATMENT_NAME
mlp_params.SHUFFLE = SHUFFLE

# %%
# Import Data
# set data file path under pathlib path for multi-system use

file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)

# set path for nomic data
nomic_df_path = pathlib.Path(
    f"../../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{mlp_params.CELL_TYPE}_cleanup4correlation.csv"
).resolve(strict=True)

df = pq.read_table(file_path).to_pandas()
nomic_df = pd.read_csv(nomic_df_path)

# %%
print(df.shape)
df = pd.merge(
    df,
    nomic_df,
    left_on=[
        "Metadata_Well",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ],
    right_on=[
        "Metadata_position_x",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ],
).drop(["Metadata_position_x"], axis=1)
print(nomic_df.shape)
print(df.shape)

# %%
# filter the oneb_Metadata_Treatment_Dose_Inhibitor_Dose column to only include the treatment and control via loc
df = df.loc[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        [mlp_params.TREATMENT_NAME, mlp_params.CONTROL_NAME]
    )
]


print("Selected Catagories are:")
print(df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique())
df_stats(df)

if mlp_params.DATA_SUBSET_OPTION == "True":
    df = df.sample(n=mlp_params.DATA_SUBSET_NUMBER)
    print("Data Subset Is On")
    print(f"Data is subset to {mlp_params.DATA_SUBSET_NUMBER}")
else:
    print("Data Subset Is Off")

# %%
np.random.seed(seed=0)
wells_to_hold = (
    df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose")
    .agg(np.random.choice)["Metadata_Well"]
    .to_list()
)
df_holdout = df[df["Metadata_Well"].isin(wells_to_hold)]
df = df[~df["Metadata_Well"].isin(wells_to_hold)]


print("Wells held out for testing:", df_holdout["Metadata_Well"].unique())
print(
    "Wells to use for training, validation, and testing", df["Metadata_Well"].unique()
)

# %%
# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)

# %%
# get all columns that contain NSU in the name
df_values_Y = df_values[df_values.columns[df_values.columns.str.contains("NSU")]]
df_values_X = df_values.drop(columns=df_values_Y.columns)
df_values_Y["Metadata_Well"] = df_descriptive["Metadata_Well"]
print(df_values.shape)
print(df_values_X.shape)
print(df_values_Y.shape)

# %% [markdown]
# #### Split Data - All Models can proceed through this point

# %%
X_train, X_test, X_val, Y_train_well, Y_test_well, Y_val_well = data_split(
    X_vals=df_values_X,
    y_vals=df_values_Y,
    train_proportion=0.8,
    val_proportion=0.1,
    test_proportion=0.1,
    seed=0,
    params=mlp_params,
)

# %%
Y_train = Y_train_well.drop(columns=["Metadata_Well"])
Y_test = Y_test_well.drop(columns=["Metadata_Well"])
Y_val = Y_val_well.drop(columns=["Metadata_Well"])

# %%
# produce data objects for train, val and test datasets
train_data = Dataset_formatter(
    torch.FloatTensor(X_train.values), torch.FloatTensor(Y_train.values)
)
val_data = Dataset_formatter(
    torch.FloatTensor(X_val.values), torch.FloatTensor(Y_val.values)
)
test_data = Dataset_formatter(
    torch.FloatTensor(X_test.values), torch.FloatTensor(Y_test.values)
)

# %%
mlp_params.IN_FEATURES = X_train.shape[1]
print("Number of in features: ", mlp_params.IN_FEATURES)
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.OUT_FEATURES = Y_train.shape[1]
else:
    mlp_params.OUT_FEATURES = len(
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
    )

print("Number of out features: ", mlp_params.OUT_FEATURES)

# %%
# convert data class into a dataloader to be compatible with pytorch
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=mlp_params.BATCH_SIZE
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=mlp_params.BATCH_SIZE
)
test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1)

# %%
df_values_X.shape
df_values_Y.shape

# %%
# no accuracy function must be loss for regression
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.METRIC = "loss"
    pass

sampler = optuna.samplers.TPESampler(seed=0)

# wrap the objective function inside of a lambda function to pass args...
objective_lambda_func = lambda trial: objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=trial,
    params=mlp_params,
    metric=mlp_params.METRIC,
    return_info=False,
)


# Study is the object for model optimization
study = optuna.create_study(direction=f"{mlp_params.DIRECTION}", sampler=sampler)
# Here I apply the optimize function of the study to the objective function
# This optimizes each parameter specified to be optimized from the defined search space
study.optimize(objective_lambda_func, n_trials=mlp_params.N_TRIALS)
# Prints out the best trial's optimized parameters
objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=study.best_trial,
    params=mlp_params,
    metric=mlp_params.METRIC,
    return_info=True,
)

# %%
fig = optuna.visualization.plot_optimization_history(study)
graph_path = pathlib.Path(f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/")
# if path doesn't exist, make path with pathlib
graph_path.mkdir(parents=True, exist_ok=True)

graph_path = f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/plot_optimization_history_graph"
fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()

# %%
fig = optuna.visualization.plot_intermediate_values(study)
graph_path = pathlib.Path(f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/")
# if path doesn't exist, make path with pathlib
graph_path.mkdir(parents=True, exist_ok=True)

graph_path = f"../../figures/{params.MODEL_TYPE}/{params.MODEL_NAME}/plot_intermediate_values_graph"
fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()

# %%
param_dict = extract_best_trial_params(
    study.best_params, params, model_name=params.MODEL_NAME
)

# %%
# call the optimized training model
train_loss, train_acc, valid_loss, valid_acc, epochs_ran, model = train_optimized_model(
    mlp_params.TRAIN_EPOCHS,
    train_loader,
    valid_loader,
    mlp_params,
    mlp_params.MODEL_NAME,
)
# get training_metrics
if params.MODEL_TYPE == "Regression":
    training_stats = pd.DataFrame(
        zip(train_loss, valid_loss, epochs_ran),
        columns=["train_loss", "valid_loss", "epochs_ran"],
    )
else:
    training_stats = pd.DataFrame(
        zip(train_loss, train_acc, valid_loss, valid_acc, epochs_ran),
        columns=["train_loss", "train_acc", "valid_loss", "valid_acc", "epochs_ran"],
    )

# %%
if params.MODEL_TYPE == "Regression":
    pass
else:
    plot_metric_vs_epoch(
        training_stats,
        x="epochs_ran",
        y1="train_acc",
        y2="valid_acc",
        title="Accuracy vs. Epochs",
        x_axis_label="Epochs",
        y_axis_label="Accuracy",
        params=params,
        model_name=params.MODEL_NAME,
    )

# %%
plot_metric_vs_epoch(
    training_stats,
    x="epochs_ran",
    y1="train_loss",
    y2="valid_loss",
    title="Loss vs. Epochs",
    x_axis_label="Epochs",
    y_axis_label="Loss",
    params=params,
    model_name=params.MODEL_NAME,
)

# %%
# calling the testing function and outputting list values of tested model
if params.MODEL_TYPE == "Multi_Class" or params.MODEL_TYPE == "Regression":
    y_pred_list = test_optimized_model(
        model, test_loader, params, model_name=params.MODEL_NAME
    )
elif params.MODEL_TYPE == "Binary_Classification":
    y_pred_list, y_pred_prob_list = test_optimized_model(
        model, test_loader, params, model_name=params.MODEL_NAME
    )
else:
    raise Exception("Model type must be specified for proper model testing")


# un-nest list if nested i.e. length of input data does not match length of output data
if len(y_pred_list) != len(Y_test):
    print("yes")
    y_pred_list = un_nest(y_pred_list)
    y_pred_prob_list = un_nest(y_pred_prob_list)
else:
    pass

# %%
# make the prediction list into a dataframe using column names from the test data
prediction_df = pd.DataFrame(y_pred_list, columns=Y_test.columns, index=Y_test.index)
prediction_df["Metadata_Well"] = Y_test_well["Metadata_Well"]
prediction_df["Metadata_Well"].unique()
prediction_df

# %%
# treatment_col = prediction_df["Metadata_Well"]
# index_dict = {}
# for i in test_col.unique():
#     index_list = []
#     for j in enumerate(test_col):
#         if i == j[1]:
#             index_list.append(j[0])
#     index_dict[i] = index_list
# new_value_dict = {}
# for i in index_dict:
#     new_value_dict[i] = [pred_col[index_dict[i]].median(),treatment_col[index_dict[i]].unique().tolist()[0]]


# df = pd.DataFrame.from_dict(new_value_dict, orient="index").reset_index().rename(columns={"index": "Actual", 0: "Average Predicted", 1: "Metadata_Well"})
# df['cytokine'] = col
# df['model_type'] = "test"

# %%
Y_test.head(2)

# %%
import matplotlib.pyplot as plt

# %%
testing_values = pd.DataFrame(columns=["Actual", "Average Predicted", "cytokine"])

# %%
from sklearn.metrics import mean_squared_error, r2_score

# get the list of columns in the test data
test_data_columns = Y_test.columns.to_list()
# loop through the columns
for col in Y_test:
    # get the column from test data
    test_col = Y_test[col]
    # get the column from prediction data
    pred_col = prediction_df[col]
    # list of treatment names
    treatment_col = prediction_df["Metadata_Well"]
    # get the mse and r2 for the columns
    mse = mean_squared_error(test_col, pred_col)
    r_square = r2_score(test_col, pred_col)
    a, b = np.polyfit(test_col, pred_col, 1)
    index_dict = {}
    for i in test_col.unique():
        index_list = []
        for j in enumerate(test_col):
            if i == j[1]:
                index_list.append(j[0])
        index_dict[i] = index_list
    new_value_dict = {}
    for i in index_dict:
        new_value_dict[i] = [
            pred_col[index_dict[i]].median(),
            treatment_col[index_dict[i]].unique().tolist()[0],
        ]

    df = (
        pd.DataFrame.from_dict(new_value_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "Actual", 0: "Average Predicted", 1: "Metadata_Well"})
    )
    df["cytokine"] = col
    df["model_type"] = "test"

    testing_values = pd.concat([testing_values, df], axis=0)
    # plt.scatter(df['Actual'], df['Average Predicted'])
    # # plt.plot([min(test_col), max(test_col)], [min(test_col), max(test_col)], color='red', linestyle='--')
    # plt.plot(
    #     test_col,
    #     a * test_col + b,
    #     color="red",
    #     label="R2={0:0.2f}".format(r_square),
    # )
    # plt.title(
    #     f"Regression Nerual Network Prediction vs. True \n {col}", fontsize=25
    # )
    # plt.ylabel("Predicted", fontsize=18)
    # plt.xlabel("Target", fontsize=18)
    # # make data continuous
    # plt.show()
    # plt.close()

# %%
# get the unique rows in a dataframe
key_df = df_descriptive[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
].drop_duplicates()
key_df

# %%
# add treatment column based on well
testing_values = pd.merge(testing_values, key_df, on="Metadata_Well", how="left")

# %%
testing_values

# %% [markdown]
# # Hold out

# %%
# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df_holdout.columns[df_holdout.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df_holdout[df_metadata]
df_values = df_holdout.drop(columns=df_metadata)

# %%
# get all columns that contain NSU in the name
df_values_Y = df_values[df_values.columns[df_values.columns.str.contains("NSU")]]
df_values_X = df_values.drop(columns=df_values_Y.columns)
df_values_Y["Metadata_Well"] = df_descriptive["Metadata_Well"]

print(df_values.shape)
print(df_values_X.shape)
print(df_values_Y.shape)

# %%
df_values_Y_well = df_values_Y
df_values_Y = df_values_Y.drop(columns=["Metadata_Well"])

# %%
test_data = Dataset_formatter(
    torch.FloatTensor(df_values_X.values), torch.FloatTensor(df_values_Y.values)
)

# convert data class into a dataloader to be compatible with pytorch
test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1)

# %%
# calling the testing function and outputting list values of tested model
if params.MODEL_TYPE == "Multi_Class" or params.MODEL_TYPE == "Regression":
    y_pred_list = test_optimized_model(
        model, test_loader, params, model_name=params.MODEL_NAME
    )
elif params.MODEL_TYPE == "Binary_Classification":
    y_pred_list, y_pred_prob_list = test_optimized_model(
        model, test_loader, params, model_name=params.MODEL_NAME
    )
else:
    raise Exception("Model type must be specified for proper model testing")


# un-nest list if nested i.e. length of input data does not match length of output data
if len(y_pred_list) != len(df_values_Y):
    y_pred_list = un_nest(y_pred_list)
    y_pred_prob_list = un_nest(y_pred_prob_list)
else:
    pass

# %%
# make the prediction list into a dataframe using column names from the test data
prediction_df = pd.DataFrame(
    y_pred_list, columns=df_values_Y.columns, index=df_values_Y.index
)
prediction_df.head(2)

# %%
print(df_values_Y.shape, prediction_df.shape)
df_values_Y.head(2)
prediction_df.head(2)

# %%
# make the prediction list into a dataframe using column names from the test data
prediction_df = pd.DataFrame(
    y_pred_list, columns=df_values_Y.columns, index=df_values_Y.index
)
prediction_df["Metadata_Well"] = df_values_Y_well["Metadata_Well"]
prediction_df["Metadata_Well"].unique()
prediction_df

# %%
df_values_Y.reset_index(drop=True, inplace=True)
prediction_df.reset_index(drop=True, inplace=True)

# %%
from sklearn.metrics import mean_squared_error, r2_score

# get the list of columns in the test data
test_data_columns = df_values_Y.columns.to_list()
# loop through the columns
for col in test_data_columns:
    # get the column from test data
    test_col = df_values_Y[col]
    # get the column from prediction data
    pred_col = prediction_df[col]
    # list of treatment names
    treatment_col = prediction_df["Metadata_Well"]
    # get the mse and r2 for the columns
    mse = mean_squared_error(test_col, pred_col)
    r_square = r2_score(test_col, pred_col)
    a, b = np.polyfit(test_col, pred_col, 1)
    index_dict = {}
    for i in test_col.unique():
        index_list = []
        for j in enumerate(test_col):
            if i == j[1]:
                index_list.append(j[0])
        index_dict[i] = index_list
    new_value_dict = {}
    for i in index_dict:
        new_value_dict[i] = [
            pred_col[index_dict[i]].median(),
            treatment_col[index_dict[i]].unique().tolist()[0],
        ]

    df = (
        pd.DataFrame.from_dict(new_value_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "Actual", 0: "Average Predicted", 1: "Metadata_Well"})
    )
    df["cytokine"] = col
    df["model_type"] = "holdout"

    testing_values = pd.concat([testing_values, df], axis=0)
    # plt.scatter(df['Actual'], df['Average Predicted'])
    # # plt.plot([min(test_col), max(test_col)], [min(test_col), max(test_col)], color='red', linestyle='--')
    # plt.plot(
    #     test_col,
    #     a * test_col + b,
    #     color="red",
    #     label="R2={0:0.2f}".format(r_square),
    # )
    # plt.title(
    #     f"Regression Nerual Network Prediction vs. True \n {col}", fontsize=25
    # )
    # plt.ylabel("Predicted", fontsize=18)
    # plt.xlabel("Target", fontsize=18)
    # # make data continuous
    # plt.show()
    # plt.close()

# %%
# get the unique rows in a dataframe
key_df = df_descriptive[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
].drop_duplicates()
key_df

# %%
# add treatment column based on well
testing_values = pd.merge(testing_values, key_df, on="Metadata_Well", how="left")

# %%
testing_values

# %%
# combine two columns into one oneb_Metadata_Treatment_Dose_Inhibitor_Dose_x and oneb_Metadata_Treatment_Dose_Inhibitor_Dose_y
testing_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = testing_values[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose_x"
].fillna(testing_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose_y"])
testing_values.drop(
    columns=[
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose_x",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose_y",
    ],
    inplace=True,
)
testing_values

# %%
