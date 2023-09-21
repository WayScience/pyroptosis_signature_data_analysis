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
mlp_params.SHUFFLE = False

# %%
# Import Data
# set data file path under pathlib path for multi-system use

# Commented out for now, using a different data set to trobleshoot
file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)


# set path for nomic data
nomic_df_path = pathlib.Path(
    f"../../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{mlp_params.CELL_TYPE}_cleanup4correlation.csv"
).resolve(strict=True)

df = pd.read_parquet(file_path)
nomic_df = pd.read_csv(nomic_df_path)

# %%
# change the nomic df to standard scaler
# select the columns that contain "NSU"
nomic_df_scaled = nomic_df.filter(regex="NSU")
# standardize the nomic data
# scaler = preprocessing.StandardScaler()
# nomic_df_scaled = pd.DataFrame(scaler.fit_transform(nomic_df_scaled), columns=nomic_df_scaled.columns)
# add the nomic data metadata back
nomic_df_scaled[
    [
        "Metadata_position_x",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ]
] = nomic_df[
    [
        "Metadata_position_x",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ]
]
nomic_df = nomic_df_scaled.copy()
del nomic_df_scaled

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
# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = df.columns[df.columns.str.contains("Metadata")].to_list()

# define which columns are data and which are descriptive
df_values = df.drop(columns=df_metadata)

# %%
df_values[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
] = df_metadata[["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]]
df = (
    df_values.groupby(["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"])
    .median()
    .reset_index()
)

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
# get random wells from each treatment group to hold out
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
df_metadata = df.columns[df.columns.str.contains("Metadata")].to_list()

# define which columns are data and which are descriptive
df_values = df.drop(columns=df_metadata)

# %%
# get all columns that contain NSU in the name
df_values_Y = df_values[df_values.columns[df_values.columns.str.contains("NSU")]]
df_values_X = df_values.drop(columns=df_values_Y.columns)
# drop all columns except for IL1B and TNFa
col = ["IL-1 beta [NSU]"]
df_values_Y = df_values_Y[col]
df_values_Y["Metadata_Well"] = df_metadata["Metadata_Well"]
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
