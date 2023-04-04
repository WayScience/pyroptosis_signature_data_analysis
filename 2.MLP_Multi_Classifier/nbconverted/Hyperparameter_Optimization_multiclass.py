#!/usr/bin/env python
# coding: utf-8

# ## Hyperparameter tuning via Optuna for Binary MLP model

# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.
# ### Here I will be predicting if a cell received a treatment or not

# In[1]:


import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import plotly
import seaborn as sns
import toml
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

sys.path.append("..")

from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    data_split,
    extract_best_trial_params,
    objective_model_optimizer,
    optimized_model_create,
    parameter_set,
    plot_metric_vs_epoch,
    results_output,
    test_optimized_model,
    train_optimized_model,
    un_nest,
)
from sklearn.metrics import (
    auc,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)

from utils.utils import df_stats

# In[2]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_fs_cellprofiler.csv.gz"
)
df = pd.read_csv(
    file_path,
    low_memory=False,
)


# In[3]:


data = Path("MLP_utils/config.toml")
config = toml.load(data)
params = Parameters()
params = parameter_set(params, config)


# In[4]:


# Combine treatment with dosage to be able to discern treatments with different doses as a different condition
# Combine treatment and dose
df = df.assign(
    Metadata_Treatment_and_Dose=lambda x: df["Metadata_treatment"]
    + "_"
    + df["Metadata_dose"]
)

# df["Metadata_treatment"] = df["Metadata_treatment"] + "_" + df["Metadata_dose"]
print(df["Metadata_Treatment_and_Dose"].unique())

# Generate df speceific to analysis and model
df = df.query(
    "Metadata_Treatment_and_Dose == 'LPS_10µg/ml'| Metadata_Treatment_and_Dose == 'Media only_0' | Metadata_Treatment_and_Dose == 'Disulfiram_2.5µM'"
)
print(df["Metadata_Treatment_and_Dose"].unique())

df_stats(df)
# Drop na and reindex accordingly
df = df.dropna()
df = df.reset_index(drop=True)

# Check for Nans again
df_stats(df)
# Understand categorical data such as treatment and dosing
df[["Metadata_Treatment_and_Dose"]].drop_duplicates()
print(params.DATA_SUBSET_OPTION)
print(params.DATA_SUBSET_NUMBER)
if params.DATA_SUBSET_OPTION == True:
    df = df.sample(n=params.DATA_SUBSET_NUMBER)
    print("yes")
else:
    pass

# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.startswith("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)
print(len(df_values))


#  ### Setting up data for network training

# In[5]:


# Creating label encoder
le = preprocessing.LabelEncoder()
# Converting strings into numbers
df_values["Metadata_Treatment_and_Dose"] = le.fit_transform(
    df_descriptive["Metadata_Treatment_and_Dose"]
)
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop("Metadata_Treatment_and_Dose", axis=1)
df_values_Y = df_values["Metadata_Treatment_and_Dose"]

df_values_X.head()

X_train, X_test, X_val, Y_train, Y_test, Y_val = data_split(
    X_vals=df_values_X,
    y_vals=df_values_Y,
    train_proportion=0.8,
    val_proportion=0.1,
    test_proportion=0.1,
    seed=1,
)


# In[6]:


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

IN_FEATURES = X_train.shape[1]
print("Number of in features: ", IN_FEATURES)
OUT_FEATURES = len(df_values["Metadata_Treatment_and_Dose"].unique())
print("Number of out features: ", OUT_FEATURES)


# In[7]:


# convert data class into a dataloader to be compatible with pytorch
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=params.BATCH_SIZE
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=params.BATCH_SIZE
)
test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1)


# In[8]:


# wrap the objective function inside of a lambda function to pass args...
objective_lambda_func = lambda trial: objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=trial,
    in_features=IN_FEATURES,
    out_features=OUT_FEATURES,
    params=params,
    metric=params.METRIC,
    return_info=False,
)
# Study is the object for model optimization
study = optuna.create_study(direction=f"{params.DIRECTION}")
# Here I apply the optimize function of the study to the objective function
# This optimizes each parameter specified to be optimized from the defined search space
study.optimize(objective_lambda_func, n_trials=params.N_TRIALS)
# Prints out the best trial's optimized parameters
objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=study.best_trial,
    in_features=IN_FEATURES,
    out_features=OUT_FEATURES,
    params=params,
    metric=params.METRIC,
    return_info=True,
)


# In[9]:


optuna.visualization.plot_optimization_history(study)


# In[10]:


optuna.visualization.plot_intermediate_values(study)


# In[11]:


# call function for best trial parameter extraction
param_dict = extract_best_trial_params(study.best_params)


# In[12]:


# call the optimized training model
train_loss, train_acc, valid_loss, valid_acc, epochs_ran, model = train_optimized_model(
    params.TRAIN_EPOCHS,
    train_loader,
    valid_loader,
    IN_FEATURES,
    OUT_FEATURES,
    param_dict,
    params,
)
# create a DataFrame of each stat
training_stats = pd.DataFrame(
    zip(train_loss, train_acc, valid_loss, valid_acc, epochs_ran),
    columns=["train_loss", "train_acc", "valid_loss", "valid_acc", "epochs_ran"],
)


# In[13]:


plot_metric_vs_epoch(
    training_stats,
    x="epochs_ran",
    y1="train_acc",
    y2="valid_acc",
    title="Accuracy vs. Epochs",
    x_axis_label="Epochs",
    y_axis_label="Accuracy",
)


# In[14]:


plot_metric_vs_epoch(
    training_stats,
    x="epochs_ran",
    y1="train_loss",
    y2="valid_loss",
    title="Loss vs. Epochs",
    x_axis_label="Epochs",
    y_axis_label="Loss",
)


# In[15]:


# calling the testing function and outputing list values of tested model
y_pred_list = test_optimized_model(
    model, test_loader, IN_FEATURES, OUT_FEATURES, param_dict, params
)
# un-nest list if nested i.e. length of input data does not match length of output data
if len(y_pred_list) != len(Y_test):
    y_pred_list = un_nest(y_pred_list)
    y_pred_prob_list = un_nest(y_pred_prob_list)
else:
    pass


# In[16]:


# Call visualization function
results_output(y_pred_list, Y_test, OUT_FEATURES)


# In[ ]:
