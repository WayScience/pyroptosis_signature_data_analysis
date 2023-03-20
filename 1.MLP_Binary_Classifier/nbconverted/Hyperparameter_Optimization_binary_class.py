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
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn import preprocessing
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

sys.path.append("..")
from MLP_utils.parameters import Parameters

# from MLP_utils.utils import build_model_custom
# from MLP_utils.utils import train_n_validate
from MLP_utils.utils import (
    Dataset,
    data_split,
    extract_best_trial_params,
    objective,
    optimized_model,
    plot_metric_vs_epoch,
    results_output,
    test_optimized_model,
    train_optimized_model,
    un_nest,
)

from utils.utils import df_stats

params = Parameters()


# In[2]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz"
)
df = pd.read_csv(
    file_path,
    low_memory=False,
)


# Combine treatment with dosage to be able to discern treatments with different doses as a different condition

# In[3]:


# Combine treatment and dose
df["Metadata_treatment"] = df["Metadata_treatment"] + "_" + df["Metadata_dose"]
print(df["Metadata_treatment"].unique())

# Generate df speceific to analysis and model
df = df.query(
    "Metadata_treatment == 'LPS_10µg/ml'| Metadata_treatment == 'Media only_0'"
)
print(df["Metadata_treatment"].unique())

df_stats(df)
# Drop na and reindex accordingly
df = df.dropna()
df.reindex()
# Check for Nans again
df_stats(df)
# Understand categorical data such as treatment and dosing
df[["Metadata_treatment", "Metadata_dose"]].drop_duplicates()
if params.SUBSET_OPTION:
    df = df.sample(n=params.SUBSET_NUMBER)
else:
    pass
# Code snipptet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.startswith("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# ### Setting up data for network training

# In[4]:


# Creating label encoder
le = preprocessing.LabelEncoder()
# Converting strings into numbers
df_values["Metadata_treatment"] = le.fit_transform(df_descriptive["Metadata_treatment"])
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop("Metadata_treatment", axis=1)
df_values_Y = df_values["Metadata_treatment"]

X_train, X_test, X_val, Y_train, Y_test, Y_val = data_split(
    df_values_X, df_values_Y, 0.8, 0.1, 0.1
)


# In[5]:


# produce data objects for train, val and test datasets
train_data = Dataset(
    torch.FloatTensor(X_train.values), torch.FloatTensor(Y_train.values)
)
val_data = Dataset(torch.FloatTensor(X_val.values), torch.FloatTensor(Y_val.values))
test_data = Dataset(torch.FloatTensor(X_test.values), torch.FloatTensor(Y_test.values))

IN_FEATURES = X_train.shape[1]
print("Number of in features: ", IN_FEATURES)
out_features = len(df_values["Metadata_treatment"].unique())
print("Number of out features: ", out_features)


# In[6]:


# convert data class into a dataloader to be compatible with pytorch
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=params.BATCH_SIZE
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=params.BATCH_SIZE
)
test_loader = torch.utils.data.DataLoader(
    dataset=test_data, batch_size=params.BATCH_SIZE
)


# In[7]:


# wrap the objective function inside of a lambda function to pass args...
objective_lambda_func = lambda trial: objective(
    trial, IN_FEATURES, train_loader, valid_loader, params, False
)
# Study is the object for model optimzation
study = optuna.create_study(direction="minimize")
# Here I appply the optimize function of the study to the objective function
# This optimizes each parameter specified to be optinmized from the defined search space
study.optimize(objective_lambda_func, n_trials=params.N_TRIALS)
# Prints out the best trial's optimized parameters
objective(study.best_trial, IN_FEATURES, train_loader, valid_loader, params, True)


# In[ ]:


# In[8]:


optuna.visualization.plot_optimization_history(study)


# In[9]:


optuna.visualization.plot_intermediate_values(study)


# In[10]:


# call function
param_dict = extract_best_trial_params(study.best_params)


# In[11]:


# call the optimized trainig model
train_loss, train_acc, valid_loss, valid_acc, epochs_ran, model = train_optimized_model(
    params.TRAIN_EPOCHS, train_loader, valid_loader, IN_FEATURES, param_dict, params
)


# In[12]:


training_stats = pd.DataFrame(
    zip(train_loss, train_acc, valid_loss, valid_acc, epochs_ran),
    columns=["train_loss", "train_acc", "valid_loss", "valid_acc", "epochs_ran"],
)


# In[13]:


plot_metric_vs_epoch(training_stats, "epochs_ran", "train_acc", "valid_acc")


# In[14]:


plot_metric_vs_epoch(training_stats, "epochs_ran", "train_loss", "valid_loss")


# In[15]:


# calling the testing function and outputing list values of tested model
y_pred_list, y_pred_prob_list = test_optimized_model(
    model, test_loader, IN_FEATURES, param_dict, params
)


# un-nest list if nested i.e. length of input data does not match length of output data
if len(y_pred_list) != len(Y_test):
    y_pred_list = un_nest(y_pred_list)
    y_pred_prob_list = un_nest(y_pred_prob_list)
else:
    pass


# In[16]:


# Call visulalization function
results_output(y_pred_list, y_pred_prob_list, Y_test)
