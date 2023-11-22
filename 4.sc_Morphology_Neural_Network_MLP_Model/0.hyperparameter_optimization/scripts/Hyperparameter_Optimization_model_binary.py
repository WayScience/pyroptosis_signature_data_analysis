#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8


# <span style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">An Exception was encountered at '<a href="#papermill-error-cell">In [3]</a>'.</span>

# ## Hyperparameter tuning via Optuna

# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.<br>
# ### Here I will be predicting if a cell received a treatment or not

# In[1]:

# In[ ]:


import pathlib
import sys

import numpy as np
import optuna
import pandas as pd
import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

# In[ ]:


# In[ ]:


sys.path.append("../..")


# In[ ]:


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

# In[ ]:


sys.path.append("../../..")
from utils.utils import df_stats

# ## Papermill is used for executing notebooks in the CLI with multiple parameters<br>
# Here the `injected-parameters` cell is used to inject parameters into the notebook via papermill.<br>
# This enables multiple notebooks to be executed with different parameters, preventing to manually update parameters or have multiple copies of the notebook.

# In[2]:

# Parameters

# In[ ]:


CELL_TYPE = "PBMC"
MODEL_NAME = "MultiClass_MLP"


# <span id="papermill-error-cell" style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">Execution using papermill encountered an exception here and stopped:</span>

# In[3]:

# In[ ]:


ml_configs_file = pathlib.Path("../../MLP_utils/binary_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)


# overwrite params via command line arguments from papermill

# In[ ]:


mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
mlp_params.CONTROL_NAME = CONTROL_NAME
mlp_params.TREATMENT_NAME = TREATMENT_NAME
mlp_params.MODEL_NAME = MODEL_NAME


# In[ ]:

# Import Data<br>
# set data file path under pathlib path for multi-system use

# In[ ]:


file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)


# In[ ]:


df = pq.read_table(file_path).to_pandas()


# In[ ]:

# In[ ]:


df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()


# #### Set up Data to be compatible with model

# ##### Classification Models:<br>
# Comment out code if using regression

# In[ ]:

# filter the oneb_Metadata_Treatment_Dose_Inhibitor_Dose column to only include the treatment and control via loc

# In[ ]:


df = df.loc[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        [mlp_params.TREATMENT_NAME, mlp_params.CONTROL_NAME]
    )
]


# In[ ]:


print("Selected Catagories are:")
print(df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique())
df_stats(df)


# In[ ]:


if mlp_params.DATA_SUBSET_OPTION == "True":
    df = df.sample(n=mlp_params.DATA_SUBSET_NUMBER)
    print("Data Subset Is On")
    print(f"Data is subset to {mlp_params.DATA_SUBSET_NUMBER}")
else:
    print("Data Subset Is Off")


# In[ ]:

# In[ ]:


np.random.seed(seed=0)
wells_to_hold = (
    df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose")
    .agg(np.random.choice)["Metadata_Well"]
    .to_list()
)
df_holdout = df[df["Metadata_Well"].isin(wells_to_hold)]
df = df[~df["Metadata_Well"].isin(wells_to_hold)]


# In[ ]:


print("Wells held out for testing:", df_holdout["Metadata_Well"].unique())
print(
    "Wells to use for training, validation, and testing", df["Metadata_Well"].unique()
)


# In[ ]:

# Code snippet for metadata extraction by Jenna Tomkinson

# In[ ]:


df_metadata = list(df.columns[df.columns.str.startswith("Metadata")])


# define which columns are data and which are descriptive

# In[ ]:


df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# In[ ]:

# Creating label encoder

# In[ ]:


le = preprocessing.LabelEncoder()
# Converting strings into numbers
df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = le.fit_transform(
    df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
)
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop(
    [
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "twob_Metadata_Treatment_Dose_Inhibitor_Dose",
        "threeb_Metadata_Treatment_Dose_Inhibitor_Dose",
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ],
    axis=1,
)
df_values_Y = df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]


# #### Split Data - All Models can proceed through this point

# In[ ]:

# In[ ]:


X_train, X_test, X_val, Y_train, Y_test, Y_val = data_split(
    X_vals=df_values_X,
    y_vals=df_values_Y,
    train_proportion=0.8,
    val_proportion=0.1,
    test_proportion=0.1,
    seed=1,
    params=params,
)


# In[ ]:

# produce data objects for train, val and test datasets

# In[ ]:


train_data = Dataset_formatter(
    torch.FloatTensor(X_train.values), torch.FloatTensor(Y_train.values)
)
val_data = Dataset_formatter(
    torch.FloatTensor(X_val.values), torch.FloatTensor(Y_val.values)
)
test_data = Dataset_formatter(
    torch.FloatTensor(X_test.values), torch.FloatTensor(Y_test.values)
)


# In[ ]:

# In[ ]:


mlp_params.IN_FEATURES = X_train.shape[1]
print("Number of in features: ", mlp_params.IN_FEATURES)
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.OUT_FEATURES = 1
else:
    mlp_params.OUT_FEATURES = len(
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
    )


# In[ ]:


print("Number of out features: ", mlp_params.OUT_FEATURES)


# In[ ]:


if mlp_params.OUT_FEATURES > 2:
    mlp_params.MODEL_TYPE = "Multi_Class"
elif mlp_params.OUT_FEATURES == 2:
    mlp_params.OUT_FEATURES = mlp_params.OUT_FEATURES - 1
    mlp_params.MODEL_TYPE = "Binary_Classification"
elif mlp_params.OUT_FEATURES == 1:
    mlp_params.MODEL_TYPE = "Regression"
else:
    pass
print(mlp_params.MODEL_TYPE)


# In[ ]:

# convert data class into a dataloader to be compatible with pytorch

# In[ ]:


train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=mlp_params.BATCH_SIZE
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=mlp_params.BATCH_SIZE
)
test_loader = torch.utils.data.DataLoader(
    dataset=test_data,
    batch_size=1,
)


# In[ ]:

# In[ ]:


print(mlp_params.DEVICE)


# In[ ]:

# no accuracy function must be loss for regression

# In[ ]:


if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.METRIC = "loss"
    pass


# wrap the objective function inside of a lambda function to pass args...

# In[ ]:


objective_lambda_func = lambda trial: objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=trial,
    params=params,
    metric=mlp_params.METRIC,
    return_info=False,
)


# Study is the object for model optimization

# In[ ]:


study = optuna.create_study(direction=f"{mlp_params.DIRECTION}")
# Here I apply the optimize function of the study to the objective function
# This optimizes each parameter specified to be optimized from the defined search space
study.optimize(objective_lambda_func, n_trials=mlp_params.N_TRIALS)
# Prints out the best trial's optimized parameters
objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=study.best_trial,
    params=params,
    metric=mlp_params.METRIC,
    return_info=True,
)


# In[ ]:

# create graph directory for this model

# In[ ]:


graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
)


# In[ ]:


pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_optimization_history(study)


# In[ ]:


graph_path = f"{graph_path}/plot_optimization_history_graph"


# In[ ]:


fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[ ]:

# create graph directory for this model

# In[ ]:


graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
).resolve(strict=True)


# In[ ]:


pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_intermediate_values(study)


# In[ ]:


graph_path = f"{graph_path}/plot_intermediate_values_graph"


# In[ ]:


fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[ ]:

# In[ ]:


param_dict = extract_best_trial_params(
    study.best_params, params, model_name=mlp_params.MODEL_NAME
)
