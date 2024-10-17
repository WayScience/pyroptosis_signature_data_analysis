#!/usr/bin/env python
# coding: utf-8

# ## Hyperparameter tuning via Optuna

# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.
# ### Here I will be predicting if a cell received a treatment or not

# In[1]:


import argparse
import json
import pathlib
import sys

import numpy as np
import optuna
import pandas as pd
import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

MLP_parent_path = pathlib.Path("../../../utils/")
sys.path.append(str(MLP_parent_path.resolve()))
MLP_path = pathlib.Path("../../../utils/MLP_utils").resolve()

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
from sklearn.model_selection import train_test_split

from utils import df_stats

# In[2]:


# set up the parser
parser = argparse.ArgumentParser(description="Run hyperparameter optimization")
parser.add_argument(
    "--cell_type",
    type=str,
    default="all",
    help="Cell type to run hyperparameter optimization for",
)
parser.add_argument(
    "--model_name",
    type=str,
    default="all",
    help="Model name to run hyperparameter optimization for",
)

parser.add_argument(
    "--channel_combination_key",
    type=str,
    default="all",
    help="key to a dictionary containing the feature types to split the data into",
)

# parse arguments
args = parser.parse_args()

CELL_TYPE = args.cell_type
MODEL_NAME = args.model_name
channel_combination_key = args.channel_combination_key


# In[3]:


# load in the channel combinations file
channel_combinations_file_path = pathlib.Path(
    f"../../0.data_splits/results/feature_combinations_{CELL_TYPE}.toml"
).resolve(strict=True)
channel_combinations = toml.load(channel_combinations_file_path)
channel_combinations = channel_combinations[channel_combination_key]


# In[4]:


ml_configs_file = pathlib.Path(MLP_path / "multi_class_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs, "Bulk_Leave_One_Channel_Out_Params")

# overwrite params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = f"{MODEL_NAME}_{channel_combination_key}"
MODEL_TYPE = mlp_params.MODEL_TYPE
HYPERPARAMETER_BATCH_SIZE = mlp_params.HYPERPARAMETER_BATCH_SIZE


# In[5]:


# Import Data
# set data file path under pathlib path for multi-system use

file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)

df1 = pd.read_parquet(file_path, columns=channel_combinations)
print(df1.shape)


# In[6]:


# get paths for toml files
ground_truth_file_path = pathlib.Path(MLP_path / "ground_truth.toml").resolve(
    strict=True
)
treatment_splits_file_path = pathlib.Path(MLP_path / "splits.toml").resolve(strict=True)
# read toml files
ground_truth = toml.load(ground_truth_file_path)
treatment_splits = toml.load(treatment_splits_file_path)


# In[7]:


# get information from toml files
apoptosis_groups_list = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_groups_list = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
healthy_groups_list = ground_truth["Healthy"]["healthy_groups_list"]
test_split_100 = treatment_splits["splits"]["data_splits_100"]
test_split_75 = treatment_splits["splits"]["data_splits_75"]


# In[8]:


np.random.seed(0)
if mlp_params.DATA_SUBSET_OPTION == "True":
    df1 = df1.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").apply(
        lambda x: x.sample(n=mlp_params.DATA_SUBSET_NUMBER, random_state=0)
    )
    print("Data Subset Is On")
    print(f"Data is subset to {mlp_params.DATA_SUBSET_NUMBER} per treatment group")
    print(df1.shape)
    df1.reset_index(drop=True, inplace=True)
else:
    print("Data Subset Is Off")


# In[9]:


# add apoptosis, pyroptosis and healthy columns to dataframe
df1["apoptosis"] = df1["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    apoptosis_groups_list
)
df1["pyroptosis"] = df1["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    pyroptosis_groups_list
)
df1["healthy"] = df1["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    healthy_groups_list
)

# merge apoptosis, pyroptosis, and healthy columns into one column
conditions = [
    (df1["apoptosis"] == True),
    (df1["pyroptosis"] == True),
    (df1["healthy"] == True),
]
choices = ["apoptosis", "pyroptosis", "healthy"]
df1["labels"] = np.select(conditions, choices, default="healthy")

# drop apoptosis, pyroptosis, and healthy columns
df1.drop(columns=["apoptosis", "pyroptosis", "healthy"], inplace=True)


# ### Split said data

# In[10]:


index_path = pathlib.Path(f"../../0.data_splits/indexes/{CELL_TYPE}/multi_class/")


# save indexes as tsv file
index_data_path = pathlib.Path(
    f"{index_path}/{params.CELL_TYPE}_data_split_indexes.tsv", sep="\t", index=False
)
indexes = pd.read_csv(index_data_path, sep="\t")

# get the labeld_data_index column from the indexes dataframe for train label
training_data_set_index = indexes.loc[
    indexes["label"] == "train", "labeled_data_index"
].values
val_data_set_index = indexes.loc[indexes["label"] == "val", "labeled_data_index"].values
testing_data_set_index = indexes.loc[
    indexes["label"] == "test", "labeled_data_index"
].values
treatment_holdout_index = indexes.loc[
    indexes["label"] == "treatment_holdout", "labeled_data_index"
].values
holdout_index = indexes.loc[indexes["label"] == "holdout", "labeled_data_index"].values


# #### Set up Data to be compatible with model

# ##### Classification Models:
# Comment out code if using regression

# In[11]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df1.columns[df1.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df1[df_metadata]
df_values = df1.drop(columns=df_metadata)


# In[12]:


# get weights
class_weights_file = pathlib.Path(
    f"../../0.data_splits/class_weights/{CELL_TYPE}/multi_class/class_weights.json"
).resolve(strict=True)
# read in the class weights file json into a dict
with open(class_weights_file, "r") as file:
    class_weights = json.load(file)

class_weights


# In[13]:


# Creating label encoder
le = preprocessing.LabelEncoder()
df_values["new_labels"] = le.fit_transform(df_values["labels"])
# make a dict encoder for the labels
encoder = dict(zip(le.classes_, le.transform(le.classes_)))
# change the keys and valules to be strings
encoder = {str(key): str(value) for key, value in encoder.items()}
# save the encoder as a json
encoder_file = pathlib.Path(
    f"../../0.data_splits/class_weights/{CELL_TYPE}/multi_class/encoder.json"
)
with open(encoder_file, "w") as file:
    json.dump(encoder, file)

# get mini dataframe that contains the decoder
decoder = df_values[["labels", "new_labels"]].drop_duplicates()
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop(
    ["new_labels", "labels"],
    axis=1,
)
df_values_Y = df_values["new_labels"]
df_values_Y.head()
df_values_Y.unique()


# In[14]:


# replace the class weights key with the encoder key's matching value
class_weights = {encoder[key]: value for key, value in class_weights.items()}
class_weights = [f for f in class_weights.values()]
class_weights


# #### Split Data - All Models can proceed through this point

# In[15]:


# split into train and test sets from indexes previously defined

X_train = df_values_X.loc[training_data_set_index]
X_val = df_values_X.loc[val_data_set_index]
X_test = df_values_X.loc[testing_data_set_index]
X_holdout = df_values_X.loc[holdout_index]
X_treatment_holdout = df_values_X.loc[treatment_holdout_index]

Y_train = df_values_Y.loc[training_data_set_index]
Y_val = df_values_Y.loc[val_data_set_index]
Y_test = df_values_Y.loc[testing_data_set_index]
Y_holdout = df_values_Y.loc[holdout_index]
Y_treatment_holdout = df_values_Y.loc[treatment_holdout_index]


# In[16]:


mlp_params.OUT_FEATURES = len(df_values_Y.unique())
mlp_params.OUT_FEATURES


# In[17]:


Y_train = torch.tensor(Y_train.values)
Y_train = torch.nn.functional.one_hot(
    Y_train, num_classes=mlp_params.OUT_FEATURES
).float()

Y_val = torch.tensor(Y_val.values)
Y_val = torch.nn.functional.one_hot(Y_val, num_classes=mlp_params.OUT_FEATURES).float()

Y_test = torch.tensor(Y_test.values)
Y_test = torch.nn.functional.one_hot(
    Y_test, num_classes=mlp_params.OUT_FEATURES
).float()

Y_holdout = torch.tensor(Y_holdout.values)
Y_holdout = torch.nn.functional.one_hot(
    Y_holdout, num_classes=mlp_params.OUT_FEATURES
).float()

Y_treatment_holdout = torch.tensor(Y_treatment_holdout.values)
Y_treatment_holdout = torch.nn.functional.one_hot(
    Y_treatment_holdout, num_classes=mlp_params.OUT_FEATURES
).float()

# convert the X dataframes to tensors
X_train = torch.tensor(X_train.values)
X_val = torch.tensor(X_val.values)
X_test = torch.tensor(X_test.values)
X_holdout = torch.tensor(X_holdout.values)
X_treatment_holdout = torch.tensor(X_treatment_holdout.values)


# In[18]:


# produce data objects for train, val and test datasets
train_data = torch.utils.data.TensorDataset(X_train, Y_train)
val_data = torch.utils.data.TensorDataset(X_val, Y_val)
test_data = torch.utils.data.TensorDataset(X_test, Y_test)


# In[19]:


mlp_params.IN_FEATURES = X_train.shape[1]
print("Number of in features: ", mlp_params.IN_FEATURES)
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.OUT_FEATURES = 1
else:
    mlp_params.OUT_FEATURES = len(df_values["labels"].unique())

print("Number of out features: ", mlp_params.OUT_FEATURES)

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


# In[20]:


# convert data class into a dataloader to be compatible with pytorch
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=mlp_params.HYPERPARAMETER_BATCH_SIZE, shuffle=True
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=mlp_params.HYPERPARAMETER_BATCH_SIZE, shuffle=False
)


# In[21]:


mlp_params.OUT_FEATURES


# In[22]:


# check device
print(mlp_params.DEVICE)


# In[23]:


# no accuracy function must be loss for regression
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.METRIC = "loss"
    pass


# wrap the objective function inside of a lambda function to pass args...
objective_lambda_func = lambda trial: objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=trial,
    params=params,
    metric=mlp_params.METRIC,
    return_info=False,
    class_weights=class_weights,
)


# Study is the object for model optimization
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
    class_weights=class_weights,
)


# In[24]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}_{channel_combination_key}/hyperparameter_optimization"
)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_optimization_history(study)


graph_path = f"{graph_path}/plot_optimization_history_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))


# In[25]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}_{channel_combination_key}/hyperparameter_optimization"
).resolve(strict=True)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_intermediate_values(study)

graph_path = f"{graph_path}/plot_intermediate_values_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))


# In[26]:


param_dict = extract_best_trial_params(
    study.best_params, params, model_name=mlp_params.MODEL_NAME
)
