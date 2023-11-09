#!/usr/bin/env python
# coding: utf-8

# ## Hyperparameter tuning via Optuna

# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.
# ### Here I will be predicting if a cell received a treatment or not

# In[1]:


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
from sklearn.model_selection import train_test_split

sys.path.append("../../..")
from utils.utils import df_stats

# ## Papermill is used for executing notebooks in the CLI with multiple parameters
# Here the `injected-parameters` cell is used to inject parameters into the notebook via papermill.
# This enables multiple notebooks to be executed with different parameters, preventing to manually update parameters or have multiple copies of the notebook.

# In[2]:


# Parameters
CELL_TYPE = "SHSY5Y"
MODEL_NAME = "MultiClass_MLP"


# In[3]:


ml_configs_file = pathlib.Path("../../MLP_utils/multi_class_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)

# overwrite params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
mlp_params.MODEL_NAME = MODEL_NAME
MODEL_TYPE = mlp_params.MODEL_TYPE


# In[4]:


# Import Data
# set data file path under pathlib path for multi-system use

file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)

df1 = pd.read_parquet(file_path)


# In[5]:


# get paths for toml files
ground_truth_file_path = pathlib.Path(f"../../MLP_utils/ground_truth.toml").resolve(
    strict=True
)
treatment_splits_file_path = pathlib.Path(f"../../MLP_utils/splits.toml").resolve(
    strict=True
)
# read toml files
ground_truth = toml.load(ground_truth_file_path)
treatment_splits = toml.load(treatment_splits_file_path)


# In[6]:


# get information from toml files
apoptosis_groups_list = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_groups_list = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
healthy_groups_list = ground_truth["Healthy"]["healthy_groups_list"]
test_split_100 = treatment_splits["splits"]["data_splits_100"]
test_split_75 = treatment_splits["splits"]["data_splits_75"]


# In[7]:


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


# In[8]:


# add apoptosis, pyroptosis and healthy columns to dataframe
df1["apoptosis"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in apoptosis_groups_list,
    axis=1,
)
df1["pyroptosis"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in pyroptosis_groups_list,
    axis=1,
)
df1["healthy"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in healthy_groups_list,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column
df1["labels"] = df1.apply(
    lambda row: "apoptosis"
    if row["apoptosis"]
    else "pyroptosis"
    if row["pyroptosis"]
    else "healthy",
    axis=1,
)
# drop apoptosis, pyroptosis, and healthy columns
df1.drop(columns=["apoptosis", "pyroptosis", "healthy"], inplace=True)


# ### Split said data

# In[9]:


# randomly select wells to hold out for testing one per treatment group
# stratified by treatment group
np.random.seed(seed=0)
wells_to_hold = (
    df1.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose")
    .agg(np.random.choice)["Metadata_Well"]
    .to_list()
)
df_holdout = df1[df1["Metadata_Well"].isin(wells_to_hold)]
df = df1[~df1["Metadata_Well"].isin(wells_to_hold)]


print("Wells held out for testing:", df_holdout["Metadata_Well"].unique())
print(
    "Wells to use for training, validation, and testing", df1["Metadata_Well"].unique()
)


# In[10]:


# variable test and train set splits
# 100% test set
# subset the following treatments for test set
test_set_all = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_100)
]
# 75% test set and 25% train set
test_set_75 = df[df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_75)]

test_100_and_75 = test_split_100 + test_split_75

# 50% test set and 50% train set
# get all treatments that are not in the_test_set_all and the test_set_75
test_set_50 = df[
    ~df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_100_and_75)
]

print(test_set_all.shape, test_set_75.shape, test_set_50.shape)


# In[11]:


# get the train test splits from each group
# 100% test set
test_set_all

# 75% test set and 25% train set
test_ratio = 0.75
training_data_set_75, testing_data_set_75 = train_test_split(
    test_set_75,
    test_size=test_ratio,
    stratify=test_set_75["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    random_state=0,
)

# 50% test set and 50% train set
test_ratio = 0.5
training_data_set_50, testing_data_set_50 = train_test_split(
    test_set_50,
    test_size=test_ratio,
    stratify=test_set_50["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    random_state=0,
)

# verify that the correct splits have been made
# 100% test set
print(f"Shape for the 100% test set: {test_set_all.shape}\n")

# 75% test set and 25% train set
print(
    f"Shape for the 75% test set: {training_data_set_75.shape};\nShape for the 75% train set: {testing_data_set_75.shape}\n"
)

# 50% test set and 50% train set
print(
    f"Shape for the 50% test set: {training_data_set_50.shape};\nShape for the 50% train set: {testing_data_set_50.shape}"
)

print(f"Shape for the holdout set: {df_holdout.shape}")


# In[12]:


# combine all testing sets together while preserving the index
testing_data_set = pd.concat(
    [test_set_all, testing_data_set_75, testing_data_set_50], axis=0
)
testing_data_set = testing_data_set.sort_index()
testing_data_set

# combine all training sets together while preserving the index
training_data_set = pd.concat([training_data_set_75, training_data_set_50], axis=0)
training_data_set = training_data_set.sort_index()
training_data_set

val_data_set, training_data_set = train_test_split(
    training_data_set,
    test_size=0.20,
    stratify=training_data_set["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
)
print(
    f"""
    Testing set length: {len(testing_data_set)}\n
    Training set length: {len(training_data_set)}\n
    Validation set length: {len(val_data_set)}\n
    Holdout set length: {len(df_holdout)}"""
)
# get the indexes for the training and testing sets
testing_data_set_index = testing_data_set.index
training_data_set_index = training_data_set.index
val_data_set_index = val_data_set.index
df_holdout_index = df_holdout.index


# In[13]:


# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in training_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in testing_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "test"})
for index in df_holdout_index:
    index_data.append({"labeled_data_index": index, "label": "holdout"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data).sort_values(["labeled_data_index"])
index_data


# In[14]:


save_path = pathlib.Path(f"../indexes/{CELL_TYPE}/multi_class/")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)


# In[15]:


# save indexes as tsv file
index_data.to_csv(f"{save_path}/multi_class_data_split_indexes.tsv", sep="\t")


# #### Set up Data to be compatible with model

# ##### Classification Models:
# Comment out code if using regression

# In[16]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df1[df_metadata]
df_values = df1.drop(columns=df_metadata)


# In[17]:


# Creating label encoder
le = preprocessing.LabelEncoder()
# Converting strings into numbers
df_values["labels"] = le.fit_transform(df_values["labels"])
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop(
    ["labels"],
    axis=1,
)
df_values_Y = df_values["labels"]
df_values_Y.head()


# #### Split Data - All Models can proceed through this point

# In[18]:


# split into train and test sets from indexes previously defined

X_train = df_values_X.loc[training_data_set_index]
X_val = df_values_X.loc[val_data_set_index]
X_test = df_values_X.loc[testing_data_set_index]
X_holdout = df_values_X.loc[df_holdout_index]

Y_train = df_values_Y.loc[training_data_set_index]
Y_val = df_values_Y.loc[val_data_set_index]
Y_test = df_values_Y.loc[testing_data_set_index]
Y_holdout = df_values_Y.loc[df_holdout_index]


# In[19]:


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


# In[20]:


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


# In[21]:


# convert data class into a dataloader to be compatible with pytorch
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


# In[22]:


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
)


# In[24]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_optimization_history(study)


graph_path = f"{graph_path}/plot_optimization_history_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[25]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
).resolve(strict=True)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_intermediate_values(study)

graph_path = f"{graph_path}/plot_intermediate_values_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[26]:


param_dict = extract_best_trial_params(
    study.best_params, params, model_name=mlp_params.MODEL_NAME
)
