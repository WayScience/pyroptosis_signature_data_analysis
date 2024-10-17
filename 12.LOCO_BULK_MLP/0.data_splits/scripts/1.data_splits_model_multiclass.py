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

# parse arguments
args = parser.parse_args()

CELL_TYPE = args.cell_type
MODEL_NAME = args.model_name


# In[3]:


ml_configs_file = pathlib.Path(MLP_path / "multi_class_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)

# overwrite params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
MODEL_TYPE = mlp_params.MODEL_TYPE
HYPERPARAMETER_BATCH_SIZE = mlp_params.HYPERPARAMETER_BATCH_SIZE


# In[4]:


# Import Data
# set data file path under pathlib path for multi-system use

file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)

df1 = pd.read_parquet(file_path)


# In[5]:


# get paths for toml files
ground_truth_file_path = pathlib.Path(MLP_path / "ground_truth.toml").resolve(
    strict=True
)
treatment_splits_file_path = pathlib.Path(MLP_path / "splits.toml").resolve(strict=True)
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


# In[8]:


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
print(df_holdout.shape, df.shape)


# In[10]:


# variable test and train set splits
# 100% test set
# subset the following treatments for test set
treatment_holdout = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_100)
]
df = df[~df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_100)]
print(treatment_holdout.shape, df.shape)


# In[11]:


training_data_set, testing_data_set = train_test_split(
    df,
    test_size=0.20,
    stratify=df["labels"],
)

print(training_data_set.shape, testing_data_set.shape)

training_data_set, val_data_set = train_test_split(
    training_data_set,
    test_size=0.20,
    stratify=training_data_set["labels"],
)
print(
    f"""
    Testing set length: {len(testing_data_set)}\n
    Training set length: {len(training_data_set)}\n
    Validation set length: {len(val_data_set)}\n
    Treatment Holdout set length: {len(treatment_holdout)}\n
    Holdout set length: {len(df_holdout)}
    Added set length: {len(testing_data_set) + len(training_data_set) + len(val_data_set) + len(treatment_holdout) + len(df_holdout)}
    Total actual set length: {len(df1)}
"""
)


# In[12]:


# get the indexes for the training and testing sets

training_data_set_index = training_data_set.index
val_data_set_index = val_data_set.index
testing_data_set_index = testing_data_set.index
treatment_holdout_index = treatment_holdout.index
df_holdout_index = df_holdout.index

assert len(training_data_set_index) + len(val_data_set_index) + len(
    testing_data_set_index
) + len(treatment_holdout_index) + len(df_holdout_index) == len(df1)


# In[13]:


print(
    training_data_set_index.shape,
    val_data_set_index.shape,
    testing_data_set_index.shape,
    treatment_holdout_index.shape,
    df_holdout_index.shape,
)
print(
    training_data_set_index.shape[0]
    + val_data_set_index.shape[0]
    + testing_data_set_index.shape[0]
    + treatment_holdout_index.shape[0]
    + df_holdout_index.shape[0]
)


# In[14]:


# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in training_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in val_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "val"})
for index in testing_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "test"})
for index in treatment_holdout_index:
    index_data.append({"labeled_data_index": index, "label": "treatment_holdout"})
for index in df_holdout_index:
    index_data.append({"labeled_data_index": index, "label": "holdout"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data)
index_data


# In[15]:


index_data["label"].unique()


# In[16]:


save_path = pathlib.Path(f"../indexes/{CELL_TYPE}/multi_class/")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)
# save indexes as tsv file
index_data.to_csv(
    f"{save_path}/{params.CELL_TYPE}_data_split_indexes.tsv", sep="\t", index=False
)


# In[17]:


# get the class weights for the loss function to account for class imbalance
# get the number of samples in each class
targets, counts = np.unique(df1["labels"], return_counts=True)
print(targets, counts)
total_counts = np.sum(counts)
# get the class weights
class_weights = []
class_targets = []
for class_name in enumerate(targets):
    class_targets.append(class_name[1])
for count in enumerate(counts):
    class_weights.append(1 - (count[1] / total_counts))
print(class_targets, class_weights)
# write the class weights to a file for use in the model
class_weights_file = pathlib.Path(f"../class_weights/{CELL_TYPE}/multi_class/")
class_weights_file.mkdir(parents=True, exist_ok=True)
class_targets_dicts = {
    class_targets[i]: class_weights[i] for i in range(len(class_targets))
}
# write the file to json
class_weights_file = class_weights_file / "class_weights.json"
with open(class_weights_file, "w") as file:
    json.dump(class_targets_dicts, file)
