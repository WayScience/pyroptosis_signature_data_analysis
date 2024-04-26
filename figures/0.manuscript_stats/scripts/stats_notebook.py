#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import pathlib

import numpy as np
import pandas as pd
import toml

# ## Morphology Feature space stats

# In[2]:


# set paths to data
norm_data_path = pathlib.Path("../../../data/PBMC_sc_norm.parquet").resolve(strict=True)

# fs data
norm_fs_data_path = pathlib.Path(
    "../../../data/PBMC_preprocessed_sc_norm.parquet"
).resolve(strict=True)


# ## Check Raw features shape

# In[3]:


# load in the normalized data
norm_data = pd.read_parquet(norm_data_path)


# In[4]:


# get columns that contain Metadata
metadata_cols = [col for col in norm_data.columns if "Metadata" in col]
metadata_df = norm_data[metadata_cols]
features_df = norm_data.drop(metadata_cols, axis="columns")
print(f"metadata_df shape: {metadata_df.shape}")
print(f"features_df shape: {features_df.shape}")
print(f"There are {metadata_df.shape[0]} cells in the dataset")
print(f"There are {features_df.shape[1]} features in the dataset")


# In[5]:


# remove the dfs from memory
del norm_data
del metadata_df
del features_df
# collect the garbage
gc.collect()


# ## Check feature selected shape

# In[6]:


# load in the feature selected data
norm_fs_df = pd.read_parquet(norm_fs_data_path)
# get columns that contain Metadata
metadata_cols = [col for col in norm_fs_df.columns if "Metadata" in col]
metadata_df = norm_fs_df[metadata_cols]
features_df = norm_fs_df.drop(metadata_cols, axis="columns")
print(f"metadata_df shape: {metadata_df.shape}")
print(f"features_df shape: {features_df.shape}")
print(f"There are {metadata_df.shape[0]} cells in the dataset")
print(f"There are {features_df.shape[1]} features in the dataset")

# remove the dfs from memory
del norm_fs_df
del metadata_df
del features_df
# collect the garbage
gc.collect()

norm_fs_df_subset = pd.read_parquet(
    norm_fs_data_path,
    columns=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
)
print(norm_fs_df_subset.shape)
norm_fs_df_subset.head()


# In[7]:


import toml

# load in the ground truth
ground_truth_file_path = pathlib.Path(
    "../../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
).resolve(strict=True)
# read in the ground truth toml file
ground_truth = toml.load(ground_truth_file_path)
# get information from toml files
apoptosis_groups_list = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_groups_list = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
healthy_groups_list = ground_truth["Healthy"]["healthy_groups_list"]
# add apoptosis, pyroptosis and healthy columns to dataframe
norm_fs_df_subset["apoptosis"] = norm_fs_df_subset.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in apoptosis_groups_list,
    axis=1,
)
norm_fs_df_subset["pyroptosis"] = norm_fs_df_subset.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in pyroptosis_groups_list,
    axis=1,
)
norm_fs_df_subset["healthy"] = norm_fs_df_subset.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in healthy_groups_list,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column
norm_fs_df_subset["labels"] = norm_fs_df_subset.apply(
    lambda row: "apoptosis"
    if row["apoptosis"]
    else "pyroptosis"
    if row["pyroptosis"]
    else "healthy",
    axis=1,
)
# drop apoptosis, pyroptosis, and healthy columns
norm_fs_df_subset.drop(columns=["apoptosis", "pyroptosis", "healthy"], inplace=True)


# In[8]:


# print the number of samples in each class
print(norm_fs_df_subset["labels"].value_counts())


# ## Stats for the Elastic Net models

# In[9]:


# set path for models performances
model_performances_path = pathlib.Path(
    "../../../6.bulk_Morphology_Elastic_Network/4.model_performance/results/regression/PBMC/all_model_performance.csv"
).resolve(strict=True)
# load in the model performances
model_performances = pd.read_csv(model_performances_path)


# In[10]:


# drop uneeded columns
columns_to_drop = [
    "feature_names",
    "coefficients",
    "cell_type",
    "alpha",
    "l1_ratio",
]
model_performances.drop(columns=columns_to_drop, inplace=True)
# drop duplicates
print(model_performances.shape)
model_performances.drop_duplicates(inplace=True)
print(model_performances.shape)
model_performances.head()


# In[11]:


# split the shuffled and final model performances
suffled_models = model_performances.loc[model_performances["shuffle"] == "shuffled"]
final_models = model_performances.loc[model_performances["shuffle"] == "final"]
print(suffled_models.shape)
print(final_models.shape)


# In[12]:


# sort the final models by r2 score
final_models.sort_values(by="r2", ascending=False, inplace=True)
final_models.head()


# In[13]:


# get the percentage of models that are above the threshold
threshold = 0.8
final_models_above_threshold = final_models.loc[final_models["r2"] >= threshold]
print(
    f"Percentage of models with r2 score above {threshold}: "
    f"{(final_models_above_threshold.shape[0] / final_models.shape[0]) * 100}",
    f"\n"
    f"The total number of models above the threshold is: {final_models_above_threshold.shape[0]}",
)


# In[14]:


# sort the shuffled models by r2 score from low to high
final_models.sort_values(by="r2", ascending=True, inplace=True)
final_models.head()
