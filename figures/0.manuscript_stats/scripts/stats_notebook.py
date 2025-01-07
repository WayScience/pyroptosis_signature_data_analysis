#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import pathlib

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
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
# norm_data = pd.read_parquet(norm_data_path)

norm_schema = pq.read_schema(norm_data_path)

# get a list of column names
norm_cols = [col.name for col in norm_schema]
print(len(norm_cols))
# get columns that contain Metadata
metadata_cols = [col for col in norm_cols if "Metadata" in col]
# remove metadata columns from the list of columns
data_cols = [col for col in norm_cols if col not in metadata_cols]

print(f"There are {len(data_cols)} data columns")
print(f"There are {len(metadata_cols)} metadata columns")


# In[ ]:


# ## Check feature selected shape

# In[4]:


norm_fs_schema = pq.read_schema(norm_fs_data_path)

# get a list of column names
norm_cols = [col.name for col in norm_schema]
print(len(norm_cols))
# get columns that contain Metadata
metadata_cols = [col for col in norm_cols if "Metadata" in col]
# remove metadata columns from the list of columns
data_cols = [col for col in norm_cols if col not in metadata_cols]

print(f"There are {len(data_cols)} data columns")
print(f"There are {len(metadata_cols)} metadata columns")

# get columns that contain Metadata

norm_fs_df_subset = pd.read_parquet(
    norm_fs_data_path,
    columns=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
)
print(norm_fs_df_subset.shape)
norm_fs_df_subset.head()


# In[5]:


# path to the ground truth file
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


# In[6]:


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


# ## LOCO ENET stats

# In[58]:


# set path for models performances
model_performances_path = pathlib.Path(
    "../../../11.bulk_Morphology_Elastic_Network_LOCO/2.test_models/results/regression/PBMC_aggregated_with_nomic/model_stats.csv"
).resolve(strict=True)

variance_r2_stats_path = pathlib.Path(
    "../../../11.bulk_Morphology_Elastic_Network_LOCO/2.test_models/results/regression/PBMC_aggregated_with_nomic/variance_r2_stats.csv"
).resolve(strict=True)

model_performances = pd.read_csv(model_performances_path)
print(model_performances.shape)
model_performances.head()


# In[59]:


variance_r2_stats = pd.read_csv(variance_r2_stats_path)
print(variance_r2_stats.shape)
variance_r2_stats.head()


# In[60]:


# get only select keys
model_performances = model_performances.loc[
    model_performances["channel_feature_combinations_key"].isin(
        [
            "All_channels",
            "CorrDNA_CorrGasdermin_CorrMito_CorrER",
            "CorrDNA_CorrPM_CorrGasdermin_CorrER",
            "CorrDNA_CorrPM_CorrGasdermin_CorrMito",
            "CorrDNA_CorrPM_CorrMito_CorrER",
            "CorrPM_CorrGasdermin_CorrMito_CorrER",
        ]
    )
]
# replace string values with more readable names
model_performances["channel_feature_combinations_key"] = model_performances[
    "channel_feature_combinations_key"
].replace(
    {
        "All_channels": "All channels",
        "CorrDNA_CorrGasdermin_CorrMito_CorrER": "PM removed",
        "CorrDNA_CorrPM_CorrGasdermin_CorrER": "Mito removed",
        "CorrDNA_CorrPM_CorrGasdermin_CorrMito": "ER removed",
        "CorrDNA_CorrPM_CorrMito_CorrER": "Gasdermin removed",
        "CorrPM_CorrGasdermin_CorrMito_CorrER": "DNA removed",
    }
)
model_performances["channel_feature_combinations_key"].unique()


# get only select keys
variance_r2_stats = variance_r2_stats.loc[
    variance_r2_stats["channel_feature_combinations_key"].isin(
        [
            "All_channels",
            "CorrDNA_CorrGasdermin_CorrMito_CorrER",
            "CorrDNA_CorrPM_CorrGasdermin_CorrER",
            "CorrDNA_CorrPM_CorrGasdermin_CorrMito",
            "CorrDNA_CorrPM_CorrMito_CorrER",
            "CorrPM_CorrGasdermin_CorrMito_CorrER",
        ]
    )
]

# replace string values with more readable names
variance_r2_stats["channel_feature_combinations_key"] = variance_r2_stats[
    "channel_feature_combinations_key"
].replace(
    {
        "All_channels": "All channels",
        "CorrDNA_CorrGasdermin_CorrMito_CorrER": "PM removed",
        "CorrDNA_CorrPM_CorrGasdermin_CorrER": "Mito removed",
        "CorrDNA_CorrPM_CorrGasdermin_CorrMito": "ER removed",
        "CorrDNA_CorrPM_CorrMito_CorrER": "Gasdermin removed",
        "CorrPM_CorrGasdermin_CorrMito_CorrER": "DNA removed",
    }
)

# drop the shuffled models
model_performances = model_performances.loc[model_performances["shuffle"] == "final"]
variance_r2_stats = variance_r2_stats.loc[variance_r2_stats["shuffle"] == "final"]
print(model_performances.shape)

print(variance_r2_stats.shape)


# In[61]:


model_performances
# get the explained variance, MSE, R2 for each cytokine, data split, channel combination
model_performances_grouped = model_performances.groupby(
    ["cytokine", "data_split", "channel_feature_combinations_key"]
).agg(
    {
        "explained_variance": "mean",
        "neg_mean_squared_error": "mean",
        "r2": "mean",
    }
)
model_performances_grouped.reset_index(inplace=True)
print(model_performances_grouped.shape)


# ## Stats for 11A-C

# In[63]:


# get the global average of neg mean squared error, explained variance, and r2 for each channel combination
channel_feature_combinations_key_global_avg = model_performances_grouped.groupby(
    "channel_feature_combinations_key"
).agg(
    {
        "explained_variance": "mean",
        "neg_mean_squared_error": "mean",
        "r2": "mean",
    }
)
channel_feature_combinations_key_global_avg.reset_index(inplace=True)
channel_feature_combinations_key_global_avg[
    "percent_change_in_negMSE_compared_to_all_channels"
] = (
    (
        channel_feature_combinations_key_global_avg["neg_mean_squared_error"]
        - channel_feature_combinations_key_global_avg.loc[
            channel_feature_combinations_key_global_avg[
                "channel_feature_combinations_key"
            ]
            == "All channels",
            "neg_mean_squared_error",
        ].values[0]
    )
    / channel_feature_combinations_key_global_avg.loc[
        channel_feature_combinations_key_global_avg["channel_feature_combinations_key"]
        == "All channels",
        "neg_mean_squared_error",
    ].values[0]
    * 100
)
channel_feature_combinations_key_global_avg[
    "percent_change_in_explained_variance_compared_to_all_channels"
] = (
    (
        channel_feature_combinations_key_global_avg["explained_variance"]
        - channel_feature_combinations_key_global_avg.loc[
            channel_feature_combinations_key_global_avg[
                "channel_feature_combinations_key"
            ]
            == "All channels",
            "explained_variance",
        ].values[0]
    )
    / channel_feature_combinations_key_global_avg.loc[
        channel_feature_combinations_key_global_avg["channel_feature_combinations_key"]
        == "All channels",
        "explained_variance",
    ].values[0]
    * 100
)
channel_feature_combinations_key_global_avg[
    "percent_change_in_r2_compared_to_all_channels"
] = (
    (
        channel_feature_combinations_key_global_avg["r2"]
        - channel_feature_combinations_key_global_avg.loc[
            channel_feature_combinations_key_global_avg[
                "channel_feature_combinations_key"
            ]
            == "All channels",
            "r2",
        ].values[0]
    )
    / channel_feature_combinations_key_global_avg.loc[
        channel_feature_combinations_key_global_avg["channel_feature_combinations_key"]
        == "All channels",
        "r2",
    ].values[0]
    * 100
)
channel_feature_combinations_key_global_avg


# In[64]:


# get the min and max r2 values for each channel combination
channel_feature_combinations_key_min_max = model_performances_grouped.groupby(
    "channel_feature_combinations_key"
).agg(
    {
        "r2": ["min", "max"],
    }
)
channel_feature_combinations_key_min_max.reset_index(inplace=True)
channel_feature_combinations_key_min_max


# In[ ]:


# subset for IL-1beta across all channel combinations
IL1beta_model_performances = model_performances_grouped.loc[
    model_performances_grouped["cytokine"] == "IL-1beta"
]
