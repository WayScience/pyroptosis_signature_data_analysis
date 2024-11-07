#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import itertools
import pathlib

import joblib
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
import tqdm
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import explained_variance_score, mean_squared_error, r2_score
from sklearn.utils import parallel_backend

# In[2]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", type=str, default="all")
argparser.add_argument("--shuffle", type=str, default=False)
argparser.add_argument("--cytokine", type=str, default="cytokine")
argparser.add_argument("--feature_combinations_key", type=str, default="all")
argparser.add_argument("--data_split", type=str)

args = argparser.parse_args()

cell_type = args.cell_type
shuffle = args.shuffle
cytokine = args.cytokine
feature_combinations_key = args.feature_combinations_key
data_split = args.data_split


# cell_type = "PBMC"
# cytokine = "IL-1 beta [NSU]"
# shuffle = "False"
# feature_combinations_key = "CorrDNA"
# data_split = "train"

if shuffle == "True":
    shuffle = True
elif shuffle == "False":
    shuffle = False
else:
    raise ValueError("shuffle must be True or False")


# In[3]:


# Parameters
aggregation = True
nomic = True


# In[4]:


# set shuffle value
if shuffle == True:
    shuffle = "shuffled_baseline"
else:
    shuffle = "final"


# In[5]:


MODEL_TYPE = "regression"


# In[6]:


# load training data from indexes and features dataframe
data_split_path = pathlib.Path(
    f"../../0.split_data/indexes/{cell_type}/regression/aggregated_sc_and_nomic_data_split_indexes.tsv"
).resolve(strict=True)

feature_combinations_file_path = pathlib.Path(
    f"../../0.split_data/results/channel_feature_combinations_{cell_type}.toml"
).resolve(strict=True)

data_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
).resolve(strict=True)

feature_combination_key_file = pathlib.Path(
    "../../0.split_data/results/channel_feature_combinations_keys.txt"
).resolve(strict=True)

# load the feature combinations file
feature_combinations = toml.load(feature_combinations_file_path)
feature_combinations_columns = feature_combinations[feature_combinations_key]

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pd.read_parquet(data_path)
data_df = data_df[feature_combinations_columns]

data_split_indexes = pd.read_csv(data_split_path, sep="\t")


# In[7]:


# select tht indexes for the training and test set
data_split_indexes = data_split_indexes.loc[data_split_indexes["label"] == data_split]
# subset data_df by indexes in data_split_indexes
data_split_data = data_df.loc[data_split_indexes["labeled_data_index"]]
# define metadata columns
# subset each column that contains metadata
metadata = data_split_data.filter(regex="Metadata")
# drop all metadata columns
data_x = data_split_data.drop(metadata.columns, axis=1)
labeled_data = data_split_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name
data_y_cols = data_x.filter(regex="NSU").columns
data_split_y = data_split_data[data_y_cols]
data_split_x = data_x.drop(data_y_cols, axis=1)
# drop the oneb_Treatment_Dose_Inhibitor_Dose column if it exists
if "oneb_Treatment_Dose_Inhibitor_Dose" in data_split_x.columns:
    data_split_x = data_split_x.drop(columns="oneb_Treatment_Dose_Inhibitor_Dose")


# In[8]:


print(data_split_x.shape, data_split_y.shape)


# In[9]:


# set model path from parameters
if (aggregation is True) and (nomic is True):
    model_path = pathlib.Path(f"models/regression/{cell_type}/aggregated_with_nomic/")
elif (aggregation is True) and (nomic is False):
    model_path = pathlib.Path(f"models/regression/{cell_type}/aggregated/")
elif (aggregation is False) and (nomic is True):
    model_path = pathlib.Path(f"models/regression/{cell_type}/sc_with_nomic/")
elif (aggregation is False) and (nomic is False):
    model_path = pathlib.Path(f"models/regression/{cell_type}/sc/")
else:
    print("Error")


# In[10]:


data_dict = {
    "test_data": {
        "data_x": data_split_x,
        "data_y": data_split_y,
        "col_names": data_y_cols,
        "metadata": metadata,
    },
}


# In[11]:


# list of metrics to use
output_metric_scores = {}


# In[12]:


# blank df for concatenated results
results_df = pd.DataFrame(
    columns=[
        "explained_variance",
        "neg_mean_absolute_error",
        "neg_mean_squared_error",
        "well",
        "treatment",
        "r2",
        "cytokine",
        "data_split",
        "shuffle",
        "predicted_value",
        "actual_value",
        "log10_neg_mean_absolute_error",
        "log10_neg_mean_squared_error",
        "log10_explained_variance",
    ]
)


# In[13]:


data_x = data_split_x
data_y = data_split_y
if shuffle == "shuffled_baseline":
    model = joblib.load(
        f"../../1.train_models/{model_path}/{cytokine}_{feature_combinations_key}_shuffled_baseline__all_nomic.joblib"
    )
elif shuffle == "final":
    model = joblib.load(
        f"../../1.train_models/{model_path}/{cytokine}_{feature_combinations_key}_final__all_nomic.joblib"
    )
else:
    print("Error")

# get the cytokine column of choice
y_selected = data_y[cytokine]

if shuffle == "shuffled_baseline":
    for column in data_x:
        np.random.shuffle(data_x[column].values)

# get predictions
predictions = model.predict(data_x)


# In[14]:


explained_variance = explained_variance_score(y_selected, predictions)
output_metric_scores["explained_variance"] = explained_variance
neg_mean_absolute_error = -mean_squared_error(y_selected, predictions)
output_metric_scores["neg_mean_absolute_error"] = neg_mean_absolute_error
neg_mean_squared_error = -mean_squared_error(y_selected, predictions)
output_metric_scores["neg_mean_squared_error"] = neg_mean_squared_error
r2 = r2_score(y_selected, predictions)
output_metric_scores["treatment"] = metadata[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
output_metric_scores["well"] = metadata["Metadata_Well"].values
df = pd.DataFrame.from_dict(output_metric_scores)
df["r2"] = r2
df["cytokine"] = cytokine
df["data_split"] = data_split
df["shuffle"] = shuffle
df["predicted_value"] = predictions
df["actual_value"] = y_selected
df["log10_neg_mean_absolute_error"] = -np.log10(-df["neg_mean_absolute_error"])
df["log10_neg_mean_squared_error"] = -np.log10(-df["neg_mean_squared_error"])
df["log10_explained_variance"] = -np.log10(df["explained_variance"])

# replace "[NSU]" with """
df["cytokine"] = df["cytokine"].replace("[ \[\]NSU]", "", regex=True)
df["cytokine"] = df["cytokine"].replace(" ", "_", regex=True)

# concat the dataframes
results_df = pd.concat([results_df, df], axis=0)
results_df["channel_feature_combinations_key"] = feature_combinations_key


# In[15]:


results_df.head()


# In[16]:


var_df = results_df.drop(
    columns=[
        "explained_variance",
        "neg_mean_absolute_error",
        "neg_mean_squared_error",
        "well",
        "treatment",
        "r2",
        "log10_neg_mean_absolute_error",
        "log10_neg_mean_squared_error",
        "log10_explained_variance",
    ]
)
# calculate the variance of the actual and predicted values per cytokine
var_df = var_df.groupby(
    ["cytokine", "data_split", "shuffle", "channel_feature_combinations_key"]
).var()
var_df = pd.merge(
    var_df,
    results_df.groupby(
        ["cytokine", "data_split", "shuffle", "channel_feature_combinations_key"]
    ).r2.unique(),
    left_index=True,
    right_index=True,
)
var_df.reset_index(inplace=True)
var_df.head()


# In[17]:


# set model path from parameters
if aggregation and nomic:
    results_path = pathlib.Path(
        f"../results/regression/{cell_type}_aggregated_with_nomic/"
    )
elif (aggregation is True) and (nomic is False):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_aggregated/")
elif (aggregation is False) and (nomic is True):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_sc_with_nomic/")
elif (aggregation is False) and (nomic is False):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_sc/")
else:
    print("Error")
pathlib.Path(results_path).mkdir(parents=True, exist_ok=True)


# In[18]:


# check if the model training metrics file exists
metrics_file = pathlib.Path(
    f"{results_path}/{cytokine}_{shuffle}_{feature_combinations_key}_{data_split}_data_model_stats.csv"
)

results_df.to_csv(metrics_file, index=False)

# do the same for the variance df
# check if the model training metrics file exists
metrics_file = pathlib.Path(
    f"{results_path}/{cytokine}_{shuffle}_{feature_combinations_key}_{data_split}_data_variance_r2_stats.csv"
)
var_df.to_csv(metrics_file, index=False)
