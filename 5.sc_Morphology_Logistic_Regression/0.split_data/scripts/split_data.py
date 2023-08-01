#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ast
import itertools
import pathlib

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.utils import parallel_backend, shuffle

# In[2]:


# Parameters
cell_type = "PBMC"
aggregation = False
nomic = False
flag = True
control = "DMSO_0.100_DMSO_0.025"
treatment = "Thapsigargin_1.000_DMSO_0.025"


# In[3]:


if flag == False:
    # read in toml file and get parameters
    toml_path = pathlib.Path("../1.train_models/single_class_config.toml")
    with open(toml_path, "r") as f:
        config = toml.load(f)
    f.close()
    control = config["logistic_regression_params"]["control"]
    treatment = config["logistic_regression_params"]["treatments"]
    aggregation = ast.literal_eval(config["logistic_regression_params"]["aggregation"])
    nomic = ast.literal_eval(config["logistic_regression_params"]["nomic"])
    cell_type = config["logistic_regression_params"]["cell_type"]
    print(aggregation, nomic, cell_type)


# In[4]:


path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet")

data_df = pq.read_table(path).to_pandas()

data_df.head()


# In[5]:


if nomic == True:
    nomic_df_path = pathlib.Path(
        f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
    )
    df_nomic = pd.read_csv(nomic_df_path)
    # drop columns that contain [pgML]
    df_nomic = df_nomic.drop(
        columns=[col for col in df_nomic.columns if "[pgML]" in col]
    )
    # drop first 25 columns
    df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
    df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])
else:
    df_nomic = None


# In[6]:


# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata["Metadata_Well"]

data = pd.merge(data, metadata_well, left_index=True, right_index=True)


# In[7]:


if (aggregation == True) and (nomic == True):

    # subset each column that contains metadata
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df, metadata, left_on="Metadata_Well", right_on="Metadata_Well"
    )
    data_df = pd.merge(
        data_df, df_nomic, left_on="Metadata_Well", right_on="position_x"
    )
    data_df = data_df.drop(columns=["position_x"])

elif (aggregation == True) and (nomic == False):
    # subset each column that contains metadata
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df, metadata, left_on="Metadata_Well", right_on="Metadata_Well"
    )
elif (aggregation == False) and (nomic == True):
    data_df = pd.merge(
        data_df, df_nomic, left_on="Metadata_Well", right_on="position_x"
    )
    data_df = data_df.drop(columns=["position_x"])
elif aggregation == False and nomic == False:
    pass
else:
    print("Error")


# In[8]:


# drop all metadata columns
data_x = data_df.drop(metadata.columns, axis=1)
labeled_data = data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]


# In[9]:


# https://github.com/WayScience/phenotypic_profiling_model/blob/main/1.split_data/split_data.ipynb


# In[10]:


# get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
data_df = data_df[
    data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin([control, treatment])
]


# In[11]:


# ratio of data to be used for testing (ex 0.15 = 15%)
test_ratio = 0.25

# get indexes of training and testing data
training_data, testing_data = train_test_split(
    data_df,
    test_size=test_ratio,
    stratify=data_df[["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]],
    random_state=1,
)
train_indexes = training_data.index.to_numpy()
test_indexes = testing_data.index.to_numpy()

print(f"Training data has shape: {training_data.shape}")
print(f"Testing data has shape: {testing_data.shape}")


# In[12]:


# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in train_indexes:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in test_indexes:
    index_data.append({"labeled_data_index": index, "label": "test"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data).sort_values(["labeled_data_index"])


# In[13]:


# set save path
if aggregation == True:
    if nomic == True:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{control}_{treatment}")
    elif nomic == False:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{control}_{treatment}")
elif aggregation == False:
    if nomic == True:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{control}_{treatment}")
    elif nomic == False:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{control}_{treatment}")
else:
    print("Error")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)


# In[14]:


# save indexes as tsv file
if aggregation == True:
    if nomic == True:
        index_data.to_csv(
            f"{save_path}/aggregated_sc_and_nomic_data_split_indexes.tsv", sep="\t"
        )
    elif nomic == False:
        index_data.to_csv(f"{save_path}/aggregated_sc_data_split_indexes.tsv", sep="\t")
elif aggregation == False:
    if nomic == True:
        index_data.to_csv(f"{save_path}/sc_and_nomic_data_split_indexes.tsv", sep="\t")
    elif nomic == False:
        index_data.to_csv(f"{save_path}/sc_split_indexes.tsv", sep="\t")
else:
    print("Error")
