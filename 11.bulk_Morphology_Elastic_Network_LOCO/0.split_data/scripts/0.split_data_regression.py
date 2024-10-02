#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

# In[2]:


# argparser = argparse.ArgumentParser()
# argparser.add_argument("--cell_type", default="all")

# args = argparser.parse_args()

# cell_type = args.cell_type

cell_type = "PBMC"


# In[3]:


# Parameters
aggregation = True
nomic = True


# In[4]:


MODEL_TYPE = "regression"


# In[5]:


# toml file path
TOML_PATH = pathlib.Path("../splits.toml")
# read toml file via toml
data_splits_by_treatments = toml.load(TOML_PATH)

# define the 100% test set data treatments
test_100_percent = data_splits_by_treatments["splits"]["data_splits_100"]
test_75_percent = data_splits_by_treatments["splits"]["data_splits_75"]


# In[6]:


path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
)

data_df = pd.read_parquet(path)

data_df.head()


# In[7]:


# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
]

data_df = pd.merge(data, metadata_well, left_index=True, right_index=True)


# This model and code is both inspired and reused from: https://github.com/WayScience/phenotypic_profiling_model/blob/main/1.split_data/split_data.ipynb
# The bulk of this work was done by **Roshan Kern** I have only made minor changes to the code to make it more modular and easier to use for my purposes.

# In[8]:


# variable test and train set splits
# 100% test set
# subset the following treatments for test set
test_set_all = data_df[
    data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_100_percent)
]
# 75% test set and 25% train set
test_set_75 = data_df[
    data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_75_percent)
]

test_100_and_75 = test_100_percent + test_75_percent

# 50% test set and 50% train set
# get all treatments that are not in the_test_set_all and the test_set_75
test_set_50 = data_df[
    ~data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_100_and_75)
]

print(test_set_all.shape, test_set_75.shape, test_set_50.shape)


# In[9]:


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


# In[10]:


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

print(
    f"Testing set length: {len(testing_data_set)}\nTraining set length: {len(training_data_set)}"
)

# get the indexes for the training and testing sets
testing_data_set_index = testing_data_set.index
training_data_set_index = training_data_set.index


# In[11]:


# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in training_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in testing_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "test"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data).sort_values(["labeled_data_index"])


# In[12]:


# set save path
if aggregation == True:
    if nomic == True:
        save_path = pathlib.Path(f"../indexes/{cell_type}/{MODEL_TYPE}")
    elif nomic == False:
        save_path = pathlib.Path(f"../indexes/{cell_type}/{MODEL_TYPE}")
elif aggregation == False:
    if nomic == True:
        save_path = pathlib.Path(f"../indexes/{cell_type}/{MODEL_TYPE}")
    elif nomic == False:
        save_path = pathlib.Path(f"../indexes/{cell_type}/{MODEL_TYPE}")
else:
    print("Error")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)


# In[13]:


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
