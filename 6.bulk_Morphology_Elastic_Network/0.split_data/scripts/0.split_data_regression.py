# ---
# jupyter:
#   jupytext:
#     formats: ipynb,../scripts//py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Interstellar
#     language: python
#     name: python3
# ---

# %% papermill={"duration": 1.020566, "end_time": "2023-07-23T19:03:53.072209", "exception": false, "start_time": "2023-07-23T19:03:52.051643", "status": "completed"}
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

# %% papermill={"duration": 0.005665, "end_time": "2023-07-23T19:03:53.081555", "exception": false, "start_time": "2023-07-23T19:03:53.075890", "status": "completed"} tags=["injected-parameters"]
# Parameters
cell_type = "SHSY5Y"
aggregation = True
nomic = True
flag = True

# %% papermill={"duration": 0.005655, "end_time": "2023-07-23T19:03:53.088993", "exception": false, "start_time": "2023-07-23T19:03:53.083338", "status": "completed"}
MODEL_TYPE = "regression"

# %%
# toml file path
TOML_PATH = pathlib.Path("../splits.toml")
# read toml file via toml
data_splits_by_treatments = toml.load(TOML_PATH)

# define the 100% test set data treatments
test_100_percent = data_splits_by_treatments["splits"]["data_splits_100"]
test_75_percent = data_splits_by_treatments["splits"]["data_splits_75"]

# %% papermill={"duration": 82.220985, "end_time": "2023-07-23T19:05:15.311705", "exception": false, "start_time": "2023-07-23T19:03:53.090720", "status": "completed"}
path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet")
# path = pathlib.Path(
# "../../data/PBMC_subset_sc_norm_DMSO_0.100_DMSO_0.025_LPS_100.000_DMSO_0.025.parquet"
# )

data_df = pq.read_table(path).to_pandas()

data_df.head()

# %% papermill={"duration": 0.007743, "end_time": "2023-07-23T19:05:15.335255", "exception": false, "start_time": "2023-07-23T19:05:15.327512", "status": "completed"}
if nomic == True:
    # import nomic data
    nomic_df_path = pathlib.Path(
        f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_cleanup4correlation.csv"
    )
    df_nomic = pd.read_csv(nomic_df_path)

    # drop columns that contain [pgML]
    df_nomic = df_nomic.drop(
        columns=[col for col in df_nomic.columns if "[pgML]" in col]
    )
    # drop first 25 columns
    # df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
    # df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])
else:
    df_nomic = None

# %% papermill={"duration": 14.650757, "end_time": "2023-07-23T19:05:29.988218", "exception": false, "start_time": "2023-07-23T19:05:15.337461", "status": "completed"}
# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
]

data_df = pd.merge(data, metadata_well, left_index=True, right_index=True)

# %% papermill={"duration": 0.007782, "end_time": "2023-07-23T19:05:30.004155", "exception": false, "start_time": "2023-07-23T19:05:29.996373", "status": "completed"}
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
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["Metadata_position_x", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["Metadata_position_x"])
    # drop all metadata columns
    labeled_data = data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    data_x = data_df.drop(metadata.columns, axis=1)

elif (aggregation == True) and (nomic == False):
    # subset each column that contains metadata
    metadata = data.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["Metadata_position_x", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    )
elif (aggregation == False) and (nomic == True):
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["Metadata_position_x", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["Metadata_position_x"])
elif aggregation == False and nomic == False:
    pass
else:
    print("Error")

# %% [markdown] papermill={"duration": 0.006981, "end_time": "2023-07-23T19:05:37.283405", "exception": false, "start_time": "2023-07-23T19:05:37.276424", "status": "completed"}
# This model and code is both inspired and reused from: https://github.com/WayScience/phenotypic_profiling_model/blob/main/1.split_data/split_data.ipynb
# The bulk of this work was done by **Roshan Kern** I have only made minor changes to the code to make it more modular and easier to use for my purposes.

# %% papermill={"duration": 0.696472, "end_time": "2023-07-23T19:05:37.981908", "exception": false, "start_time": "2023-07-23T19:05:37.285436", "status": "completed"}
# get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
# data_df = data_df[
#     data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin([control, treatment])
# ]

# %%
data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()

# %%
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

# 50% test set and 50% train set
# get all treatments that are not in the_test_set_all and the test_set_75
test_set_50 = data_df[
    ~data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        test_set_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    )
]
test_set_50 = test_set_50[
    ~test_set_50["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        test_set_75["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    )
]

print(test_set_all.shape, test_set_75.shape, test_set_50.shape)

# %%
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

# %%
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

# %% papermill={"duration": 4.305535, "end_time": "2023-07-23T19:05:42.289632", "exception": false, "start_time": "2023-07-23T19:05:37.984097", "status": "completed"}
# # ratio of data to be used for testing (ex 0.15 = 15%)
# test_ratio = 0.5

# # get indexes of training and testing data
# training_data, testing_data = train_test_split(
#     data_df,
#     test_size=test_ratio,
#     stratify=data_df[["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]],
#     random_state=0,
# )
# train_indexes = training_data.index.to_numpy()
# test_indexes = testing_data.index.to_numpy()

# print(f"Training data has shape: {training_data.shape}")
# print(f"Testing data has shape: {testing_data.shape}")

# %% papermill={"duration": 0.808538, "end_time": "2023-07-23T19:05:43.104343", "exception": false, "start_time": "2023-07-23T19:05:42.295805", "status": "completed"}
# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in training_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in testing_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "test"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data).sort_values(["labeled_data_index"])

# %% papermill={"duration": 0.008037, "end_time": "2023-07-23T19:05:43.115913", "exception": false, "start_time": "2023-07-23T19:05:43.107876", "status": "completed"}
# set save path
if aggregation == True:
    if nomic == True:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{MODEL_TYPE}")
    elif nomic == False:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{MODEL_TYPE}")
elif aggregation == False:
    if nomic == True:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{MODEL_TYPE}")
    elif nomic == False:
        save_path = pathlib.Path(f"./indexes/{cell_type}/{MODEL_TYPE}")
else:
    print("Error")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)

# %% papermill={"duration": 0.320546, "end_time": "2023-07-23T19:05:43.438785", "exception": false, "start_time": "2023-07-23T19:05:43.118239", "status": "completed"}
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
