# ---
# jupyter:
#   jupytext:
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

# %% papermill={"duration": 1.193087, "end_time": "2023-07-22T05:40:47.784792", "exception": false, "start_time": "2023-07-22T05:40:46.591705", "status": "completed"}
import itertools
import pathlib
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV, LogisticRegression, MultiTaskElasticNetCV

# import RepeatedKFold
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedKFold,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.utils import parallel_backend, shuffle

# import mse

# %% papermill={"duration": 0.005605, "end_time": "2023-07-22T05:40:47.792730", "exception": false, "start_time": "2023-07-22T05:40:47.787125", "status": "completed"} tags=["injected-parameters"]
# Parameters
cell_type = "SHSY5Y"
aggregation = True
nomic = True
flag = True
shuffle = True

# %%
# set shuffle value
if shuffle:
    shuffle = "shuffled_baseline"
else:
    shuffle = "final"

# %% papermill={"duration": 0.005391, "end_time": "2023-07-22T05:40:47.799818", "exception": false, "start_time": "2023-07-22T05:40:47.794427", "status": "completed"}
MODEL_TYPE = "regression"

# %% papermill={"duration": 214.792447, "end_time": "2023-07-22T05:44:22.593856", "exception": false, "start_time": "2023-07-22T05:40:47.801409", "status": "completed"}
# load training data from indexes and features dataframe
# data_split_path = pathlib.Path(f"../0.split_data/indexes/data_split_indexes.tsv")
data_path = pathlib.Path(f"../../../data/{cell_type}_preprocessed_sc_norm.parquet")

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pq.read_table(data_path).to_pandas()

# import nomic data
nomic_df_path = pathlib.Path(
    f"../../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_cleanup4correlation.csv"
)
df_nomic = pd.read_csv(nomic_df_path)

# clean up nomic data
df_nomic = df_nomic.drop(columns=[col for col in df_nomic.columns if "[pgML]" in col])
# drop first 25 columns (Metadata that is not needed)
# df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
# df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])

# %%
print(df_nomic["Activin A [NSU]"].describe())

# %% papermill={"duration": 0.162062, "end_time": "2023-07-22T05:44:22.812544", "exception": false, "start_time": "2023-07-22T05:44:22.650482", "status": "completed"}
if (aggregation == True) and (nomic == True):
    data_split_path = pathlib.Path(
        f"../../0.split_data/indexes/{cell_type}/{MODEL_TYPE}/aggregated_sc_and_nomic_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
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
elif (aggregation == True) and (nomic == False):
    data_split_path = pathlib.Path(
        f"../../0.split_data/indexes/{cell_type}/{MODEL_TYPE}/aggregated_sc_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
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
    data_split_path = pathlib.Path(
        f"../../0.split_data/indexes/{cell_type}/{MODEL_TYPE}/sc_and_nomic_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["Metadata_position_x", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["Metadata_position_x"])
elif aggregation == False and nomic == False:
    data_split_path = pathlib.Path(
        f"../../0.split_data/indexes/{cell_type}/{MODEL_TYPE}/sc_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
else:
    print("Error")

# %% papermill={"duration": 0.028408, "end_time": "2023-07-22T05:44:22.842747", "exception": false, "start_time": "2023-07-22T05:44:22.814339", "status": "completed"}
# select tht indexes for the training and test set
train_indexes = data_split_indexes.loc[data_split_indexes["label"] == "train"]

# %% papermill={"duration": 0.529615, "end_time": "2023-07-22T05:44:23.373993", "exception": false, "start_time": "2023-07-22T05:44:22.844378", "status": "completed"}
# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[train_indexes["labeled_data_index"]]

# %%
training_data.head()

# %% papermill={"duration": 0.352839, "end_time": "2023-07-22T05:44:23.728560", "exception": false, "start_time": "2023-07-22T05:44:23.375721", "status": "completed"}
# # get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
# training_data = training_data[
#     training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
#         [control, treatment]
#     )
# ]

# %%
# TODO holdout certain treatments and different percentages of hold out for each treatment
# where holdout is the test set

# %% papermill={"duration": 1.687722, "end_time": "2023-07-22T05:44:25.418506", "exception": false, "start_time": "2023-07-22T05:44:23.730784", "status": "completed"}
# # at random downsample the DMSO treatment to match the number of wells in the LPS treatment
# seed = 0
# # get the number of wells in the LPS treatment
# trt_wells = training_data[
#     training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatment
# ].shape[0]
# # get the number of wells in the DMSO treatment
# dmso_wells = training_data[
#     training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
# ].shape[0]
# # downsample the DMSO treatment to match the number of wells in the LPS treatment
# dmso_holdout = training_data[
#     training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
# ].sample(n=trt_wells, random_state=seed)
# # remove the downsampled DMSO wells from the data
# training_data = training_data.drop(dmso_holdout.index)

# %% papermill={"duration": 0.179094, "end_time": "2023-07-22T05:44:25.601137", "exception": false, "start_time": "2023-07-22T05:44:25.422043", "status": "completed"}
# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name
data_y_cols = data_x.filter(regex="NSU").columns
train_y = training_data[data_y_cols]
train_x = data_x.drop(data_y_cols, axis=1)

# %%
from sklearn.model_selection import LeaveOneOut

loo = LeaveOneOut()
loo.get_n_splits(train_x)
loo.get_n_splits(train_y)

# %%
for cytokine in train_y.columns:
    train_data_y = train_y[cytokine]
    model = ElasticNetCV(
        random_state=0,
        max_iter=100000,
        cv=loo,
        l1_ratio=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.99],
        alphas=[0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
        fit_intercept=True,
        selection="random",
    )
    # train model on training data on all combinations of model types, feature types, and phenotypic classes

    if shuffle == "shuffled_baseline":
        print("Shuffling data")
        for column in train_x:
            np.random.shuffle(train_x[column].values)
    else:
        print("Not shuffling data")
    # define parameters to search over
    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )
            # create a logistic regression model
            model.fit(train_x, train_data_y)
            scores = cross_val_score(
                model,
                train_x,
                train_data_y,
                scoring="neg_mean_absolute_error",
                cv=loo,
                n_jobs=-1,
            )
            print(scores)
            print(f"Mean MAE: {scores.mean()}")
            print(f"Std MAE: {scores.std()}")
            print(f"R2: {model.score(train_x, train_data_y)}")

    if (aggregation == True) and (nomic == True):
        results_dir = f"../models/regression/{cell_type}/aggregated_with_nomic/"
    elif (aggregation == True) and (nomic == False):
        results_dir = f"../models/regression/{cell_type}/aggregated/"
    elif (aggregation == False) and (nomic == True):
        results_dir = f"../models/regression/{cell_type}/sc_with_nomic/"
    elif (aggregation == False) and (nomic == False):
        results_dir = f"../models/regression/{cell_type}/sc/"
    else:
        print("Error")

    # create results directory if it doesn't exist
    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

    # save final estimator
    if shuffle == "shuffled_baseline":
        dump(
            model,
            f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.joblib",
        )
    elif shuffle == "final":
        dump(
            model,
            f"{results_dir}/{cytokine}_final__all_nomic.joblib",
        )
    else:
        print("Error")

    # save condfig copy specific to this model to the folder with the results
    # use pathlib
    if shuffle == "shuffled_baseline":
        config_copy_path = pathlib.Path(
            f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.toml"
        )
    elif shuffle == "final":
        config_copy_path = pathlib.Path(
            f"{results_dir}/{cytokine}_final__all_nomic.toml"
        )
    else:
        print("Error")

    # write toml file with parameters used from injected parameters

    with open(config_copy_path, "w") as f:
        f.write(f"model_type='{shuffle}'\n")
        f.write(f"aggregation={aggregation}\n")
        f.write(f"nomic={nomic}\n")
        f.write(f"cell_type='{cell_type}'\n")
        f.write(f"feature=all\n")
