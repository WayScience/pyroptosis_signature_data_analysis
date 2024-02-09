#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging
import pathlib
import sys
from typing import Optional

import numpy as np
import pandas as pd
import toml
from copairs.map import run_pipeline
from pycytominer import feature_select

# imports src
sys.path.append("../")
from src import utils

# setting up logger
logging.basicConfig(
    filename="map_analysis_testing.log",
    level=logging.DEBUG,
    format="%(levelname)s:%(asctime)s:%(name)s:%(message)s",
)


# ## Helper functions
# Set of helper functions to help out throughout the notebook

# In[2]:


## Helper function


def shuffle_meta_labels(
    dataset: pd.DataFrame, target_col: str, seed: Optional[int] = 0
) -> pd.DataFrame:
    """shuffles labels or values within a single selected column

    Parameters
    ----------
    dataset : pd.DataFrame
        dataframe containing the dataset

    target_col : str
        Column to select in order to conduct the shuffling

    seed : int
        setting random seed

    Returns
    -------
    pd.DataFrame
        shuffled dataset

    Raises
    ------
    TypeError
        raised if incorrect types are provided
    """
    # setting seed
    np.random.seed(seed)

    # type checking
    if not isinstance(target_col, str):
        raise TypeError("'target_col' must be a string type")
    if not isinstance(dataset, pd.DataFrame):
        raise TypeError("'dataset' must be a pandas dataframe")

    # selecting column, shuffle values within column, add to dataframe
    dataset[target_col] = np.random.permutation(dataset[target_col].values)
    return dataset


def shuffle_features(feature_vals: np.array, seed: Optional[int] = 0) -> np.array:
    """suffles all values within feature space

    Parameters
    ----------
    feature_vals : np.array
        shuffled

    seed : Optional[int]
        setting random seed

    Returns
    -------
    np.array
        Returns shuffled values within the feature space

    Raises
    ------
    TypeError
        Raised if a numpy array is not provided
    """
    # setting seed
    np.random.seed(seed)

    # shuffle given array
    if not isinstance(feature_vals, np.ndarray):
        raise TypeError("'feature_vals' must be a numpy array")
    if feature_vals.ndim != 2:
        raise TypeError("'feature_vals' must be a 2x2 matrix")

    # creating a copy for feature vales to prevent overwriting of global variables
    feature_vals = np.copy(feature_vals)

    # shuffling feature space
    n_cols = feature_vals.shape[1]
    for col_idx in range(0, n_cols):
        # selecting column, shuffle, and update:
        feature_vals[:, col_idx] = np.random.permutation(feature_vals[:, col_idx])

    return feature_vals


# ## Setting up Paths and loading data

# In[3]:


# load in the treatment groups
ground_truth = pathlib.Path(
    "../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
).resolve(strict=True)
# load in the ground truth
ground_truth = toml.load(ground_truth)
apoptosis_ground_truth = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_ground_truth = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control_ground_truth = ground_truth["Healthy"]["healthy_groups_list"]


# In[4]:


single_cell_data = pathlib.Path(
    f"../../data/PBMC_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)
df = pd.read_parquet(single_cell_data)


# In[5]:


# out paths
map_out_dir = pathlib.Path("../data/processed/mAP_scores/morphology/")
map_out_dir.mkdir(exist_ok=True, parents=True)

# regular data output
# saving to csv
regular_feat_map_path = pathlib.Path(map_out_dir / "mAP_scores_regular_class.csv")

# shuffled data output
shuffled_feat_map_path = pathlib.Path(map_out_dir / "mAP_scores_shuffled_class.csv")

# shuffled feature space output
shuffled_feat_space_map_path = pathlib.Path(
    map_out_dir / "mAP_scores_shuffled_feature_space_class.csv"
)


# ### Clean up data

# In[6]:


# add apoptosis, pyroptosis and healthy columns to dataframe
df["Apoptosis"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in apoptosis_ground_truth,
    axis=1,
)
df["Pyroptosis"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in pyroptosis_ground_truth,
    axis=1,
)
df["Control"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in control_ground_truth,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column
df["Metadata_labels"] = df.apply(
    lambda row: "Apoptosis"
    if row["Apoptosis"]
    else "Pyroptosis"
    if row["Pyroptosis"]
    else "Control",
    axis=1,
)
# # drop apoptosis, pyroptosis, and healthy columns
df.drop(columns=["Apoptosis", "Pyroptosis", "Control"], inplace=True)
df.drop(columns=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"], inplace=True)


# In[7]:


# output directories
map_out_dir = pathlib.Path("../data/processed/mAP_scores/")
map_out_dir.mkdir(parents=True, exist_ok=True)


# ### mAP Pipeline Parameters

# The null size needs to be determined for the mAP pipeline. This is the size of the null class that is used to determine the mAP score.

# In[8]:


tmp = (
    df.groupby(["Metadata_labels"])
    .count()
    .reset_index()[["Metadata_Well", "Metadata_labels"]]
)
# get the Pyroptosis number of Metadata_Well
Pyroptosis_count = tmp[tmp["Metadata_labels"] == "Pyroptosis"]["Metadata_Well"].values[
    0
]
Pyroptosis_count


# In[9]:


pos_sameby = [
    "Metadata_labels",
]
pos_diffby = ["Metadata_Well"]

neg_sameby = []
neg_diffby = ["Metadata_labels"]

null_size = Pyroptosis_count
batch_size = 1

# number of resampling
n_resamples = 10


# ### mAP analysis for non-shuffled data

# In[10]:


# This will generated 100 values [0..100] as seed values
# This will occur per phenotype

# spliting metadata and raw feature values
logging.info("splitting data set into metadata and raw feature values")
df_meta, df_feats = utils.split_data(df)
df_feats = np.array(df_feats)

# execute pipeline on negative control with training dataset with cp features
# print(negative_training_cp_meta)
# print(negative_training_cp_feats)
try:
    # execute pipeline on negative control with trianing dataset with cp features
    # print(negative_training_cp_meta)
    # print(negative_training_cp_feats)
    logging.info(f"Running pipeline on CP features using phenotype")
    result = run_pipeline(
        meta=df_meta,
        feats=df_feats,
        pos_sameby=pos_sameby,
        pos_diffby=pos_diffby,
        neg_sameby=neg_sameby,
        neg_diffby=neg_diffby,
        batch_size=batch_size,
        null_size=null_size,
    )

    # adding columns
    result["shuffled"] = "non-shuffled"


except ZeroDivisionError as e:
    logging.warning(f"{e} captured on phenotye:. Skipping")
# concatenating all datasets
result.to_csv(regular_feat_map_path, index=False)
result.head()


# ### mAP analysis for shuffled data (Phenotype)

# In[11]:


logging.info("Running mAP pipeline with shuffled phenotype labeled data")
seed = 0
# running process
# for loop selects one single phenotype
# then splits the data into metadata and raw feature values
# two different groups that contains 3 splits caused by the types of features

# This will generated 100 values [0..100] as seed values
# splitting metadata labeled shuffled data
logging.info("splitting shuffled data set into metadata and raw feature values")
df = shuffle_meta_labels(dataset=df, target_col="Metadata_labels", seed=seed)
(
    df_meta,
    df_feats,
) = utils.split_data(df)


df_feats = np.array(df_feats)

try:
    # execute pipeline on negative control with trianing dataset with cp features
    logging.info(
        f"Running pipeline on CP features using  phenotype, data is shuffled by phenoptype labels"
    )
    shuffled_result = run_pipeline(
        meta=df_meta,
        feats=df_feats,
        pos_sameby=pos_sameby,
        pos_diffby=pos_diffby,
        neg_sameby=neg_sameby,
        neg_diffby=neg_diffby,
        batch_size=batch_size,
        null_size=null_size,
    )

    # adding shuffle label column
    shuffled_result["shuffled"] = "phenotype_shuffled"


except ZeroDivisionError as e:
    logging.warning(f"{e} captured on phenotye: Skipping")

# saving to csv
shuffled_result.to_csv(shuffled_feat_map_path, index=False)
shuffled_result


# ### mAP analysis for shuffled data (Feature space)

# In[12]:


seed = 0


# split the shuffled dataset
# spliting metadata and raw feature values
logging.info("splitting shuffled data set into metadata and raw feature values")
(
    df_meta,
    df_feats,
) = utils.split_data(df)

df_feats = np.array(df_feats)


# shuffling the features, this will overwrite the generated feature space from above with the shuffled one
df_feats = shuffle_features(feature_vals=df_feats, seed=seed)


try:
    # execute pipeline on negative control with trianing dataset with cp features
    logging.info(
        f"Running pipeline on CP features using phenotype, feature space is shuffled"
    )
    shuffled_feat_map = run_pipeline(
        meta=df_meta,
        feats=df_feats,
        pos_sameby=pos_sameby,
        pos_diffby=pos_diffby,
        neg_sameby=neg_sameby,
        neg_diffby=neg_diffby,
        batch_size=batch_size,
        null_size=null_size,
    )

    # adding shuffle label column
    shuffled_feat_map["shuffled"] = "features_shuffled"


except ZeroDivisionError as e:
    logging.warning(f"{e} captured on phenotype:  Skipping")


# saving to csv
shuffled_feat_map.to_csv(shuffled_feat_space_map_path, index=False)
shuffled_feat_map
