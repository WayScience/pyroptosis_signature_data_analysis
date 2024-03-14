#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
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
    f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_PBMC_clean.parquet"
).resolve(strict=True)
df = pd.read_parquet(single_cell_data)


# In[5]:


df


# In[6]:


# rename columns
df = df.rename(
    columns={
        "position_x": "Metadata_Well",
        "oneb_Treatment_Dose_Inhibitor_Dose": "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    }
)


# In[7]:


# out paths
map_out_dir = pathlib.Path("../data/processed/mAP_scores/secretome/")
map_out_dir.mkdir(exist_ok=True, parents=True)

# regular data output
# saving to csv
regular_feat_map_path = pathlib.Path(map_out_dir / "mAP_scores_regular_treatment.csv")

# shuffled data output
shuffled_feat_map_path = pathlib.Path(map_out_dir / "mAP_scores_shuffled_treatment.csv")

# shuffled feature space output
shuffled_feat_space_map_path = pathlib.Path(
    map_out_dir / "mAP_scores_shuffled_feature_space_treatment.csv"
)


# ### Clean up data

# In[8]:


# keep columns that contain Metdata and ['NSU']
df = df.filter(regex="Metadata|NSU")
df.head()

df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
# replace values in the oneb_Metadata_Treatment_Dose_Inhibitor_Dose column
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_0.100_ug_per_ml_DMSO_0.000_%", "Flagellin_0.100_ug_per_ml_DMSO_0.025_%"
)
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_0.100_ug_per_ml_DMSO_0.0_%", "Flagellin_0.100_ug_per_ml_DMSO_0.025_%"
)
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_1.000_ug_per_ml_DMSO_0.0_%", "Flagellin_1.000_ug_per_ml_DMSO_0.025_%"
)
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("Flagellin_1.000_0_DMSO_0.025_%", "Flagellin_1.000_ug_per_ml_DMSO_0.025_%")
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_1.000_ug_per_ml_DMSO_0.000_%", "Flagellin_1.000_ug_per_ml_DMSO_0.025_%"
)
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("media_ctr_0.0_0_Media_0_0", "Media")
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("media_ctr_0.0_0_Media_ctr_0.0_0", "Media")
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_1.000_0_Disulfiram_1.000_uM",
    "Flagellin_1.000_ug_per_ml_Disulfiram_1.000_uM",
)
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("media_ctr_0.0_0_Media_0.0_0", "Media")
print(len(df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()))


# In[9]:


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


# In[10]:


# output directories
map_out_dir = pathlib.Path("../data/processed/mAP_scores/")
map_out_dir.mkdir(parents=True, exist_ok=True)


# ### mAP Pipeline Parameters

# The null size needs to be determined for the mAP pipeline. This is the size of the null class that is used to determine the mAP score.

# In[11]:


tmp = (
    df.groupby(["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"])
    .count()
    .reset_index()[["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]]
)
# get the Pyroptosis number of Metadata_Well
# get the counts of each oneb_Metadata_Treatment_Dose_Inhibitor_Dose
min_count = tmp["Metadata_Well"].min()
print(min_count)
tmp


# In[12]:


# generate the permutations of cell death labels via itertools
pos_samby_permutations = list(itertools.combinations(df["Metadata_labels"].unique(), 2))
pos_samby_permutations


# In[13]:


pos_sameby = ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
pos_diffby = ["Metadata_Well"]

neg_sameby = ["Metadata_labels"]
neg_diffby = [
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
]

null_size = min_count
batch_size = 1

# number of resampling
n_resamples = 10


# ### mAP analysis for non-shuffled data

# In[14]:


results_df = pd.DataFrame(
    columns=[
        "Metadata_Well",
        "Metadata_labels",
        "average_precision",
        "p_value",
        "n_pos_pairs",
        "n_total_pairs",
        "shuffled",
        "comparison",
    ]
)


# In[15]:


for i in df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():

    # manually get treatment
    tmp = df[
        df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.contains(i)
    ].reset_index(drop=True)
    # get the label
    label = tmp["Metadata_labels"].unique().tolist()[0]
    # add all labels to the df except for the LPS treatment label
    tmp1 = df[~df["Metadata_labels"].str.contains(label)].reset_index(drop=True)
    # concat tmp and tmp1
    tmp1 = pd.concat([tmp, tmp1]).reset_index(drop=True)

    # drop rows that contain the label
    _pos_samby_permutations = [x for x in pos_samby_permutations if label in x]

    for j in _pos_samby_permutations:
        tmp1 = tmp1[tmp1["Metadata_labels"].isin(j)].reset_index(drop=True)

        # spliting metadata and raw feature values
        logging.info("splitting data set into metadata and raw feature values")
        df_meta, df_feats = utils.split_data(tmp1)
        df_feats = np.array(df_feats)
        try:
            # execute pipeline on negative control with trianing dataset with cp features

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

            result["shuffled"] = "non-shuffled"
            comparison = i
            comparison = comparison + "_" + "_vs_".join(j)

            result["comparison"] = comparison

        except ZeroDivisionError as e:
            logging.warning(f"{e} captured on phenotye:. Skipping")

        # concatenating all datasets
        results_df = pd.concat([results_df, result], ignore_index=True)
results_df.to_csv(regular_feat_map_path, index=False)
results_df.head()


# ### mAP analysis for shuffled data (Phenotype)

# In[16]:


results_df = pd.DataFrame(
    columns=[
        "Metadata_Well",
        "Metadata_labels",
        "average_precision",
        "p_value",
        "n_pos_pairs",
        "n_total_pairs",
        "shuffled",
        "comparison",
    ]
)


# In[17]:


for i in df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():

    # manually get treatment
    tmp = df[
        df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.contains(i)
    ].reset_index(drop=True)
    # get the label
    label = tmp["Metadata_labels"].unique().tolist()[0]
    # add all labels to the df except for the LPS treatment label
    tmp1 = df[~df["Metadata_labels"].str.contains(label)].reset_index(drop=True)
    # concat tmp and tmp1
    tmp1 = pd.concat([tmp, tmp1]).reset_index(drop=True)

    # drop rows that contain the label
    _pos_samby_permutations = [x for x in pos_samby_permutations if label in x]

    for j in _pos_samby_permutations:
        tmp1 = tmp1[tmp1["Metadata_labels"].isin(j)].reset_index(drop=True)

        # spliting metadata and raw feature values
        logging.info("splitting data set into metadata and raw feature values")
        df_meta, df_feats = utils.split_data(tmp1)
        df_feats = np.array(df_feats)
        seed = np.random.randint(0, 100)

        # shuffling the features, this will overwrite the generated feature space from above with the shuffled one
        df_feats = shuffle_features(feature_vals=df_feats, seed=seed)

        try:
            # execute pipeline on negative control with trianing dataset with cp features

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

            result["shuffled"] = "shuffled"
            comparison = i
            comparison = comparison + "_" + "_vs_".join(j)

            result["comparison"] = comparison

        except ZeroDivisionError as e:
            logging.warning(f"{e} captured on phenotye:. Skipping")

        # concatenating all datasets
        results_df = pd.concat([results_df, result], ignore_index=True)

# saving to csv
results_df.to_csv(shuffled_feat_space_map_path, index=False)
results_df
