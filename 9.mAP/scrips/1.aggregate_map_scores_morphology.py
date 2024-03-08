#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
from copairs.map import aggregate

warnings.filterwarnings("ignore")


# In[2]:


# Directories
processed_data_dir = pathlib.Path("../data/processed/")
sc_ap_scores_dir = (processed_data_dir / "mAP_scores/morphology").resolve()
agg_sc_ap_scores_dir = (processed_data_dir / "aggregate_mAPs/morphology").resolve()
agg_sc_ap_scores_dir.mkdir(parents=True, exist_ok=True)


# ## Preparing the dataset
#

# In[3]:


all_files = list(sc_ap_scores_dir.glob("*.csv"))
# get the files that contain the string class
class_files = [file for file in all_files if "class" in file.stem]
mAPs = []
for file in class_files:
    df = pd.read_csv(file)
    df["file"] = file.stem
    mAPs.append(df)
# single-cell mAP scores
mAPs = pd.concat(mAPs)
mAPs.head()
mAPs["comparison"].unique()


# In[4]:


# grabbing all cp features (regular, feature shuffled and labeled shuffled)
reg_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "non-shuffled"]
shuffled_feat_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "features_shuffled"]


# In[5]:


# calculating sampling error
# grouping dataframe based on phenotype levels, feature and feature types
df_group = mAPs.groupby(by=["Metadata_labels", "shuffled", "comparison"])
df_group
sampling_error_df = []
for name, df in df_group:
    pheno, shuffled_type, comparison = name

    # caclulating sampling error
    avg_percision = df["average_precision"].values
    sampling_error = np.std(avg_percision) / np.sqrt(len(avg_percision))

    sampling_error_df.append([pheno, shuffled_type, sampling_error, comparison])
cols = ["Metadata_labels", "shuffled", "sampling_error", "comparison"]
sampling_error_df = pd.DataFrame(sampling_error_df, columns=cols)


sampling_error_df.head()


# In[6]:


# Generating aggregate scores with a threshold p-value of 0.05
mAP_dfs = []
for name, df in tuple(mAPs.groupby(by=["Metadata_labels", "shuffled", "comparison"])):
    agg_df = aggregate(df, sameby=["Metadata_labels"], threshold=0.05)
    agg_df["Metadata_labels"] = name[0]
    agg_df["shuffled"] = name[1]
    agg_df["comparison"] = name[2]
    mAP_dfs.append(agg_df)

mAP_dfs = pd.concat(mAP_dfs)
mAP_dfs.to_csv(agg_sc_ap_scores_dir / "mAP_scores_class.csv", index=False)
mAP_dfs.head()
