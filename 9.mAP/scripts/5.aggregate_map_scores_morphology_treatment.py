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
class_files = [file for file in all_files if "treatment" in file.stem]
mAPs = []
for file in class_files:
    df = pd.read_csv(file)
    df["file"] = file.stem
    mAPs.append(df)
# single-cell mAP scores
mAPs = pd.concat(mAPs)
mAPs.head()


# In[4]:


# grabbing all cp features (regular, feature shuffled and labeled shuffled)
reg_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "non-shuffled"]
shuffled_feat_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "features_shuffled"]


# In[6]:


# grouping dataframe based on phenotype levels, feature and feature types
df_group = mAPs.groupby(by=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled"])

# calculating sampling error
sampling_error_df = []
for name, df in df_group:
    pheno, shuffled_type = name

    # caclulating sampling error
    avg_percision = df["average_precision"].values
    sampling_error = np.std(avg_percision) / np.sqrt(len(avg_percision))

    sampling_error_df.append([pheno, shuffled_type, sampling_error])
cols = ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled", "sampling_error"]
sampling_error_df = pd.DataFrame(sampling_error_df, columns=cols)

sampling_error_df.head()


# In[8]:


# Generating aggregate scores with a threshold p-value of 0.05
mAP_dfs = []
for name, df in tuple(
    mAPs.groupby(by=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled"])
):
    agg_df = aggregate(
        df, sameby=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"], threshold=0.05
    )
    agg_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = name[0]
    agg_df["shuffled"] = name[1]
    mAP_dfs.append(agg_df)

mAP_dfs = pd.concat(mAP_dfs)
mAP_dfs.to_csv(agg_sc_ap_scores_dir / "mAP_scores_treatment.csv", index=False)
mAP_dfs.head()
