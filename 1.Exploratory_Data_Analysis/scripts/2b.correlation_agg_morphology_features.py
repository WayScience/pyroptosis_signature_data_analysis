#!/usr/bin/env python
# coding: utf-8

# # Correlation of Aggregated Features
# Each well is median aggregated and normalized.  
# The correlation of the aggregated features is calculated per:
# * well
# * per treatment
# * per selected treatment

# In[1]:


import pathlib

import numpy as np
import pandas as pd


# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# set import data paths
feature_df_path = pathlib.Path(f"../data/{cell_type}_sc_aggregate_norm.parquet")

# read in data
feature_df = pd.read_parquet(feature_df_path)


# ## Clean up df

# In[4]:


# remove uM in each row of the Metadata_inducer1_concentration column if it is present
if "Metadata_inducer1_concentration" in feature_df.columns:
    feature_df["Metadata_inducer1_concentration"] = feature_df[
        "Metadata_inducer1_concentration"
    ].str.replace("µM", "")
# replace nan values with 0
feature_df["Metadata_inducer1_concentration"] = feature_df[
    "Metadata_inducer1_concentration"
].fillna(0)
feature_df["Metadata_inducer2_concentration"] = feature_df[
    "Metadata_inducer2_concentration"
].fillna(0)
feature_df["Metadata_inhibitor_concentration"] = feature_df[
    "Metadata_inhibitor_concentration"
].fillna(0)
# treatment column merge
conditions = [
    (feature_df["Metadata_inducer2"].isnull()),
    feature_df["Metadata_inducer2"].notnull(),
]
results = [
    (feature_df["Metadata_inducer1"]).astype(str),
    (feature_df["Metadata_inducer1"] + "_" + feature_df["Metadata_inducer2"]).astype(
        str
    ),
]
feature_df["Metadata_Treatment"] = np.select(condlist=conditions, choicelist=results)

# dose column merge
conditions = [
    (feature_df["Metadata_inducer2"].isnull()),
    feature_df["Metadata_inducer2"].notnull(),
]

results = [
    (feature_df["Metadata_inducer1_concentration"].astype(str)).astype(str),
    (
        feature_df["Metadata_inducer1_concentration"].astype(str)
        + "_"
        + feature_df["Metadata_inducer2_concentration"].astype(str)
    ).astype(str),
]
feature_df["Metadata_Dose"] = np.select(condlist=conditions, choicelist=results)
feature_df["Metadata_inducer1_concentration"] = pd.to_numeric(
    feature_df["Metadata_inducer1_concentration"]
)
# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_Dose"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)


# ## Correlation of Aggregated Features per Well

# In[5]:


# drop metadata columns if the column contains 'Metadata'
cell_df = feature_df
feature_df = cell_df.loc[:, ~cell_df.columns.str.contains("Metadata")]
metadata_df = cell_df.loc[:, cell_df.columns.str.contains("Metadata")]


# In[6]:


feature_df.loc[:, "Metadata_Well"] = metadata_df["Metadata_Well"]
feature_df.index = feature_df["Metadata_Well"]
feature_df = feature_df.drop(columns=["Metadata_Well"])


# In[7]:


feature_df = (
    feature_df.T.corr()
    .reset_index()
    .rename(columns={"Metadata_Well": "Wells"})
    .melt(
        id_vars="Wells",
        var_name="Metadata_Well",
        value_name="correlation",
    )
)


# In[8]:


save_path = pathlib.Path(f"./results/correlation/{cell_type}/aggregated_morphology")
save_path.mkdir(parents=True, exist_ok=True)
feature_df.to_csv(f"{save_path}/wells.csv")


# # Correlation of Aggregated Features per Treatment

# In[9]:


feature_df = cell_df.loc[:, ~cell_df.columns.str.contains("Metadata")]
metadata_df = cell_df.loc[:, cell_df.columns.str.contains("Metadata")]
feature_df.loc[:, "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = metadata_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
# groupby oneb_Metadata_Treatment_Dose_Inhibitor_Dose column
feature_df = feature_df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").median()


# In[10]:


feature_df = (
    feature_df.T.corr()
    .reset_index()
    .rename(columns={"oneb_Metadata_Treatment_Dose_Inhibitor_Dose": "Treatments"})
    .melt(
        id_vars="Treatments",
        var_name="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        value_name="correlation",
    )
)


# In[11]:


save_path = pathlib.Path(f"./results/correlation/{cell_type}/aggregated_morphology")
save_path.mkdir(parents=True, exist_ok=True)
feature_df.to_csv(f"{save_path}/treatments.csv")

