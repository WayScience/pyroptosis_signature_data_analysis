#!/usr/bin/env python
# coding: utf-8

# # Correlation of single cell morhologies and nELISA Cytokine/Chemokine Panel
# Each well is median aggregated and normalized.
# The correlation of the aggregate morphology features and nELISA features is calculated per:
# * well
# * per treatment
# * per selected treatment

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as gg
import scipy.stats as stats
import seaborn as sns
import statsmodels.api as sm

# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# set import data paths
PBMC_df_path = pathlib.Path(f"../data/{cell_type}_sc_aggregate_norm.parquet")
nomic_df_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)

# read in data
feature_df = pd.read_parquet(PBMC_df_path)
nomic_df = pd.read_csv(nomic_df_path)


# ## Clean up morphology df

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


# ## Clean up Nomic df

# In[5]:


# remove column if colname has pgml in it
nomic_df = nomic_df.loc[:, ~nomic_df.columns.str.contains("pgml", case=False, na=False)]
# if column does not contain [NSU] then prefix with Metadata_
for col in nomic_df.columns:
    if not any(x in col for x in ["NSU"]):
        nomic_df = nomic_df.rename(columns={col: "Metadata_" + col})


# In[6]:


# drop first column of metadata
nomic_df.columns[3:25]
nomic_df = nomic_df.drop(nomic_df.columns[3:25], axis=1)
nomic_df = nomic_df.drop(nomic_df.columns[0:2], axis=1)


# ### Merge Nomic and Aggregated Morphology dfs

# In[7]:


# merge nomic and morphology data by metadata_well and position_x
merged_df = pd.merge(
    feature_df,
    nomic_df,
    left_on=["Metadata_Well"],
    right_on=["Metadata_position_x"],
)
merged_df = merged_df.drop(["Metadata_position_x"], axis=1)


# ## Correlation of Aggregated Features per Well + Nomic Data

# In[8]:


# drop metadata columns if the column contains 'Metadata'
cell_df = merged_df
feature_df = cell_df.loc[:, ~cell_df.columns.str.contains("Metadata")]
metadata_df = cell_df.loc[:, cell_df.columns.str.contains("Metadata")]


# In[9]:


feature_df.loc[:, "Metadata_Well"] = metadata_df["Metadata_Well"]
feature_df.index = feature_df["Metadata_Well"]
feature_df = feature_df.drop(columns=["Metadata_Well"])


# In[10]:


well_corr_df = feature_df.T.corr()
save_path = pathlib.Path(
    f"./results/correlation/{cell_type}/aggregated_morphology_and_nomic"
)
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/wells_corr.csv")


# # All Treatment Correlation

# In[11]:


feature_df.reset_index(inplace=True)
feature_df.drop(columns=["Metadata_Well"], inplace=True)
feature_df.loc[:, "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = metadata_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
# group by oneb_Metadata_Treatment_Dose_Inhibitor_Dose
feature_df = feature_df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").mean()


# In[12]:


well_corr_df = feature_df.T.corr()
save_path = pathlib.Path(
    f"./results/correlation/{cell_type}/aggregated_morphology_and_nomic"
)
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/treatments_corr.csv")


# # Selected Treatment correlation

# In[13]:


feature_df.reset_index(inplace=True)
treatments_agg = feature_df


# In[14]:


list_of_treatments = [
    "LPS_0.010_DMSO_0.025",
    "LPS_0.100_DMSO_0.025",
    "LPS_1.000_DMSO_0.025",
    "LPS_10.000_DMSO_0.025",
    "LPS_100.000_DMSO_0.025",
    "DMSO_0.100_DMSO_0.025",
    "Thapsigargin_1.000_DMSO_0.025",
    "Thapsigargin_10.000_DMSO_0.025",
]


# In[15]:


# subset the data to only include the treatments of interest from list_of_treatments
treatments_agg = treatments_agg[
    treatments_agg["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        list_of_treatments
    )
]
# set the index to the treatment column
treatments_agg.set_index("oneb_Metadata_Treatment_Dose_Inhibitor_Dose", inplace=True)


# In[16]:


well_corr_df = treatments_agg.T.corr()
save_path = pathlib.Path(
    f"./results/correlation/{cell_type}/aggregated_morphology_and_nomic"
)
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/selected_treatments_corr.csv")
