#!/usr/bin/env python
# coding: utf-8

# # Correlation of Nomic Features
# Each nELISA features in used for the correlation of the wells.
# The correlation of the aggregated features is calculated per:
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
nomic_df_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)

# read in data
nomic_df_raw = pd.read_csv(nomic_df_path)


# In[4]:


# remove column if colname has pgml in it
nomic_df = nomic_df_raw.loc[
    :, ~nomic_df_raw.columns.str.contains("pgml", case=False, na=False)
]
# if column does not contain [NSU] then prefix with Metadata_
for col in nomic_df.columns:
    if not any(x in col for x in ["NSU"]):
        nomic_df = nomic_df.rename(columns={col: "Metadata_" + col})


# In[5]:


## Clean up df
# remove uM in each row of the Metadata_inducer1_concentration column if it is present
# if "inducer1_concentration_value" in nomic_df.columns:
#     nomic_df["inducer1_concentration_value"] = nomic_df[
#         "inducer1_concentration_value"
#     ].str.replace("µM", "")
# replace nan values with 0
nomic_df["Metadata_inducer1_concentration_value"] = nomic_df[
    "Metadata_inducer1_concentration_value"
].fillna(0)
nomic_df["Metadata_inducer2_concentration_value"] = nomic_df[
    "Metadata_inducer2_concentration_value"
].fillna(0)
nomic_df["Metadata_inhibitor_concentration_value"] = nomic_df[
    "Metadata_inhibitor_concentration_value"
].fillna(0)
# treatment column merge
conditions = [
    (nomic_df["Metadata_inducer2"].isnull()),
    nomic_df["Metadata_inducer2"].notnull(),
]
results = [
    (nomic_df["Metadata_inducer1"]).astype(str),
    (nomic_df["Metadata_inducer1"] + "_" + nomic_df["Metadata_inducer2"]).astype(str),
]
nomic_df["Metadata_Treatment"] = np.select(condlist=conditions, choicelist=results)

# dose column merge
conditions = [
    (nomic_df["Metadata_inducer2"].isnull()),
    nomic_df["Metadata_inducer2"].notnull(),
]

results = [
    (nomic_df["Metadata_inducer1_concentration_value"].astype(str)).astype(str),
    (
        nomic_df["Metadata_inducer1_concentration_value"].astype(str)
        + "_"
        + nomic_df["Metadata_inducer2_concentration_value"].astype(str)
    ).astype(str),
]
nomic_df["Metadata_Dose"] = np.select(condlist=conditions, choicelist=results)
nomic_df["Metadata_inducer1_concentration_value"] = pd.to_numeric(
    nomic_df["Metadata_inducer1_concentration_value"]
)
# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
nomic_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    nomic_df["Metadata_Treatment"]
    + "_"
    + nomic_df["Metadata_Dose"].astype(str)
    + "_"
    + nomic_df["Metadata_inhibitor"].astype(str)
    + "_"
    + nomic_df["Metadata_inhibitor_concentration_value"].astype(str)
).astype(str)


# In[6]:


nomic_cleaned = nomic_df.copy()
# drop first column of metadata
nomic_df.columns[3:25]
nomic_df = nomic_df.drop(nomic_df.columns[3:25], axis=1)
nomic_df = nomic_df.drop(nomic_df.columns[0:2], axis=1)
nomic_df.drop(nomic_df.columns[0], axis=1, inplace=True)
# drop Metadata_Dose column
nomic_df = nomic_df.drop(["Metadata_Dose"], axis=1)
nomic_df = nomic_df.drop(["Metadata_Treatment"], axis=1)
nomic_df = nomic_df.drop(["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"], axis=1)


# # Normalization of Values

# In[7]:


# min-max normalization of nomic data from scipy
from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler()
nomic_df = pd.DataFrame(scaler.fit_transform(nomic_df), columns=nomic_df.columns)


# In[8]:


# summary statistics of df to check min-max normalization
nomic_df.describe()


# In[9]:


# add position_x back to df
nomic_df.loc[:, "Metadata_position_x"] = nomic_df_raw["position_x"]


# # Correlation of Wells

# In[10]:


# set index to Metadata_Well
nomic_df = nomic_df.set_index("Metadata_position_x")


# In[11]:


well_corr_df = nomic_df.T.corr()
save_path = pathlib.Path(f"./results/correlation/{cell_type}/nomic/")
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/wells_corr.csv")


# # All Treatment correlation

# In[12]:


nomic_df.reset_index(inplace=True)
nomic_df.drop(["Metadata_position_x"], axis=1, inplace=True)
nomic_df.loc[:, "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = nomic_cleaned[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
# groupby Metadata_Treatment_Dose_Inhibitor_Dose
nomic_df = nomic_df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").mean()


# In[13]:


well_corr_df = nomic_df.T.corr()
save_path = pathlib.Path(f"./results/correlation/{cell_type}/nomic/")
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/treatments_corr.csv")


# # Treatment correlation for selected treatments

# In[14]:


list_of_treatments = [
    "LPS_0.01_DMSO_0.025",
    "LPS_0.1_DMSO_0.025",
    "LPS_1.0_DMSO_0.025",
    "LPS_10.0_DMSO_0.025",
    "LPS_100.0_DMSO_0.025",
    "DMSO_0.1_DMSO_0.025",
    "Thapsigargin_1.0_DMSO_0.025",
    "Thapsigargin_10.0_DMSO_0.025",
]


# In[15]:


nomic_df = nomic_df.reset_index()
# subset the data to only include the treatments of interest from list_of_treatments
nomic_df = nomic_df[
    nomic_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(list_of_treatments)
]
# aggregate by treatment and dose
nomic_df = nomic_df.groupby(["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]).mean()


# In[16]:


well_corr_df = nomic_df.T.corr()
save_path = pathlib.Path(f"./results/correlation/{cell_type}/nomic/")
save_path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{save_path}/selected_treatments_corr.csv")
