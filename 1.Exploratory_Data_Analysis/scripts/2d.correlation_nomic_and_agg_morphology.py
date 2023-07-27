#!/usr/bin/env python
# coding: utf-8

# # Correlation of single cell morhologies and nELISA Cytokine/Chemokine Panel
# Each well is median aggregated and normalized.
# The correlation of the aggregate morphology features and nELISA features is calculated per:
# * well
# * per treatment
# * per selected treatment

# In[2]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as gg
import scipy.stats as stats
import seaborn as sns
import statsmodels.api as sm

# In[3]:


cell_type = "PBMC"
data_type = "norm"


# In[4]:


# set import data paths
PBMC_df_path = pathlib.Path(f"../data/{cell_type}_sc_aggregate_{data_type}.parquet")
nomic_df_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)

# read in data
feature_df = pd.read_parquet(PBMC_df_path)
nomic_df = pd.read_csv(nomic_df_path)


# ## Clean up df

# In[5]:


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


# In[6]:


# remove column if colname has pgml in it
nomic_df = nomic_df.loc[:, ~nomic_df.columns.str.contains("pgml", case=False, na=False)]
# if column does not contain [NSU] then prefix with Metadata_
for col in nomic_df.columns:
    if not any(x in col for x in ["NSU"]):
        nomic_df = nomic_df.rename(columns={col: "Metadata_" + col})


# In[7]:


# drop first column of metadata
nomic_df.columns[3:25]
nomic_df = nomic_df.drop(nomic_df.columns[3:25], axis=1)
nomic_df = nomic_df.drop(nomic_df.columns[0:2], axis=1)


# In[8]:


# merge nomic and morphology data by metadata_well and position_x
merged_df = pd.merge(
    feature_df,
    nomic_df,
    left_on=["Metadata_Well"],
    right_on=["Metadata_position_x"],
)
merged_df = merged_df.drop(["Metadata_position_x"], axis=1)
merged_df


# ## Correlation of Aggregated Features per Well + Nomic Data

# In[9]:


# drop metadata columns if the column contains 'Metadata'
PBMC_df = merged_df
feature_df = PBMC_df.loc[:, ~PBMC_df.columns.str.contains("Metadata")]
metadata_df = PBMC_df.loc[:, PBMC_df.columns.str.contains("Metadata")]


# In[10]:


feature_df.loc[:, "Metadata_Well"] = metadata_df["Metadata_Well"]
feature_df.index = feature_df["Metadata_Well"]
feature_df = feature_df.drop(columns=["Metadata_Well"])


# In[11]:


feature_df


# In[12]:


well_corr_df = feature_df.T.corr()
pathlib.Path("./results/correlation/{cell_type}/nomic_and_{data_type}")
well_corr_df.to_csv("./results/well_corr_df_.csv")


# # Selected Treatment correlation

# In[13]:


feature_df.reset_index(inplace=True)
feature_df.loc[:, "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = metadata_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
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
# drop metadata_well
treatments_agg = treatments_agg.drop(columns=["Metadata_Well"])
# aggregate by treatment and dose
treatments_agg = treatments_agg.groupby(
    ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
).mean()


# In[16]:


well_corr_df = treatments_agg.T.corr()
path = pathlib.Path(f"./results/correlation/{cell_type}/nomic_and_{data_type}")
path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{path}/selected_treatments.csv")


# In[17]:


## Correlation of Well per nomic feature


# In[18]:


# drop metadata columns if the column contains 'Metadata'

feature_df = nomic_df.loc[:, ~nomic_df.columns.str.contains("Metadata")]
metadata_df = nomic_df.loc[:, nomic_df.columns.str.contains("Metadata")]


# In[19]:


# from sklearn.preprocessing import MinMaxScaler

# scaler = MinMaxScaler()
# scaler.fit(feature_df)
# feature_df = pd.DataFrame(scaler.transform(feature_df), columns=feature_df.columns)

# standardize the nomic data via standard scaling
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()
scaler.fit(feature_df)
feature_df = pd.DataFrame(scaler.transform(feature_df), columns=feature_df.columns)


feature_df.loc[:, "Metadata_position_x"] = metadata_df["Metadata_position_x"]
feature_df.index = feature_df["Metadata_position_x"]
feature_df = feature_df.drop(columns=["Metadata_position_x"])
# normalize the nomic data via min max scaling


# In[20]:


well_corr_df = feature_df.T.corr()
path = pathlib.Path(f"./results/correlation/{cell_type}")
path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{path}/nomic.csv")


# In[21]:


well_corr_df


# In[27]:


treatments_agg


# In[28]:


# drop columns that do not contain '[NUS]
nomic_df = treatments_agg.loc[:, treatments_agg.columns.str.contains("NSU")]


# In[30]:


well_corr_df = nomic_df.T.corr()
path = pathlib.Path(f"./results/correlation/{cell_type}")
path.mkdir(parents=True, exist_ok=True)
well_corr_df.to_csv(f"{path}/nomic_selected_features.csv")


# In[ ]:
