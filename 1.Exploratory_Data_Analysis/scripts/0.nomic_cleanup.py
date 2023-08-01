#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# set import data paths
nomic_df_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)

nomic_df_filtered_out_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_cleanup4correlation.csv"
)


# In[4]:


# read in data
nomic_df_raw = pd.read_csv(nomic_df_path)


# In[5]:


# get the dimensions of the df
print(f"{nomic_df_raw.shape}_before_filtering")
# remove column if colname has pgml in it as we are using the Normalised values ('NSU') columns
nomic_df = nomic_df_raw.loc[
    :, ~nomic_df_raw.columns.str.contains("pgml", case=False, na=False)
]
print(f"{nomic_df.shape}_after_filtering")
# if column does not contain [NSU] then prefix with Metadata_
for col in nomic_df.columns:
    if not any(x in col for x in ["NSU"]):
        nomic_df = nomic_df.rename(columns={col: "Metadata_" + col})


# In[6]:


# show all columns
pd.set_option("display.max_columns", None)
nomic_df.head(3)


def add_trailing_zeros(x):
    return "{:.3f}".format(x)


# Apply the function to the 'column_name' column
nomic_df["Metadata_inducer1_concentration_value"] = nomic_df[
    "Metadata_inducer1_concentration_value"
].apply(add_trailing_zeros)
nomic_df["Metadata_inducer2_concentration_value"] = nomic_df[
    "Metadata_inducer2_concentration_value"
].apply(add_trailing_zeros)
nomic_df["Metadata_inhibitor_concentration_value"] = nomic_df[
    "Metadata_inhibitor_concentration_value"
].apply(add_trailing_zeros)


# In[7]:


## Clean up df
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


# In[8]:


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


# In[9]:


scaler = MinMaxScaler()
nomic_df = pd.DataFrame(scaler.fit_transform(nomic_df), columns=nomic_df.columns)


# In[10]:


# summary statistics of df to check min-max normalization
nomic_df.describe()


# In[11]:


# add position_x back to df
nomic_df.loc[:, "Metadata_position_x"] = nomic_df_raw["position_x"]


# In[12]:


nomic_df = nomic_df.assign(
    oneb_Metadata_Treatment_Dose_Inhibitor_Dose=nomic_cleaned[
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    ]
)


# In[13]:


nomic_df.head(3)


# In[14]:


nomic_df.to_csv(nomic_df_filtered_out_path, index=False)
