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

import numpy as np
import pandas as pd

# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# set import data paths
nomic_df_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_cleanup4correlation.csv"
)

# read in data
nomic_df_raw = pd.read_csv(nomic_df_path)


# # Correlation of Wells

# In[4]:


# set index to Metadata_Well
nomic_df = nomic_df_raw.set_index("Metadata_position_x")


# In[5]:


# drop oneb_Metadata_Treatment_Dose_Inhibitor_Dose
nomic_df = nomic_df.drop(columns=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"])


# In[6]:


nomic_df = (
    nomic_df.T.corr()
    .reset_index()
    .rename(columns={"Metadata_position_x": "Wells"})
    .melt(
        id_vars="Wells",
        var_name="Metadata_position_x",
        value_name="correlation",
    )
    .rename(columns={"Metadata_position_x": "Metadata_Well"})
)


# In[7]:


save_path = pathlib.Path(f"./results/correlation/{cell_type}/nomic/")
save_path.mkdir(parents=True, exist_ok=True)
nomic_df.to_csv(f"{save_path}/wells.csv")


# # All Treatment correlation

# In[8]:


nomic_df = nomic_df_raw.copy()
nomic_df.reset_index(inplace=True)
nomic_df.drop(["Metadata_position_x"], axis=1, inplace=True)
nomic_df.loc[:, "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = nomic_df_raw[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
# groupby Metadata_Treatment_Dose_Inhibitor_Dose
nomic_df = nomic_df.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").mean()


# In[9]:


nomic_df = (
    nomic_df.T.corr()
    .reset_index()
    .rename(columns={"oneb_Metadata_Treatment_Dose_Inhibitor_Dose": "Treatments"})
    .melt(
        id_vars="Treatments",
        var_name="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        value_name="correlation",
    )
    .rename(columns={"Metadata_position_x": "Metadata_Well"})
)


# In[10]:


save_path = pathlib.Path(f"./results/correlation/{cell_type}/nomic/")
save_path.mkdir(parents=True, exist_ok=True)
nomic_df.to_csv(f"{save_path}/treatments.csv")
