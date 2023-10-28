#!/usr/bin/env python
# coding: utf-8

# <span style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">An Exception was encountered at '<a href="#papermill-error-cell">In [8]</a>'.</span>

# # This notebook looks into the cell heterogeneity in the control treatments

# In[1]:


import itertools
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
import toml
import umap
from sklearn.cluster import KMeans

# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# read in toml file

# set up the path
toml_path = pathlib.Path("./utils/params.toml")
# read in the toml file
params = toml.load(toml_path)
list_of_treatments = params["list_of_treatments"]["treatments"]
print(len(list_of_treatments))
print(list_of_treatments)


# In[4]:


# Set path to parquet file
path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_norm.parquet")
# Read in parquet file
df = pq.read_table(path).to_pandas()
df


# In[5]:


df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# subset the df for the control
df = df[df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == "DMSO_0.100_DMSO_0.025"]
df


# In[6]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# In[7]:


# set umap parameters
umap_params = umap.UMAP(
    n_components=2,
    spread=1.1,
    random_state=0,
    n_neighbors=6,
    min_dist=0.8,
    metric="cosine",
)


# <span id="papermill-error-cell" style="color:red; font-family:Helvetica Neue, Helvetica, Arial, sans-serif; font-size:2em;">Execution using papermill encountered an exception here and stopped:</span>

# In[8]:


# fit and transform data for umap
proj_2d = umap_params.fit_transform(df_values)

# add umap coordinates to dataframe of metadata and raw data
df_values["umap_1"] = proj_2d[:, 0]
df_values["umap_2"] = proj_2d[:, 1]


# In[ ]:


df_values["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df_descriptive[
    "fourb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
df_values["Metadata_Well"] = df_descriptive["Metadata_Well"]


# In[ ]:


# Figure Showing UMAP of Clusters vs Treatment
sns.scatterplot(
    data=df_values,
    x="umap_1",
    y="umap_2",
    hue="Metadata_Well",
    legend="full",
    # make points smaller,
    s=10,
    alpha=0.3,
)
# add contour lines to the plot

plt.title(f"Visualized {cell_type} on umap")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

# if path does not exist create it
plt.savefig(
    f"Figures/umap_plate2/cell_heterogeneity_{cell_type}_umap.png", bbox_inches="tight"
)
