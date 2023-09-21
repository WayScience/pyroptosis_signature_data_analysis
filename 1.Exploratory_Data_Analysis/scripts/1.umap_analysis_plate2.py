#!/usr/bin/env python
# coding: utf-8

# ## Plate 2

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

get_ipython().run_line_magic("matplotlib", "inline")
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

# In[2]:


# Parameters
celltype = "SHSY5Y"


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
path = pathlib.Path(f"../data/{celltype}_preprocessed_sc_norm.parquet")
# Read in parquet file
df1 = pq.read_table(path).to_pandas()
# subset data frame to 1000 samples too much data results in poor clustering

df = df1
# save memory by deleting df1
del df1


# In[5]:


# get rows that have values in column fourb_Metadata_Treatment_Dose_Inhibitor_Dose that match treatment list
df = df[df["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(list_of_treatments)]

df = (
    df.groupby("fourb_Metadata_Treatment_Dose_Inhibitor_Dose")
    .apply(lambda x: x.sample(n=100, random_state=0))
    .droplevel(0)
)


# In[6]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# In[8]:


# set umap parameters
umap_params = umap.UMAP(
    n_components=2,
    spread=1.1,
    init="random",
    random_state=0,
)


# In[9]:


# fit and transform data for umap
proj_2d = umap_params.fit_transform(df_values)

# add umap coordinates to dataframe of metadata and raw data
df_values["umap_1"] = proj_2d[:, 0]
df_values["umap_2"] = proj_2d[:, 1]


# In[ ]:


df_values["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df_descriptive[
    "fourb_Metadata_Treatment_Dose_Inhibitor_Dose"
]


# In[10]:


# Figure Showing UMAP of Clusters vs Treatment
sns.scatterplot(
    data=df_values,
    x="umap_1",
    y="umap_2",
    hue="fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
    legend="full",
    alpha=0.7,
)
plt.title("Visualized on umap")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)

# if path does not exist create it
plt.savefig(f"Figures/umap_plate2/{celltype}_umap.png", bbox_inches="tight")
