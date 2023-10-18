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
df = pq.read_table(path).to_pandas()


# In[5]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])
# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# In[6]:


# set umap parameters
umap_params = umap.UMAP(
    n_components=2,
    spread=1.1,
    min_dist=0.8,
    init="random",
    metric="cosine",
    random_state=0,
    n_jobs=-1,
)


# In[7]:


# fit and transform data for umap
proj_2d = umap_params.fit_transform(df_values)
# add umap coordinates to dataframe of metadata and raw data
df_values["umap_1"] = proj_2d[:, 0]
df_values["umap_2"] = proj_2d[:, 1]


# In[8]:


df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = df_descriptive[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]


# In[9]:


# randomize the rows of the dataframe to plot the order of the data evenly
df_values = df_values.sample(frac=1, random_state=0)

df_values_path = pathlib.Path(
    f"./results/{celltype}_umap_values_morphology_all_cells.parquet"
)
# if path does not exist create it
df_values_path.parent.mkdir(parents=True, exist_ok=True)
# save the dataframe as a parquet file
df_values.to_parquet(df_values_path)


# In[10]:


# Figure Showing UMAP of Clusters vs Treatment
# figure size set
plt.figure(figsize=(10, 10))
sns.scatterplot(
    data=df_values,
    x="umap_1",
    y="umap_2",
    hue="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    legend="full",
    alpha=0.3,
)
plt.title("Visualized on umap")
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
# if path does not exist create it
plt.savefig(f"Figures/umap_plate2/{celltype}_umap.png", bbox_inches="tight")
