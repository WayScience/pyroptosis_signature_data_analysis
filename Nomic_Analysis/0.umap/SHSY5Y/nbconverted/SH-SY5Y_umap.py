#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly_express as px
import seaborn as sns
import umap

# In[2]:


nELISA_plate_430420_SH_SY5Y_path = pathlib.Path(
    "../../Data/clean/Plate2/nELISA_plate_430420_SH_SY5Y.csv"
)
manual_cluster_1_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_1.csv"
)

manual_cluster_2_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_2.csv"
)


nELISA_plate_430420_SH_SY5Y = pd.read_csv(nELISA_plate_430420_SH_SY5Y_path)
manual_clusters_1 = pd.read_csv(manual_cluster_1_path)
manual_clusters_2 = pd.read_csv(manual_cluster_2_path)

nELISA_orgingal_plate = nELISA_plate_430420_SH_SY5Y.copy()


# In[3]:


# select data only columns and make floats
nELISA_data_values = nELISA_orgingal_plate.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")
nELISA_data_values.head()


# In[4]:


# normalize data via max value in each column
max_values = nELISA_data_values.max()  # find max value in each column
nELISA_data_values_sensor_max_norm = nELISA_data_values.div(
    max_values
)  # divide each value in each column by max value in that column
nELISA_data_values_sensor_max_norm.head()


# In[5]:


# print mean and stdev of first data column before and after normalization to check normalization
print(f"NSU nELISA mean of Activin A: {nELISA_data_values['Activin A [NSU]'].mean()}")
print(f"NSU nELISA STDEV of Activin A: {nELISA_data_values['Activin A [NSU]'].std()}")

print(
    f"NSU sensor normalized nELISA mean of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].mean()}"
)
print(
    f"NSU sensor normalized nELISA STDEV of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].std()}"
)


# In[6]:


# display max columns pandas
pd.set_option("display.max_columns", None)
nELISA_orgingal_plate.head()


# In[7]:


# rename columns to remove special character "/"
nELISA_orgingal_plate.columns = nELISA_orgingal_plate.columns.str.replace("/", "_")

# set umap parameters
umap_params = umap.UMAP(
    n_neighbors=6,
    min_dist=0.8,
    n_components=2,
    metric="cosine",
    spread=1.1,
    init="random",
    random_state=0,
)

# fit and transform data for umap
proj_2d = umap_params.fit_transform(nELISA_data_values_sensor_max_norm)

# add umap coordinates to dataframe of metadata and raw data
nELISA_orgingal_plate["umap_1"] = proj_2d[:, 0]
nELISA_orgingal_plate["umap_2"] = proj_2d[:, 1]

# add manual clusters columns to dataframe
nELISA_plate_430420 = pd.merge(
    nELISA_orgingal_plate,
    manual_clusters_2,
    on=(
        "inducer1",
        "inducer1_concentration_value",
        "inhibitor",
        "inhibitor_concentration_value",
        "inducer2",
        "inducer2_concentration_value",
    ),
    how="inner",
)

nELISA_plate_430420["inducer1_plus_concentration"] = (
    nELISA_plate_430420["inducer1"]
    + "_"
    + nELISA_plate_430420["inducer1_concentration"]
)

nELISA_plate_430420["inducer1_plus_inhibitor"] = (
    nELISA_plate_430420["inducer1"] + "_" + nELISA_plate_430420["inhibitor"]
)

nELISA_plate_430420["inducer1_plus_concentration_plus_inhibitor"] = (
    nELISA_plate_430420["inducer1"]
    + "_"
    + nELISA_plate_430420["inducer1_concentration"]
    + "_"
    + nELISA_plate_430420["inhibitor"]
)

# define output paths
nELISA_plate_430420_out_path = pathlib.Path(
    "./results/nELISA_plate_430420_umap_SH-SY5Y.csv"
)
# write to csv
nELISA_plate_430420.to_csv(nELISA_plate_430420_out_path, index=False)
