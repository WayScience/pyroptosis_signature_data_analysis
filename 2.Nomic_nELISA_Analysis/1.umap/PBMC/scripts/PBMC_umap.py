#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import toml
import umap

# In[2]:


nELISA_plate_430420_PBMC_path = pathlib.Path(
    "../../Data/clean/Plate2/nELISA_plate_430420_PBMC_cleanup4correlation.csv"
)
manual_cluster_1_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_1.csv"
)

manual_cluster_2_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_2.csv"
)

treatment_clusters_path = pathlib.Path(
    "../../../1.Exploratory_Data_Analysis/utils/params.toml"
)


nELISA_plate_430420_PBMC = pd.read_csv(nELISA_plate_430420_PBMC_path)
manual_clusters_1 = pd.read_csv(manual_cluster_1_path)
manual_clusters_2 = pd.read_csv(manual_cluster_2_path)
treatments = toml.load(treatment_clusters_path)["list_of_treatments"]["treatments"]

nELISA_orgingal_plate = nELISA_plate_430420_PBMC.copy()


# In[3]:


# select data only columns and make floats
nELISA_data_values = nELISA_orgingal_plate.filter(like="NSU", axis=1).astype("float")
nELISA_data_values.head()


# In[4]:


# print mean and stdev of first data column before and after normalization to check normalization
print(f"NSU nELISA mean of Activin A: {nELISA_data_values['Activin A [NSU]'].mean()}")
print(f"NSU nELISA STDEV of Activin A: {nELISA_data_values['Activin A [NSU]'].std()}")
print(f"NSU nELISA min of Activin A: {nELISA_data_values['Activin A [NSU]'].min()}")
print(f"NSU nELISA max of Activin A: {nELISA_data_values['Activin A [NSU]'].max()}")


# In[5]:


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


# In[6]:


# fit and transform data for umap
proj_2d = umap_params.fit_transform(nELISA_data_values)

# add umap coordinates to dataframe of metadata and raw data
nELISA_orgingal_plate["umap_1"] = proj_2d[:, 0]
nELISA_orgingal_plate["umap_2"] = proj_2d[:, 1]


# In[7]:


# define output paths
nELISA_plate_430420_out_path = pathlib.Path(
    "./results/nELISA_plate_430420_umap_PBMC.csv"
)
# write to csv
nELISA_orgingal_plate.to_csv(nELISA_plate_430420_out_path, index=False)


# ### Selected Treatments

# In[8]:


# select treatments from the lsit of treatments from the df
nELISA_plate_430420_PBMC_treatments = nELISA_plate_430420_PBMC[
    nELISA_plate_430420_PBMC["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        treatments
    )
]
# select data only columns and make floats
nELISA_plate_430420_PBMC_treatments_values = nELISA_plate_430420_PBMC_treatments.filter(
    like="NSU", axis=1
)
nELISA_plate_430420_PBMC_treatments_values = (
    nELISA_plate_430420_PBMC_treatments_values.astype("float")
)
# fit and transform data for umap
proj_2d = umap_params.fit_transform(nELISA_plate_430420_PBMC_treatments_values)

# add umap coordinates to dataframe of metadata and raw data
nELISA_plate_430420_PBMC_treatments["umap_1"] = proj_2d[:, 0]
nELISA_plate_430420_PBMC_treatments["umap_2"] = proj_2d[:, 1]

# define output paths
nELISA_plate_430420_selected_treatments_out_path = pathlib.Path(
    "./results/nELISA_plate_430420_umap_PBMC_selected_treatments.csv"
)
# write to csv
nELISA_plate_430420_PBMC_treatments.to_csv(
    nELISA_plate_430420_selected_treatments_out_path, index=False
)


# In[ ]:
