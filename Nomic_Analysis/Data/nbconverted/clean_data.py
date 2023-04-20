#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly
import plotly_express as px
import seaborn as sns
import umap

# ## Clean Data from Nomic

# In[2]:


nELISA_data_input_path = pathlib.Path("./raw/nELISA_Data_UCA1_2023.04.11.csv")
MetaData_input_path = pathlib.Path("./raw/Metadata_UCA1_2023.04.11.csv")

# import data
nELISA_data_all = pd.read_csv(nELISA_data_input_path)
MetaData = pd.read_csv(MetaData_input_path)


# In[3]:


# Change the 'A1' cell format to 'A01' format
position = []
for i in MetaData["position"].astype(str):
    if len(i) < 3:
        position.append(i[:1] + f"{0}" + i[1:])
    else:
        position.append(i)
MetaData["position"] = position
MetaData.head()


# In[4]:


# Change column names
nELISA_data_all = nELISA_data_all.rename({"user_well_loc": "position"}, axis=1)
nELISA_data_all.head()
# nELISA_data

MetaData["plate_position"] = (
    MetaData["plate_barcode"].astype(str) + "_" + MetaData["position"].astype(str)
)
MetaData["plate_position"]

nELISA_data_all["plate_position"] = (
    nELISA_data_all["user_plate_id"].astype(str)
    + "_"
    + nELISA_data_all["position"].astype(str)
)

# Fix plate naming
nELISA_data_all.replace(regex=[" and "], value="_", inplace=True)


# In[5]:


# Seperate df out by plate
MetaData_plate_430418_430419 = MetaData.loc[
    MetaData["plate_barcode"] == "430418_430419"
]

MetaData_plate_430420 = MetaData.loc[MetaData["plate_barcode"] == "430420"]


# In[6]:


# seperate out by plate
nELISA_data_all_plate_430418_430419 = nELISA_data_all.loc[
    nELISA_data_all["user_plate_id"] == "430418_430419"
]

nELISA_data_all_plate_430420 = nELISA_data_all.loc[
    nELISA_data_all["user_plate_id"] == "430420"
]


# In[7]:


# Merge the two dataframes for plate 430420 via concat because pandas is being a pain
plate_430420 = pd.merge(
    MetaData_plate_430420,
    nELISA_data_all_plate_430420,
    on="plate_position",
    how="inner",
)

# remove empty wells and qc_fail
plate_430420 = plate_430420[
    ~plate_430420.nelisa_sample_comments.str.contains("empty_well", na=False)
]
plate_430420 = plate_430420[
    ~plate_430420.nelisa_sample_comments.str.contains("qc_fail", na=False)
]

plate_430420.shape


# In[8]:


# Merge the two dataframes for plate 430418_430419
plate_430418_430419 = pd.merge(
    MetaData_plate_430418_430419,
    nELISA_data_all_plate_430418_430419,
    on="plate_position",
    how="inner",
)

# remove empty wells and qc_fail
plate_430418_430419 = plate_430418_430419[
    ~plate_430418_430419.nelisa_sample_comments.str.contains("empty_well", na=False)
]
plate_430418_430419 = plate_430418_430419[
    ~plate_430418_430419.nelisa_sample_comments.str.contains("qc_fail", na=False)
]
plate_430418_430419.shape


# In[9]:


# seperate out by cell type
plate_430420_PBMC = plate_430420.loc[plate_430420["cell_type"] == "PBMC"]
plate_430420_PBMC.shape

plate_430420_SH_SY5Y = plate_430420.loc[plate_430420["cell_type"] == "SH-SY5Y"]
plate_430420_SH_SY5Y.shape


# In[10]:


# define output paths
nELISA_plate_430420_path = pathlib.Path("./clean/Plate2/nELISA_plate_430420.csv")
nELISA_plate_430420_PBMC_path = pathlib.Path(
    "./clean/Plate2/nELISA_plate_430420_PBMC.csv"
)
nELISA_plate_430420_SH_SY5Y_path = pathlib.Path(
    "./clean/Plate2/nELISA_plate_430420_SH_SY5Y.csv"
)
nELISA_430418_430419_path = pathlib.Path(
    "./clean/Plate1/nELISA_plate_430418_430419.csv"
)
# write to csv
plate_430420.to_csv(nELISA_plate_430420_path, index=False)
plate_430420.to_csv(nELISA_plate_430420_PBMC_path, index=False)
plate_430420.to_csv(nELISA_plate_430420_SH_SY5Y_path, index=False)
plate_430418_430419.to_csv(nELISA_430418_430419_path, index=False)
