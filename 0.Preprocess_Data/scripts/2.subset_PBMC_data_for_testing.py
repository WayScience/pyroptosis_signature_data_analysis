#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# In[3]:


import pathlib

import numpy as np
import pandas as pd
import papermill as pm
import pyarrow as pa
import pyarrow.parquet as pq

# In[4]:


# Parameters
celltype = "SHSY5Y"
CONTROL_NAME = "DMSO_0.100_DMSO_0.025"
TREATMENT_NAME = "LPS_100.000_DMSO_0.025"


# In[5]:


# Define inputs
feature_file = pathlib.Path("../data/PBMC_preprocessed_sc_norm.parquet")
df = pq.read_table(feature_file).to_pandas()
output_dir = pathlib.Path(
    f"../data/PBMC_subset_sc_norm_{CONTROL_NAME}_{TREATMENT_NAME}.parquet"
)


# In[6]:


# filter the oneb_Metadata_Treatment_Dose_Inhibitor_Dose column to only include the treatment and control via loc
df = df.loc[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        [TREATMENT_NAME, CONTROL_NAME]
    )
]


# In[7]:


df.to_parquet(output_dir)


# In[ ]:
