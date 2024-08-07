#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import umap

# In[2]:


cell_type = "PBMC"


# In[3]:


bulk_profile_file_path = pathlib.Path(
    f"../../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)
pathlib.Path("../results").mkdir(parents=True, exist_ok=True)
umap_output_file_path = pathlib.Path(
    f"../results/{cell_type}_umap_bulk_profile.parquet"
).resolve()


bulk_profile = pd.read_parquet(bulk_profile_file_path)
print(bulk_profile.shape)
bulk_profile.head()


# In[4]:


# get the Metadata columns
metadata_columns = bulk_profile.columns[bulk_profile.columns.str.contains("Metadata_")]
metadata_df = bulk_profile[metadata_columns]
feature_df = bulk_profile.drop(metadata_columns, axis=1)

# UMAP
# set umap parameters
umap_params = umap.UMAP(
    n_components=2,
    spread=1.1,
    min_dist=0.8,
    init="random",
    metric="cosine",
    random_state=0,
)

# fit umap
umap_output = umap_params.fit_transform(feature_df)
umap_output_df = pd.DataFrame(umap_output, columns=["UMAP1", "UMAP2"])
umap_output_df = pd.concat([metadata_df, umap_output_df], axis=1)

# save umap output
umap_output_df.to_parquet(umap_output_file_path)
