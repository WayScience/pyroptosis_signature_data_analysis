#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import pathlib

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

# In[2]:


cell_type = "PBMC"


# In[3]:


final_df = pd.DataFrame(
    columns=[
        "group1",
        "group2",
        "meandiff",
        "lower",
        "upper",
        "p-adj",
        "reject",
        "features",
    ]
)


# In[4]:


# directory to get files from
data_dir = pathlib.Path(f"../results/{cell_type}/").resolve(strict=True)

# save directory
output_file_path = pathlib.Path(f"../results/{cell_type}_combined.parquet")

# get list of files in directory
file_list = [x for x in data_dir.iterdir() if x.is_file()]

list_of_dfs = []
# loop through files
for file in file_list:
    tmp_df = pd.read_parquet(file)
    list_of_dfs.append(tmp_df)
final_df = pd.concat(list_of_dfs, ignore_index=True)

final_df.shape


# In[5]:


# correct for multiple testing
final_df["reject"], final_df["p-adj_fdr_bh"], _, _ = multipletests(
    final_df["p-adj"], alpha=0.001, method="fdr_bh"
)
final_df.head()


# In[6]:


# save the final_df
final_df.to_parquet(output_file_path)
