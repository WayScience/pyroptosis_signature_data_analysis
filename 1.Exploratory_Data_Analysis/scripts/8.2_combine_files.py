#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import pathlib

import numpy as np
import pandas as pd

# In[2]:


cell_type = "SHSY5Y"


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

# loop through files
for file in file_list:
    tmp_df = pd.read_parquet(file)
    final_df = pd.concat([final_df, tmp_df])
    # del the tmp_df to save memory
    del tmp_df
    # run garbage collection to free up memory
    gc.collect()

final_df.shape


# In[5]:


# save the final_df
final_df.to_parquet(output_file_path)
