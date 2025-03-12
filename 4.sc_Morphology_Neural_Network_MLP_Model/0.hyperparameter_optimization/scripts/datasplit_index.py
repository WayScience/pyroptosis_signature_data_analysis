#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
import pyarrow.parquet as pq

# In[2]:


cell_type = "PBMC"


# In[3]:


path_to_indexes = pathlib.Path(
    f"../indexes/{cell_type}/MultiClass_MLP_data_split_indexes.tsv"
).resolve()
file_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm.parquet"
).resolve(strict=True)


# In[4]:


indexes_df = pd.read_csv(path_to_indexes, sep="\t")
print(indexes_df.shape)
indexes_df.head()


# In[5]:


# replace the index with the labeled_data index
indexes_df = indexes_df.set_index("labeled_data_index")
# sort the index
indexes_df = indexes_df.sort_index()
indexes_df.head()


# In[6]:


df = pd.read_parquet(
    file_path,
    columns=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
)
print(df.shape)
df.head()


# In[7]:


# add data split label to the dataframe via index
df["data_split"] = indexes_df["label"]
df.head()


# In[8]:


# get counts of each data split and treatment
grouped_df = pd.DataFrame(
    df.groupby(
        ["data_split", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    ).value_counts()
)
# melt the dataframe to get datasplits as columns
grouped_df = grouped_df.unstack(level=0)
# sort by treatment
grouped_df = grouped_df.sort_index()
grouped_df
