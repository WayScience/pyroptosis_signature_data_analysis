#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pathlib

import pandas as pd
import pyarrow.parquet as pq

# In[12]:


cell_type = "PBMC"


# In[13]:


path_to_indexes = pathlib.Path(
    f"../indexes/{cell_type}/MultiClass_MLP_data_split_indexes.tsv"
).resolve()
file_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm.parquet"
).resolve(strict=True)


# In[14]:


indexes_df = pd.read_csv(path_to_indexes, sep="\t")
print(indexes_df.shape)
indexes_df.head()


# In[15]:


# replace the index with the labeled_data index
indexes_df = indexes_df.set_index("labeled_data_index")
# sort the index
indexes_df = indexes_df.sort_index()
indexes_df.head()


# In[27]:


df = pd.read_parquet(
    file_path,
    columns=[
        # "Metadata_Well",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    ],
)
print(df.shape)
df.head()


# In[28]:


# add data split label to the dataframe via index
df["data_split"] = indexes_df["label"]
df.head()


# In[30]:


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


# In[ ]:
