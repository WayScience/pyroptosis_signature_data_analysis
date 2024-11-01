#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gc
import pathlib

import pandas as pd
import tqdm

# In[2]:


# import jump data # absolute paths
jump_data_path = pathlib.Path("/home/lippincm/Desktop/18TB/normalized_sc_data").resolve(
    strict=True
)

# get the list of parquet files in the directory that are not aggregated
jump_data_files = list(jump_data_path.glob("*.parquet"))
jump_data_files = [x for x in jump_data_files if "agg" not in str(x)]


# In[3]:


# loop through the files and aggregate the data
for file in tqdm.tqdm(jump_data_files):
    # read file
    # read file
    jump_df = pd.read_parquet(file)
    # extract the file name and path
    file_name = file.stem
    file_path = file.parent
    # define save path
    agg_file_name = file_path / f"{file_name}_agg.parquet"

    # separate the data into the different types metadata, data
    metadata = jump_df[jump_df.columns[jump_df.columns.str.contains("Metadata")]]
    features = jump_df[jump_df.columns[~jump_df.columns.str.contains("Metadata")]]
    features = features.copy()
    features.loc[:, "Metadata_Well"] = metadata["Metadata_Well"]
    # aggregate the data
    agg_df = features.groupby("Metadata_Well").agg("mean")
    metadata = metadata.drop_duplicates(subset="Metadata_Well")
    # add the metadata back to the aggregated data
    agg_df = agg_df.merge(metadata, on="Metadata_Well", how="left")
    # save the aggregated data
    agg_df.to_parquet(agg_file_name)
    print(agg_df.shape)
    agg_df.head()
    del jump_df, agg_df, metadata, features
    gc.collect()
