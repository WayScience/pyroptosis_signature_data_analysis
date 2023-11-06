#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import itertools
import pathlib

import numpy as np
import pandas as pd
import toml

# In[2]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", default="all")

args = argparser.parse_args()

cell_type = args.cell_type


# In[4]:


# Parameters
aggregation = True
nomic = True


# In[5]:


MODEL_TYPE = "regression"


# In[15]:


# toml file path
TOML_PATH = pathlib.Path("../splits.toml")
# read toml file via toml
data_splits_by_treatments = toml.load(TOML_PATH)

# define the 100% test set data treatments
test_100_percent = data_splits_by_treatments["splits"]["data_splits_100"]
test_75_percent = data_splits_by_treatments["splits"]["data_splits_75"]


# In[60]:


aggregate_and_nomic_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
).resolve(strict=True)
aggregate_path = pathlib.Path(
    f"../../../data/SHSY5Y_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)
data_df = pd.read_parquet(aggregate_and_nomic_path)

data_df.head()

morphology_df = pd.read_parquet(aggregate_path)


# In[67]:


# get the NSU columns
nsu_cols = [col for col in data_df.columns if "NSU" in col]
nomic_df = data_df[nsu_cols]
nomic_df.loc["Metadata_Well"] = data_df["Metadata_Well"]
nomic_df.loc["oneb_Treatment_Dose_Inhibitor_Dose"] = data_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
]


# In[43]:


# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
]

data_df = pd.merge(data, metadata_well, left_index=True, right_index=True)


# In[44]:


# drop morphology metadata
morphology_df = morphology_df.drop(
    morphology_df.filter(regex="Metadata").columns, axis=1
)
morphology_df.head()


# In[45]:


# define the list of the channels
channel_list = ["DNA", "Gasdermin", "ER", "Mito", "PM"]


# In[46]:


# combiantions of channels
channel_combinations = []
for i in range(1, len(channel_list) + 1):
    tmp_list = list(itertools.combinations(channel_list, i))
    channel_combinations += tmp_list


# In[47]:


# set up the LOO channel with recursion for dropping multiple channels
def channel_drop(df, channel):
    df = df.drop(df.filter(regex=channel).columns, axis=1)
    return df


# In[48]:


# dictionary for each df to go into
results_dict = {}


# In[49]:


# get all of the the channel combinations
for i in channel_combinations:
    # keep all channels and drop all channels
    if len(i) == 5:
        # keep all channels
        tmp_df = morphology_df
        # get the remaining channels for indexing purposes
        channel_list_index = [x for x in channel_list if x in i]
        channel_list_index = "_".join(channel_list_index)
        results_dict[channel_list_index] = tmp_df

        # drop all channels
        tmp = channel_drop(morphology_df, i[0])
        for j in range(1, len(i)):
            tmp = channel_drop(tmp, i[j])
        tmp_df = tmp
        # get the remaining channels for indexing purposes
        channel_list_index = "No Channels"
        results_dict[channel_list_index] = tmp_df
    # drop 4 channels
    elif len(i) == 4:
        tmp = channel_drop(morphology_df, i[0])
        for j in range(1, len(i)):
            tmp = channel_drop(tmp, i[j])
        tmp_df = tmp
        # get the remaining channels for indexing purposes
        channel_list_index = [x for x in channel_list if x not in i]
        channel_list_index = "_".join(channel_list_index)
        results_dict[channel_list_index] = tmp_df
    # drop 3 channels
    elif len(i) == 3:
        tmp = channel_drop(morphology_df, i[0])
        for j in range(1, len(i)):
            tmp = channel_drop(tmp, i[j])
        tmp_df = tmp
        # get the remaining channels for indexing purposes
        channel_list_index = [x for x in channel_list if x not in i]
        channel_list_index = "_".join(channel_list_index)
        results_dict[channel_list_index] = tmp_df
    # drop 2 channels
    elif len(i) == 2:
        tmp = channel_drop(morphology_df, i[0])
        for j in range(1, len(i)):
            tmp = channel_drop(tmp, i[j])
        tmp_df = tmp
        # get the remaining channels for indexing purposes
        channel_list_index = [x for x in channel_list if x not in i]
        channel_list_index = "_".join(channel_list_index)
        results_dict[channel_list_index] = tmp_df
    # drop 1 channel
    elif len(i) == 1:
        tmp = channel_drop(morphology_df, i[0])
        tmp_df = tmp
        # get the remaining channels for indexing purposes
        channel_list_index = [x for x in channel_list if x not in i]
        channel_list_index = "_".join(channel_list_index)
        results_dict[channel_list_index] = tmp_df
    else:
        print("channel length error")


# In[71]:


# set path to save
pathlib.Path(f"../indexes/{cell_type}/regression/channels").mkdir(
    parents=True, exist_ok=True
)

# loop through the dictionary and save each dataframe
for i in results_dict:
    print(i)
    print(results_dict[i].shape)
    # rename the dictionary keys
    # combine the metadata and morphology dataframes
    new_df = pd.merge(results_dict[i], metadata_well, left_index=True, right_index=True)
    # combine the cytokine dataframes
    new_df = pd.merge(new_df, nomic_df, left_index=True, right_index=True)
    # set file path
    file_path = pathlib.Path(f"../indexes/{cell_type}/regression/channels/{i}.parquet")
    # save the dataframe
    new_df.to_parquet(file_path)


# In[73]:


# get the list of the dictionary keys
index_list = list(results_dict.keys())
index_list_new = []
for i in index_list:
    index_list_new.append(i + ".parquet")
# write the list to a text file
# file path
file_write_path = pathlib.Path(f"../cytokine_list/channel_splits.txt")
with open(file_write_path, "w") as f:
    for i in index_list_new:
        f.write("%s\n" % i)
