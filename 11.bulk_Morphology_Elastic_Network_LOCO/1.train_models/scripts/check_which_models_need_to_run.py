#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas


# In[2]:


# set path to array
SHSY5Y_path = pathlib.Path("../models/SHSY5Y.txt").resolve()
PBMC_path = pathlib.Path("../models/PBMC.txt").resolve()
# read in the array
SHSY5Y_df = pandas.read_csv(SHSY5Y_path, sep="\t", header=None)
PBMC_df = pandas.read_csv(PBMC_path, sep="\t", header=None)
SHSY5Y_df["cell_type"] = "SHSY5Y"
PBMC_df["cell_type"] = "PBMC"
# combine the two arrays
df_both_cell_types = pandas.concat([SHSY5Y_df, PBMC_df])


# In[3]:


# split the columns into two on [NSU]
df = df_both_cell_types[0].str.split(" \[NSU]_", expand=True)
df.rename(columns={0: "Cytokine"}, inplace=True)
df1 = df[1].str.split("__", expand=True)
df1.rename(columns={0: "channels", 1: "file"}, inplace=True)
# add df1 to df
df = pandas.concat([df, df1], axis=1)
df.drop(columns=1, inplace=True)
df["cell_type"] = df_both_cell_types["cell_type"]
# drop rows that contain .toml
df = df[~df["file"].str.contains(".toml")]
df.reset_index(drop=True, inplace=True)
df.head()


# ### Channel possibilities

# In[4]:


channel_possibilities = [
    "All_channels_final",
    "All_channels_shuffled_baseline",
    "CorrDNA_CorrER_final",
    "CorrDNA_CorrER_shuffled_baseline",
    "CorrDNA_CorrGasdermin_CorrER_final",
    "CorrDNA_CorrGasdermin_CorrER_shuffled_baseline",
    "CorrDNA_CorrGasdermin_CorrMito_CorrER_final",
    "CorrDNA_CorrGasdermin_CorrMito_CorrER_shuffled_baseline",
    "CorrDNA_CorrGasdermin_CorrMito_final",
    "CorrDNA_CorrGasdermin_CorrMito_shuffled_baseline",
    "CorrDNA_CorrGasdermin_final",
    "CorrDNA_CorrGasdermin_shuffled_baseline",
    "CorrDNA_CorrMito_CorrER_final",
    "CorrDNA_CorrMito_CorrER_shuffled_baseline",
    "CorrDNA_CorrMito_final",
    "CorrDNA_CorrMito_shuffled_baseline",
    "CorrDNA_CorrPM_CorrER_final",
    "CorrDNA_CorrPM_CorrER_shuffled_baseline",
    "CorrDNA_CorrPM_CorrGasdermin_CorrER_final",
    "CorrDNA_CorrPM_CorrGasdermin_CorrER_shuffled_baseline",
    "CorrDNA_CorrPM_CorrGasdermin_CorrMito_final",
    "CorrDNA_CorrPM_CorrGasdermin_CorrMito_shuffled_baseline",
    "CorrDNA_CorrPM_CorrGasdermin_final",
    "CorrDNA_CorrPM_CorrGasdermin_shuffled_baseline",
    "CorrDNA_CorrPM_CorrMito_CorrER_final",
    "CorrDNA_CorrPM_CorrMito_CorrER_shuffled_baseline",
    "CorrDNA_CorrPM_CorrMito_final",
    "CorrDNA_CorrPM_CorrMito_shuffled_baseline",
    "CorrDNA_CorrPM_final",
    "CorrDNA_CorrPM_shuffled_baseline",
    "CorrDNA_final",
    "CorrDNA_shuffled_baseline",
    "CorrER_final",
    "CorrER_shuffled_baseline",
    "CorrGasdermin_CorrER_final",
    "CorrGasdermin_CorrER_shuffled_baseline",
    "CorrGasdermin_CorrMito_CorrER_final",
    "CorrGasdermin_CorrMito_CorrER_shuffled_baseline",
    "CorrGasdermin_CorrMito_final",
    "CorrGasdermin_CorrMito_shuffled_baseline",
    "CorrGasdermin_final",
    "CorrGasdermin_shuffled_baseline",
    "CorrMito_CorrER_final",
    "CorrMito_CorrER_shuffled_baseline",
    "CorrMito_final",
    "CorrMito_shuffled_baseline",
    "CorrPM_CorrER_final",
    "CorrPM_CorrER_shuffled_baseline",
    "CorrPM_CorrGasdermin_CorrER_final",
    "CorrPM_CorrGasdermin_CorrER_shuffled_baseline",
    "CorrPM_CorrGasdermin_CorrMito_CorrER_final",
    "CorrPM_CorrGasdermin_CorrMito_CorrER_shuffled_baseline",
    "CorrPM_CorrGasdermin_CorrMito_final",
    "CorrPM_CorrGasdermin_CorrMito_shuffled_baseline",
    "CorrPM_CorrGasdermin_final",
    "CorrPM_CorrGasdermin_shuffled_baseline",
    "CorrPM_CorrMito_CorrER_final",
    "CorrPM_CorrMito_CorrER_shuffled_baseline",
    "CorrPM_CorrMito_final",
    "CorrPM_CorrMito_shuffled_baseline",
    "CorrPM_final",
    "CorrPM_shuffled_baseline",
    "No_channels_final",
    "No_channels_shuffled_baseline",
]


# ### Find the missing channel combinations

# In[5]:


# get value counts per cytokine
counts = df["Cytokine"].value_counts()
dict_of_reruns = {
    "cytokine": [],
    "channel": [],
    "cell_type": [],
}


for cytokine in counts.index:
    if counts[cytokine] < 128:
        for channel in channel_possibilities:
            for cell_type in df["cell_type"].unique():
                if (
                    df[
                        (df["Cytokine"] == cytokine)
                        & (df["channels"] == channel)
                        & (df["cell_type"] == cell_type)
                    ].shape[0]
                    == 0
                ):
                    dict_of_reruns["cytokine"].append(cytokine)
                    dict_of_reruns["channel"].append(channel)
                    dict_of_reruns["cell_type"].append(cell_type)
print(len(dict_of_reruns["cytokine"]))


# In[6]:


dict_of_reruns_df = pandas.DataFrame(dict_of_reruns)
# make a new column for shuffle if final is in the channel
# if not then it is shuffled
dict_of_reruns_df["shuffle"] = dict_of_reruns_df["channel"].str.contains(
    "shuffled_baseline"
)
# replace values in the shuffle column
# remove _final and _shuffled_baseline from the channels column
dict_of_reruns_df["channel"] = (
    dict_of_reruns_df["channel"]
    .str.replace("_final", "")
    .str.replace("_shuffled_baseline", "")
)
# add [NSU] to the cytokine column
dict_of_reruns_df["cytokine"] = dict_of_reruns_df["cytokine"] + " [NSU]"
dict_of_reruns_df.head()


# In[7]:


# iterate over the rows in the rerun dataframe and write to a file
file_path = pathlib.Path("../manual_jobs.txt").resolve()
# remove file if it exists
if file_path.exists():
    file_path.unlink()
for i, row in dict_of_reruns_df.iterrows():
    with open(file_path, "a") as f:
        # write each item in the line in ''
        f.write(
            f"'{row['cell_type']}' '{row['shuffle']}' '{row['channel']}' '{row['cytokine']}'\n"
        )

