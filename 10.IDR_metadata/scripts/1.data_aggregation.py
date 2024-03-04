#!/usr/bin/env python
# coding: utf-8

# This noteboook aggregates the data from the previous notebooks and creates the final dataset for the analysis barring the data are aggregated in the analysis.

# In[1]:


import pathlib

import numpy as np
import pandas as pd


# In[2]:


# Parameters
cell_type = "PBMC"
aggregation = True
nomic = True


# In[3]:


if aggregation and nomic:
    aggregated_data_path = pathlib.Path(
        f"../../data/{cell_type}_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
    )
elif not aggregation and nomic:
    aggregated_data_path = pathlib.Path(
        f"../../data/{cell_type}_preprocess_sc_norm_no_fs_nomic.parquet"
    )
elif aggregation and not nomic:
    aggregated_data_path = pathlib.Path(
        f"../../data/{cell_type}_preprocess_sc_norm_no_fs_aggregated.parquet"
    )
elif not aggregation and not nomic:
    pass
else:
    raise ValueError("Wrong parameters")


# In[4]:


path = pathlib.Path(f"../../data/{cell_type}_preprocess_sc_norm_no_fs.parquet")

data_df = pd.read_parquet(path)

data_df.head()

if nomic:
    # import nomic data
    nomic_df_path = pathlib.Path(
        f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_clean.parquet"
    )
    df_nomic = pd.read_parquet(nomic_df_path)

    # drop columns that contain [pgML]
    df_nomic = df_nomic.drop(
        columns=[col for col in df_nomic.columns if "[pgML]" in col]
    )
    # drop first 25 columns (metadata that does not contain metadata in the title)
    df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
    df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])
elif not nomic:
    pass
else:
    raise ValueError("Nomic data not imported")


# In[5]:


# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
]

data_df = pd.merge(data, metadata_well, left_index=True, right_index=True)


# In[6]:


if nomic:
    df_nomic.drop(
        columns=[
            "Treatment",
            "Dose",
            "twob_Treatment_Dose_Inhibitor_Dose",
            "threeb_Treatment_Dose_Inhibitor_Dose",
            "fourb_Treatment_Dose_Inhibitor_Dose",
        ],
        inplace=True,
    )


# In[7]:


if aggregation and nomic:

    # subset each column that contains metadata
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df, metadata, left_on="Metadata_Well", right_on="Metadata_Well"
    )
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["position_x", "oneb_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["position_x"])
    # drop all metadata columns
    labeled_data = data_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    data_x = data_df.drop(metadata.columns, axis=1)


elif aggregation and not nomic:
    # get metadata columns
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    metadata
    data_df = pd.concat([data_df, metadata], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby(
        ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    ).mean()
    # # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # # get the metadata for each well
    # # set path to save the data
    # reset the index
    data_df = data_df.reset_index()

elif not aggregation and nomic:
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["position_x", "oneb_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["position_x"])
elif aggregation == False and nomic == False:
    pass
else:
    raise ValueError("Wrong parameters nomic and/or aggregation not defined")


# In[8]:


# save the data
data_df.to_parquet(aggregated_data_path)

