#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd

# In[2]:


# Parameters
cell_type = "SHSY5Y"
aggregation = True
nomic = True


# In[ ]:


if aggregation == True & nomic == True:
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
    )
elif aggregation == False & nomic == True:
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_nomic.parquet"
    )
elif aggregation and not nomic:
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
    )
elif aggregation == False & nomic == False:
    pass
else:
    raise ValueError("Wrong parameters")


# In[10]:


path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_norm.parquet")

data_df = pd.read_parquet(path)

data_df.head()

if nomic == True:
    # import nomic data
    nomic_df_path = pathlib.Path(
        f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_clean.csv"
    )
    df_nomic = pd.read_csv(nomic_df_path)

    # drop columns that contain [pgML]
    df_nomic = df_nomic.drop(
        columns=[col for col in df_nomic.columns if "[pgML]" in col]
    )
    # drop first 25 columns
    df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
    df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])
else:
    df_nomic = None


# In[11]:


df_nomic


# In[12]:


# subset each column that contains metadata
metadata = data_df.filter(regex="Metadata")

# get all columns that are not metadata except for metadata_Well
data = data_df.drop(metadata.columns, axis=1)

# get the metadata_Well column
metadata_well = metadata[
    ["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
]

data_df = pd.merge(data, metadata_well, left_index=True, right_index=True)


# In[21]:


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


# In[24]:


if (aggregation == True) and (nomic == True):

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
    # set path to save the data
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
    )


elif (aggregation == True) and (nomic == False):
    # subset each column that contains metadata
    metadata = data.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["position_x", "oneb_Treatment_Dose_Inhibitor_Dose"],
    )
    # set path to save the data
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
    )
elif (aggregation == False) and (nomic == True):
    data_df = pd.merge(
        data_df,
        df_nomic,
        left_on=["Metadata_Well", "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        right_on=["position_x", "oneb_Treatment_Dose_Inhibitor_Dose"],
    )
    data_df = data_df.drop(columns=["position_x"])
    # set path to save the data
    aggregated_data_path = pathlib.Path(
        f"../data/{cell_type}_preprocessed_sc_norm_with_nomic.parquet"
    )
elif aggregation == False and nomic == False:
    pass
else:
    print("Error")


# In[ ]:


# set path to save the data

# save the data
data_df.to_parquet(aggregated_data_path)
