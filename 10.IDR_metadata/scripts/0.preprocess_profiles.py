#!/usr/bin/env python
# coding: utf-8

# This notebooks curates and preprocesses the data for IDR.
# Specifically, it preprocesses the normalized, non-feature selected data.

# In[1]:


# Parameters
cell_type = "PBMC"


# In[2]:


import pathlib

import numpy as np
import pandas as pd
import papermill as pm
import pyarrow as pa
import pyarrow.parquet as pq

# In[3]:


# Define inputs
feature_file = pathlib.Path(f"../../data/{cell_type}_sc_norm.parquet")
feature_df = pd.read_parquet(feature_file)


# In[4]:


len(feature_df["Metadata_Well"].unique())


# In[5]:


# replace all " " with "_" in all values of the dataframe
feature_df = feature_df.replace(to_replace=" ", value="_", regex=True)


# In[6]:


# remove uM in each row of the Metadata_inducer1_concentration column
feature_df["Metadata_inducer1_concentration"] = feature_df[
    "Metadata_inducer1_concentration"
].str.replace("µM", "")


# In[7]:


feature_df["Metadata_inducer1_concentration"].unique()


# In[8]:


# define output file path
feature_df_out_path = pathlib.Path(
    f"../../data/{cell_type}_preprocess_sc_norm_no_fs.parquet"
)


# In[9]:


print(feature_df.shape)
feature_df.head()


# In[10]:


# removing costes features as they behave with great variance across all data
feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df.head()


# In[11]:


# replacing '/' in treatment dosage column to avoid errors in file interpolation including such strings
feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)


# In[12]:


# replace nan values with 0

columns_to_fill = [
    "Metadata_inducer1_concentration",
    "Metadata_inducer2_concentration",
    "Metadata_inhibitor_concentration",
]
feature_df[columns_to_fill].fillna(0, inplace=True)


# In[13]:


# replace all None values with 0
feature_df["Metadata_inducer1_concentration"].fillna(0, inplace=True)


# In[14]:


# create a list of columns to be converted to float
col_list = [
    "Metadata_inducer1_concentration",
    "Metadata_inducer2_concentration",
    "Metadata_inhibitor_concentration",
]
# loop through the list and convert each column to float
for i in col_list:
    feature_df[i] = feature_df[i].apply(
        lambda x: f"{float(x):.3f}" if float(x) != 0 else float(x)
    )


# In[15]:


len(feature_df["Metadata_Well"].unique())


# #### Combine Inducer1 and Inducer2 into one column

# In[16]:


# treatment column merge
conditions = [
    (feature_df["Metadata_inducer2"].isnull()),
    feature_df["Metadata_inducer2"].notnull(),
]

results = [
    (feature_df["Metadata_inducer1"]).astype(str),
    (
        feature_df["Metadata_inducer1"]
        + "_"
        + feature_df["Metadata_inducer2"].astype(str)
    ),
]
feature_df["Metadata_Treatment"] = np.select(condlist=conditions, choicelist=results)


# dose column merge
results = [
    (
        feature_df["Metadata_inducer1_concentration"].astype(str)
        + "_"
        + feature_df["Metadata_inducer1_concentration_unit"].astype(str)
    ),
    (
        feature_df["Metadata_inducer1_concentration"].astype(str)
        + "_"
        + feature_df["Metadata_inducer1_concentration_unit"].astype(str)
        + "_"
        + feature_df["Metadata_inducer2_concentration"].astype(str)
        + "_"
        + feature_df["Metadata_inducer2_concentration_unit"].astype(str)
    ),
]
feature_df["Metadata_Dose"] = np.select(condlist=conditions, choicelist=results)


# ## Create a unique treatment column
# This columnd will be used later to annotate data into cell death groups

# In[17]:


# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_Dose"].astype(str)
    # + "_"
    # + feature_df['Metadata_inducer1_concentration_unit'].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration_unit"].astype(str)
).astype(str)


# In[18]:


replacement_dict = {
    "None": "0",
    "µ": "u",
    "nan": "0",
}
for pattern, replacement in replacement_dict.items():
    print(pattern, replacement)
    feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    ].replace(to_replace=str(pattern), value=str(replacement), regex=True)


# In[19]:


feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("media_ctr_0.0_0_Media_ctr_0_0", "media_ctr_0.0_0_Media_ctr_0.0_0")


# In[20]:


# need to convert to strings to save as parquet
# if the column is an object then convert it to a string
for column in feature_df.columns:
    if feature_df[column].dtype == "object":
        feature_df[column] = feature_df[column].astype(str)


# In[21]:


print(cell_type, len(feature_df["Metadata_Well"].unique()))


# In[22]:


# write to parquet file
feature_df.to_parquet(feature_df_out_path)
