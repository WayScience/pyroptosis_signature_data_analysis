#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import papermill as pm
import pyarrow as pa
import pyarrow.parquet as pq

# In[2]:


# Parameters
celltype = "PBMC"


# In[3]:


# Define inputs
feature_file = pathlib.Path(f"../data/{celltype}_sc_norm_fs.parquet")
feature_df = pq.read_table(feature_file).to_pandas()


# In[4]:


# remove uM in each row of the Metadata_inducer1_concentration column if it is present
if "Metadata_inducer1_concentration" in feature_df.columns:
    feature_df["Metadata_inducer1_concentration"] = feature_df[
        "Metadata_inducer1_concentration"
    ].str.replace("µM", "")


# In[5]:


feature_df["Metadata_inducer1_concentration"].unique()


# In[6]:


# define output file path
feature_df_out_path = pathlib.Path(f"../data/{celltype}_preprocessed_sc_norm.parquet")


# In[7]:


print(feature_df.shape)
feature_df.head()


# In[8]:


# removing costes features as they behave with great variance across all data
feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df.head()


# In[9]:


# replacing '/' in treatment dosage column to avoid errors in file interpolation including such strings
feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)


# In[10]:


# replace nan values with 0
feature_df["Metadata_inducer1_concentration"] = feature_df[
    "Metadata_inducer1_concentration"
].fillna(0)
feature_df["Metadata_inducer2_concentration"] = feature_df[
    "Metadata_inducer2_concentration"
].fillna(0)
feature_df["Metadata_inhibitor_concentration"] = feature_df[
    "Metadata_inhibitor_concentration"
].fillna(0)


# #### Combine Inducer1 and Inducer2 into one column

# In[11]:


# treatment column merge
conditions = [
    (feature_df["Metadata_inducer2"].isnull()),
    feature_df["Metadata_inducer2"].notnull(),
]

results = [
    (feature_df["Metadata_inducer1"]).astype(str),
    (feature_df["Metadata_inducer1"] + "_" + feature_df["Metadata_inducer2"]).astype(
        str
    ),
]
feature_df["Metadata_Treatment"] = np.select(condlist=conditions, choicelist=results)

# dose column merge
conditions = [
    (feature_df["Metadata_inducer2"].isnull()),
    feature_df["Metadata_inducer2"].notnull(),
]

results = [
    (feature_df["Metadata_inducer1_concentration"].astype(str)).astype(str),
    (
        feature_df["Metadata_inducer1_concentration"].astype(str)
        + "_"
        + feature_df["Metadata_inducer2_concentration"].astype(str)
    ).astype(str),
]
feature_df["Metadata_Dose"] = np.select(condlist=conditions, choicelist=results)


# In[12]:


feature_df["Metadata_inducer1_concentration"] = pd.to_numeric(
    feature_df["Metadata_inducer1_concentration"]
)


# ## N Beta Column condition generation
# columns generated to used for linear modeling where terms separated by '__' will be a beta coefficient

# In[13]:


# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_Dose"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)


# two beta of inducer1, inhibitor, and inhibitor concentration all as 1 beta term + inducer1 concentration as 2nd beta term
feature_df["twob_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
).astype(str)

# three beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor and inhibitor concentration as 3rd beta term
feature_df["threeb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)

# four beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor as 3rd beta term, and inhibitor concentration as 4th beta term
feature_df["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)


# In[14]:


feature_df_table = pa.Table.from_pandas(feature_df)
pq.write_table(feature_df_table, feature_df_out_path)
