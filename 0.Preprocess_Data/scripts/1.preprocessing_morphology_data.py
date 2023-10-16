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
cell_type = "SHSY5Y"


# In[3]:


# Define inputs
feature_file = pathlib.Path(f"../data/{cell_type}_sc_norm_fs.parquet")
feature_df = pq.read_table(feature_file).to_pandas()


# In[4]:


# replace all " " with "_" in all values of the dataframe
feature_df = feature_df.replace(to_replace=" ", value="_", regex=True)


# In[5]:


# remove uM in each row of the Metadata_inducer1_concentration column
feature_df["Metadata_inducer1_concentration"] = feature_df[
    "Metadata_inducer1_concentration"
].str.replace("µM", "")


# In[6]:


feature_df["Metadata_inducer1_concentration"].unique()


# In[7]:


# define output file path
feature_df_out_path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_norm.parquet")


# In[8]:


print(feature_df.shape)
feature_df.head()


# In[9]:


# removing costes features as they behave with great variance across all data
feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df.head()


# In[10]:


# replacing '/' in treatment dosage column to avoid errors in file interpolation including such strings
feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)


# In[11]:


# replace nan values with 0

columns_to_fill = [
    "Metadata_inducer1_concentration",
    "Metadata_inducer2_concentration",
    "Metadata_inhibitor_concentration",
]
feature_df[columns_to_fill].fillna(0, inplace=True)


# #### Combine Inducer1 and Inducer2 into one column

# In[12]:


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


# In[13]:


feature_df["Metadata_inducer1_concentration"] = pd.to_numeric(
    feature_df["Metadata_inducer1_concentration"]
)


# ## N Beta Column condition generation
# columns generated to used for linear modeling where terms separated by '__' will be a beta coefficient

# In[14]:


# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_Dose"].astype(str)
    + "_"
    # + feature_df['Metadata_inducer1_concentration_unit'].astype(str)
    # + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration_unit"].astype(str)
).astype(str)


# two beta of inducer1, inhibitor, and inhibitor concentration all as 1 beta term + inducer1 concentration as 2nd beta term
feature_df["twob_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    # + feature_df["Metadata_inhibitor_concentration"].astype(str)
    # + "__"
    + feature_df["Metadata_Dose"].astype(str)
).astype(str)

# three beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor and inhibitor concentration as 3rd beta term
feature_df["threeb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
    + "__"
    # + feature_df['Metadata_inducer1_concentration_unit'].astype(str)
    # + "_"
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
    # + feature_df['Metadata_inducer1_concentration_unit'].astype(str)
    # + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)


# In[15]:


# fix strings in Metadata_Treatment column
# replaceall "None" with "0"
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("None", "0")
# replace all mu with u
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("µ", "u")
# replace all "nan" with "0"
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("nan", "0")


# In[16]:


# _0.03 to _0.025
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.03", "_0.025")
# _0.10_ to _0.100_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.10_", "_0.100_")
# _0.10_ to _0.100_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.10_", "_0.100_")
# _1_ to _1.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_1.0_", "_1.000_")
# _10_ to _10.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_10.0_", "_10.000_")
# _100_ to _100.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_100.0_", "_100.000_")
# _5_ to _5.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_5.0_", "_5.000_")
# _20_ to _20.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_20.0_", "_20.000_")
# _3.0_ to _3.000_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_3.0_", "_3.000_")
# _0.1_ to _0.100_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.1_", "_0.100_")
# _2.5_ to _2.500_
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("_2.5_", "_2.500_")

# mix the media treatment
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("media_ctr_0_ug_per_ml_Media_ctr_0_0", "media_ctr_0_0_Media_ctr_0.0_0")
feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = feature_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].str.replace("media_ctr_0_0_Media_ctr_0_0", "media_ctr_0_0_Media_ctr_0.0_0")


# In[17]:


feature_df_table = pa.Table.from_pandas(feature_df)
pq.write_table(feature_df_table, feature_df_out_path)
