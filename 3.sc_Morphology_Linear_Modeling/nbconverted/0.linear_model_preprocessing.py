#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import numpy as np
import pandas as pd
import pdfkit
import plotly.express as px
import pyarrow as pa
import pyarrow.parquet as pq
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

# import Union


sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

# In[2]:


# Define inputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/SHSY5Y_run_sc_norm.parquet"
)
feature_df = pq.read_table(feature_file).to_pandas()
# feature_df = pd.read_csv(feature_file, engine="pyarrow")


# In[3]:


# define output file path
feature_df_out_path = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/feature_df_sc_norm.parquet"
)


# In[4]:


print(feature_df.shape)
feature_df.head()


# In[5]:


# removing costes features as they behave with great variance across all data
feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df


# In[6]:


# replacing '/' in treatment dosage column to avoid errors in file interpolation including such strings
feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)


# In[7]:


# Recycled code from: https://github.com/WayScience/NF1_SchwannCell_data/blob/main/5_analyze_data/notebooks/linear_model/fit_linear_model.ipynb
cell_count_df = (
    feature_df.groupby("Metadata_Well")["Metadata_Plate"]
    .count()
    .reset_index()
    .rename(columns={"Metadata_Plate": "Metadata_number_of_singlecells"})
)


# In[8]:


feature_df.head()


# In[9]:


# show max column in pandas df
pd.set_option("display.max_columns", 100)


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
feature_df["Metadata_Treatment"] = np.select(conditions, results)

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
feature_df["Metadata_Dose"] = np.select(conditions, results)


# ## 1 Beta Column condition generation

# In[12]:


feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_Dose"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)

feature_df["twob_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "_"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
).astype(str)

feature_df["threeb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "_"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)
feature_df["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"] = (
    feature_df["Metadata_Treatment"]
    + "__"
    + feature_df["Metadata_Dose"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor"].astype(str)
    + "__"
    + feature_df["Metadata_inhibitor_concentration"].astype(str)
).astype(str)


# In[ ]:


feature_df["Metadata_inducer1_concentration"] = pd.to_numeric(
    feature_df["Metadata_inducer1_concentration"]
)


# In[ ]:


# feature_df.to_csv(feature_df_out_path, index=False)
feature_df_table = pa.Table.from_pandas(feature_df)
pq.write_table(feature_df_table, feature_df_out_path)
