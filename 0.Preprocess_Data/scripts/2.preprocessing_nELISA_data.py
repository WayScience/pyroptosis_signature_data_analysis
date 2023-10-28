#!/usr/bin/env python
# coding: utf-8

# This noteboook pre-processes the nELISA data to be ready for exploratory analysis and machine learning.

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing

# In[2]:


# Parameters
cell_type = "PBMC"


# In[3]:


# set path to data
data_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
).resolve(strict=True)

preprocessing_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_clean.parquet"
).resolve(strict=True)

# read in data
nomic_df = pd.read_csv(data_path)


# In[4]:


# select data only columns and make floats
nELISA_data_values = nomic_df.filter(like="NSU", axis=1).astype("float")


# In[5]:


# normalize data via max value in each column
max_values = nELISA_data_values.max()  # find max value in each column
nELISA_data_values_sensor_max_norm = nELISA_data_values.div(
    max_values
)  # divide each value in each column by max value in that column
nELISA_data_values_sensor_max_norm.head()
# min max normalization via sklearn

# normalize data via min max normalization
min_max_scaler = preprocessing.MinMaxScaler()
nELISA_data_values_min_max_norm = min_max_scaler.fit_transform(nELISA_data_values)
nELISA_data_values_min_max_norm = pd.DataFrame(
    nELISA_data_values_min_max_norm, columns=nELISA_data_values.columns
)


# In[6]:


# drop columns that are named with NSU
Metadata = nomic_df.drop(nomic_df.filter(like="NSU", axis=1), axis=1).drop(
    nomic_df.filter(like="pgML", axis=1), axis=1
)


# In[7]:


# merge metadata and normalized data values
analysis_df = pd.concat([Metadata, nELISA_data_values_min_max_norm], axis=1)


# In[8]:


# get rid of spaces
analysis_df.replace(" ", "_", regex=True, inplace=True)
# replace all "/" with "_per_"
analysis_df.replace("/", "_per_", regex=True, inplace=True)


# In[9]:


# replace nans with 0 in this case this is okay because the real nans were already removed on the basis of treatment name
analysis_df["inducer1_concentration"].replace(np.nan, 0, inplace=True)
analysis_df["inducer2_concentration"].replace(np.nan, 0, inplace=True)
analysis_df["inhibitor_concentration"].replace(np.nan, 0, inplace=True)


# In[10]:


def perform_replacements(text: str) -> str:
    """
    Function to replace special characters in text.
    Replaces `%`, `_µM`, `_nM`, `_µg_per_ml` with empty string.

    Parameters
    ----------
    text : str
        Text to be modified.

    Returns
    -------
    str
        Modified text.
    """
    replacements = {
        "%": "",
        "_µM": "",
        "_nM": "",
        "_µg_per_ml": "",
    }
    for key, value in replacements.items():
        text = str(text).replace(key, value)
    return text


# Columns to which you want to apply the changes
columns_to_apply = [
    "inducer1_concentration",
    "inducer2_concentration",
    "inhibitor_concentration",
]

# Applying the custom function to selected columns using apply
analysis_df[columns_to_apply] = analysis_df[columns_to_apply].apply(
    lambda x: x.apply(perform_replacements)
)


# In[11]:


# using an f string make "inducer1_concentration" have 3 decimal places
analysis_df["inducer1_concentration"] = analysis_df["inducer1_concentration"].apply(
    lambda x: f"{float(x):.3f}" if float(x) != 0 else float(x)
)
analysis_df["inducer2_concentration"] = analysis_df["inducer2_concentration"].apply(
    lambda x: f"{float(x):.3f}" if float(x) != 0 else float(x)
)
analysis_df["inhibitor_concentration"] = analysis_df["inhibitor_concentration"].apply(
    lambda x: f"{float(x):.3f}" if float(x) != 0 else float(x)
)


# In[12]:


# treatment column merge
conditions = [
    (analysis_df["inducer2"].isnull()),
    analysis_df["inducer2"].notnull(),
]

results = [
    (analysis_df["inducer1"]).astype(str),
    (analysis_df["inducer1"] + "_" + analysis_df["inducer2"].astype(str)),
]
analysis_df["Treatment"] = np.select(condlist=conditions, choicelist=results)


# dose column merge
results = [
    (
        analysis_df["inducer1_concentration"].astype(str)
        + "_"
        + analysis_df["inducer1_concentration_unit"].astype(str)
    ),
    (
        analysis_df["inducer1_concentration"].astype(str)
        + "_"
        + analysis_df["inducer1_concentration_unit"].astype(str)
        + "_"
        + analysis_df["inducer2_concentration"].astype(str)
        + "_"
        + analysis_df["inducer2_concentration_unit"].astype(str)
    ),
]
analysis_df["Dose"] = np.select(condlist=conditions, choicelist=results)


# In[13]:


# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "_"
    + analysis_df["Dose"].astype(str)
    + "_"
    + analysis_df["inhibitor"].astype(str)
    + "_"
    + analysis_df["inhibitor_concentration"].astype(str)
    + "_"
    + analysis_df["inhibitor_concentration_unit"].astype(str)
).astype(str)


# two beta of inducer1, inhibitor, and inhibitor concentration all as 1 beta term + inducer1 concentration as 2nd beta term
analysis_df["twob_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "_"
    + analysis_df["inhibitor"].astype(str)
    + "_"
    + analysis_df["Dose"].astype(str)
).astype(str)

# three beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor and inhibitor concentration as 3rd beta term
analysis_df["threeb_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "__"
    + analysis_df["Dose"].astype(str)
    + "__"
    + analysis_df["inhibitor"].astype(str)
    + "_"
    + analysis_df["inhibitor_concentration"].astype(str)
).astype(str)

# four beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor as 3rd beta term, and inhibitor concentration as 4th beta term
analysis_df["fourb_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "__"
    + analysis_df["Dose"].astype(str)
    + "__"
    + analysis_df["inhibitor"].astype(str)
    + "__"
    + analysis_df["inhibitor_concentration"].astype(str)
).astype(str)


# In[14]:


replacement_dict = {
    "None": "0",
    "µ": "u",
    "nan": "0",
}
for pattern, replacement in replacement_dict.items():
    print(pattern, replacement)
    analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
        "oneb_Treatment_Dose_Inhibitor_Dose"
    ].replace(to_replace=str(pattern), value=str(replacement), regex=True)


# In[15]:


# _0.03 to _0.025 for the DMSO concentration
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.03", "_0.025", regex=False)

# _0.0250 to _0.025 for the DMSO concentration
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.0250", "_0.025", regex=False)


# In[16]:


# need to convert to strings to save as parquet
# if the column is an object then convert it to a string
for column in analysis_df.columns:
    if analysis_df[column].dtype == "object":
        analysis_df[column] = analysis_df[column].astype(str)


# In[17]:


analysis_df.to_parquet(preprocessing_path)
