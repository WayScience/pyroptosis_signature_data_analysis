#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing

# In[ ]:


# Parameters
cell_type = "SHSY5Y"


# In[ ]:


# set path to data
data_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)

preprocessing_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_clean.csv"
)

# read in data
nomic_df = pd.read_csv(data_path)


# In[ ]:


# select data only columns and make floats
nELISA_data_values = nomic_df.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")


# In[ ]:


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


# In[ ]:


# drop columns that are named with NSU
Metadata = nomic_df.drop(nomic_df.filter(like="NSU", axis=1), axis=1)
Metadata = Metadata.drop(nomic_df.filter(like="pgML", axis=1), axis=1)


# In[ ]:


# merge metadata and normalized data values
analysis_df = pd.concat([Metadata, nELISA_data_values_min_max_norm], axis=1)


# In[ ]:


# get rid of spaces
analysis_df.replace(" ", "_", regex=True, inplace=True)
# replace all "/" with "_per_"
analysis_df.replace("/", "_per_", regex=True, inplace=True)


# In[ ]:


analysis_df["inducer1_concentration"].replace(np.nan, 0, inplace=True)
analysis_df["inducer2_concentration"].replace(np.nan, 0, inplace=True)
analysis_df["inhibitor_concentration"].replace(np.nan, 0, inplace=True)


# In[ ]:


# replace %, µg/ml, µM, nM, and nan with ""
analysis_df["inducer1_concentration"].replace("%", "", regex=True, inplace=True)
analysis_df["inducer1_concentration"].replace(
    "_µg_per_ml", "", regex=True, inplace=True
)
analysis_df["inducer1_concentration"].replace("_µM", "", regex=True, inplace=True)
analysis_df["inducer1_concentration"].replace("_nM", "", regex=True, inplace=True)


# In[ ]:


# replace %, µg/ml, µM, nM, and nan with ""
analysis_df["inducer2_concentration"].replace("%", "", regex=True, inplace=True)
analysis_df["inducer2_concentration"].replace(
    "_µg_per_ml", "", regex=True, inplace=True
)
analysis_df["inducer2_concentration"].replace("_µM", "", regex=True, inplace=True)
analysis_df["inducer2_concentration"].replace("_nM", "", regex=True, inplace=True)


# In[ ]:


# replace %, µg/ml, µM, nM, and nan with ""
analysis_df["inhibitor_concentration"].replace("%", "", regex=True, inplace=True)
analysis_df["inhibitor_concentration"].replace(
    "_µg_per_ml", "", regex=True, inplace=True
)
analysis_df["inhibitor_concentration"].replace("_µM", "", regex=True, inplace=True)
analysis_df["inhibitor_concentration"].replace("_nM", "", regex=True, inplace=True)


# In[ ]:


# analysis_df['inhibitor_concentration'] to numeric
analysis_df["inhibitor_concentration"] = pd.to_numeric(
    analysis_df["inhibitor_concentration"]
)

# make each 3 decimal places
analysis_df["inhibitor_concentration"] = analysis_df["inhibitor_concentration"].round(3)


# In[ ]:


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


# In[ ]:


# one beta of inudcer1, inducer1 concentration, inhibitor, and inhibitor concentration all as 1 beta term
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "_"
    + analysis_df["Dose"].astype(str)
    + "_"
    # + analysis_df['inducer1_concentration_unit'].astype(str)
    # + "_"
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
    # + analysis_df["inhibitor_concentration"].astype(str)
    # + "__"
    + analysis_df["Dose"].astype(str)
).astype(str)

# three beta of inducer 1 as 1 beta term, inducer1 concentration as 2nd beta term, inhibitor and inhibitor concentration as 3rd beta term
analysis_df["threeb_Treatment_Dose_Inhibitor_Dose"] = (
    analysis_df["Treatment"]
    + "__"
    + analysis_df["Dose"].astype(str)
    + "__"
    # + analysis_df['inducer1_concentration_unit'].astype(str)
    # + "_"
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
    # + analysis_df['inducer1_concentration_unit'].astype(str)
    # + "_"
    + analysis_df["inhibitor"].astype(str)
    + "__"
    + analysis_df["inhibitor_concentration"].astype(str)
).astype(str)


# In[ ]:


# fix strings in Metadata_Treatment column
# replaceall "None" with "0"
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("None", "0")
# replace all mu with u
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("µ", "u")
# replace all "nan" with "0"
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("nan", "0")


# In[ ]:


# _0.03 to _0.025
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.03", "_0.025", regex=False)
# _0.010_ to _0.100_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.01_", "_0.010_", regex=False)
# _0.10_ to _0.100_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.10_", "_0.100_", regex=False)
# _0.10_ to _0.100_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_0.1_", "_0.100_", regex=False)
# _1_ to _1.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_1_", "_1.000_", regex=False)
# _1.0_ to _1.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_1.0_", "_1.000_", regex=False)
# _10_ to _10.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_10_", "_10.000_", regex=False)
# _100_ to _100.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_100_", "_100.000_", regex=False)
# _100.0_ to _100.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_100.0_", "_100.000_", regex=False)
# _5_ to _5.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_5_", "_5.000_", regex=False)
# _20_ to _20.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_20_", "_20.000_", regex=False)
# _30_ to _30.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_30_", "_30.000_", regex=False)
# _2.5 to _2.500_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_2.5_", "_2.500_", regex=False)
# _3.0 to _3.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_3.0_", "_3.000_", regex=False)
# _3 to _3.000_
analysis_df["oneb_Treatment_Dose_Inhibitor_Dose"] = analysis_df[
    "oneb_Treatment_Dose_Inhibitor_Dose"
].str.replace("_3_", "_3.000_", regex=False)


# In[ ]:


# save nELISA_plate_430420_no_inhibitors dataframe to csv file
analysis_df.to_csv(preprocessing_path, index=False)
