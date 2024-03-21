#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
import toml

# In[2]:


output_path = pathlib.Path("../screenA/idr0000-screenA-processed.txt")


# In[3]:


# path to the data (processed)
PBMC_path = pathlib.Path(
    "../../data/PBMC_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
)
SHSY5Y_path = pathlib.Path(
    "../../data/SHSY5Y_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
)
ground_truth_path = pathlib.Path(
    "../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
)
# read the data
PBMC = pd.read_parquet(PBMC_path)
# read the ground truth
ground_truth = toml.load(ground_truth_path)
apoptosis = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control = ground_truth["Healthy"]["healthy_groups_list"]
# # make a column for the ground truth
PBMC
PBMC["ground_truth"] = "other"
PBMC["apoptosis"] = PBMC["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(apoptosis)
PBMC["pyroptosis"] = PBMC["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    pyroptosis
)
PBMC["control"] = PBMC["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(control)

# merge the apoptosis, pyroptosis and control columns into the ground_truth column
PBMC.loc[PBMC["apoptosis"], "ground_truth"] = "apoptosis"
PBMC.loc[PBMC["pyroptosis"], "ground_truth"] = "pyroptosis"
PBMC.loc[PBMC["control"], "ground_truth"] = "control"

# drop the apoptosis, pyroptosis and control columns
PBMC = PBMC.drop(
    columns=[
        "apoptosis",
        "pyroptosis",
        "control",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ]
)


SHSY5Y = pd.read_parquet(SHSY5Y_path)
# read the ground truth
ground_truth = toml.load(ground_truth_path)

apoptosis = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control = ground_truth["Healthy"]["healthy_groups_list"]
# make a column for the ground truth
SHSY5Y["ground_truth"] = "other"
SHSY5Y["apoptosis"] = SHSY5Y["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    apoptosis
)
SHSY5Y["pyroptosis"] = SHSY5Y["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
    pyroptosis
)
SHSY5Y["control"] = SHSY5Y["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(control)

# merge the apoptosis, pyroptosis and control columns into the ground_truth column
SHSY5Y.loc[SHSY5Y["apoptosis"], "ground_truth"] = "apoptosis"
SHSY5Y.loc[SHSY5Y["pyroptosis"], "ground_truth"] = "pyroptosis"
SHSY5Y.loc[SHSY5Y["control"], "ground_truth"] = "control"

# drop the apoptosis, pyroptosis and control columns
SHSY5Y = SHSY5Y.drop(
    columns=[
        "apoptosis",
        "pyroptosis",
        "control",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    ]
)

# concat the two dataframes
data = pd.concat([PBMC, SHSY5Y], axis=0)
data.head()
data.reset_index(drop=True, inplace=True)
data


# In[4]:


# read in the nELISA data
nELISA_PBMC_path = pathlib.Path(
    "/home/lippincm/Documents/ML/Interstellar_Analysis/2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_PBMC_clean.parquet"
)
nELISA_SHSY5Y_path = pathlib.Path(
    "/home/lippincm/Documents/ML/Interstellar_Analysis/2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_SHSY5Y_clean.parquet"
)
nELISA_PBMC = pd.read_parquet(nELISA_PBMC_path)
nELISA_SHSY5Y = pd.read_parquet(nELISA_SHSY5Y_path)
nELISA_PBMC
# drop all all columns that are not position x or do not conatin [NSU]
nELISA_PBMC = nELISA_PBMC[
    nELISA_PBMC.columns[nELISA_PBMC.columns.str.contains("position|NSU")]
]
nELISA_SHSY5Y = nELISA_SHSY5Y[
    nELISA_SHSY5Y.columns[nELISA_SHSY5Y.columns.str.contains("position|NSU")]
]
nELISA_PBMC.drop(columns=["plate_position", "position_y"], inplace=True)
nELISA_SHSY5Y.drop(columns=["plate_position", "position_y"], inplace=True)
# merge the two dataframes
nELISA = pd.concat([nELISA_PBMC, nELISA_SHSY5Y], axis=0)
nELISA


# In[5]:


# add the nELISA data to the data dataframe merge on the position_x column
data = pd.merge(
    data, nELISA, left_on="Metadata_Well", right_on="position_x", how="left"
)


# In[6]:


print(PBMC.shape, SHSY5Y.shape)


# In[7]:


path = pathlib.Path("../screenA/idr0000-screenA-library.txt")
df = pd.read_csv(path, sep="\t")
df.head()
print(df.columns.to_list())
# keep only the columns we need
df = df[
    [
        "Well",
        "Characteristics[cell type]",
        "inhibitor",
        "inhibitor_concentration",
        "inhibitor_concentration_unit",
        "inducer1",
        "inducer1_concentration",
        "inducer1_concentration_unit",
        "inducer2",
        "inducer2_concentration",
        "inducer2_concentration_unit",
        "Plate",
    ]
]


# In[8]:


# add df and data
data = data.merge(df, left_on="Metadata_Well", right_on="Well")
# replace NA with ""
data.fillna("", inplace=True)
data


# In[9]:


# write the data to a txt file
data.to_csv(output_path, sep="\t", index=False)
