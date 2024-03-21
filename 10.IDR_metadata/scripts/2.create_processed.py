#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd
import toml

# In[2]:


output_path = pathlib.Path("../IDR_metadata/screenA/")
output_path.mkdir(exist_ok=True, parents=True)
output_path = output_path / "idr0000-screenA-processed.txt"


# In[3]:


ground_truth_path = pathlib.Path(
    "../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
)

# read the ground truth
ground_truth = toml.load(ground_truth_path)
apoptosis = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control = ground_truth["Healthy"]["healthy_groups_list"]


# In[4]:


# path to the data (processed)
PBMC_path = pathlib.Path(
    "../../data/PBMC_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
)
SHSY5Y_path = pathlib.Path(
    "../../data/SHSY5Y_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
)
# read the data
PBMC = pd.read_parquet(PBMC_path)
SHSY5Y = pd.read_parquet(SHSY5Y_path)


# In[5]:


print(SHSY5Y.shape, PBMC.shape)
SHSY5Y["Metadata_Well"].unique()


# In[6]:


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


# In[7]:


# write the data to a txt file
data.to_csv(output_path, sep="\t", index=False)
