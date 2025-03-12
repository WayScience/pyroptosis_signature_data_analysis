#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pathlib
import random

import numpy as np
import pandas as pd
import toml
from copairs import map
from copairs.matching import assign_reference_index

# check if in a jupyter notebook
try:
    cfg = get_ipython().config
    in_notebook = True
except NameError:
    in_notebook = False


# In[2]:


if not in_notebook:
    parser = argparse.ArgumentParser(description="Match pairs of samples")
    parser.add_argument("--shuffle", action="store_true", help="Shuffle the data")

    args = parser.parse_args()
    shuffle = args.shuffle
else:
    shuffle = True


# In[3]:


# load in the treatment groups
ground_truth = pathlib.Path(
    "../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
).resolve(strict=True)
# load in the ground truth
ground_truth = toml.load(ground_truth)
apoptosis_ground_truth = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_ground_truth = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control_ground_truth = ground_truth["Healthy"]["healthy_groups_list"]

map_out_dir = pathlib.Path("../data/processed/mAP_scores/secretome/")
map_out_dir.mkdir(exist_ok=True, parents=True)


# In[4]:


agg_data = pathlib.Path(
    "../../data/PBMC_preprocessed_sc_norm_aggregated_nomic.parquet"
).resolve(strict=True)
df = pd.read_parquet(agg_data)
# rename oneb_Metadata_Treatment_Dose_Inhibitor_Dose to Metadata_Treatment
df = df.rename(
    columns={"oneb_Metadata_Treatment_Dose_Inhibitor_Dose": "Metadata_Treatment"}
)
df = df.filter(regex="Metadata|NSU")
df.head()


# In[5]:


# add apoptosis, pyroptosis and healthy columns to dataframe
df["Apoptosis"] = df.apply(
    lambda row: row["Metadata_Treatment"] in apoptosis_ground_truth,
    axis=1,
)
df["Pyroptosis"] = df.apply(
    lambda row: row["Metadata_Treatment"] in pyroptosis_ground_truth,
    axis=1,
)
df["Control"] = df.apply(
    lambda row: row["Metadata_Treatment"] in control_ground_truth,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column
df["Metadata_labels"] = df.apply(
    lambda row: "Apoptosis"
    if row["Apoptosis"]
    else "Pyroptosis"
    if row["Pyroptosis"]
    else "Control",
    axis=1,
)
metadata_labels = df.pop("Metadata_labels")
df.insert(1, "Metadata_labels", metadata_labels)
# # drop apoptosis, pyroptosis, and healthy columns
df.drop(columns=["Apoptosis", "Pyroptosis", "Control"], inplace=True)


# In[6]:


if shuffle:
    random.seed(0)
    # permutate the data
    for col in df.columns:
        df[col] = np.random.permutation(df[col])


# In[7]:


reference_col = "Metadata_reference_index"
df_activity = assign_reference_index(
    df,
    "Metadata_Treatment == 'DMSO_0.100_%_DMSO_0.025_%'",
    reference_col=reference_col,
    default_value=-1,
)
df_activity.head()


# In[8]:


pos_sameby = ["Metadata_Treatment", "Metadata_labels", reference_col]
pos_diffby = []
neg_sameby = []
neg_diffby = ["Metadata_Treatment", reference_col]
metadata = df_activity.filter(regex="Metadata")
profiles = df_activity.filter(regex="^(?!Metadata)").values

activity_ap = map.average_precision(
    metadata, profiles, pos_sameby, pos_diffby, neg_sameby, neg_diffby
)

activity_ap = activity_ap.query("Metadata_Treatment != 'DMSO_0.100_%_DMSO_0.025_%'")
activity_ap.head()


# In[9]:


activity_map = map.mean_average_precision(
    activity_ap, pos_sameby, null_size=1000000, threshold=0.05, seed=0
)
activity_map["-log10(p-value)"] = -activity_map["corrected_p_value"].apply(np.log10)
# flatten the multi-index columns to make it easier to work with
activity_map.reset_index(inplace=True)
activity_map.head()


# In[10]:


if shuffle:
    activity_map.to_parquet(map_out_dir / "activity_map_shuffled.parquet")
else:
    activity_map.to_parquet(map_out_dir / "activity_map.parquet")
