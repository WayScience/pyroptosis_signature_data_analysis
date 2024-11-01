#!/usr/bin/env python
# coding: utf-8

# This notebook performs data splits on the wells on the bulk data.

# In[1]:


import argparse
import pathlib

import pandas as pd
import toml

# In[2]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", default="all")

args = argparser.parse_args()

cell_type = args.cell_type


# In[3]:


# set path to the import data
data_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
)
# results path
results_path = pathlib.Path("../results/").resolve()
results_path.mkdir(exist_ok=True)

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pd.read_parquet(data_path)
data_df.head()


# In[4]:


# get all of the features
features = data_df.columns.to_list()

nuclei_features = []
pm_features = []
gsdmd_features = []
mito_features = []
er_features = []
other_features = []
metadata_features = []
correlation_features = []

# Separate the features into distinct types
# code loop modified from:
# https://github.com/WayScience/nuclear_speckles_analysis/blob/main/1.regression_modelling/1.train_models.ipynb
# Thank you Jenna Tomkinson for the code!
for feature in features:
    # check if the feature contains "Metadata"
    if "Metadata" in feature:
        metadata_features.append(feature)
    else:
        parts = feature.split("_")

        if "Correlation" in parts:  # Check if it's a correlation feature
            correlation_features.append(feature)

        else:  # Non-correlation features
            if "CorrDNA" in parts:
                nuclei_features.append(feature)
            elif "CorrPM" in parts:
                pm_features.append(feature)
            elif "CorrGasdermin" in parts:
                gsdmd_features.append(feature)
            elif "CorrMito" in parts:
                mito_features.append(feature)
            elif "CorrER" in parts:
                er_features.append(feature)
            else:
                other_features.append(feature)


# In[5]:


print(f"nuclei_features: {len(nuclei_features)}")
print(f"pm_features: {len(pm_features)}")
print(f"gsdmd_features: {len(gsdmd_features)}")
print(f"mito_features: {len(mito_features)}")
print(f"er_features: {len(er_features)}")
print(f"other_features: {len(other_features)}")
print(f"correlation_features: {len(correlation_features)}")
print(f"metadata_features: {len(metadata_features)}")
print(f"total features: {len(features)}")
# check if all featues are accounted for
assert len(features) == len(nuclei_features) + len(pm_features) + len(
    gsdmd_features
) + len(mito_features) + len(er_features) + len(other_features) + len(
    correlation_features
) + len(
    metadata_features
)


# In[6]:


# create the combinations of features
# out of 5 channels, we can have 0, 1, 2, 3, 4, or 5 channels
# even with 0 channels, we still have the non channel features (object-based)
# these are areshape features
# set the feature combination lists
dict_of_feature_combinations = {
    "No_channels": other_features,
    "CorrDNA": nuclei_features,
    "CorrPM": pm_features,
    "CorrGasdermin": gsdmd_features,
    "CorrMito": mito_features,
    "CorrER": er_features,
    "CorrDNA_CorrPM": nuclei_features + pm_features,
    "CorrDNA_CorrGasdermin": nuclei_features + gsdmd_features,
    "CorrDNA_CorrMito": nuclei_features + mito_features,
    "CorrDNA_CorrER": nuclei_features + er_features,
    "CorrPM_CorrGasdermin": pm_features + gsdmd_features,
    "CorrPM_CorrMito": pm_features + mito_features,
    "CorrPM_CorrER": pm_features + er_features,
    "CorrGasdermin_CorrMito": gsdmd_features + mito_features,
    "CorrGasdermin_CorrER": gsdmd_features + er_features,
    "CorrMito_CorrER": mito_features + er_features,
    "CorrDNA_CorrPM_CorrGasdermin": nuclei_features + pm_features + gsdmd_features,
    "CorrDNA_CorrPM_CorrMito": nuclei_features + pm_features + mito_features,
    "CorrDNA_CorrPM_CorrER": nuclei_features + pm_features + er_features,
    "CorrDNA_CorrGasdermin_CorrMito": nuclei_features + gsdmd_features + mito_features,
    "CorrDNA_CorrGasdermin_CorrER": nuclei_features + gsdmd_features + er_features,
    "CorrDNA_CorrMito_CorrER": nuclei_features + mito_features + er_features,
    "CorrPM_CorrGasdermin_CorrMito": pm_features + gsdmd_features + mito_features,
    "CorrPM_CorrGasdermin_CorrER": pm_features + gsdmd_features + er_features,
    "CorrPM_CorrMito_CorrER": pm_features + mito_features + er_features,
    "CorrGasdermin_CorrMito_CorrER": gsdmd_features + mito_features + er_features,
    "CorrDNA_CorrPM_CorrGasdermin_CorrMito": nuclei_features
    + pm_features
    + gsdmd_features
    + mito_features,
    "CorrDNA_CorrPM_CorrGasdermin_CorrER": nuclei_features
    + pm_features
    + gsdmd_features
    + er_features,
    "CorrDNA_CorrPM_CorrMito_CorrER": nuclei_features
    + pm_features
    + mito_features
    + er_features,
    "CorrDNA_CorrGasdermin_CorrMito_CorrER": nuclei_features
    + gsdmd_features
    + mito_features
    + er_features,
    "CorrPM_CorrGasdermin_CorrMito_CorrER": pm_features
    + gsdmd_features
    + mito_features
    + er_features,
    "All_channels": nuclei_features
    + pm_features
    + gsdmd_features
    + mito_features
    + er_features
    + other_features
    + correlation_features,
}
# loop through each feature combination and add the metadata features
for combination in dict_of_feature_combinations:
    if combination == "No_channels":
        temp_correlation_features = other_features
    elif "_" in combination:
        channels = combination.split("_")
        temp_correlation_features = []
        for feature in correlation_features:
            if all(channel not in feature for channel in channels):
                temp_correlation_features.append(feature)
    else:
        temp_correlation_features = []
        for feature in correlation_features:
            if combination not in feature:
                temp_correlation_features.append(feature)

    num_featues = len(dict_of_feature_combinations[combination])
    dict_of_feature_combinations[combination] += temp_correlation_features

    print(
        f"{len(dict_of_feature_combinations[combination]) - num_featues} correlation features added to {combination}"
    )
    dict_of_feature_combinations[combination] += metadata_features


# In[7]:


# save the dict to a toml file

toml_path = pathlib.Path(f"../results/feature_combinations_{cell_type}.toml")
with open(toml_path, "w") as f:
    toml.dump(dict_of_feature_combinations, f)

# write the keys to a txt file with each key on a new line
# this is for easy retrieval of the keys in bash
txt_path = pathlib.Path("../results/feature_combinations_keys.txt")
with open(txt_path, "w") as f:
    for key in dict_of_feature_combinations:
        f.write(f"{key}\n")
