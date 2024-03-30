#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import numpy as np
import pandas as pd
import toml
import umap
from tqdm import tqdm

# In[2]:


# Parameters
cell_type = "SHSY5Y"


# In[3]:


# read in toml file

# set up the path
toml_path = pathlib.Path("../utils/params.toml")
# read in the toml file
params = toml.load(toml_path)
list_of_treatments = params["list_of_treatments"]["treatments"]
print(len(list_of_treatments))
print(list_of_treatments)


# In[4]:


# Set path to parquet file
path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet").resolve(
    strict=True
)
# Read in parquet file
df = pd.read_parquet(path)


# In[5]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])
# define which columns are data and which are descriptive
df_descriptive = df[df_metadata]
df_values = df.drop(columns=df_metadata)


# In[6]:


anova_path = pathlib.Path(f"../results/{cell_type}_combined.parquet")

anova_results = pd.read_parquet(anova_path)
anova_results.head()


# Where
# a_h = apoptosis vs healthy
# a_p = apoptosis vs pyroptosis
# h_p = healthy vs pyroptosis
#

# In[ ]:


# create a column that adds group1 and group2 together
anova_results["group"] = anova_results["group1"] + "_" + anova_results["group2"]
print(anova_results.shape)

# filter out rows that have p-adj_abs > 0.05
anova_results = anova_results[anova_results["p-adj_abs"] < 0.05]
print(anova_results.shape)

# change the group names to replace healthy with control
anova_results["group"] = anova_results["group"].str.replace("healthy", "control")

# create the three df sets for a venn diagram
a_h = anova_results[anova_results["group"] == "apoptosis_control"]["features"]
a_p = anova_results[anova_results["group"] == "apoptosis_pyroptosis"]["features"]
h_p = anova_results[anova_results["group"] == "control_pyroptosis"]["features"]

# create a list of the three df sets
a_h_list = a_h.tolist()
a_p_list = a_p.tolist()
h_p_list = h_p.tolist()

# add sets together
a_h__a_p = np.union1d(a_h_list, a_p_list)
a_h__h_p = np.union1d(a_h_list, h_p_list)
a_p__h_p = np.union1d(a_p_list, h_p_list)


# In[ ]:


# get the unique features for each set
a_h_unique = np.setdiff1d(a_h_list, a_p__h_p)
print(len(a_h_unique))

a_p_unique = np.setdiff1d(a_p_list, a_h__h_p)
print(len(a_p_unique))

h_p_unique = np.setdiff1d(h_p_list, a_h__a_p)
print(len(h_p_unique))

# get the common features for each set
a_h__a_p_common = np.intersect1d(a_h_list, a_p_list)
a_h__a_p_common = np.setdiff1d(a_h__a_p_common, h_p_list)
print(len(a_h__a_p_common))

a_h__h_p_common = np.intersect1d(a_h_list, h_p_list)
a_h__h_p_common = np.setdiff1d(a_h__h_p_common, a_p_list)
print(len(a_h__h_p_common))

a_p__h_p_common = np.intersect1d(a_p_list, h_p_list)
a_p__h_p_common = np.setdiff1d(a_p__h_p_common, a_h_list)
print(len(a_p__h_p_common))

# all three set intersection
a_h__a_p__h_p_common = np.intersect1d(a_h_list, a_p_list)
a_h__a_p__h_p_common = np.intersect1d(a_h__a_p__h_p_common, h_p_list)
print(len(a_h__a_p__h_p_common))


# In[ ]:


# create a list of each list of features
dict_of_feature_lists = {}
dict_of_feature_lists["a_h_unique"] = list(a_h_unique)
dict_of_feature_lists["a_p_unique"] = list(a_p_unique)
dict_of_feature_lists["h_p_unique"] = list(h_p_unique)
dict_of_feature_lists["a_h__a_p_common"] = list(a_h__a_p_common)
dict_of_feature_lists["a_h__h_p_common"] = list(a_h__h_p_common)
dict_of_feature_lists["a_p__h_p_common"] = list(a_p__h_p_common)
dict_of_feature_lists["a_h__a_p__h_p_common"] = list(a_h__a_p__h_p_common)


# In[ ]:


# set umap parameters
umap_params = umap.UMAP(
    n_components=2,
    spread=1.1,
    min_dist=0.8,
    init="random",
    metric="cosine",
    # random_state=0,
    n_jobs=-1,
)


# In[ ]:


final_df_dict = {}
for key, value in tqdm(dict_of_feature_lists.items()):
    print(key)
    print(len(value))
    df = df_values[df_values.columns[df_values.columns.isin(value)]]
    umap_results = umap_params.fit_transform(df)
    results_df = pd.DataFrame(umap_results, columns=["UMAP1", "UMAP2"])
    results_df.loc[:, "Metadata_Treatment_Dose_Inhibitor_Dose"] = df_descriptive[
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    ]
    results_df.loc[:, "Dataset_comparison"] = key
    final_df_dict[key] = results_df
final_df = pd.concat(final_df_dict.values(), ignore_index=True)


# In[ ]:


# write out the results
out_path = pathlib.Path(f"../results/{cell_type}_combined_subset_UMAP_results.parquet")
final_df.to_parquet(out_path)
final_df.head()
