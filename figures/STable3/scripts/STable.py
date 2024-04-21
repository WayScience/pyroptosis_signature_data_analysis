#!/usr/bin/env python
# coding: utf-8

# This file generates a csv and markdown table for the features of the dataset.

# In[1]:


import pathlib

import pandas as pd
import toml

# In[2]:


# set the path to the data
path = pathlib.Path(
    "../../../1.Exploratory_Data_Analysis/results/PBMC_combined.parquet"
)

# load in the data
df = pd.read_parquet(path)


# In[3]:


# path to the ground truth
ground_truth_path = pathlib.Path(
    "../../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
)
# load in the ground truth
ground_truth = toml.load(ground_truth_path)
apoptosis_ground_truth = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_ground_truth = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control_ground_truth = ground_truth["Healthy"]["healthy_groups_list"]


# In[4]:


# change the p-adj into absolute values
df["p-adj"] = df["p-adj"].abs()
df.head()
# select row that have p-adj < 0.05
df = df[df["p-adj"] < 0.05]


# In[5]:


# add the group1 and group2 columns into 1 column
df["group"] = df["group1"] + "_" + df["group2"]


# In[6]:


# set theory for apoptosis control and pyroptosis
apoptosis_vs_healthy = df[df["group"] == "apoptosis_healthy"]
pyroptosis_vs_healthy = df[df["group"] == "healthy_pyroptosis"]
pyroptosis_vs_apoptosis = df[df["group"] == "apoptosis_pyroptosis"]

# get thee list of genes that are significant from each comparision
# define the sets
A = set(apoptosis_vs_healthy["features"].tolist())  # apoptosis_vs_healthy_list
B = set(pyroptosis_vs_apoptosis["features"].tolist())  # pyroptosis_vs_apoptosis_list
C = set(pyroptosis_vs_healthy["features"].tolist())  # pyroptosis_vs_healthy_list


# In[7]:


# get the the intersections and union of the genes
U = set(df["features"].tolist())

# get the union of the genes
# Apoptosis vs control
A_int_b_un = U.difference(B.union(C))
# Pyroptosis vs apoptosis
B_int_c_un = U.difference(A.union(C))
# Pyroptosis vs control
C_int_a_un = U.difference(A.union(B))

# get the inersection of each of the groups


print(len(A_int_b_un), len(B_int_c_un), len(C_int_a_un))


# In[8]:


# get the features that are in A_int_b_un
A_int_b_un_df = df[df["features"].isin(A_int_b_un)]
B_int_c_un_df = df[df["features"].isin(B_int_c_un)]
C_int_a_un_df = df[df["features"].isin(C_int_a_un)]

print(A_int_b_un_df.shape, B_int_c_un_df.shape, C_int_a_un_df.shape)
# concat all the dataframes
all_selected_features_df = pd.concat([A_int_b_un_df, B_int_c_un_df, C_int_a_un_df])


# In[9]:


# drop columns from the df
all_selected_features_df = all_selected_features_df.drop(
    columns=[
        "group1",
        "group2",
        "meandiff",
        "lower",
        "upper",
        "reject",
        "p-adj",
        "pos_neg",
    ]
)
all_selected_features_df.rename(columns={"p-adj_abs": "p.adj.value"}, inplace=True)
all_selected_features_df.reset_index(drop=True, inplace=True)
all_selected_features_df.head()


# In[10]:


# drop columns from the df
df = df.drop(
    columns=[
        "group1",
        "group2",
        "meandiff",
        "lower",
        "upper",
        "reject",
        "p-adj",
        "pos_neg",
    ]
)
df.head()


# In[11]:


# set the output file path
output_file_path = pathlib.Path("../results/")
output_file_path.mkdir(exist_ok=True, parents=True)


# In[12]:


# print the table to a csv
df.to_csv("../results/all_features.csv", index=False)


# In[13]:


all_selected_features_df.to_markdown("../results/all_features.md", index=False)


# In[14]:


# print the table to a csv
all_selected_features_df.to_csv("../results/STable_features.csv", index=False)


# In[15]:


# csv to markdown
all_selected_features_df.to_markdown("../results/STable_features.md", index=False)
