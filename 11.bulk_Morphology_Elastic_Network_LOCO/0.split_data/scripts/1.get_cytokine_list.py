#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


# only need one since all the files are in the same directory
cell_type = "PBMC"


# In[3]:


nomic_df_path = pathlib.Path(
    f"../../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_cleanup4correlation.csv"
)
df_nomic = pd.read_csv(nomic_df_path)

# drop columns that contain [pgML]
df_nomic = df_nomic.drop(columns=[col for col in df_nomic.columns if "[pgML]" in col])


# In[4]:


# drop columns that contain "metadata"
df_nomic = df_nomic.drop(
    columns=[col for col in df_nomic.columns if "Metadata" in col]
).columns.to_list()


# In[5]:


# write the list to a file
# set the path to save the file
cytokine_list_path = pathlib.Path(f"../cytokine_list/cytokine_list.txt")
# if the directory does not exist, create it
cytokine_list_path.parent.mkdir(parents=True, exist_ok=True)

# write the list to a file
with open(cytokine_list_path, "w") as f:
    for item in df_nomic:
        f.write(f"{item}\n")
