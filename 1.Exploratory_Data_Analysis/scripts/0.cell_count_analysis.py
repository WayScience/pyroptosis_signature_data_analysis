#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt

get_ipython().run_line_magic("matplotlib", "inline")
import pathlib

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns

# In[3]:


# Parameters
celltype = "PBMC"


# In[4]:


# Set Path
path = pathlib.Path(f"../data/{celltype}_preprocessed_sc_norm.parquet")


# In[5]:


# import data
df = pq.read_table(path).to_pandas()


# In[6]:


# barplot of number of single cells per treatment
sns.barplot(
    x="Metadata_Treatment",
    y="Metadata_number_of_singlecells",
    data=df,
    estimator=np.mean,
    errorbar=("sd"),
)

plt.ylabel("Number of single cells")
plt.xlabel("Treatment")
plt.xticks(rotation=45)
plt.title("Number of single cells per treatment")

# if path does not exist, create it
pathlib.Path(f"Figures/cell_counts_plate2/{celltype}").mkdir(
    parents=True, exist_ok=True
)
# save figure
plt.savefig(
    f"Figures/cell_counts_plate2/{celltype}/Number_of_single_cells_per_treatment.png",
    bbox_inches="tight",
)
plt.show()
plt.close()


# In[7]:


# if path does not exist, create it
pathlib.Path(f"Figures/cell_counts_plate2/{celltype}").mkdir(
    parents=True, exist_ok=True
)

# Number of single cells per treatment and dose level
for i in df["Metadata_Treatment"].unique():
    tmp_df = df[df["Metadata_Treatment"] == i]
    sns.barplot(
        x="Metadata_Dose",
        y="Metadata_number_of_singlecells",
        # hue="Metadata_Treatment",
        estimator=np.median,
        data=tmp_df,
        errorbar=("sd"),
    )
    plt.xlabel(f"{i}_doasage")
    plt.ylabel("Number of single cells")
    plt.xticks(rotation=45)
    plt.title(f"Number of single cells per {i}")

    plt.savefig(
        f"Figures/cell_counts_plate2/{celltype}/Number_of_single_cells_per_{i}.png",
        bbox_inches="tight",
    )
    plt.show()
    plt.close()
