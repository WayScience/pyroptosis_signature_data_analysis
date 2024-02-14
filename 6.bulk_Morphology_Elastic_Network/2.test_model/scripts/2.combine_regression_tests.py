#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


cell_type = "PBMC"


# In[3]:


# define paths to the results
results_dir_path = pathlib.Path(
    f"../results/regression/{cell_type}/aggregated_with_nomic/"
).resolve(strict=True)

model_stats_final_output_path = pathlib.Path(
    f"../results/regression/{cell_type}/combined/model_stats.csv"
)
# make the directory if it does not exist
model_stats_final_output_path.parent.mkdir(parents=True, exist_ok=True)

variance_r2_stats_final_output_path = pathlib.Path(
    f"../results/regression/{cell_type}/combined/variance_r2_stats.csv"
)

# get a list of all the files that contain "model_stats" in the name
model_stats_files = list(results_dir_path.glob("*model_stats*"))


# get a list of all the files that contain "variance_r2_stats" in the name
variance_r2_stats_files = list(results_dir_path.glob("*variance_r2_stats*"))


# In[4]:


# # concate all the model_stats files
# model_stats_df = pd.concat([pd.read_csv(f) for f in model_stats_files])

# # concate all the variance_r2_stats files
# variance_r2_stats_df = pd.concat([pd.read_csv(f) for f in variance_r2_stats_files])


# In[5]:


# # save the final output to csv
# model_stats_df.to_csv(model_stats_final_output_path, index=False)
# variance_r2_stats_df.to_csv(variance_r2_stats_final_output_path, index=False)


# In[6]:


list_of_dfs = []
for f in model_stats_files:
    df = pd.read_csv(f)
    f = str(f.stem)
    f = (
        f.split("[NSU]_")[1]
        .split("_model_stats")[0]
        .strip("final_")
        .strip("shuffled_baseline_")
    )
    df["dataset"] = f
    list_of_dfs.append(df)

model_stats_df = pd.concat(list_of_dfs)
print(model_stats_df.shape)
# save the final output to csv
model_stats_df.to_csv(model_stats_final_output_path, index=False)


# In[7]:


list_of_dfs = []
for f in variance_r2_stats_files:
    df = pd.read_csv(f)
    f = str(f.stem)
    f = (
        f.split("[NSU]_")[1]
        .split("_variance_r2_stats")[0]
        .strip("final_")
        .strip("shuffled_baseline_")
    )
    df["dataset"] = f
    list_of_dfs.append(df)

variance_r2_stats_df = pd.concat(list_of_dfs)
print(variance_r2_stats_df.shape)
# save the final output to csv
variance_r2_stats_df.to_csv(variance_r2_stats_final_output_path, index=False)

variance_r2_stats_df.head()
