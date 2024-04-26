#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import pathlib

import pandas as pd

# In[2]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", type=str, default="cell_type")

args = argparser.parse_args()

cell_type = args.cell_type


# In[ ]:


results_dir_path = pathlib.Path(
    f"../results/regression/{cell_type}/aggregated_with_nomic/"
).resolve(strict=True)

model_stats_final_output_path = pathlib.Path(
    f"../results/regression/{cell_type}/aggregated_with_nomic/model_stats.csv"
)

variance_r2_stats_final_output_path = pathlib.Path(
    f"../results/regression/{cell_type}/aggregated_with_nomic/variance_r2_stats.csv"
)

# get a list of all the files that contain "model_stats" in the name
model_stats_files = list(results_dir_path.glob("*model_stats*"))


# get a list of all the files that contain "variance_r2_stats" in the name
variance_r2_stats_files = list(results_dir_path.glob("*variance_r2_stats*"))


# concate all the model_stats files
model_stats_df = pd.concat([pd.read_csv(f) for f in model_stats_files])

# concate all the variance_r2_stats files
variance_r2_stats_df = pd.concat([pd.read_csv(f) for f in variance_r2_stats_files])

# save the final output to csv
model_stats_df.to_csv(model_stats_final_output_path, index=False)
variance_r2_stats_df.to_csv(variance_r2_stats_final_output_path, index=False)


# In[ ]:


print("Completed!")
