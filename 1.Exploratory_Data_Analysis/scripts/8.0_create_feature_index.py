#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import pathlib

import pyarrow.parquet as pq

# In[ ]:


# parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--cell_type", help="cell type to analyze", type=str, default="all")

args = parser.parse_args()
cell_type = args.cell_type


# In[3]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet")
columns_list = pq.read_schema(file_path).names
columns_list = [x for x in columns_list if "Metadata" not in x]
columns_list = [x for x in columns_list if "__index_level_0__" not in x]


# In[4]:


# index output for features
output_file = pathlib.Path(f"../features/{cell_type}_feature_index.txt")
# create output directory if it doesn't exist
output_file.parent.mkdir(parents=True, exist_ok=True)


# In[5]:


# write each feature to a file
with open(output_file, "w") as f:
    for item in columns_list:
        f.write("%s\n" % item)
