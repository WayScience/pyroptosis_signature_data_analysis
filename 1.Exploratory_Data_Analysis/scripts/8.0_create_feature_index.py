#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import toml
from matplotlib import rcParams
from tqdm import tqdm

rcParams.update({"figure.autolayout": True})

# create a venn diagram of the features that are significant in all conditions

warnings.filterwarnings("ignore")
from pycytominer.cyto_utils import infer_cp_features
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# In[ ]:


cell_type = "SHSY5Y"


# In[ ]:


# # parse command line arguments
# parser = argparse.ArgumentParser()
# parser.add_argument("--cell_type", help="cell type to analyze", type=str, default="all")

# args = parser.parse_args()
# cell_type = args.cell_type


# In[4]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet")
df = pd.read_parquet(file_path)

# index output for features
output_file = pathlib.Path(f"../features/{cell_type}_feature_index.txt")
# create output directory if it doesn't exist
output_file.parent.mkdir(parents=True, exist_ok=True)


# In[5]:


df_metadata = df.filter(regex="Metadata")
df_data = df.drop(df_metadata.columns, axis=1)
df_data["Metadata_number_of_singlecells"] = df_metadata[
    "Metadata_number_of_singlecells"
]
cp_features = infer_cp_features(df)

# write each feature to a file
with open(output_file, "w") as f:
    for item in cp_features:
        f.write("%s\n" % item)
