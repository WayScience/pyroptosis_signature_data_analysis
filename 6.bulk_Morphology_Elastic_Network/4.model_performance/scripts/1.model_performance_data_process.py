#!/usr/bin/env python
# coding: utf-8

# ## This notebook pre-processes the results of the model performance evaluation to be used in the model performance visualization notebook.

# In[1]:


import pathlib

import numpy as np
import pandas as pd

# In[2]:


cell_type = "PBMC"


# In[3]:


# set path of models
model_path = pathlib.Path(f"../../3.model_coefficients/results/regression/{cell_type}/")
# output all models path
output_path = pathlib.Path(f"../results/regression/{cell_type}/")
output_path.mkdir(parents=True, exist_ok=True)


# In[4]:


# declare a blank dataframe to store all the model performance
all_model_performance = pd.DataFrame()


# In[5]:


for model_file in model_path.glob("*.csv"):
    all_model_performance = pd.concat(
        # concat the model performance
        [all_model_performance, pd.read_csv(model_file)]
    )


# In[6]:


# all_model_performance write to csv
all_model_performance.to_csv(
    output_path / f"all_model_performance.csv",
    index=False,
)
