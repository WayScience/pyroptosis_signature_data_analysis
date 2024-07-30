#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


# set path to data
valid_data_file_path = pathlib.Path(
    "./raw/nELISA - Luminex comp in LPS stimulated PBMCs.xlsx"
).resolve()

# get sheet names
xls = pd.ExcelFile(valid_data_file_path)
sheet_names = xls.sheet_names
print(sheet_names)
nELISA = pd.read_excel(valid_data_file_path, sheet_name="nELISA_pgmL")
Luminex = pd.read_excel(valid_data_file_path, sheet_name="xMAP_pgmL")


# In[3]:


print(nELISA.columns)
print(Luminex.columns)
Luminex.rename(columns={"IFNgamma": "IFN gamma", "TNFalpha": "TNF alpha"}, inplace=True)


# In[4]:


# add a column to each dataframe to indicate the type of data
nELISA["data_type"] = "nELISA"
Luminex["data_type"] = "Luminex"

# combine the two dataframes
validation_data = pd.concat([nELISA, Luminex])
print(validation_data.shape)
# rename the columns to remove spaces
validation_data.columns = [col.replace(" ", "_") for col in validation_data.columns]
# cas
validation_data.head()


# In[5]:


# min - max scaling
validation_data_min_max = validation_data.copy()

for col in validation_data_min_max.columns:
    if col not in ["LPS_concentration", "data_type"]:
        validation_data_min_max[col] = (
            validation_data_min_max[col] - validation_data_min_max[col].min()
        ) / (validation_data_min_max[col].max() - validation_data_min_max[col].min())


# In[6]:


# convert both dataframes to long format
validation_data_long = pd.melt(
    validation_data,
    id_vars=["LPS_concentration", "data_type"],
    var_name="cytokine",
    value_name="concentration",
)
validation_data_min_max_long = pd.melt(
    validation_data_min_max,
    id_vars=["LPS_concentration", "data_type"],
    var_name="cytokine",
    value_name="concentration",
)


# In[7]:


# write the data to a new file
output_file_path = pathlib.Path("./clean/validation/").resolve()
output_file_path.mkdir(parents=True, exist_ok=True)
output_file_path = pathlib.Path(
    "./clean/validation/nELISA_luminex_validation_data.parquet"
).resolve()
validation_data_long.to_parquet(output_file_path, index=False)

output_file_path = pathlib.Path("./clean/validation/").resolve()
output_file_path.mkdir(parents=True, exist_ok=True)
output_file_path = pathlib.Path(
    "./clean/validation/nELISA_luminex_validation_data_min_max.parquet"
).resolve()
validation_data_min_max_long.to_parquet(output_file_path, index=False)
