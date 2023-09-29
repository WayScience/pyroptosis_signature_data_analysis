#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import pandas as pd

# In[2]:


cell_type = "SHSY5Y"


# In[3]:


# set import paths
model_path = pathlib.Path(
    f"../../1.train_models/models/regression/{cell_type}/aggregated_with_nomic/IL-1 beta [NSU]_final__all_nomic.joblib"
)

# set and create output paths
output_path = pathlib.Path(f"../results/regression/{cell_type}/")
output_path.mkdir(parents=True, exist_ok=True)


# In[4]:


# set path of models
model_path = pathlib.Path(
    f"../../1.train_models/models/regression/{cell_type}/aggregated_with_nomic/"
)


# In[5]:


for model_file in list(model_path.glob("*.joblib")):
    # import model
    model = joblib.load(model_file)
    # create a df with the coefficients
    df_coefficients = pd.DataFrame(
        model.coef_, index=model.feature_names_in_, columns=["coefficients"]
    )
    # sort by absolute value of coefficients
    df_coefficients = df_coefficients.reindex(
        df_coefficients["coefficients"].abs().sort_values(ascending=False).index
    )
    df_coefficients["secreted_proteins"] = "IL-1 beta [NSU]"
    df_coefficients["shuffle"] = "final"
    df_coefficients["cell_type"] = cell_type
    df_coefficients["alpha"] = model.alpha_
    df_coefficients["l1_ratio"] = model.l1_ratio_
    df_coefficients = df_coefficients.reset_index(inplace=True)
    df_coefficients = df_coefficients.rename(
        columns={"index": "feature_names"}, inplace=True
    )

    # get the file name of the model
    # replace the joblib ending with csv
    file_name = pathlib.Path(model_file).name.replace(".joblib", ".csv")
    # write the output path
    df_coefficients.to_csv(output_path / file_name, index=False)
