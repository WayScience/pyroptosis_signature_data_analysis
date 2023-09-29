#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import joblib
import numpy as np
import pandas as pd

# In[2]:


cell_type = "PBMC"


# In[3]:


# set and create output paths
output_path = pathlib.Path(f"../results/regression/{cell_type}/")
output_path.mkdir(parents=True, exist_ok=True)


# In[4]:


# set path of models
model_path = pathlib.Path(
    f"../../1.train_models/models/regression/{cell_type}/aggregated_with_nomic/"
)
model_stats_path = pathlib.Path(
    f"../../2.test_model/results/regression/{cell_type}/aggregated_with_nomic/"
)


# In[15]:


for model_file in list(model_path.glob("*.joblib")):
    # get the basename of the model to load the test results too
    basename = pathlib.Path(model_file).name
    basename = basename.split("__")[0] + "_variance_r2_stats.csv"
    test_results_file = pathlib.Path(model_stats_path / basename)
    # load the test results
    df_test_results = pd.read_csv(test_results_file)
    # get the r2 score
    r2 = df_test_results.loc[(df_test_results["data_split"] == "test_data")][
        "r2"
    ].unique()[0]
    r2 = np.float64(r2.strip("[]"))
    assert r2.dtype == np.float64
    # import model
    model = joblib.load(model_file)
    # create a df with the coefficients
    df_coefficients = pd.DataFrame(
        model.coef_, index=model.feature_names_in_, columns=["coefficients"]
    )

    # print(basename.split("_")[0], "__",basename.split("_")[1])

    # sort by absolute value of coefficients
    df_coefficients = df_coefficients.reindex(
        df_coefficients["coefficients"].abs().sort_values(ascending=False).index
    )
    df_coefficients["secreted_proteins"] = basename.split("_")[0]
    if basename.split("_")[1] != "shuffled" and basename.split("_")[1] != "final":
        df_coefficients["shuffle"] = basename.split("_")[2]
    else:
        df_coefficients["shuffle"] = basename.split("_")[1]
    df_coefficients["cell_type"] = cell_type
    df_coefficients["alpha"] = model.alpha_
    df_coefficients["l1_ratio"] = model.l1_ratio_
    df_coefficients = df_coefficients.reset_index()
    df_coefficients = df_coefficients.rename(
        columns={"index": "feature_names"},
    )
    df_coefficients["r2"] = r2
    # get the file name of the model
    # replace the joblib ending with csv
    file_name = pathlib.Path(model_file).name.replace(".joblib", ".csv")
    # write the output path
    df_coefficients.to_csv(output_path / file_name, index=False)
