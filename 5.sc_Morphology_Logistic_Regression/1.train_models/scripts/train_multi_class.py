#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.utils import parallel_backend, shuffle

# In[2]:


# load training data from indexes and features dataframe
data_split_path = pathlib.Path(f"../0.split_data/indexes/data_split_indexes.tsv")
data_path = pathlib.Path("../../data/SHSY5Y_preprocessed_sc_norm.parquet")
# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pq.read_table(data_path).to_pandas()
data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)

# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[data_split_indexes["labeled_data_index"]]


# In[3]:


# read in toml file and get parameters
import toml

toml_path = pathlib.Path("multi_class_config.toml")
with open(toml_path, "r") as f:
    config = toml.load(f)
f.close()
control = config["logistic_regression_params"]["control"]
treatments = config["logistic_regression_params"]["treatments"]


# In[4]:


# get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
# set tmp_lst tp select on
tmp_lst = [control] + treatments
training_data = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(tmp_lst)
]


# In[5]:


# at random downsample the DMSO treatment to match the number of wells in the LPS treatment
# get the number of wells in the LPS treatment
trt_wells = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatments[0]
].shape[0]
# get the number of wells in the DMSO treatment
dmso_wells = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
].shape[0]
# downsample the DMSO treatment to match the number of wells in the LPS treatment
dmso_holdout = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
].sample(n=trt_wells)
# remove the downsampled DMSO wells from the data
training_data = training_data.drop(dmso_holdout.index)


# In[6]:


# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]


# In[7]:


# define the model save name
treatments_str = ""
for i in tmp_lst:
    # if not last in list
    if i != tmp_lst[-1]:
        treatments_str += i + "__"
    else:
        treatments_str += i
treatments_str


# In[8]:


# specify model types, feature types, and phenotypic classes
model_types = ["final", "shuffled_baseline"]
feature_types = ["CP"]
# create stratified data sets for k-fold cross validation
straified_k_folds = StratifiedKFold(n_splits=2, shuffle=False)

# create logistic regression model with following parameters
log_reg_model = LogisticRegression(
    penalty="elasticnet",
    solver="saga",
    max_iter=10,
    n_jobs=-1,
    random_state=0,
    class_weight="balanced",
)

# specify parameters to tune for
parameters = {"C": np.logspace(-3, 3, 7), "l1_ratio": np.linspace(0, 1, 11)}
print(f"Parameters being tested during grid search: {parameters}\n")

# create grid search with cross validation with hypertuning params
grid_search_cv = GridSearchCV(
    log_reg_model, parameters, cv=straified_k_folds, n_jobs=-1, scoring="f1_weighted"
)

# train model on training data on all combinations of model types, feature types, and phenotypic classes
for model_type, feature_type in itertools.product(model_types, feature_types):

    if model_type == "shuffled_baseline":
        for column in data_x:
            np.random.shuffle(data_x[column].values)

    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )
            grid_search_cv = grid_search_cv.fit(data_x, labeled_data)

    # print info for best estimator
    print(f"Best parameters: {grid_search_cv.best_params_}")
    print(f"Score of best estimator: {grid_search_cv.best_score_}\n")

    results_dir = f"./models/multi_class/{treatments_str}"
    # create results directory if it doesn't exist
    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
    # save final estimator
    dump(
        grid_search_cv.best_estimator_,
        f"{results_dir}/{model_type}__{feature_type}.joblib",
    )


# In[9]:


# save condfig copy specific to this model to the folder with the results
# use pathlib
config_copy_path = pathlib.Path(f"{results_dir}/{model_type}__{feature_type}.toml")
with open(config_copy_path, "w") as f:
    toml.dump(config, f)
f.close()
