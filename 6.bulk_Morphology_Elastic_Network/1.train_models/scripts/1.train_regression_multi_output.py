#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import itertools
import pathlib
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV, LogisticRegression, MultiTaskElasticNetCV

# import RepeatedKFold
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedKFold,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.utils import parallel_backend

# In[ ]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", type=str, default="all")
argparser.add_argument("--shuffle", type=bool, default=False)
argparser.add_argument("--cytokine", type=str, default="cytokine")

args = argparser.parse_args()

cell_type = args.cell_type
cytokine = args.cytokine
shuffle = args.shuffle
shuffle = bool(shuffle == "True")
print(cell_type, shuffle, cytokine)


# In[ ]:


aggregation = True
nomic = True


# In[ ]:


# set shuffle value
if shuffle:
    shuffle = "shuffled_baseline"
else:
    shuffle = "final"


# In[ ]:


MODEL_TYPE = "regression"


# In[ ]:


# load training data from indexes and features dataframe
data_split_path = pathlib.Path(
    f"../../0.split_data/indexes/{cell_type}/regression/aggregated_sc_and_nomic_data_split_indexes.tsv"
)
data_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
)

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pd.read_parquet(data_path)

data_split_indexes = pd.read_csv(data_split_path, sep="\t")


# In[ ]:


# select tht indexes for the training and test set
train_indexes = data_split_indexes.loc[data_split_indexes["label"] == "train"]


# In[ ]:


# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[train_indexes["labeled_data_index"]]


# In[ ]:


# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name
data_y_cols = data_x.filter(regex="NSU").columns
train_y = training_data[data_y_cols]
train_x = data_x.drop(data_y_cols, axis=1)


# In[ ]:


from sklearn.model_selection import LeaveOneOut

loo = LeaveOneOut()
loo.get_n_splits(train_x)
loo.get_n_splits(train_y)


# In[ ]:


train_data_y = train_y[cytokine]
model = ElasticNetCV(
    random_state=0,
    max_iter=100000,
    cv=loo,
    l1_ratio=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.99],
    alphas=[0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
    fit_intercept=True,
    selection="random",
)
# train model on training data on all combinations of model types, feature types, and phenotypic classes

if shuffle == "shuffled_baseline":
    print("Shuffling data")
    for column in train_x:
        np.random.shuffle(train_x[column].values)
else:
    print("Not shuffling data")
# define parameters to search over
with parallel_backend("multiprocessing"):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn")
        # create a logistic regression model
        model.fit(train_x, train_data_y)
        scores = cross_val_score(
            model,
            train_x,
            train_data_y,
            scoring="neg_mean_absolute_error",
            cv=loo,
            n_jobs=-1,
        )
        print(scores)
        print(f"Mean MAE: {scores.mean()}")
        print(f"Std MAE: {scores.std()}")
        print(f"R2: {model.score(train_x, train_data_y)}")

if (aggregation == True) and (nomic == True):
    results_dir = f"../models/regression/{cell_type}/aggregated_with_nomic/"
elif (aggregation == True) and (nomic == False):
    results_dir = f"../models/regression/{cell_type}/aggregated/"
elif (aggregation == False) and (nomic == True):
    results_dir = f"../models/regression/{cell_type}/sc_with_nomic/"
elif (aggregation == False) and (nomic == False):
    results_dir = f"../models/regression/{cell_type}/sc/"
else:
    print("Error")

# create results directory if it doesn't exist
pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

# save final estimator
if shuffle == "shuffled_baseline":
    dump(
        model,
        f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.joblib",
    )
elif shuffle == "final":
    dump(
        model,
        f"{results_dir}/{cytokine}_final__all_nomic.joblib",
    )
else:
    print("Error")

# save condfig copy specific to this model to the folder with the results
# use pathlib
if shuffle == "shuffled_baseline":
    config_copy_path = pathlib.Path(
        f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.toml"
    )
elif shuffle == "final":
    config_copy_path = pathlib.Path(f"{results_dir}/{cytokine}_final__all_nomic.toml")
else:
    print("Error")

# write toml file with parameters used from injected parameters

with open(config_copy_path, "w") as f:
    f.write(f"model_type='{shuffle}'\n")
    f.write(f"aggregation={aggregation}\n")
    f.write(f"nomic={nomic}\n")
    f.write(f"cell_type='{cell_type}'\n")
    f.write(f"feature=all\n")
