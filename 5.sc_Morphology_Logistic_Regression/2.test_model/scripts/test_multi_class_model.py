#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
from joblib import dump, load
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
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

toml_path = pathlib.Path("../1.train_models/multi_class_config.toml")
with open(toml_path, "r") as f:
    config = toml.load(f)
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
training_data


# In[6]:


# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
labeled_data


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


# set model path from parameters
model_path = pathlib.Path(f"./models/multi_class/{treatments_str}")


# In[9]:


# set path for figures
figure_path = pathlib.Path(f"./figures/multi_class/{treatments_str}")
figure_path.mkdir(parents=True, exist_ok=True)


# In[10]:


model_types = ["final", "shuffled_baseline"]
feature_types = ["CP"]

# test model on testing data
for model_type, feature_type in itertools.product(model_types, feature_types):
    print(model_type, feature_type)
    # load model
    print(f"../1.train_models/{model_path}/{model_type}__{feature_type}.joblib")
    print(
        "../1.train_models/models/multi_class/DMSO_0.100_DMSO_0.025__LPS_1.000_DMSO_0.025__LPS_100.000_DMSO_0.025/final__CP.joblib"
    )
    model = load(f"../1.train_models/{model_path}/{model_type}__{feature_type}.joblib")

    # get predictions
    predictions = model.predict(data_x)

    # get probabilities

    probabilities = model.predict_proba(data_x)

    # get accuracy
    accuracy = accuracy_score(labeled_data, predictions)

    # get f1 score
    f1 = f1_score(labeled_data, predictions, average="weighted")

    # plot confusion matrix heatmap
    sns.heatmap(confusion_matrix(labeled_data, predictions), annot=True, fmt="g")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title(
        f"Confusion matrix of {model_type} model on {feature_type} features for all phenotypic classes"
    )
    plt.savefig(f"{figure_path}/{model_type}__{feature_type}__confusion_matrix.png")
    plt.show()


# In[ ]:
