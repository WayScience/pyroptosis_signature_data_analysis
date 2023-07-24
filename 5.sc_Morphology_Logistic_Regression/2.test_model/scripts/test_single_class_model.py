#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ast
import itertools
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns
import toml
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


# Parameters
cell_type = "PBMC"
aggregation = False
nomic = False
flag = True
control = "DMSO_0.100_DMSO_0.025"
treatment = "Thapsigargin_1.000_DMSO_0.025"


# In[3]:


if flag == False:
    # read in toml file and get parameters
    toml_path = pathlib.Path("../1.train_models/single_class_config.toml")
    with open(toml_path, "r") as f:
        config = toml.load(f)
    control = config["logistic_regression_params"]["control"]
    treatment = config["logistic_regression_params"]["treatments"]
    aggregation = ast.literal_eval(config["logistic_regression_params"]["aggregation"])
    nomic = ast.literal_eval(config["logistic_regression_params"]["nomic"])
    cell_type = config["logistic_regression_params"]["cell_type"]
    print(aggregation, nomic, cell_type)


# In[4]:


if flag == False:
    # read in toml file and get parameters
    toml_path = pathlib.Path("../1.train_models/single_class_config.toml")
    with open(toml_path, "r") as f:
        config = toml.load(f)
    f.close()
    control = config["logistic_regression_params"]["control"]
    treatment = config["logistic_regression_params"]["treatments"]
    aggregation = ast.literal_eval(config["logistic_regression_params"]["aggregation"])
    nomic = ast.literal_eval(config["logistic_regression_params"]["nomic"])
    cell_type = config["logistic_regression_params"]["cell_type"]
    print(aggregation, nomic, cell_type)


# In[5]:


# load training data from indexes and features dataframe
# data_split_path = pathlib.Path(f"../0.split_data/indexes/data_split_indexes.tsv")
data_path = pathlib.Path(f"../../data/{cell_type}_preprocessed_sc_norm.parquet")

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pq.read_table(data_path).to_pandas()

# import nomic data
nomic_df_path = pathlib.Path(
    f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)
df_nomic = pd.read_csv(nomic_df_path)

# clean up nomic data
df_nomic = df_nomic.drop(columns=[col for col in df_nomic.columns if "[pgML]" in col])
# drop first 25 columns (Metadata that is not needed)
df_nomic = df_nomic.drop(columns=df_nomic.columns[3:25])
df_nomic = df_nomic.drop(columns=df_nomic.columns[0:2])


# In[6]:


if (aggregation == True) and (nomic == True):

    data_split_path = pathlib.Path(
        f"../0.split_data/indexes/{cell_type}/{control}_{treatment}/aggregated_sc_and_nomic_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
    # subset each column that contains metadata
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df, metadata, left_on="Metadata_Well", right_on="Metadata_Well"
    )
    data_df = pd.merge(
        data_df, df_nomic, left_on="Metadata_Well", right_on="position_x"
    )
    data_df = data_df.drop(columns=["position_x"])

elif (aggregation == True) and (nomic == False):
    data_split_path = pathlib.Path(
        f"../0.split_data/indexes/{cell_type}/{control}_{treatment}/aggregated_sc_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
    # subset each column that contains metadata
    metadata = data_df.filter(regex="Metadata")
    data_df = data_df.drop(metadata.columns, axis=1)
    data_df = pd.concat([data_df, metadata["Metadata_Well"]], axis=1)
    # groupby well and take mean of each well
    data_df = data_df.groupby("Metadata_Well").mean()
    # drop duplicate rows in the metadata_well column
    metadata = metadata.drop_duplicates(subset=["Metadata_Well"])
    # get the metadata for each well
    data_df = pd.merge(
        data_df, metadata, left_on="Metadata_Well", right_on="Metadata_Well"
    )
elif (aggregation == False) and (nomic == True):
    data_split_path = pathlib.Path(
        f"../0.split_data/indexes/{cell_type}/{control}_{treatment}/sc_and_nomic_data_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
    data_df = pd.merge(
        data_df, df_nomic, left_on="Metadata_Well", right_on="position_x"
    )
    data_df = data_df.drop(columns=["position_x"])
elif aggregation == False and nomic == False:
    data_split_path = pathlib.Path(
        f"../0.split_data/indexes/{cell_type}/{control}_{treatment}/sc_split_indexes.tsv"
    )
    data_split_indexes = pd.read_csv(data_split_path, sep="\t", index_col=0)
else:
    print("Error")
data_df


# In[7]:


data_split_indexes.index = data_split_indexes["labeled_data_index"]


# In[8]:


# subset data_df by indexes in data_split_indexes
data_all = data_df.loc[data_split_indexes["labeled_data_index"]]
data_all["label"] = data_split_indexes["label"]


# In[9]:


# get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
data_all = data_all[
    data_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin([control, treatment])
]


# In[10]:


# at random downsample the DMSO treatment to match the number of wells in the LPS treatment
seed = 0
# get the number of wells in the LPS treatment
trt_wells = data_all[
    data_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatment
].shape[0]
# get the number of wells in the DMSO treatment
dmso_wells = data_all[
    data_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
].shape[0]
if dmso_wells > trt_wells:
    # downsample the DMSO treatment to match the number of wells in the LPS treatment
    dmso_holdout = data_all[
        data_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
    ].sample(n=trt_wells, random_state=seed)
    # remove the downsampled DMSO wells from the data
    data_all = data_all


# In[11]:


# set model path from parameters
if (aggregation == True) and (nomic == True):
    model_path = pathlib.Path(
        f"models/single_class/{cell_type}/aggregated_with_nomic/{control}__{treatment}"
    )
elif (aggregation == True) and (nomic == False):
    model_path = pathlib.Path(
        f"models/single_class/{cell_type}/aggregated/{control}__{treatment}"
    )
elif (aggregation == False) and (nomic == True):
    model_path = pathlib.Path(
        f"models/single_class/{cell_type}/sc_with_nomic/{control}__{treatment}"
    )
elif (aggregation == False) and (nomic == False):
    model_path = pathlib.Path(
        f"models/single_class/{cell_type}/sc/{control}__{treatment}"
    )
else:
    print("Error")


# In[12]:


model_types = ["final", "shuffled_baseline"]
feature_types = ["CP"]
phenotypic_classes = [treatment]


# In[13]:


# define metadata columns
# subset each column that contains metadata
metadata = data_all.filter(regex="Metadata")
# drop all metadata columns
data_x = data_all.drop(metadata.columns, axis=1)
labeled_data = data_all["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
evaluation_types = ["train", "test"]


train_labeled_data = data_all.loc[data_all["label"] == "train"][
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]
test_labeled_data = data_all.loc[data_all["label"] == "test"][
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
]


# In[14]:


# set path for figures
if (aggregation == True) and (nomic == True):
    figure_path = pathlib.Path(
        f"./figures/single_class/{cell_type}/aggregated_with_nomic/{control}__{treatment}"
    )
    results_path = pathlib.Path(
        f"./results/single_class/{cell_type}/aggregated_with_nomic/{control}__{treatment}"
    )
elif (aggregation == True) and (nomic == False):
    figure_path = pathlib.Path(
        f"./figures/single_class/{cell_type}/aggregated/{control}__{treatment}"
    )
    results_path = pathlib.Path(
        f"./results/single_class/{cell_type}/aggregated/{control}__{treatment}"
    )
elif (aggregation == False) and (nomic == True):
    figure_path = pathlib.Path(
        f"./figures/single_class/{cell_type}/sc_with_nomic/{control}__{treatment}"
    )
    results_path = pathlib.Path(
        f"./results/single_class/{cell_type}/sc_with_nomic/{control}__{treatment}"
    )
elif (aggregation == False) and (nomic == False):
    figure_path = pathlib.Path(
        f"./figures/single_class/{cell_type}/sc/{control}__{treatment}"
    )
    results_path = pathlib.Path(
        f"./results/single_class/{cell_type}/sc/{control}__{treatment}"
    )
else:
    print("Error")
figure_path.mkdir(parents=True, exist_ok=True)
results_path.mkdir(parents=True, exist_ok=True)


# In[15]:


data_x.reset_index(drop=True, inplace=True)
# create empty dataframe to store predictions
compiled_predictions = pd.DataFrame(
    columns=[
        "Phenotypic_Class_Predicted",
        "Phenotypic_Class_True",
        "data_split",
        "shuffled",
        "feature_type",
    ],
)

# test model on testing data
for model_type, feature_type, phenotypic_class, evaluation_type in itertools.product(
    model_types, feature_types, phenotypic_classes, evaluation_types
):
    print(model_type, feature_type, phenotypic_class, evaluation_type)
    # load model
    model = load(f"../1.train_models/{model_path}/{model_type}__{feature_type}.joblib")
    print(model)

    if evaluation_type == "train":
        # get row that are labeled train in label column
        train_data_x = data_x.loc[data_x["label"] == "train"]
        train_data_x = train_data_x.drop("label", axis=1)

        predictions = model.predict(train_data_x)
        # get probabilities
        probabilities = model.predict_proba(train_data_x)
        # get accuracy
        accuracy = accuracy_score(train_labeled_data, predictions)
        f1 = f1_score(train_labeled_data, predictions, average="weighted")
        train_predictions_df = pd.DataFrame(
            {
                "Phenotypic_Class_Predicted": predictions,
                "Phenotypic_Class_True": train_labeled_data,
                "data_split": evaluation_type,
                "shuffled": "shuffled" in model_type,
                "feature_type": feature_type,
            }
        )

        compiled_predictions = pd.concat(
            [compiled_predictions, train_predictions_df], axis=0, ignore_index=True
        )
    elif evaluation_type == "test":
        # get row that are labeled test in label column
        test_data_x = data_x.loc[data_x["label"] == "test"]
        test_data_x = test_data_x.drop("label", axis=1)
        predictions = model.predict(test_data_x)
        # get probabilities
        probabilities = model.predict_proba(test_data_x)
        # get accuracy
        accuracy = accuracy_score(test_labeled_data, predictions)
        # get f1 score
        f1 = f1_score(test_labeled_data, predictions, average="weighted")
        test_predictions_df = pd.DataFrame(
            {
                "Phenotypic_Class_Predicted": predictions,
                "Phenotypic_Class_True": test_labeled_data,
                "data_split": evaluation_type,
                "shuffled": "shuffled" in model_type,
                "feature_type": feature_type,
            }
        )
        compiled_predictions = pd.concat(
            [compiled_predictions, test_predictions_df], axis=0, ignore_index=True
        )


# In[16]:


# write compiled predictions to csv file in results folder
compiled_predictions.to_csv(f"{results_path}/compiled_predictions.csv", index=False)


# In[17]:


compiled_predictions
