#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ast
import itertools
import pathlib
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, train_test_split
from sklearn.utils import parallel_backend, shuffle

# In[2]:


# Parameters
cell_type = "SHSY5Y"
aggregation = False
nomic = False
flag = True
control = "DMSO_0.100_DMSO_0.025"
treatment = "Thapsigargin_1.000_DMSO_0.025"


# In[3]:


if flag == False:
    # read in toml file and get parameters
    toml_path = pathlib.Path("single_class_config.toml")
    with open(toml_path, "r") as f:
        config = toml.load(f)
    f.close()
    control = config["logistic_regression_params"]["control"]
    treatment = config["logistic_regression_params"]["treatments"]
    aggregation = ast.literal_eval(config["logistic_regression_params"]["aggregation"])
    nomic = ast.literal_eval(config["logistic_regression_params"]["nomic"])
    cell_type = config["logistic_regression_params"]["cell_type"]


# In[4]:


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


# In[5]:


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


# In[6]:


# select tht indexes for the training and test set
train_indexes = data_split_indexes.loc[data_split_indexes["label"] == "train"]


# In[7]:


# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[train_indexes["labeled_data_index"]]


# In[8]:


# get oneb_Metadata_Treatment_Dose_Inhibitor_Dose  =='DMSO_0.100_DMSO_0.025' and 'LPS_100.000_DMSO_0.025 and Thapsigargin_10.000_DMSO_0.025'
training_data = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        [control, treatment]
    )
]


# In[9]:


# at random downsample the DMSO treatment to match the number of wells in the LPS treatment
seed = 0
# get the number of wells in the LPS treatment
trt_wells = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatment
].shape[0]
# get the number of wells in the DMSO treatment
dmso_wells = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
].shape[0]
# downsample the DMSO treatment to match the number of wells in the LPS treatment
dmso_holdout = training_data[
    training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == control
].sample(n=trt_wells, random_state=seed)
# remove the downsampled DMSO wells from the data
training_data = training_data.drop(dmso_holdout.index)


# In[10]:


# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]


# In[11]:


# specify model types, feature types, and phenotypic classes
model_types = ["final", "shuffled_baseline"]
feature_types = ["CP"]
phenotypic_classes = [f"{treatment}"]
# create stratified data sets forLPS k-fold cross validation
straified_k_folds = StratifiedKFold(n_splits=3, shuffle=False)

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
for model_type, feature_type, phenotypic_class in itertools.product(
    model_types, feature_types, phenotypic_classes
):
    phenotypic_class_counts = training_data.loc[
        training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == phenotypic_class
    ].shape[0]
    print(
        f"Training {model_type} model on {feature_type} features for {phenotypic_class} with {phenotypic_class_counts} samples"
    )

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

    if (aggregation == True) and (nomic == True):
        results_dir = f"./models/single_class/{cell_type}/aggregated_with_nomic/{control}__{treatment}"
    elif (aggregation == True) and (nomic == False):
        results_dir = (
            f"./models/single_class/{cell_type}/aggregated/{control}__{treatment}"
        )
    elif (aggregation == False) and (nomic == True):
        results_dir = (
            f"./models/single_class/{cell_type}/sc_with_nomic/{control}__{treatment}"
        )
    elif (aggregation == False) and (nomic == False):
        results_dir = f"./models/single_class/{cell_type}/sc/{control}__{treatment}"
    else:
        print("Error")

    # create results directory if it doesn't exist
    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

    # save final estimator
    dump(
        grid_search_cv.best_estimator_,
        f"{results_dir}/{model_type}__{feature_type}.joblib",
    )


# In[12]:


# save condfig copy specific to this model to the folder with the results
# use pathlib
config_copy_path = pathlib.Path(f"{results_dir}/{model_type}__{feature_type}.toml")
# write toml file with parameters used from injected parameters
with open(config_copy_path, "a") as f:
    f.write(f"model_type='{model_type}'\n")
    f.write(f"control='{control}'\n")
    f.write(f"treatments='{treatment}'\n")
    f.write(f"aggregation={aggregation}\n")
    f.write(f"nomic={nomic}\n")
    f.write(f"cell_type='{cell_type}'\n")
f.close()
