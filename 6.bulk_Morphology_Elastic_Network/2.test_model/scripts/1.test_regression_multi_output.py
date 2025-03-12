#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import ast
import pathlib

import joblib
import numpy as np
import pandas as pd

# import mse
from sklearn.metrics import explained_variance_score, mean_squared_error, r2_score
from sklearn.utils import parallel_backend

# In[2]:


argparser = argparse.ArgumentParser()
argparser.add_argument("--cell_type", default="all")
argparser.add_argument("--shuffle", default="False")
argparser.add_argument("--cytokine", default="all")

args = argparser.parse_args()

cell_type = args.cell_type
shuffle = ast.literal_eval(args.shuffle)
cytokine = args.cytokine
print(cell_type, shuffle, cytokine)

# cell_type = "PBMC"
# shuffle = False
# cytokine = "XCL1 (Lymphotactin) [NSU]"


# In[3]:


# Parameters
aggregation = True
nomic = True


# In[4]:


# set shuffle value
if shuffle == True:
    shuffle = "shuffled_baseline"
else:
    shuffle = "final"


# In[5]:


MODEL_TYPE = "regression"


# In[6]:


# load training data from indexes and features dataframe
data_split_path = pathlib.Path(
    f"../../0.split_data/indexes/{cell_type}/regression/aggregated_sc_and_nomic_data_split_indexes.tsv"
)
data_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
)

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pd.read_parquet(data_path)

data_split_indexes = pd.read_csv(data_split_path, sep="\t")


# In[7]:


# rename column that contain the treatment dose to be a metadata column
data_df.rename(
    columns={
        "oneb_Treatment_Dose_Inhibitor_Dose": "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    },
    inplace=True,
)


# In[8]:


# remove duplicate columns
data_df = data_df.loc[:, ~data_df.columns.duplicated()]


# In[9]:


# select tht indexes for the training and test set
train_indexes = data_split_indexes.loc[data_split_indexes["label"] == "train"]
test_indexes = data_split_indexes.loc[data_split_indexes["label"] == "test"]


# In[10]:


# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[train_indexes["labeled_data_index"]]
testing_data = data_df.loc[test_indexes["labeled_data_index"]]


# In[11]:


# define metadata columns
# subset each column that contains metadata
metadata_train = training_data.filter(regex="Metadata")
# drop all metadata columns
train_data_x = training_data.drop(metadata_train.columns, axis=1)
train_treatments = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name where NSU = normalized signal units
train_data_y_cols = train_data_x.filter(regex="NSU").columns
train_data_y = training_data[train_data_y_cols]
train_data_x = train_data_x.drop(train_data_y_cols, axis=1)


# define metadata columns
# subset each column that contains metadata
metadata_test = testing_data.filter(regex="Metadata")
# drop all metadata columns
test_data_x = testing_data.drop(metadata_test.columns, axis=1)
test_treatments = testing_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name
test_data_y_cols = test_data_x.filter(regex="NSU").columns
test_data_y = testing_data[test_data_y_cols]
test_data_x = test_data_x.drop(test_data_y_cols, axis=1)


# In[12]:


print(train_data_x.shape, train_data_y.shape, test_data_x.shape, test_data_y.shape)


# In[13]:


# set model path from parameters
if (aggregation == True) and (nomic == True):
    model_path = pathlib.Path(f"models/regression/{cell_type}_aggregated_with_nomic/")
elif (aggregation == True) and (nomic == False):
    model_path = pathlib.Path(f"models/regression/{cell_type}_aggregated/")
elif (aggregation == False) and (nomic == True):
    model_path = pathlib.Path(f"models/regression/{cell_type}_sc_with_nomic/")
elif (aggregation == False) and (nomic == False):
    model_path = pathlib.Path(f"models/regression/{cell_type}_sc/")
else:
    print("Error")


# In[14]:


data_dict = {
    "train_data": {
        "data_x": train_data_x,
        "data_y": train_data_y,
        "col_names": train_data_y_cols,
        "metadata": metadata_train,
    },
    "test_data": {
        "data_x": test_data_x,
        "data_y": test_data_y,
        "col_names": test_data_y_cols,
        "metadata": metadata_test,
    },
}


# In[15]:


# list of metrics to use
output_metric_scores = {}


# In[16]:


# blank df for concatenated results
results_df = pd.DataFrame(
    columns=[
        "explained_variance",
        "neg_mean_absolute_error",
        "neg_mean_squared_error",
        "well",
        "treatment",
        "r2",
        "cytokine",
        "data_split",
        "shuffle",
        "predicted_value",
        "actual_value",
        "log10_neg_mean_absolute_error",
        "log10_neg_mean_squared_error",
        "log10_explained_variance",
    ]
)


# In[17]:


list_of_dfs = []
for data_split in data_dict.keys():
    data_x = data_dict[data_split]["data_x"]
    data_y = data_dict[data_split]["data_y"]
    col_names = data_dict[data_split]["col_names"]
    metadata = data_dict[data_split]["metadata"]

    if shuffle == "shuffled_baseline":
        model = joblib.load(
            f"../../1.train_models/{model_path}/{cytokine}_shuffled_baseline__all_nomic.joblib"
        )
    elif shuffle == "final":
        model = joblib.load(
            f"../../1.train_models/{model_path}/{cytokine}_final__all_nomic.joblib"
        )
    else:
        print("Error")

    # get the cytokine column of choice
    y_selected = data_y[cytokine]

    if shuffle == "shuffled_baseline":
        for column in data_x:
            np.random.shuffle(data_x[column].values)

    # get predictions
    predictions = model.predict(data_x)

    explained_variance = explained_variance_score(y_selected, predictions)
    output_metric_scores["explained_variance"] = explained_variance
    neg_mean_absolute_error = -mean_squared_error(y_selected, predictions)
    output_metric_scores["neg_mean_absolute_error"] = neg_mean_absolute_error
    neg_mean_squared_error = -mean_squared_error(y_selected, predictions)
    output_metric_scores["neg_mean_squared_error"] = neg_mean_squared_error
    r2 = r2_score(y_selected, predictions)
    output_metric_scores["treatment"] = metadata[
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
    ]
    output_metric_scores["well"] = metadata["Metadata_Well"].values
    df = pd.DataFrame.from_dict(output_metric_scores)
    df["r2"] = r2
    df["cytokine"] = cytokine
    df["data_split"] = data_split
    df["shuffle"] = shuffle
    df["predicted_value"] = predictions
    df["actual_value"] = y_selected
    df["log10_neg_mean_absolute_error"] = -np.log10(-df["neg_mean_absolute_error"])
    df["log10_neg_mean_squared_error"] = -np.log10(-df["neg_mean_squared_error"])
    df["log10_explained_variance"] = -np.log10(df["explained_variance"])

    # replace "[NSU]" with """
    df["cytokine"] = df["cytokine"].replace("[ \[\]NSU]", "", regex=True)
    df["cytokine"] = df["cytokine"].replace(" ", "_", regex=True)
    list_of_dfs.append(df)
    # concat the dataframes
results_df = pd.concat(list_of_dfs, axis=0)


# In[18]:


results_df.head()


# In[19]:


var_df = results_df.drop(
    columns=[
        "explained_variance",
        "neg_mean_absolute_error",
        "neg_mean_squared_error",
        "well",
        "treatment",
        "r2",
        "log10_neg_mean_absolute_error",
        "log10_neg_mean_squared_error",
        "log10_explained_variance",
    ]
)
# calculate the variance of the actual and predicted values per cytokine
var_df = var_df.groupby(["cytokine", "data_split", "shuffle"]).var()
var_df = pd.merge(
    var_df,
    results_df.groupby(["cytokine", "data_split", "shuffle"]).r2.unique(),
    left_index=True,
    right_index=True,
)
var_df.reset_index(inplace=True)
var_df.head()


# In[20]:


# set model path from parameters
if (aggregation == True) and (nomic == True):
    results_path = pathlib.Path(
        f"../results/regression/{cell_type}_aggregated_with_nomic/"
    )
elif (aggregation == True) and (nomic == False):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_aggregated/")
elif (aggregation == False) and (nomic == True):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_sc_with_nomic/")
elif (aggregation == False) and (nomic == False):
    results_path = pathlib.Path(f"../results/regression/{cell_type}_sc/")
else:
    print("Error")
pathlib.Path(results_path).mkdir(parents=True, exist_ok=True)


# In[21]:


# check if the model training metrics file exists
metrics_file = pathlib.Path(f"{results_path}/{cytokine}_{shuffle}_model_stats.csv")

results_df.to_csv(metrics_file, index=False)

# do the same for the variance df
# check if the model training metrics file exists
metrics_file = pathlib.Path(
    f"{results_path}/{cytokine}_{shuffle}_variance_r2_stats.csv"
)
var_df.to_csv(metrics_file, index=False)
