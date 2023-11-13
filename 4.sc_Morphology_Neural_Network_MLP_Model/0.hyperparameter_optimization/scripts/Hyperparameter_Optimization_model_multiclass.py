#!/usr/bin/env python
# coding: utf-8

# ## Hyperparameter tuning via Optuna

# ### Being a binary model this notebook will be limited to predicting one class 1 or 0, yes or no.
# ### Here I will be predicting if a cell received a treatment or not

# In[1]:


import pathlib
import sys

import numpy as np
import optuna
import pandas as pd
import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

sys.path.append("../..")

from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    data_split,
    extract_best_trial_params,
    objective_model_optimizer,
    parameter_set,
    plot_metric_vs_epoch,
    results_output,
    test_optimized_model,
    train_optimized_model,
    un_nest,
)
from sklearn.model_selection import train_test_split

sys.path.append("../../..")
from utils.utils import df_stats

# ## Papermill is used for executing notebooks in the CLI with multiple parameters
# Here the `injected-parameters` cell is used to inject parameters into the notebook via papermill.
# This enables multiple notebooks to be executed with different parameters, preventing to manually update parameters or have multiple copies of the notebook.

# In[2]:


# Parameters
CELL_TYPE = "PBMC"
MODEL_NAME = "MultiClass_MLP"


# In[3]:


ml_configs_file = pathlib.Path("../../MLP_utils/multi_class_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)

# overwrite params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
mlp_params.MODEL_NAME = MODEL_NAME
MODEL_TYPE = mlp_params.MODEL_TYPE
HYPERPARAMETER_BATCH_SIZE = mlp_params.HYPERPARAMETER_BATCH_SIZE


# In[4]:


# Import Data
# set data file path under pathlib path for multi-system use

file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)

df1 = pd.read_parquet(file_path)


# In[5]:


if params.MODEL_NAME == "MultiClass_MLP_h202_remove":
    # drop H2O2_100.000_uM_DMSO_0.025_% and H2O2_100.000_nM_DMSO_0.025_% while keeping the index order
    df1 = df1[
        ~df1["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
            ["H2O2_100.000_uM_DMSO_0.025_%", "H2O2_100.000_nM_DMSO_0.025_%"]
        )
    ]


# In[6]:


# get paths for toml files
ground_truth_file_path = pathlib.Path(f"../../MLP_utils/ground_truth.toml").resolve(
    strict=True
)
treatment_splits_file_path = pathlib.Path(f"../../MLP_utils/splits.toml").resolve(
    strict=True
)
# read toml files
ground_truth = toml.load(ground_truth_file_path)
treatment_splits = toml.load(treatment_splits_file_path)


# In[7]:


# get information from toml files
apoptosis_groups_list = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_groups_list = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
healthy_groups_list = ground_truth["Healthy"]["healthy_groups_list"]
test_split_100 = treatment_splits["splits"]["data_splits_100"]
test_split_75 = treatment_splits["splits"]["data_splits_75"]


# In[8]:


np.random.seed(0)
if mlp_params.DATA_SUBSET_OPTION == "True":
    df1 = df1.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose").apply(
        lambda x: x.sample(n=mlp_params.DATA_SUBSET_NUMBER, random_state=0)
    )
    print("Data Subset Is On")
    print(f"Data is subset to {mlp_params.DATA_SUBSET_NUMBER} per treatment group")
    print(df1.shape)
    df1.reset_index(drop=True, inplace=True)
else:
    print("Data Subset Is Off")


# In[9]:


# add apoptosis, pyroptosis and healthy columns to dataframe
df1["apoptosis"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in apoptosis_groups_list,
    axis=1,
)
df1["pyroptosis"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in pyroptosis_groups_list,
    axis=1,
)
df1["healthy"] = df1.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    in healthy_groups_list,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column
df1["labels"] = df1.apply(
    lambda row: "apoptosis"
    if row["apoptosis"]
    else "pyroptosis"
    if row["pyroptosis"]
    else "healthy",
    axis=1,
)
# drop apoptosis, pyroptosis, and healthy columns
df1.drop(columns=["apoptosis", "pyroptosis", "healthy"], inplace=True)


# ### Split said data

# In[10]:


# randomly select wells to hold out for testing one per treatment group
# stratified by treatment group
np.random.seed(seed=0)
wells_to_hold = (
    df1.groupby("oneb_Metadata_Treatment_Dose_Inhibitor_Dose")
    .agg(np.random.choice)["Metadata_Well"]
    .to_list()
)
df_holdout = df1[df1["Metadata_Well"].isin(wells_to_hold)]
df = df1[~df1["Metadata_Well"].isin(wells_to_hold)]


print("Wells held out for testing:", df_holdout["Metadata_Well"].unique())
print(
    "Wells to use for training, validation, and testing", df1["Metadata_Well"].unique()
)


# In[11]:


# variable test and train set splits
# 100% test set
# subset the following treatments for test set
test_set_all = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_100)
]
# 75% test set and 25% train set
test_set_75 = df[df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_split_75)]

test_100_and_75 = test_split_100 + test_split_75

# 50% test set and 50% train set
# get all treatments that are not in the_test_set_all and the test_set_75
test_set_50 = df[
    ~df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(test_100_and_75)
]

print(test_set_all.shape, test_set_75.shape, test_set_50.shape)


# In[12]:


# get the train test splits from each group
# 100% test set
test_set_all

# 75% test set and 25% train set
test_ratio = 0.75
training_data_set_75, testing_data_set_75 = train_test_split(
    test_set_75,
    test_size=test_ratio,
    stratify=test_set_75["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    random_state=0,
)

# 50% test set and 50% train set
test_ratio = 0.5
training_data_set_50, testing_data_set_50 = train_test_split(
    test_set_50,
    test_size=test_ratio,
    stratify=test_set_50["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    random_state=0,
)

# verify that the correct splits have been made
# 100% test set
print(f"Shape for the 100% test set: {test_set_all.shape}\n")

# 75% test set and 25% train set
print(
    f"Shape for the 75% test set: {training_data_set_75.shape};\nShape for the 75% train set: {testing_data_set_75.shape}\n"
)

# 50% test set and 50% train set
print(
    f"Shape for the 50% test set: {training_data_set_50.shape};\nShape for the 50% train set: {testing_data_set_50.shape}"
)

print(f"Shape for the holdout set: {df_holdout.shape}")


# In[13]:


# combine all testing sets together while preserving the index
testing_data_set = pd.concat(
    [test_set_all, testing_data_set_75, testing_data_set_50], axis=0
)
testing_data_set = testing_data_set.sort_index()
testing_data_set

# combine all training sets together while preserving the index
training_data_set = pd.concat([training_data_set_75, training_data_set_50], axis=0)
training_data_set = training_data_set.sort_index()
training_data_set

training_data_set, val_data_set = train_test_split(
    training_data_set,
    test_size=0.20,
    stratify=training_data_set["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
)
print(
    f"""
    Testing set length: {len(testing_data_set)}\n
    Training set length: {len(training_data_set)}\n
    Validation set length: {len(val_data_set)}\n
    Holdout set length: {len(df_holdout)}"""
)


# In[14]:


# # train
# # downsample healthy and pyroptosis to match number of apoptosis
# # to balance classes
# df_healthy_train = training_data_set[training_data_set["labels"] == "healthy"]
# df_pyroptosis_train = training_data_set[training_data_set["labels"] == "pyroptosis"]
# df_apoptosis_train = training_data_set[training_data_set["labels"] == "apoptosis"]
# print(df_healthy_train.shape, df_pyroptosis_train.shape, df_apoptosis_train.shape)

# # downsample healthy and pyroptosis to match number of apoptosis
# df_healthy_train = df_healthy_train.sample(n=df_apoptosis_train.shape[0], random_state=0)
# df_pyroptosis_train = df_pyroptosis_train.sample(n=df_apoptosis_train.shape[0], random_state=0)
# print(df_healthy_train.shape, df_pyroptosis_train.shape, df_apoptosis_train.shape)
# training_data_set = pd.concat([df_healthy_train, df_pyroptosis_train, df_apoptosis_train])
# # show that the df was downsampled and recombined correctly
# assert (df_healthy_train + df_pyroptosis_train + df_apoptosis_train).shape[0] == training_data_set.shape[0]


# # validation
# # downsample healthy and pyroptosis to match number of apoptosis
# # to balance classes
# df_healthy_val = val_data_set[val_data_set["labels"] == "healthy"]
# df_pyroptosis_val = val_data_set[val_data_set["labels"] == "pyroptosis"]
# df_apoptosis_val = val_data_set[val_data_set["labels"] == "apoptosis"]
# print(df_healthy_val.shape, df_pyroptosis_val.shape, df_apoptosis_val.shape)

# # downsample healthy and pyroptosis to match number of apoptosis
# df_healthy_val = df_healthy_val.sample(n=df_apoptosis_val.shape[0], random_state=0)
# df_pyroptosis_val = df_pyroptosis_val.sample(n=df_apoptosis_val.shape[0], random_state=0)
# print(df_healthy_val.shape, df_pyroptosis_val.shape, df_apoptosis_val.shape)
# val_data_set = pd.concat([df_healthy_val, df_pyroptosis_val, df_apoptosis_val])
# # show that the df was downsampled and recombined correctly
# assert (df_healthy_val + df_pyroptosis_val + df_apoptosis_val).shape[0] == val_data_set.shape[0]


# # test
# # downsample healthy and pyroptosis to match number of apoptosis
# # to balance classes
# df_healthy_test = testing_data_set[testing_data_set["labels"] == "healthy"]
# df_pyroptosis_test = testing_data_set[testing_data_set["labels"] == "pyroptosis"]
# df_apoptosis_test = testing_data_set[testing_data_set["labels"] == "apoptosis"]
# print(df_healthy_test.shape, df_pyroptosis_test.shape, df_apoptosis_test.shape)

# # downsample healthy and pyroptosis to match number of apoptosis
# df_healthy_test = df_healthy_test.sample(n=df_apoptosis_test.shape[0], random_state=0)
# df_pyroptosis_test = df_pyroptosis_test.sample(n=df_apoptosis_test.shape[0], random_state=0)
# print(df_healthy_test.shape, df_pyroptosis_test.shape, df_apoptosis_test.shape)
# testing_data_set = pd.concat([df_healthy_test, df_pyroptosis_test, df_apoptosis_test])
# # show that the df was downsampled and recombined correctly
# assert (df_healthy_test + df_pyroptosis_test + df_apoptosis_test).shape[0] == testing_data_set.shape[0]


# print(len(training_data_set), len(val_data_set), len(testing_data_set), len(df_holdout))


# In[15]:


# get the indexes for the training and testing sets

training_data_set_index = training_data_set.index
val_data_set_index = val_data_set.index
testing_data_set_index = testing_data_set.index
df_holdout_index = df_holdout.index


# In[16]:


print(
    training_data_set_index.shape,
    val_data_set_index.shape,
    testing_data_set_index.shape,
    df_holdout_index.shape,
)
print(
    training_data_set_index.shape[0]
    + val_data_set_index.shape[0]
    + testing_data_set_index.shape[0]
    + df_holdout_index.shape[0]
)


# In[17]:


# create pandas dataframe with all indexes and their respective labels, stratified by phenotypic class
index_data = []
for index in training_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "train"})
for index in val_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "val"})
for index in testing_data_set_index:
    index_data.append({"labeled_data_index": index, "label": "test"})
for index in df_holdout_index:
    index_data.append({"labeled_data_index": index, "label": "holdout"})

# make index data a dataframe and sort it by labeled data index
index_data = pd.DataFrame(index_data)
index_data


# In[18]:


save_path = pathlib.Path(f"../indexes/{CELL_TYPE}/multi_class/")

print(save_path)
# create save path if it doesn't exist
save_path.mkdir(parents=True, exist_ok=True)


# In[19]:


# save indexes as tsv file
index_data.to_csv(
    f"{save_path}/{params.MODEL_NAME}_data_split_indexes.tsv", sep="\t", index=False
)


# #### Set up Data to be compatible with model

# ##### Classification Models:
# Comment out code if using regression

# In[20]:


# Code snippet for metadata extraction by Jenna Tomkinson
df_metadata = list(df.columns[df.columns.str.contains("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df1[df_metadata]
df_values = df1.drop(columns=df_metadata)


# In[21]:


# get the class weights for the loss function to account for class imbalance
# get the number of samples in each class
targets, counts = np.unique(df1["labels"], return_counts=True)
print(targets, counts)
total_counts = np.sum(counts)
# get the class weights
class_weights = []
class_targets = []
for class_name in enumerate(targets):
    class_targets.append(class_name[0])
for count in enumerate(counts):
    class_weights.append(1 - (count[1] / total_counts))
print(class_targets, class_weights)
# write the class weights to a file for use in the model
class_weights_file = pathlib.Path(f"../class_weights/{CELL_TYPE}/multi_class/")
class_weights_file.mkdir(parents=True, exist_ok=True)
with open(f"{class_weights_file}/class_weights.txt", "w") as filehandle:
    for listitem in class_weights:
        filehandle.write("%s\n" % listitem)


# In[22]:


# Creating label encoder
le = preprocessing.LabelEncoder()
df_values["new_labels"] = le.fit_transform(df_values["labels"])
# get mini dataframe that contains the decoder
decoder = df_values[["labels", "new_labels"]].drop_duplicates()
# split into X and Y where Y are the predictive column and x are the observable data
df_values_X = df_values.drop(
    ["new_labels", "labels"],
    axis=1,
)
df_values_Y = df_values["new_labels"]
df_values_Y.head()
df_values_Y.unique()


# #### Split Data - All Models can proceed through this point

# In[23]:


# split into train and test sets from indexes previously defined

X_train = df_values_X.loc[training_data_set_index]
X_val = df_values_X.loc[val_data_set_index]
X_test = df_values_X.loc[testing_data_set_index]
X_holdout = df_values_X.loc[df_holdout_index]

Y_train = df_values_Y.loc[training_data_set_index]
Y_val = df_values_Y.loc[val_data_set_index]
Y_test = df_values_Y.loc[testing_data_set_index]
Y_holdout = df_values_Y.loc[df_holdout_index]


# In[24]:


# produce data objects for train, val and test datasets
train_data = Dataset_formatter(
    torch.FloatTensor(X_train.values), torch.FloatTensor(Y_train.values)
)
val_data = Dataset_formatter(
    torch.FloatTensor(X_val.values), torch.FloatTensor(Y_val.values)
)
test_data = Dataset_formatter(
    torch.FloatTensor(X_test.values), torch.FloatTensor(Y_test.values)
)


# In[25]:


mlp_params.IN_FEATURES = X_train.shape[1]
print("Number of in features: ", mlp_params.IN_FEATURES)
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.OUT_FEATURES = 1
else:
    mlp_params.OUT_FEATURES = len(df_values["labels"].unique())

print("Number of out features: ", mlp_params.OUT_FEATURES)

if mlp_params.OUT_FEATURES > 2:
    mlp_params.MODEL_TYPE = "Multi_Class"
elif mlp_params.OUT_FEATURES == 2:
    mlp_params.OUT_FEATURES = mlp_params.OUT_FEATURES - 1
    mlp_params.MODEL_TYPE = "Binary_Classification"
elif mlp_params.OUT_FEATURES == 1:
    mlp_params.MODEL_TYPE = "Regression"
else:
    pass
print(mlp_params.MODEL_TYPE)


# In[26]:


# convert data class into a dataloader to be compatible with pytorch
train_loader = torch.utils.data.DataLoader(
    dataset=train_data, batch_size=mlp_params.HYPERPARAMETER_BATCH_SIZE
)
valid_loader = torch.utils.data.DataLoader(
    dataset=val_data, batch_size=mlp_params.HYPERPARAMETER_BATCH_SIZE
)


# In[27]:


# check device
print(mlp_params.DEVICE)


# In[28]:


# no accuracy function must be loss for regression
if mlp_params.MODEL_TYPE == "Regression":
    mlp_params.METRIC = "loss"
    pass


# wrap the objective function inside of a lambda function to pass args...
objective_lambda_func = lambda trial: objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=trial,
    params=params,
    metric=mlp_params.METRIC,
    return_info=False,
    class_weights=class_weights,
)


# Study is the object for model optimization
study = optuna.create_study(direction=f"{mlp_params.DIRECTION}")
# Here I apply the optimize function of the study to the objective function
# This optimizes each parameter specified to be optimized from the defined search space
study.optimize(objective_lambda_func, n_trials=mlp_params.N_TRIALS)
# Prints out the best trial's optimized parameters
objective_model_optimizer(
    train_loader,
    valid_loader,
    trial=study.best_trial,
    params=params,
    metric=mlp_params.METRIC,
    return_info=True,
    class_weights=class_weights,
)


# In[29]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_optimization_history(study)


graph_path = f"{graph_path}/plot_optimization_history_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[30]:


# create graph directory for this model
graph_path = pathlib.Path(
    f"../../figures/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}/hyperparameter_optimization"
).resolve(strict=True)

pathlib.Path(graph_path).mkdir(parents=True, exist_ok=True)
fig = optuna.visualization.plot_intermediate_values(study)

graph_path = f"{graph_path}/plot_intermediate_values_graph"

fig.write_image(pathlib.Path(f"{graph_path}.png"))
fig.show()


# In[31]:


param_dict = extract_best_trial_params(
    study.best_params, params, model_name=mlp_params.MODEL_NAME
)
