#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

sys.path.append("../..")
from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    optimized_model_create,
    parameter_set,
    results_output,
    test_optimized_model,
    un_nest,
)

sys.path.append("../../..")


# In[2]:


# Parameters
SHUFFLE_DATA = False
CELL_TYPE = "PBMC"
CONTROL_NAME = "DMSO_0.100_DMSO_0.025"
TREATMENT_NAME = "Thapsigargin_1.000_DMSO_0.025"
MODEL_NAME = "DMSO_0.025_vs_Thapsigargin_1"


# In[3]:


ml_configs_file = pathlib.path("../../MLP_utils/binary_config.toml").resolve(
    strict=True
)
ml_configs = toml.load(ml_configs_file)
params = Parameters()
mlp_params = parameter_set(params, ml_configs)

# overwrite mlp_params via command line arguments from papermill
mlp_params.CELL_TYPE = CELL_TYPE
mlp_params.MODEL_NAME = MODEL_NAME
mlp_params.CONTROL_NAME = CONTROL_NAME
mlp_params.TREATMENT_NAME = TREATMENT_NAME
mlp_params.MODEL_NAME = MODEL_NAME


# In[4]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)

df = pq.read_table(file_path).to_pandas()


# In[5]:


def test_loop(df, output_name, title):
    # Code snippet for metadata extraction by Jenna Tomkinson
    df_metadata = list(df.columns[df.columns.str.startswith("Metadata")])

    # define which columns are data and which are descriptive
    df_descriptive = df[df_metadata]
    df_values = df.drop(columns=df_metadata)
    # Creating label encoder
    le = preprocessing.LabelEncoder()
    # Converting strings into numbers
    df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = le.fit_transform(
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    )
    # split into X and Y where Y are the predictive column and x are the observable data
    df_values_X = df_values.drop(
        [
            "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
            "twob_Metadata_Treatment_Dose_Inhibitor_Dose",
            "threeb_Metadata_Treatment_Dose_Inhibitor_Dose",
            "fourb_Metadata_Treatment_Dose_Inhibitor_Dose",
        ],
        axis=1,
    )
    df_values_Y = df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    test_data = Dataset_formatter(
        torch.FloatTensor(df_values_X.values), torch.FloatTensor(df_values_Y.values)
    )

    mlp_params.IN_FEATURES = df_values_X.shape[1]
    print("Number of in features: ", mlp_params.IN_FEATURES)
    if mlp_params.MODEL_TYPE == "Regression":
        mlp_params.OUT_FEATURES = 1
    else:
        mlp_params.OUT_FEATURES = len(
            df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
        )

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

    # convert data class into a dataloader to be compatible with pytorch
    test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1)
    model = optimized_model_create(mlp_params, mlp_params.MODEL_NAME)
    # calling the testing function and outputting list values of tested model
    if mlp_params.MODEL_TYPE == "Multi_Class" or mlp_params.MODEL_TYPE == "Regression":
        y_pred_list = test_optimized_model(
            model, test_loader, mlp_params, model_name=mlp_params.MODEL_NAME
        )
    elif mlp_params.MODEL_TYPE == "Binary_Classification":
        y_pred_list, y_pred_prob_list = test_optimized_model(
            model, test_loader, mlp_params, model_name=mlp_params.MODEL_NAME
        )
    else:
        raise Exception("Model type must be specified for proper model testing")

    # un-nest list if nested i.e. length of input data does not match length of output data
    if len(y_pred_list) != len(df_values_Y):
        y_pred_list = un_nest(y_pred_list)
        y_pred_prob_list = un_nest(y_pred_prob_list)
    else:
        pass
    # Call visualization function
    # calling the testing function and outputting list values of tested model
    if mlp_params.MODEL_TYPE == "Multi_Class" or mlp_params.MODEL_TYPE == "Regression":
        confusion_matrix_df = results_output(
            y_pred_list,
            df_values_Y,
            mlp_params,
            test_name=f"{output_name}_all_testing",
            model_name=mlp_params.MODEL_NAME,
        )
    elif mlp_params.MODEL_TYPE == "Binary_Classification":
        results_output(
            y_pred_list,
            df_values_Y,
            mlp_params,
            y_pred_prob_list,
            test_name=f"{output_name}_all_testing",
            model_name=mlp_params.MODEL_NAME,
            title=title,
        )
    else:
        raise Exception("Model type must be specified for proper model testing")


# In[6]:


df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()


# In[7]:


paired_treatment_list = [
    ["DMSO_0.100_DMSO_0.025", "LPS_100.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_0.100_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_1.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_10.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_100.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "Flagellin_0.100_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "Flagellin_1.000_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "Flagellin_1.000_Disulfiram_1.0"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_100.000_1.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_100.000_3.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_100.000_10.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_1.000_1.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_1.000_3.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "LPS_Nigericin_1.000_10.0_DMSO_0.025"],
    ["DMSO_0.100_DMSO_0.025", "H2O2_100.000_Z-VAD-FMK_100.0"],
    ["DMSO_0.100_DMSO_0.025", "H2O2_100.000_DMSO_0.025"],
    ["LPS_100.000_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["LPS_100.000_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
    ["LPS_10.000_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["LPS_10.000_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
    ["LPS_1.000_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["LPS_1.000_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
    ["LPS_0.100_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["LPS_0.100_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
    ["LPS_0.010_DMSO_0.025", "Thapsigargin_1.000_DMSO_0.025"],
    ["LPS_0.010_DMSO_0.025", "Thapsigargin_10.000_DMSO_0.025"],
]


# In[8]:


for i, j in paired_treatment_list:
    test_df = df.query(
        f"oneb_Metadata_Treatment_Dose_Inhibitor_Dose == '{j}' | oneb_Metadata_Treatment_Dose_Inhibitor_Dose == '{i}'"
    )
    output_name = (" ").join(
        test_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
    )

    print(output_name)

    title = f'{output_name.split(" ")[0].split("_")[0]} vs {(" ").join(output_name.split(" ")[1].split("_")[:2])}'
    test_loop(test_df, output_name, title)


# In[ ]:
