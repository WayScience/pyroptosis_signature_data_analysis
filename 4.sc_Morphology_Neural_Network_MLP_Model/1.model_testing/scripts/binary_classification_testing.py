# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Interstellar
#     language: python
#     name: python3
# ---

# %% papermill={"duration": 3.760048, "end_time": "2023-08-05T20:11:33.826236", "exception": false, "start_time": "2023-08-05T20:11:30.066188", "status": "completed"} tags=[]
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly
import pyarrow.parquet as pq
import seaborn as sns
import toml
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn import preprocessing
from sklearn.metrics import (
    auc,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import train_test_split

sys.path.append("../..")
from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    data_split,
    extract_best_trial_params,
    objective_model_optimizer,
    optimized_model_create,
    parameter_set,
    plot_metric_vs_epoch,
    results_output,
    test_optimized_model,
    train_optimized_model,
    un_nest,
)

sys.path.append("../../..")
from utils.utils import df_stats

# %% papermill={"duration": 0.006762, "end_time": "2023-08-05T20:11:33.835879", "exception": false, "start_time": "2023-08-05T20:11:33.829117", "status": "completed"} tags=["injected-parameters"]
# Parameters
SHUFFLE_DATA = False
CELL_TYPE = "PBMC"
CONTROL_NAME = "DMSO_0.100_DMSO_0.025"
TREATMENT_NAME = "Thapsigargin_1.000_DMSO_0.025"
MODEL_NAME = "DMSO_0.025_vs_Thapsigargin_1"

# %% papermill={"duration": 0.006503, "end_time": "2023-08-05T20:11:33.844096", "exception": false, "start_time": "2023-08-05T20:11:33.837593", "status": "completed"} tags=[]
data = Path("../../MLP_utils/binary_config.toml")
config = toml.load(data)
params = Parameters()
params = parameter_set(params, config)

# overwrite params via command line arguments from papermill
params.CELL_TYPE = CELL_TYPE
params.MODEL_NAME = MODEL_NAME
params.CONTROL_NAME = CONTROL_NAME
params.TREATMENT_NAME = TREATMENT_NAME
params.MODEL_NAME = MODEL_NAME

# %% papermill={"duration": 311.682826, "end_time": "2023-08-05T20:16:45.528630", "exception": false, "start_time": "2023-08-05T20:11:33.845804", "status": "completed"} tags=[]
# Import Data
# set data file path under pathlib path for multi-system use
file_path = Path(f"../../../data/{params.CELL_TYPE}_preprocessed_sc_norm.parquet")

df = pq.read_table(file_path).to_pandas()


# %% papermill={"duration": 0.01145, "end_time": "2023-08-05T20:16:45.592116", "exception": false, "start_time": "2023-08-05T20:16:45.580666", "status": "completed"} tags=[]
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

    params.IN_FEATURES = df_values_X.shape[1]
    print("Number of in features: ", params.IN_FEATURES)
    if params.MODEL_TYPE == "Regression":
        params.OUT_FEATURES = 1
    else:
        params.OUT_FEATURES = len(
            df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
        )

    print("Number of out features: ", params.OUT_FEATURES)

    if params.OUT_FEATURES > 2:
        params.MODEL_TYPE = "Multi_Class"
    elif params.OUT_FEATURES == 2:
        params.OUT_FEATURES = params.OUT_FEATURES - 1
        params.MODEL_TYPE = "Binary_Classification"
    elif params.OUT_FEATURES == 1:
        params.MODEL_TYPE = "Regression"
    else:
        pass
    print(params.MODEL_TYPE)

    # convert data class into a dataloader to be compatible with pytorch
    test_loader = torch.utils.data.DataLoader(dataset=test_data, batch_size=1)
    model = optimized_model_create(params, params.MODEL_NAME)
    # calling the testing function and outputting list values of tested model
    if params.MODEL_TYPE == "Multi_Class" or params.MODEL_TYPE == "Regression":
        y_pred_list = test_optimized_model(
            model, test_loader, params, model_name=params.MODEL_NAME
        )
    elif params.MODEL_TYPE == "Binary_Classification":
        y_pred_list, y_pred_prob_list = test_optimized_model(
            model, test_loader, params, model_name=params.MODEL_NAME
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
    if params.MODEL_TYPE == "Multi_Class" or params.MODEL_TYPE == "Regression":
        confusion_matrix_df = results_output(
            y_pred_list,
            df_values_Y,
            params,
            test_name=f"{output_name}_all_testing",
            model_name=params.MODEL_NAME,
        )
    elif params.MODEL_TYPE == "Binary_Classification":
        results_output(
            y_pred_list,
            df_values_Y,
            params,
            y_pred_prob_list,
            test_name=f"{output_name}_all_testing",
            model_name=params.MODEL_NAME,
            title=title,
        )
    else:
        raise Exception("Model type must be specified for proper model testing")


# %% papermill={"duration": 0.284484, "end_time": "2023-08-05T20:16:45.878644", "exception": false, "start_time": "2023-08-05T20:16:45.594160", "status": "completed"} tags=[]
df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()

# %% papermill={"duration": 0.007489, "end_time": "2023-08-05T20:16:45.888088", "exception": false, "start_time": "2023-08-05T20:16:45.880599", "status": "completed"} tags=[]
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

# %% papermill={"duration": 1430.863697, "end_time": "2023-08-05T20:40:36.753561", "exception": false, "start_time": "2023-08-05T20:16:45.889864", "status": "completed"} tags=[]
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

# %% papermill={"duration": 0.029857, "end_time": "2023-08-05T20:40:36.829055", "exception": false, "start_time": "2023-08-05T20:40:36.799198", "status": "completed"} tags=[]
