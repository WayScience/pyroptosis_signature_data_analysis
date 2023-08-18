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

# %% papermill={"duration": 3.641528, "end_time": "2023-08-06T18:56:40.949436", "exception": false, "start_time": "2023-08-06T18:56:37.307908", "status": "completed"} tags=[]
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
import torch
from sklearn import preprocessing

sys.path.append("../..")
from MLP_utils.parameters import Parameters
from MLP_utils.utils import (
    Dataset_formatter,
    optimized_model_create,
    output_stats,
    parameter_set,
    results_output,
    test_optimized_model,
    un_nest,
)
from sklearn.metrics import (
    accuracy_score,
    auc,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)

sys.path.append("../../..")

# %% papermill={"duration": 0.00606, "end_time": "2023-08-06T18:56:40.960224", "exception": false, "start_time": "2023-08-06T18:56:40.954164", "status": "completed"} tags=["injected-parameters"]
# Parameters
CELL_TYPE = "SHSY5Y"
CONTROL_NAME = "DMSO_0.100_DMSO_0.025"
TREATMENT_NAME = "LPS_100.000_DMSO_0.025"
MODEL_NAME = "DMSO_0.025_vs_LPS_100"
SHUFFLE = False

# %% papermill={"duration": 0.00696, "end_time": "2023-08-06T18:56:40.968861", "exception": false, "start_time": "2023-08-06T18:56:40.961901", "status": "completed"} tags=[]
ml_configs_file = pathlib.Path("../../MLP_utils/binary_config.toml").resolve(
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
mlp_params.SHUFFLE = SHUFFLE

# %% papermill={"duration": 56.196209, "end_time": "2023-08-06T18:57:37.166710", "exception": false, "start_time": "2023-08-06T18:56:40.970501", "status": "completed"} tags=[]
# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.Path(
    f"../../../data/{mlp_params.CELL_TYPE}_preprocessed_sc_norm.parquet"
).resolve(strict=True)

df = pq.read_table(file_path).to_pandas()


# %% papermill={"duration": 0.010934, "end_time": "2023-08-06T18:57:37.182482", "exception": false, "start_time": "2023-08-06T18:57:37.171548", "status": "completed"} tags=[]
def test_loop(df, output_name, title):
    # Code snippet for metadata extraction by Jenna Tomkinson
    df_metadata = list(df.columns[df.columns.str.startswith("Metadata")])

    # define which columns are data and which are descriptive
    df_descriptive = df[df_metadata]
    df_values = df.drop(columns=df_metadata)
    # Creating label encoder
    le = preprocessing.LabelEncoder()
    # Converting strings into numbers
    print(df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique().tolist())
    lst_of_treatments = (
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique().tolist()
    )

    df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = le.fit_transform(
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    )
    print(df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique().tolist())
    lst_of_coded_treatments = (
        df_values["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique().tolist()
    )
    # make a dictionary of the treatments and their corresponding codes to decode later
    dict_of_treatments = {}
    for i, j in zip(
        lst_of_coded_treatments,
        lst_of_treatments,
    ):
        dict_of_treatments[i] = j
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
    # convert data class into a dataloader to be compatible with pytorch
    test_loader = torch.utils.data.DataLoader(
        dataset=test_data, batch_size=1, shuffle=mlp_params.SHUFFLE
    )
    model, _ = optimized_model_create(mlp_params, mlp_params.MODEL_NAME)
    # calling the testing function and outputting list values of tested model
    if mlp_params.MODEL_TYPE == "Multi_Class" or mlp_params.MODEL_TYPE == "Regression":
        y_pred_list = test_optimized_model(
            model,
            test_loader,
            mlp_params,
            model_name=mlp_params.MODEL_NAME,
            shuffle=mlp_params.SHUFFLE,
        )
    elif mlp_params.MODEL_TYPE == "Binary_Classification":
        y_pred_list, y_pred_prob_list = test_optimized_model(
            model,
            test_loader,
            mlp_params,
            model_name=mlp_params.MODEL_NAME,
            shuffle=mlp_params.SHUFFLE,
        )
    else:
        raise Exception("Model type must be specified for proper model testing")

    # un-nest list if nested i.e. length of input data does not match length of output data
    if len(y_pred_list) != len(df_values_Y):
        y_pred_list = un_nest(y_pred_list)
        y_pred_prob_list = un_nest(y_pred_prob_list)
    else:
        pass

    stats, recall, precision, f1, precision_, recall_, threshold_ = output_stats(
        y_pred_list,
        df_values_Y,
        mlp_params,
        y_pred_prob_list,
        test_name=f"{output_name}_all_testing",
        model_name=mlp_params.MODEL_NAME,
        title=title,
        shuffle=mlp_params.SHUFFLE,
    )
    return (
        stats,
        recall,
        precision,
        f1,
        precision_,
        recall_,
        threshold_,
        dict_of_treatments,
    )


# %% papermill={"duration": 0.284008, "end_time": "2023-08-06T18:57:37.468142", "exception": false, "start_time": "2023-08-06T18:57:37.184134", "status": "completed"} tags=[]
print(df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique().tolist())

# %% papermill={"duration": 0.007431, "end_time": "2023-08-06T18:57:37.477993", "exception": false, "start_time": "2023-08-06T18:57:37.470562", "status": "completed"} tags=[]
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

# %%
# create a dataframe to store the model stats
model_stats_df = pd.DataFrame(
    columns=[
        "treatments_tested",
        "model",
        "group",
        "shuffled_data",
        "PR_Threshold",
        "Precision",
        "Recall",
    ]
)
model_stats_df

# %% papermill={"duration": 1448.376356, "end_time": "2023-08-06T19:21:45.856103", "exception": false, "start_time": "2023-08-06T18:57:37.479747", "status": "completed"} tags=[]
for i in paired_treatment_list:
    # filter df to only include the two treatments to test
    test_df = df.query(
        f"oneb_Metadata_Treatment_Dose_Inhibitor_Dose == '{i[0]}' | oneb_Metadata_Treatment_Dose_Inhibitor_Dose == '{i[1]}'"
    )
    output_name = ("__").join(
        test_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()
    )

    print(output_name)

    title = f'{output_name.split("__")[0].split("_")[0]} vs {("__").join(output_name.split("__")[1].split("_")[:2])}'
    print(title)
    (
        stats,
        recall,
        precision,
        f1,
        precision_,
        recall_,
        threshold_,
        dict_of_treatments,
    ) = test_loop(test_df, output_name, title)
    print(recall, precision, f1)

    threshold_ = np.append(threshold_, None)
    stats_df = pd.DataFrame(
        {
            "PR_Threshold": threshold_,
            "Precision": precision_,
            "Recall": recall_,
        }
    )

    stats_df["treatments_tested"] = "0 vs 1"
    # make it so that the second treatment is always the one that is being tested as the positive label
    stats_df["treatments_tested"] = stats_df["treatments_tested"].replace(
        "0 vs 1", f"{dict_of_treatments[0]} vs {dict_of_treatments[1]}"
    )
    stats_df["model"] = mlp_params.MODEL_NAME
    stats_df["group"] = "test"
    stats_df["shuffled_data"] = mlp_params.SHUFFLE
    stats_df
    model_stats_df = pd.concat([model_stats_df, stats_df], axis=0)

# %%
model_stats_df

# %%
# set path for the model training metrics
metrics_path = pathlib.Path(
    f"../../results/{mlp_params.MODEL_TYPE}/{mlp_params.MODEL_NAME}/{mlp_params.CELL_TYPE}"
)
metrics_path.mkdir(parents=True, exist_ok=True)
# check if the model training metrics file exists
metrics_file = pathlib.Path(f"{metrics_path}/testing_metrics.csv")
if metrics_file.exists():
    metrics_df = pd.read_csv(metrics_file)
    if len(metrics_df["shuffled_data"].unique()) > 1:
        pass
    elif metrics_df["shuffled_data"].unique() == mlp_params.SHUFFLE:
        pass
    else:
        metrics_df = pd.concat([metrics_df, model_stats_df], axis=0)
        metrics_df.to_csv(metrics_file, index=False)
else:
    model_stats_df.to_csv(metrics_file, index=False)

# %%
