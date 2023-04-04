#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

import pandas as pd
import plotly.express as px
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression

# In[2]:


sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

# In[3]:


# Define inputs and outputs
df = pd.read_csv(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz",
    low_memory=False,
)
# output_dir = pathlib.Path("results")
# output_cp_file = pathlib.Path(output_dir, "linear_model_cp_features.tsv")


# In[4]:


# Recycled code from: https://github.com/WayScience/NF1_SchwannCell_data/blob/main/5_analyze_data/notebooks/linear_model/fit_linear_model.ipynb
cell_count_df = (
    df.groupby("Metadata_Well")["Metadata_Plate"]
    .count()
    .reset_index()
    .rename(columns={"Metadata_Plate": "Metadata_number_of_singlecells"})
)

df = df.merge(cell_count_df, on="Metadata_Well")


# Define CellProfiler features
cp_features = infer_cp_features(df)

print(f"We are testing {len(cp_features)} CellProfiler features")
# Drop na and reindex accordingly
df = df.dropna()
df = df.reindex()
df = df.assign(
    Metadata_Treatment_and_Dose=lambda x: df["Metadata_treatment"]
    + "_"
    + df["Metadata_dose"]
)
df["Metadata_Treatment_and_Dose"].unique()


# In[5]:


treatment = [
    "LPS_10µg/ml",
    "Disulfiram_2.5µM",
    "LPS_1µg/ml",
    "Disulfiram_0.1µM",
    "H2O2_500µM",
    "Thapsi_10µM",
    "H2O2_50µM",
    "Thapsi_1µM",
    "ATP_1mM",
    "LPS + Nigericin_1µg/ml + 10µM",
    "ATP_0.1mM",
    "LPS + Nigericin_1µg/ml + 1µM",
    "Flagellin_1µg/ml",
    "Flagellin_0.1µg/ml",
]
control = ["DMSO 0.1%_0", "Media only_0"]


def feature_importance_linear_model(
    df: pd.DataFrame, treatment_list: list, control_list: list
):
    dosage_treatments_list = []
    for i in treatment_list:
        dosage_treatments_list.append(i)
    for i in control_list:
        dosage_treatments_list.append(i)
    df = df.query("Metadata_Treatment_and_Dose == @dosage_treatments_list")
    # Add dummy matrix of categorical genotypes
    treatment = pd.get_dummies(data=df.Metadata_Treatment_and_Dose)
    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]
    X = df.loc[:, variables]

    X = pd.concat([X, treatment], axis=1)

    columns_list = []
    columns_list.append("feature")
    columns_list.append("r2score")

    for i in X:
        columns_list.append(i)

    # Fit linear model for each feature
    lm_results = []
    for cp_feature in cp_features:
        # Subset CP data to each individual feature (univariate test)
        cp_subset_df = df.loc[:, cp_feature]

        # Fit linear model
        lm = LinearRegression(fit_intercept=True)
        lm_result = lm.fit(X=X, y=cp_subset_df)

        # Extract Beta coefficients
        # (contribution of feature to X covariates)
        coef = lm_result.coef_

        # Estimate fit (R^2)
        r2_score = lm.score(X=X, y=cp_subset_df)

        # Add results to a growing list
        lm_results.append([cp_feature, r2_score] + list(coef))

        # Convert results to a pandas DataFrame
    lm_results = pd.DataFrame(lm_results, columns=columns_list)

    # Output file
    # lm_results.to_csv(output_cp_file, sep="\t", index=False)

    print(lm_results.shape)
    lm_results.head()
    return lm_results


lm_results = feature_importance_linear_model(df, treatment, control)


# In[6]:


def plot_lm(lm_df: pd.DataFrame, x: str, y: str, fill: str):
    """Plot linear model beta coefficients

    Parameters
    ----------
    lm_df : pd.DataFrame
        Data Frame of outputed linear model values
    x : str
        x-axis column name
    y : str
        y-axis column name
    fill : str
        data-point fill column name
    """
    fig = px.scatter(lm_df, x=x, y=y, color=fill, title=f"Linear Model of {y}")
    # fig.update_yaxes(range = [-0.01,0.5])
    # fig.update_xaxes(range = [-2,2])
    fig.show()


# ##### Here I plot the beta coeifccents for each treatment against the number of cells per well. Data points the drift heavily in the Y axis are features that are affected the most by the y-axis treatment while data points that drift more in the x-axis are features that are most affected by the number of cells in a well.

# In[7]:


for i in [
    "LPS_10µg/ml",
    "Disulfiram_2.5µM",
    "LPS_1µg/ml",
    "H2O2_500µM",
    "DMSO 0.1%_0",
    "Media only_0",
]:
    plot_lm(lm_results, "Metadata_number_of_singlecells", i, "feature")
