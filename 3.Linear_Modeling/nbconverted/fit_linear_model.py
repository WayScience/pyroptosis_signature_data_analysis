#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import pandas as pd
import plotly.express as px
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression

sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

# In[2]:


def linear_model_function(control, treatment, df, cp_features):
    print(control[0])
    print(treatment[0])
    dosage_treatments_list = treatment + control

    df = df.query("Metadata_Treatment_and_Dose == @dosage_treatments_list")
    # Add dummy matrix of categorical genotypes
    treatment_df = pd.get_dummies(data=df.Metadata_Treatment_and_Dose)
    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]
    X = df.loc[:, variables]

    X = pd.concat([X, treatment_df], axis=1)

    columns_list = []
    columns_list.append("feature")
    columns_list.append("r2_score")

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
        coef = list(lm_result.coef_)
        # Estimate fit (R^2)
        r2_score = lm.score(X=X, y=cp_subset_df)

        # Add results to a growing list
        lm_results.append([cp_feature, r2_score] + coef)

    # Convert results to a pandas DataFrame
    lm_results = pd.DataFrame(lm_results, columns=columns_list)

    # Output file
    output_dir = pathlib.Path("results")
    output_cp_file = pathlib.Path(
        output_dir,
        f"Linear_Model_cp_features_{treatment}_beta_coefficient_against_{control}_beta_coefficient.tsv",
    )
    print(output_cp_file)
    lm_results.to_csv(output_cp_file, sep="\t", index=False)

    # lm_results.head()
    return lm_results


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
    figure_file_path = pathlib.Path(
        f"./figures/Linear_Model_{y}_beta_coefficient_against_{x}_beta_coefficient.html"
    )
    fig.write_html(figure_file_path)
    fig.show()


# In[3]:


# Define inputs and outputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz"
)
feature_df = pd.read_csv(feature_file, low_memory=False)
output_dir = pathlib.Path("results")

output_cp_file = pathlib.Path(output_dir, "linear_model_cp_features.tsv")


# In[4]:


pd.set_option("display.max_columns", 5000)


# In[5]:


print(feature_df.shape)
feature_df.head()


# In[6]:


feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df.head()


# In[7]:


feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)
feature_df.head()


# In[8]:


# Recycled code from: https://github.com/WayScience/NF1_SchwannCell_data/blob/main/5_analyze_data/notebooks/linear_model/fit_linear_model.ipynb
cell_count_df = (
    feature_df.groupby("Metadata_Well")["Metadata_Plate"]
    .count()
    .reset_index()
    .rename(columns={"Metadata_Plate": "Metadata_number_of_singlecells"})
)

feature_df = feature_df.merge(cell_count_df, on="Metadata_Well")

# Drop na and reindex accordingly
feature_df = (
    feature_df.dropna()
    .reset_index()
    .assign(
        Metadata_Treatment_and_Dose=lambda x: feature_df["Metadata_treatment"]
        + "_"
        + feature_df["Metadata_dose"]
    )
)

cp_features = infer_cp_features(feature_df)
print(f"We are testing {len(cp_features)} CellProfiler features")

feature_df["Metadata_Treatment_and_Dose"].unique()


# ##### Here I plot the beta coeifccents for each treatment against the number of cells per well. Data points the drift heavily in the Y axis are features that are affected the most by the y-axis treatment while data points that drift more in the x-axis are features that are most affected by the number of cells in a well.

# In[9]:


feature_df["Metadata_Treatment_and_Dose"].unique()


# In[10]:


for i in feature_df["Metadata_Treatment_and_Dose"].unique():
    treatment = []
    control = ["DMSO 0.1%_0"]
    treatment.append(i)
    lm_results = linear_model_function(control, treatment, feature_df, cp_features)
    plot_lm(lm_results, "Metadata_number_of_singlecells", i, "feature")


# In[ ]:
