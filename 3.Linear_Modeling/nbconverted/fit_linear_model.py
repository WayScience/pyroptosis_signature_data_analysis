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


def plot_lm(
    lm_df: pd.DataFrame,
    x: str,
    y: str,
    fill: str,
    title: str,
    file_name: str,
    x_label: str = False,
    y_label: str = False,
):
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
    title: str
        title of graph
    file_name : str
        substring of title of file
    x_label : str
        x-axis label
    y_label : str
        y-axis label
    """

    fig = px.scatter(lm_df, x=x, y=y, color=fill, title=title)
    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title=y_label,
    )
    # fig.update_yaxes(range = [-0.01,0.5])
    # fig.update_xaxes(range = [-2,2])
    figure_file_path = f"./figures/{file_name}"
    fig.write_html(pathlib.Path(f"{figure_file_path}.html"))
    fig.write_image(pathlib.Path(f"{figure_file_path}.png"))
    fig.show()


# In[3]:


# Define inputs and outputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_fs_cellprofiler.csv.gz"
)
feature_df = pd.read_csv(feature_file, engine="pyarrow")
output_dir = pathlib.Path("results")

output_cp_file = pathlib.Path(output_dir, "linear_model_cp_features.tsv")


# In[4]:


print(feature_df.shape)
feature_df.head()


# In[5]:


# removing costes features as they behave with great variance across all data
feature_df = feature_df.drop(feature_df.filter(regex="Costes").columns, axis=1)
print(feature_df.shape)
feature_df.head()


# In[6]:


# replacing '/' in treatment dosage column to avoid errors in file interpolation including such strings
feature_df = feature_df.replace(to_replace="/", value="_per_", regex=True)
feature_df.head()


# In[7]:


# Recycled code from: https://github.com/WayScience/NF1_SchwannCell_data/blob/main/5_analyze_data/notebooks/linear_model/fit_linear_model.ipynb
cell_count_df = (
    feature_df.groupby("Metadata_Well")["Metadata_Plate"]
    .count()
    .reset_index()
    .rename(columns={"Metadata_Plate": "Metadata_number_of_singlecells"})
)


# In[8]:


feature_df = feature_df.merge(cell_count_df, on="Metadata_Well")

# Drop na and reindex accordingly
feature_df = (
    feature_df.assign(
        Metadata_Treatment_and_Dose=lambda x: feature_df["Metadata_treatment"]
        + "_"
        + feature_df["Metadata_dose"]
    )
    .dropna()
    .reset_index()
)

cp_features = infer_cp_features(feature_df)
print(f"We are testing {len(cp_features)} CellProfiler features")

feature_df["Metadata_Treatment_and_Dose"].unique()


# ##### Here I plot the beta coeifccents for each treatment against the number of cells per well. Data points the drift heavily in the Y axis are features that are affected the most by the y-axis treatment while data points that drift more in the x-axis are features that are most affected by the number of cells in a well.

# #### Simple Linear Modeling
# Here I merged the treatment and dosage and used DMSO 0.1% as the control simply comparing one dosage/treatment at a time and outputting each graph for each treatment for all features. All features and treatments will be exported into 1 file which will be the `simple model file`

# In[9]:


control = ["DMSO 0.1%_0"]
lm_results_df_all = pd.DataFrame(cp_features, columns=["feature"])
lm_results_df_all
for i in feature_df["Metadata_Treatment_and_Dose"].unique():
    treatment = []
    treatment.append(i)

    dosage_treatments_list = treatment + control

    df = feature_df.query("Metadata_Treatment_and_Dose == @dosage_treatments_list")

    treatment_df = pd.get_dummies(data=df.Metadata_Treatment_and_Dose)
    treatment_df
    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]
    X = df.loc[:, variables]
    X = pd.concat([X, treatment_df], axis=1)

    columns_list = ["feature", "r2_score"] + X.columns.tolist()

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
    lm_results_df = pd.DataFrame(lm_results, columns=columns_list)
    lm_results_df
    plot_lm(
        lm_df=lm_results_df,
        x=f"{''.join(variables)}",
        y=f"{i}",
        fill="feature",
        title=f"Linear Model of {i}",
        file_name=f"lm_of_{i}_beta_coefficient_against_{variables}_beta_coefficient",
        x_label=f"{''.join(variables)}",
        y_label=f"{i}",
    )
    x = []
    for i in lm_results_df.columns:
        x.append(f'{i}_{"_+_".join(dosage_treatments_list)}')
    lm_results_df.columns = x
    lm_results_df_all = pd.concat([lm_results_df_all, lm_results_df], axis=1)
simple_model_output_file_path = pathlib.Path(
    f'./results/lm_cp_features_all_treatments_against_{"_".join(variables)}.tsv'
)
lm_results_df_all.to_csv(simple_model_output_file_path, sep="\t", index=False)


# #### Complex Linear Model
# Here I run the same analysis as above but with dosage of a treatment being a factor in the linear model

# In[10]:


lm_results_df_all = pd.DataFrame(cp_features, columns=["feature"])
control = ["DMSO 0.1%"]
for i in feature_df["Metadata_treatment"].unique():

    treatment = []
    treatment.append(i)

    dosage_treatments_list = treatment + control
    print(dosage_treatments_list)

    df = feature_df.query("Metadata_treatment == @dosage_treatments_list")
    # Add dummy matrix of categorical genotypes
    treatment_df = feature_df[["Metadata_treatment", "Metadata_dose"]]
    tmp_df = df[["Metadata_treatment", "Metadata_dose"]]
    treatment_df = pd.get_dummies(
        data=tmp_df, columns=["Metadata_treatment", "Metadata_dose"]
    )

    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]

    X = pd.concat([df.loc[:, variables], treatment_df], axis=1)

    columns_list = ["feature", "r2_score"] + X.columns.tolist()

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
    lm_results_df = pd.DataFrame(lm_results, columns=columns_list)
    for i in lm_results_df.iloc[:, 3:]:
        print(i)
        graph_title = f"lm_of_{i}_beta_coefficient_against_{variables}_beta_coefficient"
        plot_lm(
            lm_df=lm_results_df,
            x=variables,
            y=f"{i}",
            fill="feature",
            title=f"lm of all doses of {'_'.join(dosage_treatments_list)}",
            file_name=f"lm_of_all_doses_{'_'.join(dosage_treatments_list)}_and_{i}_beta_coefficient_against_{'_'.join(variables)}_beta_coefficient",
            x_label=f"{'_'.join(variables)}",
            y_label=f"{'_'.join(treatment)} {i}",
        )
    x = []
    for i in lm_results_df.columns:
        x.append(f'{i}_{"_+_".join(dosage_treatments_list)}')
    lm_results_df.columns = x
    lm_results_df_all = pd.concat([lm_results_df_all, lm_results_df], axis=1)
simple_model_output_file_path = pathlib.Path(
    f'./results/LM_cp_features_all_doses_and_treatments__against_{"_".join(variables)}.tsv'
)
lm_results_df_all.to_csv(simple_model_output_file_path, sep="\t", index=False)
