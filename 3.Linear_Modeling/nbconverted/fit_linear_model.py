#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

# import Union
from typing import Union

import pandas as pd
import pdfkit
import plotly.express as px
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

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
    file_name: Union[bool, str] = False,
    x_label: str = False,
    y_label: str = False,
):
    """Plot linear model beta coefficients

    Parameters
    ----------
    lm_df : pd.DataFrame
        Data Frame of outputted linear model values
    x : str
        x-axis column name
    y : str
        y-axis column name
    fill : str
        data-point fill column name
    title : str
        title of graph
    file_name : Union[bool, str], optional
        substring of title of file, by default False
    x_label : str, optional
         x-axis label, by default False thus no export of file be done
    y_label : str, optional
        y-axis label, by default False

    Returns
    -------
    plotly.graph_objects.Figure
    """

    fig = px.scatter(lm_df, x=x, y=y, color=fill, title=title)
    fig.update_layout(
        xaxis_title=x_label,
        yaxis_title=y_label,
    )
    if file_name == False:
        pass
    else:
        figure_file_path = f"./figures/{file_name}"
        fig.write_html(pathlib.Path(f"{figure_file_path}.html"))
        fig.write_image(pathlib.Path(f"{figure_file_path}.png"))
        # fig.show()
    return fig


# In[3]:


# Define inputs and outputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz"
)
feature_df = pd.read_csv(feature_file, engine="pyarrow")


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


# ##### Here I plot the beta coefficients for each treatment against the number of cells per well. Data points the drift heavily in the Y axis are features that are affected the most by the y-axis treatment while data points that drift more in the x-axis are features that are most affected by the number of cells in a well.

# #### Simple Linear Modeling (cell count beta + 1 beta approach)
# Here I merged the treatment and dosage and used DMSO 0.1% as the control simply comparing one dosage/treatment at a time and outputting each graph for each treatment for all features. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count,
# $\beta _{1}$ is the beta coefficient attributed to treatment & dose combined.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[9]:


control = ["DMSO 0.1%_0"]
# lm_results_df_all = pd.DataFrame(cp_features, columns=["feature"])
# lm_results_df_all
lm_results = []
for i in feature_df["Metadata_Treatment_and_Dose"].unique():
    treatment = []
    treatment.append(i)

    dosage_treatments_list = treatment + control
    print(dosage_treatments_list)
    df = feature_df.query("Metadata_Treatment_and_Dose == @dosage_treatments_list")

    df["Metadata_Treatment_and_Dose"] = LabelEncoder().fit_transform(
        df["Metadata_Treatment_and_Dose"]
    )

    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]
    X = df.loc[:, variables]
    X = pd.concat([X, df["Metadata_Treatment_and_Dose"]], axis=1)

    columns_list = (
        ["feature", "r2_score"] + X.columns.tolist() + ["dosage_treatments_list"]
    )
    # Fit linear model for each feature
    for cp_feature in cp_features:
        # Subset CP data to each individual feature (univariate test)
        cp_subset_df = df.loc[:, cp_feature]

        # Fit linear model
        lm = LinearRegression(fit_intercept=True)
        lm_result = lm.fit(X=X, y=cp_subset_df)

        # Extract Beta coefficients(contribution of feature to X covariates)
        coef = list(lm_result.coef_)
        # Estimate fit (R^2)
        r2_score = lm.score(X=X, y=cp_subset_df)

        # Add results to a growing list
        lm_results.append(
            [cp_feature, r2_score] + coef + [f"{'_'.join(dosage_treatments_list)}"]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)
lm_results_df.head(100)


# define output file path
simple_model_output_file_path = pathlib.Path(
    f'./results/lm_cp_features_all_treatments_against_{"_".join(variables)}.tsv'
)

# write output to file
lm_results_df.to_csv(simple_model_output_file_path, sep="\t", index=False)


# In[10]:


# for loop to graph all the lm results
control = ["DMSO 0.1%_0"]
for i in feature_df["Metadata_Treatment_and_Dose"].unique():
    treatment = []
    treatment.append(i)
    dosage_treatments_list = treatment + control
    df = lm_results_df.query("Metadata_Treatment_and_Dose == @dosage_treatments_list")
    # plot results of lm
    plot_lm(
        lm_df=lm_results_df,
        x=f"{''.join(variables)}",
        y=f"Metadata_Treatment_and_Dose",
        fill="feature",
        title=f"Linear Model of {i}",
        file_name=f"lm_of_{i}_beta_coefficient_against_{''.join(variables)}_beta_coefficient",
        x_label=f"{''.join(variables)}",
        y_label=f"{i}",
    )


# #### Complex Linear Modeling (cell count btea + 2 beta approach)
# Here I run the same analysis as above but with dosage of a treatment being a factor in the linear model. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count.
# $\beta _{1}$ is the beta coefficient attributed to treatment.
# $\beta _{2}$ is the beta coefficient attributed to dose.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[11]:


# Loop for each treatment then each feature

# define the control and treatment
lm_results = []
control = ["DMSO 0.1%"]

# loop through each treatment
for i in feature_df["Metadata_treatment"].unique():
    treatment = []
    treatment.append(i)

    dosage_treatments_list = treatment + control

    # query for the treatment
    df = feature_df.query("Metadata_treatment == @dosage_treatments_list")
    # Add dummy matrix of categorical genotypes
    treatment_df = feature_df[["Metadata_treatment", "Metadata_dose"]]
    tmp_df = df.loc[:, ("Metadata_treatment", "Metadata_dose")]
    tmp_df["Metadata_treatment"] = LabelEncoder().fit_transform(
        tmp_df["Metadata_treatment"]
    )
    tmp_df["Metadata_dose"] = LabelEncoder().fit_transform(tmp_df["Metadata_dose"])

    # Setup linear modeling framework
    variables = ["Metadata_number_of_singlecells"]

    X = pd.concat([df.loc[:, variables], tmp_df], axis=1)
    columns_list = ["feature", "r2_score"] + X.columns.tolist() + ["treatment", "dose"]

    # Fit linear model for each feature
    # lm_results = []
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
        lm_results.append(
            [cp_feature, r2_score]
            + coef
            + [treatment[0], "_".join(df["Metadata_dose"].unique())]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# define output file path
complex_model_output_file_path = pathlib.Path(
    f'./results/lm_cp_features_all_treatments_and_doses_against_{"_".join(variables)}.tsv'
)
# write output to file
lm_results_df.to_csv(complex_model_output_file_path, sep="\t", index=False)


# In[12]:


# for loop to graph all the lm results
for i in feature_df["Metadata_treatment"].unique():
    treatment = []
    treatment.append(i)
    dose = lm_results_df["dose"].loc[lm_results_df.index[0]]
    df = lm_results_df.query("treatment == @treatment")

    # plot treatment lm
    fig1 = plot_lm(
        lm_df=df,
        x=variables,
        y="Metadata_treatment",
        fill="feature",
        title=f"lm of all doses of {'_'.join(treatment)}",
        x_label=f"{'_'.join(variables)}",
        y_label=f"{'_'.join(treatment)}",
    )
    fig1.show()

    # plot dose lm
    fig2 = plot_lm(
        lm_df=df,
        x=variables,
        y="Metadata_dose",
        fill="feature",
        title=f"lm of {'_'.join(treatment)} for all {'/'.join(dose)}",
        x_label=f"{' '.join(variables)}",
        y_label=f"Dose of {treatment[0]}",
    )
    fig2.show()

    # write both treatment and dose lm to html file
    filename = f"figures/lm_two_beta_{'_'.join(treatment)}"
    with open(f"{filename}.html", "a") as f:
        f.write(fig1.to_html(full_html=False, include_plotlyjs="cdn"))
        f.write(fig2.to_html(full_html=False, include_plotlyjs="cdn"))
    f.close()

    # convert HTML file to PDF
    pdfkit.from_file(f"{filename}.html", f"{filename}.pdf")
