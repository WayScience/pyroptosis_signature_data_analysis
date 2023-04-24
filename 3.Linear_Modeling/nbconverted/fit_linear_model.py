#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import pandas as pd
import pdfkit
import plotly.express as px
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

# import Union


sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

# In[2]:


# Define inputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz"
)
feature_df = pd.read_csv(feature_file, engine="pyarrow")


# In[3]:


# define output file path
simple_model_output_file_path = pathlib.Path("./results/lm_one_beta")

complex_model_output_file_path = pathlib.Path("./results/lm_two_beta")

feature_df_out_path = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/feature_df_sc_norm_fs.csv"
)

filename = pathlib.Path("figures/lm_one_beta")
filename = pathlib.Path("figures/lm_two_beta")


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

new_line = "\n"
print(
    f"The unique Treatment-Dosages are: {f', {new_line}'.join((feature_df['Metadata_Treatment_and_Dose'].unique()))}"
)


# In[9]:


feature_df.to_csv(feature_df_out_path, index=False)


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

# In[ ]:


pd.set_option("mode.chained_assignment", None)

model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO 0.1%_0"
lm_results = []
for treatment in feature_df["Metadata_Treatment_and_Dose"].unique():
    dosage_treatments_list = [treatment, control]
    # filter df for treatment and dose
    df = feature_df.query("Metadata_Treatment_and_Dose in @dosage_treatments_list")
    # encode treatment and dose as integers
    df["Metadata_Treatment_and_Dose"] = LabelEncoder().fit_transform(
        df["Metadata_Treatment_and_Dose"]
    )

    # Setup linear modeling framework

    X = df.loc[:, model_covariates]
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
            [cp_feature, r2_score] + coef + [f"{'-'.join(dosage_treatments_list)}"]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

simple_model_output_file_path = (
    f'{simple_model_output_file_path}_{"_".join(model_covariates)}.tsv'
)

# write output to file
lm_results_df.to_csv(simple_model_output_file_path, sep="\t", index=False)


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

# In[ ]:


# Loop for each treatment then each feature

# define the control and treatment
# Setup linear modeling framework
model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO 0.1%"
lm_results = []
for treatment in feature_df["Metadata_treatment"].unique():
    dosage_treatments_list = [treatment, control]
    df = feature_df.query("Metadata_treatment in @dosage_treatments_list")
    # Add dummy matrix of categorical genotypes
    treatment_df = feature_df[["Metadata_treatment", "Metadata_dose"]]
    tmp_df = df.loc[:, ("Metadata_treatment", "Metadata_dose")]
    tmp_df["Metadata_treatment"] = LabelEncoder().fit_transform(
        tmp_df["Metadata_treatment"]
    )
    tmp_df["Metadata_dose"] = LabelEncoder().fit_transform(tmp_df["Metadata_dose"])

    X = pd.concat([df.loc[:, model_covariates], tmp_df], axis=1)
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
            + [treatment, "_".join(df["Metadata_dose"].unique())]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# define output file path
complex_model_output_file_path = (
    f'{complex_model_output_file_path}_{"_".join(model_covariates)}.tsv'
)

# write output to file
lm_results_df.to_csv(complex_model_output_file_path, sep="\t", index=False)
