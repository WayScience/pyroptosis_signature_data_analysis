#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pathlib
import sys

import numpy as np
import pandas as pd
import pdfkit
import plotly.express as px
import pyarrow as pa
import pyarrow.parquet as pq
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

# import Union


sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

pd.set_option("mode.chained_assignment", None)


# In[ ]:


# Define inputs
feature_file = pathlib.Path(
    "../../Extracted_Features_(CSV_files)/feature_df_sc_norm.parquet"
)
feature_df = pq.read_table(feature_file).to_pandas()


# In[ ]:


# define output file path
one_beta_output_file_path = pathlib.Path("./results/lm_one_beta.tsv")
two_beta_output_file_path = pathlib.Path("./results/lm_two_beta.tsv")
three_beta_output_file_path = pathlib.Path("./results/lm_three_beta.tsv")
four_beta_output_file_path = pathlib.Path("./results/lm_four_beta.tsv")


# In[ ]:


cp_features = infer_cp_features(feature_df)
print(f"We are testing {len(cp_features)} CellProfiler features")

new_line = "\n"
print(
    f"The unique Treatment-Dosages are: {f', {new_line}'.join((feature_df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()))}"
)


# ##### Here I plot the beta coefficients for each treatment against the number of cells per well. Data points the drift heavily in the Y axis are features that are affected the most by the y-axis treatment while data points that drift more in the x-axis are features that are most affected by the number of cells in a well.

# #### Simple Linear Modeling (cell count beta + 1 beta approach)
# Here I merged the treatment and dosage and used DMSO 0.1% as the control simply comparing one dosage/treatment at a time and outputting each graph for each treatment for all features. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count.
# $\beta _{1}$ is the beta coefficient attributed to Inducer, Inhibitor,Inhibitor Dose and, Inducer dose.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[ ]:


model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO_0.100_DMSO_0.025"
lm_results = []
for treatment in feature_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():
    dosage_treatments_list = [treatment, control]
    print(dosage_treatments_list)
    # filter df for treatment and dose
    df = feature_df.query(
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose in @dosage_treatments_list"
    )
    # encode treatment and dose as integers

    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = LabelEncoder().fit_transform(
        df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
    )

    # Setup linear modeling framework

    X = df.loc[:, model_covariates]
    X = pd.concat([X, df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]], axis=1)

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
lm_results
# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# write output to file
lm_results_df.to_csv(one_beta_output_file_path, sep="\t", index=False)


# #### Complex Linear Modeling (cell count btea + 2 beta approach)
# Here I run the same analysis as above but with dosage of a treatment being a factor in the linear model. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count.
# $\beta _{1}$ is the beta coefficient attributed to Inducer, Inhibitor, and Inhibitor Dose.
# $\beta _{2}$ is the beta coefficient attributed to Inducer dose.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[ ]:


# Loop for each treatment then each feature

# define the control and treatment
# Setup linear modeling framework
model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO_DMSO_0.025__0.100"
lm_results = []
for treatment in feature_df["twob_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():
    dosage_treatments_list = [treatment, control]
    print(dosage_treatments_list)
    df = feature_df.query(
        "twob_Metadata_Treatment_Dose_Inhibitor_Dose in @dosage_treatments_list"
    )
    # Add dummy matrix of categorical genotypes
    # treatment_df = feature_df[["Metadata_inducer1", "Metadata_inducer1_concentration"]]
    df[["twob_Metadata_Treatment_Inhibitor_Dose", "Treatment_Dose"]] = df[
        "twob_Metadata_Treatment_Dose_Inhibitor_Dose"
    ].str.split("__", expand=True)
    tmp_df = df.loc[
        :,
        (
            "twob_Metadata_Treatment_Inhibitor_Dose",
            "Treatment_Dose",
        ),
    ]

    tmp_df["twob_Metadata_Treatment_Inhibitor_Dose"] = LabelEncoder().fit_transform(
        tmp_df["twob_Metadata_Treatment_Inhibitor_Dose"]
    )
    tmp_df["Treatment_Dose"] = LabelEncoder().fit_transform(tmp_df["Treatment_Dose"])

    X = pd.concat([df.loc[:, model_covariates], tmp_df], axis=1)
    columns_list = (
        ["feature", "r2_score"]
        + X.columns.tolist()
        + [
            "inducer1",
            "inducer1_dose",
            "inducer2",
            "inducer2_dose",
            "inhibitor",
            "inhibitor_dose",
        ]
    )

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
            + [
                treatment,
                "_".join(df["Metadata_inducer1_concentration"].astype(str).unique()),
            ]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# write output to file
lm_results_df.to_csv(two_beta_output_file_path, sep="\t", index=False)


# #### Complex Linear Modeling (cell count beta + 3 beta approach)
# Here I run the same analysis as above but with dosage of a treatment being a factor in the linear model. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \beta _{3}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count.
# $\beta _{1}$ is the beta coefficient attributed to Inducer.
# $\beta _{2}$ is the beta coefficient attributed to Inducer dose.
# $\beta _{3}$ is the beta coefficient attributed to Inhibitor, and Inhibitor Dose.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[ ]:


# Loop for each treatment then each feature

# define the control and treatment
# Setup linear modeling framework
model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO__0.100__DMSO_0.025"
lm_results = []
for treatment in feature_df["threeb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():
    dosage_treatments_list = [treatment, control]
    print(dosage_treatments_list)
    df = feature_df.query(
        "threeb_Metadata_Treatment_Dose_Inhibitor_Dose in @dosage_treatments_list"
    )
    # Add dummy matrix of categorical genotypes
    df[
        [
            "threeb_Treatment",
            "threeb_Treatment_Dose",
            "threeb_Inhibitor_and_Dose",
        ]
    ] = df["threeb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.split("__", expand=True)
    tmp_df = df.loc[
        :, ("threeb_Treatment", "threeb_Treatment_Dose", "threeb_Inhibitor_and_Dose")
    ]

    tmp_df["threeb_Treatment"] = LabelEncoder().fit_transform(
        tmp_df["threeb_Treatment"]
    )
    tmp_df["threeb_Treatment_Dose"] = LabelEncoder().fit_transform(
        tmp_df["threeb_Treatment_Dose"]
    )
    tmp_df["threeb_Inhibitor_and_Dose"] = LabelEncoder().fit_transform(
        tmp_df["threeb_Inhibitor_and_Dose"]
    )

    X = pd.concat([df.loc[:, model_covariates], tmp_df], axis=1)
    columns_list = (
        ["feature", "r2_score"]
        + X.columns.tolist()
        + [
            "inducer1",
            "inducer1_dose",
            "inducer2",
            "inducer2_dose",
            "inhibitor",
            "inhibitor_dose",
        ]
    )

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
            + [
                treatment,
                "_".join(df["Metadata_inducer1_concentration"].astype(str).unique()),
            ]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# write output to file
lm_results_df.to_csv(three_beta_output_file_path, sep="\t", index=False)


# #### Complex Linear Modeling (cell count beta + 4 beta approach)
# Here I run the same analysis as above but with dosage of a treatment being a factor in the linear model. All features and treatments will be exported into 1 file.
#
# Linear Model:
# $y = \beta _{0}x+ \beta _{1}x+ \beta _{2}x+ \beta _{3}x+ \beta _{4}x+ \epsilon$ where;
# $y$ is each feature
# $x$ is the inputed variables
# $\beta _{0}$ is the beta coefficient attributed to cell count.
# $\beta _{1}$ is the beta coefficient attributed to Inducer.
# $\beta _{2}$ is the beta coefficient attributed to Inducer dose.
# $\beta _{3}$ is the beta coefficient attributed to Inhibitor.
# $\beta _{4}$ is the beta coefficient attributed to Inhibitor Dose.
# $\epsilon$ is the residual variance not explained by factors in the model

# In[ ]:


# Loop for each treatment then each feature

# define the control and treatment
# Setup linear modeling framework
model_covariates = ["Metadata_number_of_singlecells"]
control = "DMSO__0.100__DMSO__0.025"
lm_results = []
for treatment in feature_df["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique():
    dosage_treatments_list = [treatment, control]
    print(dosage_treatments_list)
    df = feature_df.query(
        "fourb_Metadata_Treatment_Dose_Inhibitor_Dose in @dosage_treatments_list"
    )
    # Add dummy matrix of categorical genotypes
    df[
        [
            "fourb_Treatment",
            "fourb_Treatment_Dose",
            "fourb_Inhibitor",
            "fourb_Inhibitor_Dose",
        ]
    ] = df["fourb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.split("__", expand=True)
    tmp_df = df.loc[
        :,
        (
            "fourb_Treatment",
            "fourb_Treatment_Dose",
            "fourb_Inhibitor",
            "fourb_Inhibitor_Dose",
        ),
    ]

    tmp_df["fourb_Treatment"] = LabelEncoder().fit_transform(tmp_df["fourb_Treatment"])
    tmp_df["fourb_Treatment_Dose"] = LabelEncoder().fit_transform(
        tmp_df["fourb_Treatment_Dose"]
    )
    tmp_df["fourb_Inhibitor"] = LabelEncoder().fit_transform(tmp_df["fourb_Inhibitor"])
    tmp_df["fourb_Inhibitor_Dose"] = LabelEncoder().fit_transform(
        tmp_df["fourb_Inhibitor_Dose"]
    )

    X = pd.concat([df.loc[:, model_covariates], tmp_df], axis=1)
    columns_list = (
        ["feature", "r2_score"]
        + X.columns.tolist()
        + [
            "inducer1",
            "inducer1_dose",
            "inducer2",
            "inducer2_dose",
            "inhibitor",
            "inhibitor_dose",
        ]
    )

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
            + [
                treatment,
                "_".join(df["Metadata_inducer1_concentration"].astype(str).unique()),
            ]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# write output to file
lm_results_df.to_csv(four_beta_output_file_path, sep="\t", index=False)


# In[ ]:
