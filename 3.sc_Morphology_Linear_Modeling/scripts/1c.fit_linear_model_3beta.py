#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import sys

import numpy as np
import pandas as pd
import plotly.express as px
import pyarrow as pa
import pyarrow.parquet as pq
import seaborn as sns
from matplotlib import pyplot as plt
from pycytominer.cyto_utils import infer_cp_features
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder

sys.path.append("..")
# from ..utils.utils import df_stats
import matplotlib.pyplot as plt

pd.set_option("mode.chained_assignment", None)


# In[2]:


# Parameters
celltype = "PBMC"


# In[3]:


# Define inputs
feature_file = pathlib.Path(f"../data/{celltype}_preprocessed_sc_norm.parquet")
feature_df = pq.read_table(feature_file).to_pandas()


# In[4]:


# if path does not exist, create one
pathlib.Path(f"./results/{celltype}").mkdir(parents=True, exist_ok=True)

# define output file path
one_beta_output_file_path = pathlib.Path(f"./results/{celltype}/lm_one_beta.tsv")
two_beta_output_file_path = pathlib.Path(f"./results/{celltype}/lm_two_beta.tsv")
three_beta_output_file_path = pathlib.Path(f"./results/{celltype}/lm_three_beta.tsv")
four_beta_output_file_path = pathlib.Path(f"./results/{celltype}/lm_four_beta.tsv")


# In[5]:


cp_features = infer_cp_features(feature_df)
print(f"We are testing {len(cp_features)} CellProfiler features")

new_line = "\n"
print(
    f"The unique Treatment-Dosages are: {f', {new_line}'.join((feature_df['threeb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()))}"
)


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

# In[6]:


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
            "inducer1__inducer1_dose__inhibitor_inhibitor_dose",
        ]
    )

    # Fit linear model for each feature
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
            ]
        )

# Convert results to a pandas DataFrame
lm_results_df = pd.DataFrame(lm_results, columns=columns_list)

# write output to file
lm_results_df.to_csv(three_beta_output_file_path, sep="\t", index=False)
