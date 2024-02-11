#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import plotly.express as px
from copairs.map import aggregate

warnings.filterwarnings("ignore")


# In[2]:


# Directories
processed_data_dir = pathlib.Path("../data/processed/")
sc_ap_scores_dir = (processed_data_dir / "mAP_scores/secretome").resolve()
agg_sc_ap_scores_dir = (processed_data_dir / "aggregate_mAPs/secretome").resolve()
agg_sc_ap_scores_dir.mkdir(parents=True, exist_ok=True)


# ## Preparing the dataset
#

# In[3]:


all_files = list(sc_ap_scores_dir.glob("*.csv"))
# get the files that contain the string class
class_files = [file for file in all_files if "treatment" in file.stem]
mAPs = []
for file in class_files:
    df = pd.read_csv(file)
    df["file"] = file.stem
    mAPs.append(df)
# single-cell mAP scores
mAPs = pd.concat(mAPs)
mAPs.head()


# In[4]:


# grabbing all cp features (regular, feature shuffled and labeled shuffled)
reg_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "non-shuffled"]
shuffled_feat_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "features_shuffled"]
shuffled_pheno_sc_mAPs = mAPs.loc[mAPs["shuffled"] == "phenotype_shuffled"]


# In[5]:


# Generating sampling_error df
# This table will be used to merge with the aggregate table to get the sampling error a specific category.
merged_sc_ap_scores_df = pd.concat(
    [
        reg_sc_mAPs,
        shuffled_feat_sc_mAPs,
        shuffled_pheno_sc_mAPs,
    ]
)

# grouping dataframe based on phenotype levels, feature and feature types
df_group = merged_sc_ap_scores_df.groupby(
    by=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled"]
)

# calculating sampling error
sampling_error_df = []
for name, df in df_group:
    pheno, shuffled_type = name

    # caclulating sampling error
    avg_percision = df["average_precision"].values
    sampling_error = np.std(avg_percision) / np.sqrt(len(avg_percision))

    sampling_error_df.append([pheno, shuffled_type, sampling_error])
cols = ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled", "sampling_error"]
sampling_error_df = pd.DataFrame(sampling_error_df, columns=cols)

# updating name:
sampling_error_df.loc[
    sampling_error_df["shuffled"] == "phenotype_shuffled"
] = "phenotypes_shuffled"

sampling_error_df.head()


# In[6]:


# aggregate single cells scores with cell labels
data = tuple(
    merged_sc_ap_scores_df.groupby(by=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"])
)
columns = merged_sc_ap_scores_df.columns
agg_sc_ap_scores_df = []
for cell_id, df1 in data:
    for shuffle_type, df2 in df1.groupby(by="shuffled"):
        aggregated_ap_score = df2["average_precision"].mean()

        # select a single row since all the metadata is the same
        selected_row = df2.iloc[0]

        # update the average precision score of the single row
        selected_row["average_precision"] = aggregated_ap_score
        agg_sc_ap_scores_df.append(selected_row.values.tolist())

# saving into the results repo
agg_sc_ap_scores_df = pd.DataFrame(data=agg_sc_ap_scores_df, columns=columns)
agg_sc_ap_scores_df.to_csv(
    sc_ap_scores_dir / "merged_sc_agg_ap_scores_treatment.csv", index=False
)
agg_sc_ap_scores_df.head()


# In[7]:


# Generating aggregate scores with a threshold p-value of 0.05
mAP_dfs = []
for name, df in tuple(
    agg_sc_ap_scores_df.groupby(
        by=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "shuffled"]
    )
):
    agg_df = aggregate(
        df, sameby=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"], threshold=0.05
    )
    agg_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = name[0]
    agg_df["shuffled"] = name[1]
    mAP_dfs.append(agg_df)

mAP_dfs = pd.concat(mAP_dfs)
mAP_dfs.to_csv(agg_sc_ap_scores_dir / "mAP_scores_treatment.csv", index=False)
mAP_dfs.head()


# ## Forming bar plots
#

# ### Forming bar plots with CP Features
#

# In[8]:


# selecting dataset to plot
agg_reg_sc_mAPs = mAP_dfs.loc[(mAP_dfs["shuffled"] == "non-shuffled")]
agg_shuffled_feat_sc_mAPs = mAP_dfs.loc[(mAP_dfs["shuffled"] == "features_shuffled")]
agg_shuffled_pheno_sc_mAPs = mAP_dfs.loc[(mAP_dfs["shuffled"] == "phenotype_shuffled")]

# phenotypes
df = (
    pd.concat(
        [
            agg_reg_sc_mAPs,
            agg_shuffled_feat_sc_mAPs,
            agg_shuffled_pheno_sc_mAPs,
        ]
    )
    .reset_index()
    .drop("index", axis=1)
)[["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "mean_average_precision", "shuffled"]]


fig = px.bar(
    df,
    x="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    y="mean_average_precision",
    color="shuffled",
    barmode="group",
    title="Mean Average Precision for each Cell Death Phenotype",
    labels={
        "mean_average_precision": "Mean Average Precision",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose": "Cell Death Phenotypes",
    },
)


# ## Generating box plots of single cell ap scores

# In[9]:


all_df = pd.concat(
    [
        reg_sc_mAPs,
        shuffled_feat_sc_mAPs,
        shuffled_pheno_sc_mAPs,
    ]
)


# In[10]:


# Assuming all_cp_df, all_dp_df, and all_cp_dp_df are your DataFrames
categories_order = all_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].unique()

# Create individual figures with the same category order
fig1 = px.box(
    all_df,
    x="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
    y="average_precision",
    color="shuffled",
    title="Single Well Average Percision Scores",
    category_orders={"oneb_Metadata_Treatment_Dose_Inhibitor_Dose": categories_order},
    labels={
        "average_precision": "Average Precision Scores",
        "oneb_Metadata_Treatment_Dose_Inhibitor_Dose": "Treatment",
    },
)
