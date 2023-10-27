#!/usr/bin/env python
# coding: utf-8

# The goal of this notebook is to determine which cytokines and chemokines are found at high levels in pyroptotic inducing agents.
# Doing this will allow us to determine ground truth of pyroptosis occurance.

# ### Imports

# In[ ]:


import pathlib

# umap analysis of treatment groups
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import toml
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
from scipy.cluster.hierarchy import linkage
from scipy.stats import f_oneway

# post hoc test for 'VEGF-C [NSU]' column using Tukey's HSD test
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# anova test on each group


warnings.filterwarnings("ignore")
warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
warnings.simplefilter("ignore", category=NumbaPendingDeprecationWarning)
import umap

# In[ ]:


# set path

df_path = pathlib.Path(
    f"../../Data/clean/Plate2/nELISA_plate_430420_PBMC_clean.parquet"
)


# read in the data
df = pd.read_parquet(df_path)


# In[ ]:


# import selected treatmenets
# set path
toml_path = pathlib.Path("../../../1.Exploratory_Data_Analysis/utils/params.toml")

# read in toml file
params = toml.load(toml_path)
list_of_treatments = params["list_of_treatments"]["treatments"]


# In[ ]:


# get the treatments in fourb_Metadata_Treatment_Dose_Inhibitor coulumn for each treatment in the list of treatments
df = df.drop(
    columns=[
        "Dose",
        "Treatment",
        "twob_Treatment_Dose_Inhibitor_Dose",
        "threeb_Treatment_Dose_Inhibitor_Dose",
        "fourb_Treatment_Dose_Inhibitor_Dose",
    ]
)
# if column name does not contain [NSU], add Metadata_ to the beginning of the column name
df.columns = ["Metadata_" + col if "[NSU]" not in col else col for col in df.columns]

df_metadata = df[df.columns[df.columns.str.contains("Metadata")]]
# non_metadata_cols
df = df.drop(columns=df_metadata.columns)
df["oneb_Treatment_Dose_Inhibitor_Dose"] = df_metadata[
    "Metadata_oneb_Treatment_Dose_Inhibitor_Dose"
]
df["Metadata_position_x"] = df_metadata["Metadata_position_x"]


# In[ ]:


# set output path
all_cytokines_path = pathlib.Path(
    f"./results/PBMC_all_cytokine_values_per_treatment_per_well.parquet"
)
all_cytokines_path_melted = pathlib.Path(
    f"./results/PBMC_all_cytokine_values_per_treatment_per_well_melted.parquet"
)
df.to_parquet(all_cytokines_path)

df_melted = df.melt(
    id_vars=["Metadata_position_x", "oneb_Treatment_Dose_Inhibitor_Dose"],
    var_name="cytokine",
    value_name="cytokine_value",
)

df_melted.to_parquet(all_cytokines_path_melted)


# ## Anova and Post-Hoc Analysis
# Anova of all treatments and post-hoc analysis of all treatments for each cytokine and chemokine.
# This will determine the cytokines and chemokines that are found at high levels in pyroptotic inducing agents.

# In[ ]:


# define blank df
final_df_tukey = pd.DataFrame(
    {
        "group1": [""],
        "group2": [""],
        "meandiff": [""],
        "lower": [""],
        "upper": [""],
        "reject": [""],
        "p-adj": [""],
        "cytokine": [""],
    }
)


# In[ ]:


# perform anova on each column of the data frame with oneb_meta as the groupby
num = 0
alpha = 0.05
alpha_adj = alpha / (len(df.columns) - 1)
for i in df.columns:
    for treatment in list_of_treatments:
        if i == "oneb_Treatment_Dose_Inhibitor_Dose":
            continue
        one_way_anova = stats.f_oneway(
            df[i][df["oneb_Treatment_Dose_Inhibitor_Dose"] == treatment],
            df[i][df["oneb_Treatment_Dose_Inhibitor_Dose"] != treatment],
        )
        if one_way_anova.pvalue < alpha:
            num += 1
            tukey = pairwise_tukeyhsd(
                endog=df[i],
                groups=df["oneb_Treatment_Dose_Inhibitor_Dose"],
                alpha=alpha_adj,
            )
            # send the results to a dataframe
            tukey_results = pd.DataFrame(
                data=tukey._results_table.data[1:], columns=tukey._results_table.data[0]
            )
            tukey_results["cytokine"] = f"{i}"
            # concat the results to the blank df
            final_df_tukey = pd.concat([final_df_tukey, tukey_results], axis=0)
        else:
            pass
print(
    f"Out of the {len(df.columns ) - 1} cytokines tested, {num} were significantly different between groups (p < {alpha})"
)


# In[ ]:


# check for blank first row...
final_df_tukey.head(3)


# In[ ]:


# remove first row as it is blank fro some reason
final_df_tukey = final_df_tukey.iloc[1:]
final_df_tukey.head(3)


# Clean up the data and filter out tests that are not significant.

# In[ ]:


# drop rows in pvalue column that are over 0.05
final_df_tukey = final_df_tukey[final_df_tukey["p-adj"] < 0.05]


# In[ ]:


# sort the df by p-adj
final_df_tukey = final_df_tukey.sort_values(by=["p-adj"], ascending=[True])


# # filter the data for significanct post hoc tests
# If we see two high dose groups of pyroptotic treatments in this p-adj value < 0.05 data then we can toss it.
# This implies a variable treatment.
# We are primarily interested in which cytokines best differentiate between control, apoptosis, and pyroptosis

# In[ ]:


final_df_tukey["cytokine"].unique()
# create output path for the df
output_path = pathlib.Path(f"./results/tukey_filtered_nomic_results.csv")
# save the df
final_df_tukey.to_csv(output_path)


# In[ ]:


# graph each cytokine
for col in final_df_tukey["cytokine"].unique():
    sns.barplot(
        x="oneb_Treatment_Dose_Inhibitor_Dose",
        y=col,
        data=df,
        capsize=0.2,
        order=list_of_treatments,
    )
    plt.title(col)
    plt.xticks(rotation=90)
    plt.show()
# feature pick
cytokines = [
    "Activin A [NSU]",
    "IL-1 alpha [NSU]",
    "IL-1 beta [NSU]",
    "Oncostatin M (OSM) [NSU]",
    "IFN gamma [NSU]",
    "Osteopontin (OPN) [NSU]",
    "TNF alpha [NSU]",
    "EMMPRIN [NSU]",
    "G-CSF [NSU]",
    "MMP-9 [NSU]",
    "IL-6 [NSU]",
    "MIF [NSU]",
    "IL-16 [NSU]",
    "IL-22 [NSU]",
    "IL-18 [NSU]",
    "CCL24 [NSU]",
    "CCL20 [NSU]",
    "CXCL11 [NSU]",
    "CXCL1 [NSU]",
]


# In[ ]:


# drop all columns that are not in cytokines list
selected_cytokines = df[cytokines]

# plot the results of the tukey test for each cytokine
a = len(selected_cytokines.columns)
b = 6
plt.figure(figsize=(50, 100))
plt.suptitle("Cytokine Levels Across Treatments", fontsize=18)
plt.subplots_adjust(top=0.975, bottom=0.01, hspace=1, wspace=0.3)
for col in enumerate(selected_cytokines.columns):
    plt.subplot(a, b, col[0] + 1)
    sns.barplot(
        x="oneb_Treatment_Dose_Inhibitor_Dose",
        y=col[1],
        data=df,
        capsize=0.2,
        order=list_of_treatments,
    )
    # # title
    plt.title(col[1])
    # rotate xticks 90 degrees
    plt.xticks(rotation=90)
# set path for saving plot
pathlib.Path(f"./figures/").mkdir(parents=True, exist_ok=True)
# save plot
plt.savefig(f"./figures/selected_cytokines.png", bbox_inches="tight")
# # show plot
plt.show()


# In[ ]:


# save the final_df_tukey df to a csv file
final_df_tukey.to_csv("results/tukey_test_results.csv", index=False)

# write the cytokines column to a csv file
cytokines
with open("results/cytokines.csv", "w") as f:
    f.write("cytokine\n")
    for item in cytokines:
        f.write(f"{item}\n")
    f.close()


# ## Heatmaps of cytokine levels in each treatment

# In[ ]:


df_cytokines = df[cytokines]
df_cytokines = pd.concat(
    [df["oneb_Treatment_Dose_Inhibitor_Dose"], df_cytokines], axis=1
)
df_cytokines = df_cytokines.set_index("oneb_Treatment_Dose_Inhibitor_Dose")


# In[ ]:


cytokines


# In[ ]:


# aggregate the data by treatment group via mean
data_agg = df_cytokines.groupby("oneb_Treatment_Dose_Inhibitor_Dose").mean()
# heatmap of umap_clusters_with_cytokine_data_agg
# subset the columns to plot
column_list = [col for col in data_agg.columns if "[NSU]" in col]
# subset the rows to plot and label the rows with treatment groups
row_list = data_agg.index
# subset the data to plot
data = data_agg[column_list]


# In[ ]:


# order the rows by treatment group
data_agg = data_agg.reindex(list_of_treatments, axis=0)


# In[ ]:


data_agg


# In[ ]:


# create the heatmap with dendrogram and cluster the rows and columns with the euclidean distance metric
# order the rows and columns by the linkage matrix generated by the clustering algorithm
# import linkage from scipy.cluster.hierarchy to cluster the rows and columns
# define the linkage matrix
linkage_df = linkage(
    data_agg.T, metric="euclidean", method="ward", optimal_ordering=True
)
g = sns.clustermap(
    data_agg.T,
    cmap="viridis",
    metric="euclidean",
    method="ward",
    row_cluster=True,
    col_cluster=False,
    row_linkage=linkage_df,
    col_linkage=linkage_df,
    xticklabels=True,
    yticklabels=True,
    vmin=0,
    vmax=1,
)
# save the heatmap
plt.savefig("./figures/heatmap_PBMC.png", bbox_inches="tight")
# show the heatmap
plt.show()
