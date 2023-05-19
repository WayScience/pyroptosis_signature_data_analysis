#!/usr/bin/env python
# coding: utf-8

# # Exploratory Data Analysis on Wave Data

# ### Wave1 dilated 50 data
# Here we understand more about the data and how it is distrubuted.
# Wave1 data are extracted features from raw images.
# These images were processed via Cellprofiler pipelines
#
# Specifically wave1 is looking at Gasdermin-D and Nuclei Staining from a cell painting experiment.
#
# Further, nuclei were dilated using multiple values of pixel dilation. Here we use data for the 50 pixel dialation

#
# ### Wave3
# Here we look into how the data are distributed in wave3, the same plate data as wave1 but with full cell painting dataset including the Gasdermin D channel for a total of 6 channels to analyze.
# These data sets contain extracted features from raw image dataset. See README for information as to how raw image data was processed.

# In[1]:


import itertools
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

get_ipython().run_line_magic("matplotlib", "inline")

from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

sys.path.append("..")
from ..utils.utils import *

# In[2]:


# function for EDA of each input feature extraction file:
def EDA_run(df, wave_label):
    """Function runs EDA on imported data frame

    Parameters:
    df: pandas dataframe
    wave_label: string of label for outputed graphs

    Return:
    graphs, and df stats

    """
    # Call function to display df shape and # of replicates present
    df_stats(df)
    # Drop na and reindex accordingly
    df = df.dropna()
    df.reindex()
    # Check for Nans again
    df_stats(df)
    # Understand categorical data such as treatment and dosing
    df[["Metadata_treatment", "Metadata_dose"]].drop_duplicates()
    # create a list with only columns from the data frame that start with "Metadata"

    treatments = df["Metadata_treatment"].unique()
    print(treatments)

    n_clusters = len(treatments)
    # print(n_clusters)
    samples = str(treatments).strip("[]")
    # print(samples)
    new_wave_label = f"{wave_label}_treatments:{samples}"

    subset_number = 1500
    df_subset = df.sample(n=subset_number)
    # Code snipptet for metadata extraction by Jenna Tomkinson
    df_metadata = list(df_subset.columns[df_subset.columns.str.startswith("Metadata")])

    # define which columns are data and which are descriptive
    df_descriptive = df_subset[df_metadata]
    df_values = df_subset.drop(columns=df_metadata)

    treatment_ids = df_descriptive["Metadata_treatment"]
    # Cluster data
    # clustering code adapted from https://www.kaggle.com/code/aussie84/clustering-with-kmeans-pca-tsne

    kmeans = KMeans(n_clusters=n_clusters)
    clustering_ori = kmeans.fit_predict(df_values)

    X = df_values
    # n_components is the number of dimensions to reduce to
    Xtsne = TSNE(n_components=2).fit_transform(X)
    dftsneFull = pd.DataFrame(Xtsne)

    dftsneFull["cluster"] = clustering_ori
    dftsneFull.columns = ["x1", "x2", "cluster"]
    dftsneFull["Treatment"] = (
        df_descriptive["Metadata_treatment"].reset_index().drop("index", axis=1)
    )
    # Figure Showing tSNE of Clusters vs Treatment
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    plot = sns.scatterplot(
        data=dftsneFull,
        x="x1",
        y="x2",
        hue="cluster",
        legend="full",
        alpha=0.7,
        ax=ax[0],
    )
    ax[0].set_title("Visualized on TSNE")
    plot = sns.scatterplot(
        data=dftsneFull,
        x="x1",
        y="x2",
        hue="Treatment",
        legend="full",
        alpha=0.7,
        ax=ax[1],
    )
    ax[1].set_title("Visualized on TSNE")
    fig.suptitle(
        f"Comparing Clusters vs Treatment tSNE for {subset_number} subset samples for {samples} treatments"
    )
    plt.savefig(f"Figures/tSNE/{new_wave_label}_tSNE_cluster.png")
    df_values["cluster"] = clustering_ori
    # callable function for graphing features that contribute most to each cluster's grouping
    plot_features_all_cluster(
        df=df_values,
        label_col="cluster",
        n_clusters=n_clusters,
        sensitivity=0.2,
        file_label=new_wave_label,
    )


# ## Wave1 Dialte50 Analysis

# In[3]:


# Import data with low memory arg as the data are large
df = pd.read_csv(
    "../../Extracted Features (CSV files)/interstellar_wave1_dilate50_sc.csv.gz",
    low_memory=False,
)
EDA_run(df, "wave1_dialte50")


# Above tSNE shows that based on dimensionality reduction, there is no observable difference in treated cells. More sensitive methods such as machine learning models will need to be employed to achieve such.

# Each Cluster has a similar distrubution in amount of features affecting its grouping

# ## Wave3 Single Cell Feature Extraction

# In[4]:


# Import data with low memory arg as the data are large
df = pd.read_csv(
    "../../Extracted Features (CSV files)/interstellar_wave3_sc.csv.gz",
    low_memory=False,
)
EDA_run(df, "Wave3_sc")


# ## Wave3 Single Cell Normalized Features

# In[5]:


# Import data with low memory arg as the data are large
df = pd.read_csv(
    "../../Extracted Features (CSV files)/interstellar_wave3_sc_norm_cellprofiler.csv.gz",
    low_memory=False,
)
EDA_run(df, "Wave3_sc_normalized")


# ## Wave3 Single Cell Normalized Selected Features
# #### Feature Selection was performed in the data processeing repo. Check README for more information

# In[6]:


# Import data with low memory arg as the data are large
df = pd.read_csv(
    "../../Extracted Features (CSV files)/interstellar_wave3_sc_norm_fs_cellprofiler.csv.gz",
    low_memory=False,
)
EDA_run(df, "Wave3_sc_normalized_selected_features")
