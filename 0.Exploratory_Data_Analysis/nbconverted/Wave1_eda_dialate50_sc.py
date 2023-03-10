#!/usr/bin/env python
# coding: utf-8

# # Wave1 Exploratory Data Analysis on dilated 50 data

# Here we understand more about the data and how it is distrubuted.
# Wave1 data are extracted features from raw images.
# These images were processed via Cellprofiler pipelines
#
# Specifically wave1 is looking at Gasdermin-D and Nuclei Staining from a cell painting experiment.
#
# Further, nuclei were dilated using multiple values of pixel dilation. Here we use data for the 50 pixel dialation

# In[2]:


import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

get_ipython().run_line_magic("matplotlib", "inline")
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE

sys.path.append("..")
from ..utils.utils import *

# In[3]:


# Import data with low memory arg as the data are large
df = pd.read_csv(
    "../../Extracted Features (CSV files)/interstellar_wave1_dilate50_sc.csv.gz",
    low_memory=False,
)


# In[4]:


# Call function to display df shape and # of replicates present
df_stats(df)


# In[5]:


# Drop na and reindex accordingly
df = df.dropna()
df.reindex()
# Check for Nans again
df_stats(df)


# In[6]:


# Understand categorical data such as treatment and dosing
df[["Metadata_treatment", "Metadata_dose"]].drop_duplicates()


# In[7]:


# create a list with only columns from the data frame that start with "Metadata"
# Code by Jenna Tomkinson
df_subset = df.sample(n=1500)

df_metadata = list(df_subset.columns[df_subset.columns.str.startswith("Metadata")])

# define which columns are data and which are descriptive
df_descriptive = df_subset[df_metadata]
df_values = df_subset.drop(columns=df_metadata)

treatment_ids = df_descriptive["Metadata_treatment"]


# In[8]:


# Cluster data
# clustering code adapted from https://www.kaggle.com/code/aussie84/clustering-with-kmeans-pca-tsne

kmeans = KMeans(n_clusters=9)
clustering_ori = kmeans.fit_predict(df_values)

X = df_values
Xtsne = TSNE(n_components=2).fit_transform(X)
dftsneFull = pd.DataFrame(Xtsne)

dftsneFull["cluster"] = clustering_ori
dftsneFull.columns = ["x1", "x2", "cluster"]
dftsneFull["Treatment"] = (
    df_descriptive["Metadata_treatment"].reset_index().drop("index", axis=1)
)


# In[9]:


# Figure Showing tSNE of Clusters vs Treatment
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
plot = sns.scatterplot(
    data=dftsneFull, x="x1", y="x2", hue="cluster", legend="full", alpha=0.7, ax=ax[0]
)
ax[0].set_title("Visualized on TSNE")
plot = sns.scatterplot(
    data=dftsneFull, x="x1", y="x2", hue="Treatment", legend="full", alpha=0.7, ax=ax[1]
)
ax[1].set_title("Visualized on TSNE")
fig.suptitle("Comparing Clusters vs Treatment tSNE")
df_values["cluster"] = clustering_ori


# Above tSNE shows that based on dimensionality reduction, there is no observable difference in treated cells. More sensitive methods such as machine learning models will need to be employed to achieve such.

# In[10]:


# callable function for graphing features that contribute most to each cluster's grouping
plot_features_all_cluster(
    df=df_values, label_col="cluster", n_clusters=6, sensitivity=0.2
)


# Each Cluster has a similar distrubution in amount of features affecting its grouping
