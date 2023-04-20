#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import iplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly_express as px
import seaborn as sns
import umap
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot, plot
from sklearn import preprocessing

# In[2]:


def umap_graph(df: pd.DataFrame, color: str, save_name: str):
    """plot umap graph

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe with umap coordinates
    color : str
        column to color points by
    save_name : str
        name of file to be saved
    """

    # plot UMAP
    fig_2d = px.scatter(
        df,
        x="umap_1",
        y="umap_2",
        color=color,
        labels={"color": "Cell Type"},
        title="UMAP projection of the nELISA data",
    ).update_layout(xaxis_title="UMAP_1", yaxis_title="UMAP_2")

    fig_2d.update_traces(marker={"size": 12})
    fig_2d.show()

    # save UMAP
    # fig_2d.write_html(f"{save_name}.html")
    fig_2d.write_image(f"{save_name}.png")


# In[3]:


nELISA_plate_430418_430419_path = pathlib.Path(
    "../Data/clean/nELISA_plate_430418_430419.csv"
)
nELISA_plate_430420_path = pathlib.Path("../Data/clean/nELISA_plate_430420.csv")
manual_clusters_path = pathlib.Path("../Data/Manual_Treatment_Clusters.csv")

nELISA_plate_430418_430419 = pd.read_csv(nELISA_plate_430418_430419_path)
nELISA_plate_430420 = pd.read_csv(nELISA_plate_430420_path)
manual_clusters = pd.read_csv(manual_clusters_path)


# In[4]:


# select data only columns and make floats
nELISA_data_values = nELISA_plate_430420.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")
nELISA_data_values.head()


# In[5]:


# normalize data via max value in each column
max_values = nELISA_data_values.max()  # find max value in each column
nELISA_data_values_sensor_max_norm = nELISA_data_values.div(
    max_values
)  # divide each value in each column by max value in that column
nELISA_data_values_sensor_max_norm.head()


# In[6]:


# print mean and stdev of first data column before and after normalization to check normalization
print(f"NSU nELISA mean of Activin A: {nELISA_data_values['Activin A [NSU]'].mean()}")
print(f"NSU nELISA STDEV of Activin A: {nELISA_data_values['Activin A [NSU]'].std()}")

print(
    f"NSU sensor normalized nELISA mean of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].mean()}"
)
print(
    f"NSU sensor normalized nELISA STDEV of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].std()}"
)


# In[7]:


# rename columns to remove special character "/"
nELISA_plate_430420.columns = nELISA_plate_430420.columns.str.replace("/", "_")

# set umap parameters
umap_params = umap.UMAP(
    n_neighbors=6,
    min_dist=0.8,
    n_components=2,
    metric="cosine",
    spread=1.1,
    init="random",
    random_state=0,
)

# fit and transform data for umap
proj_2d = umap_params.fit_transform(nELISA_data_values_sensor_max_norm)

# add umap coordinates to dataframe of metadata and raw data
nELISA_plate_430420["umap_1"] = proj_2d[:, 0]
nELISA_plate_430420["umap_2"] = proj_2d[:, 1]

# add manual clusters columns to dataframe
nELISA_plate_430420 = pd.merge(
    nELISA_plate_430420, manual_clusters, on=("inducer1", "inhibitor"), how="inner"
)

# define output paths
nELISA_plate_430420_out_path = pathlib.Path("./results/nELISA_plate_430420_umap.csv")
# write to csv
nELISA_plate_430420.to_csv(nELISA_plate_430420_out_path, index=False)


# In[8]:


nELISA_plate_430420 = nELISA_plate_430420[
    nELISA_plate_430420.columns.drop(list(nELISA_plate_430420.filter(regex="pgML")))
]
for col in nELISA_plate_430420.columns:
    umap_graph(
        nELISA_plate_430420,
        f"{col}",
        f"./figures/nELISA_Normalized_UMAP_{col}",
    )
