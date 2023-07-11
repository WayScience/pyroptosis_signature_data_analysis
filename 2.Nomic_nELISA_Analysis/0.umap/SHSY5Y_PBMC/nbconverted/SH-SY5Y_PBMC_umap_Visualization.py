#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly_express as px
import seaborn as sns

# In[2]:


def umap_graph(df: pd.DataFrame, color: str, save_name: pathlib.Path, title: str):
    """plot umap graph

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe with umap coordinates
    color : str
        column to color points by
    save_name : pathlib.Path
        name of file to be saved
    title : str
        title of plot
    """

    # plot UMAP
    fig_2d = px.scatter(
        df,
        x="umap_1",
        y="umap_2",
        color=color,
        labels={"color": "Cell Type"},
        title=f"{title}",
    ).update_layout(xaxis_title="UMAP_1", yaxis_title="UMAP_2")

    fig_2d.update_traces(marker={"size": 12})
    fig_2d.show()

    # save UMAP
    # fig_2d.write_html(f"{save_name}.html")
    fig_2d.write_image(f"{save_name}.png")


# In[3]:


# import nELISA umap data
nELISA_plate_430420_SH_SY5Y_path = pathlib.Path(
    "./results/nELISA_plate_430420_umap.csv"
)

nELISA_plate_430420 = pd.read_csv(nELISA_plate_430420_SH_SY5Y_path)


# In[4]:


nELISA_plate_430420 = nELISA_plate_430420[
    nELISA_plate_430420.columns.drop(list(nELISA_plate_430420.filter(regex="pgML")))
]

col_list = [
    "cell_type",
    "inducer1",
    "inhibitor",
    "Function",
    "inducer1_plus_concentration",
    "inducer1_plus_inhibitor",
    "inducer1_plus_concentration_plus_inhibitor",
]

for col in col_list:
    umap_graph(
        df=nELISA_plate_430420,
        color=f"{col}",
        save_name=f"./figures/nELISA_Normalized_UMAP_{col}",
        title=f"Umap projection of nELISA data colored by {col.replace('_', ' ')}",
    )
