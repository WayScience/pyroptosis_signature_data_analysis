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


# import nELISA umap data
nELISA_plate_430420_SH_SY5Y_path = pathlib.Path(
    "./results/nELISA_plate_430420_umap_SH-SY5Y.csv"
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
    "Manual_Cluster",
    "Function",
]

for col in col_list:
    umap_graph(
        nELISA_plate_430420,
        f"{col}",
        f"./figures/nELISA_Normalized_UMAP_{col}",
    )
