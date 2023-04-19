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


nELISA_data_input_path = pathlib.Path("../Data/clean/nELISA_Data_UCA1_2023.04.11.csv")
MetaData_input_path = pathlib.Path("../Data/clean/Metadata_UCA1_2023.04.11.csv")

nELISA_data_all = pd.read_csv(nELISA_data_input_path)
MetaData = pd.read_csv(MetaData_input_path)


# In[3]:


# nELISA_data = pd.merge(MetaData, nELISA_data_all, left_on='plate_position',right_on='plate_position', how='inner')
nELISA_data = pd.merge(MetaData, nELISA_data_all, on="plate_position", how="right")
nELISA_data_values = nELISA_data.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")
print(nELISA_data_values.shape)
# pd.merge(MetaData, nELISA_data_all, on='plate_position',how='left').shape
nELISA_data_values.head()


# In[4]:


max_values = nELISA_data_values.max()
nELISA_data_values_sensor_max_norm = nELISA_data_values.div(max_values)
nELISA_data_values_sensor_max_norm.head()


# In[5]:


print(f"NSU nELISA mean of Activin A: {nELISA_data_values['Activin A [NSU]'].mean()}")
print(f"NSU nELISA STDEV of Activin A: {nELISA_data_values['Activin A [NSU]'].std()}")

print(
    f"NSU sensor normalized nELISA mean of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].mean()}"
)
print(
    f"NSU sensor normalized nELISA STDEV of Activin A: {nELISA_data_values_sensor_max_norm['Activin A [NSU]'].std()}"
)

# normalize the data min max
scaler = preprocessing.MinMaxScaler().fit(nELISA_data_values)
nELISA_data_min_max = scaler.transform(nELISA_data_values)

nELISA_data_min_max = pd.DataFrame(
    nELISA_data_min_max, columns=nELISA_data_values.columns
)
print(
    f"NSU Min Max normalized nELISA mean of Activin A: {nELISA_data_min_max['Activin A [NSU]'].mean()}"
)
print(
    f"NSU Min Max normalized nELISA STDEV of Activin A: {nELISA_data_min_max['Activin A [NSU]'].std()}"
)

# normalize the data standard scaler
scaler = preprocessing.StandardScaler().fit(nELISA_data_values)
nELISA_data_standard = scaler.transform(nELISA_data_values)
nELISA_data_standard = pd.DataFrame(
    nELISA_data_standard, columns=nELISA_data_values.columns
)

print(
    f"NSU Standard Scalar normalized nELISA mean of Activin A: {nELISA_data_standard['Activin A [NSU]'].mean()}"
)
print(
    f"NSU Standard Scalar normalized nELISA STDEV of Activin A: {nELISA_data_standard['Activin A [NSU]'].std()}"
)


# In[6]:


nELISA_data_values.head()


# In[7]:


nELISA_data_values_sensor_max_norm.head()


# In[8]:


nELISA_data_min_max.head()


# In[9]:


nELISA_data_standard.head()


# In[10]:


def umap_graph(df, df1, save_name):
    # UMAP parameters
    umap_df = umap.UMAP(
        n_neighbors=6,
        min_dist=0.8,
        n_components=2,
        metric="cosine",
        spread=1.1,
        init="random",
        random_state=0,
    )

    # fit and transform data
    proj_2d = umap_df.fit_transform(df)

    # plot UMAP
    fig_2d = px.scatter(
        proj_2d,
        x=0,
        y=1,
        color=df1.cell_type,
        labels={"color": "Cell Type"},
        title="UMAP projection of the nELISA data",
    ).update_layout(xaxis_title="UMAP_1", yaxis_title="UMAP_2")

    fig_2d.update_traces(marker={"size": 15})
    fig_2d.show()

    # save UMAP
    fig_2d.write_html(f"{save_name}.html")
    fig_2d.write_image(f"{save_name}.png")


umap_graph(nELISA_data_values, nELISA_data, "./figures/nELISA_data_values_raw")
umap_graph(
    nELISA_data_values_sensor_max_norm,
    nELISA_data,
    "./figures/nELISA_data_values_max_sensor_noramalized",
)
umap_graph(nELISA_data_min_max, nELISA_data, "./figures/nELISA_data_min_max_normalized")
umap_graph(
    nELISA_data_standard, nELISA_data, "./figures/nELISA_data_standard_normalized"
)


# In[ ]:
