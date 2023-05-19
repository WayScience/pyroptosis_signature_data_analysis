#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn import preprocessing

# In[2]:


# sourcery skip: my-custom-rule
SH_SY5Y_PATH = pathlib.Path("../../Data/clean/Plate2/nELISA_plate_430420_SH_SY5Y.csv")

manual_cluster_1_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_1.csv"
)

manual_cluster_2_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_2.csv"
)

SHSY5Y_df = pd.read_csv(SH_SY5Y_PATH)

manual_clusters_1 = pd.read_csv(manual_cluster_1_path)
manual_clusters_2 = pd.read_csv(manual_cluster_2_path)


# In[3]:


# select data only columns and make floats
nELISA_data_values = SHSY5Y_df.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")
nELISA_data_values.head()


# In[4]:


# normalize data via max value in each column
max_values = nELISA_data_values.max()  # find max value in each column
nELISA_data_values_sensor_max_norm = nELISA_data_values.div(
    max_values
)  # divide each value in each column by max value in that column
nELISA_data_values_sensor_max_norm.head()
# min max normalization via sklearn


# normalize data via min max normalization
min_max_scaler = preprocessing.MinMaxScaler()
nELISA_data_values_min_max_norm = min_max_scaler.fit_transform(nELISA_data_values)
nELISA_data_values_min_max_norm = pd.DataFrame(
    nELISA_data_values_min_max_norm, columns=nELISA_data_values.columns
)
nELISA_data_values_min_max_norm.head()


# In[5]:


# drop columns that are named with NSU
Metadata = SHSY5Y_df.drop(SHSY5Y_df.filter(like="NSU", axis=1), axis=1)
Metadata = Metadata.drop(SHSY5Y_df.filter(like="pgML", axis=1), axis=1)


# In[6]:


analysis_df = pd.concat([Metadata, nELISA_data_values_min_max_norm], axis=1)


# In[7]:


# add manual clusters columns to dataframe
nELISA_plate_430420 = pd.merge(
    analysis_df,
    manual_clusters_2,
    on=(
        "inducer1",
        "inducer1_concentration_value",
        "inhibitor",
        "inhibitor_concentration_value",
        "inducer2",
        "inducer2_concentration_value",
    ),
    how="inner",
)


# In[8]:


# dose column merge
conditions = [
    (nELISA_plate_430420["inducer2"].isnull()),
    nELISA_plate_430420["inducer2"].notnull(),
]

results = [
    (
        nELISA_plate_430420["inducer1"]
        + "_"
        + nELISA_plate_430420["inducer1_concentration"].astype(str)
        + "_"
        + nELISA_plate_430420["inhibitor"].astype(str)
        + "_"
        + nELISA_plate_430420["inhibitor_concentration"].astype(str)
    ).astype(str),
    (
        nELISA_plate_430420["inducer1"]
        + "_"
        + nELISA_plate_430420["inducer1_concentration"].astype(str)
        + "_"
        + nELISA_plate_430420["inducer2"]
        + "_"
        + nELISA_plate_430420["inducer2_concentration"].astype(str)
        + "_"
        + nELISA_plate_430420["inhibitor"].astype(str)
        + "_"
        + nELISA_plate_430420["inhibitor_concentration"].astype(str)
    ).astype(str),
]
nELISA_plate_430420["Treatment_and_Dose"] = np.select(conditions, results)


results = [
    (nELISA_plate_430420["inducer1"]).astype(str),
    (nELISA_plate_430420["inducer1"] + "_" + nELISA_plate_430420["inducer2"]).astype(
        str
    ),
]
nELISA_plate_430420["Treatments"] = np.select(conditions, results)


# In[9]:


# select rows that contain 'Thapsigargin_10 µM_DMSO_0.03%' from Treatment_and_Dose column
nELISA_plate_430420 = nELISA_plate_430420[
    nELISA_plate_430420["Treatments"].isin(
        [
            "Thapsigargin",
            "LPS",
            "DMSO",
        ]
    )
]

# select rows that contain 'Thapsigargin_10 µM_DMSO_0.03%' from Treatment_and_Dose column
nELISA_plate_430420 = nELISA_plate_430420[
    nELISA_plate_430420["inhibitor"].isin(["DMSO"])
]


# In[10]:


def plot_cytokines(df, cytokine, cell_type):
    plt.figure()
    plt.tight_layout()
    sns.set(rc={"figure.figsize": (8, 5)})
    # plot a bar chart
    sns.barplot(
        y=cytokine,
        x="Treatments",
        hue="cell_type",
        data=df,
        estimator=np.mean,
        errorbar=("sd"),
    )
    plt.xticks(rotation=90)
    plt.ylim(0, 1)
    plt.ylabel(cytokine)
    plt.savefig(f"figures/{cytokine}_{cell_type}.png", bbox_inches="tight")
    plt.show()
    plt.close()


cell_type = nELISA_plate_430420["cell_type"].unique()[0]
# plot all cytokines
plot_cytokines(nELISA_plate_430420, "IL-1 beta [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "IL-18 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "IL-6 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "IL-8 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "IL-10 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "CCL1 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "CCL2 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "CCL5 [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "TNF alpha [NSU]", cell_type)
plot_cytokines(nELISA_plate_430420, "TRAIL [NSU]", cell_type)


# In[11]:


# open pdf file
with PdfPages(f"figures/inducers_{cell_type}.pdf") as pdf:
    # plot all cytokines and selected inducer and plot them in a pdf
    for i in nELISA_plate_430420.filter(like="NSU", axis=1).columns:
        plt.figure()
        plt.tight_layout()
        sns.set(rc={"figure.figsize": (8, 5)})
        # plot a bar chart
        sns.barplot(
            y=nELISA_plate_430420[i],
            x="inducer1",
            data=nELISA_plate_430420,
            hue="cell_type",
            estimator=np.mean,
            # standard deviation errorbars
            errorbar=("sd"),
        )
        plt.xticks(rotation=45)
        pdf.savefig(bbox_inches="tight")
        plt.close()


# In[ ]:


# open pdf file
with PdfPages(f"figures/death_type_{cell_type}.pdf") as pdf:
    # plot all cytokines and selected inducer and plot them in a pdf
    for i in nELISA_plate_430420.filter(like="NSU", axis=1).columns:
        plt.figure()
        plt.tight_layout()
        sns.set(rc={"figure.figsize": (8, 5)})
        # plot a bar chart
        sns.barplot(
            y=nELISA_plate_430420[i],
            x="Function",
            data=nELISA_plate_430420,
            hue="cell_type",
            estimator=np.mean,
            errorbar=("sd"),
        )
        pdf.savefig(bbox_inches="tight")
        plt.close()
