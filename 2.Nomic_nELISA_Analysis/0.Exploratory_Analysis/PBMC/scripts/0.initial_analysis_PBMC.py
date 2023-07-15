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

# In[ ]:


# set path to data
PBMC_path = pathlib.Path("../../Data/clean/Plate2/nELISA_plate_430420_PBMC.csv")


manual_cluster_1_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_1.csv"
)
manual_cluster_2_path = pathlib.Path(
    "../../Data/clean/Plate2/Manual_Treatment_Clusters_2.csv"
)

# read in data
PBMC_df = pd.read_csv(PBMC_path)

manual_clusters_1 = pd.read_csv(manual_cluster_1_path)
manual_clusters_2 = pd.read_csv(manual_cluster_2_path)


# In[ ]:


# select data only columns and make floats
nELISA_data_values = PBMC_df.filter(like="NSU", axis=1)
nELISA_data_values = nELISA_data_values.astype("float")


# In[ ]:


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


# In[ ]:


# drop columns that are named with NSU
Metadata = PBMC_df.drop(PBMC_df.filter(like="NSU", axis=1), axis=1)
Metadata = Metadata.drop(PBMC_df.filter(like="pgML", axis=1), axis=1)


# In[ ]:


# merge metadata and normalized data values
analysis_df = pd.concat([Metadata, nELISA_data_values_min_max_norm], axis=1)


# In[ ]:


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


# In[ ]:


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


# In[ ]:


# inducer and dose column merge
nELISA_plate_430420["Inducer1_and_dose"] = (nELISA_plate_430420["inducer1"]).astype(
    str
) + ("_" + nELISA_plate_430420["inducer1_concentration"].astype(str))
# select inhibitors that are 'DMSO'
nELISA_plate_430420_no_inhibitors = nELISA_plate_430420[
    nELISA_plate_430420["inhibitor"] == "DMSO"
]


# In[ ]:


# select rows that contain 'Thapsigargin_10 µM_DMSO_0.03%' from Treatment_and_Dose column
nELISA_plate_430420_Thapsi_sub = nELISA_plate_430420[
    nELISA_plate_430420["Treatments"].isin(
        [
            "Thapsigargin",
            "LPS",
            "DMSO",
        ]
    )
]

# select rows that contain 'Thapsigargin_10 µM_DMSO_0.03%' from Treatment_and_Dose column
nELISA_plate_430420_Thapsi_sub = nELISA_plate_430420_Thapsi_sub[
    nELISA_plate_430420_Thapsi_sub["inhibitor"].isin(["DMSO"])
]


# In[ ]:


# resort treatment list
treatments = [
    "DMSO_0.10%",
    "Topotecan_5 nM",
    "Topotecan_10 nM",
    "Topotecan_20 nM",
    "Disulfiram_0.1 µM",
    "Disulfiram_1 µM",
    "Disulfiram_2.5 µM",
    "Flagellin_0.1 µg/ml",
    "Flagellin_1 µg/ml",
    "LPS_0.01 µg/ml",
    "LPS_0.1 µg/ml",
    "LPS_1 µg/ml",
    "LPS_10 µg/ml",
    "LPS_100 µg/ml",
    "Thapsigargin_1 µM",
    "Thapsigargin_10 µM",
    "H2O2_100 nM",
    "H2O2_100 µM",
]


# In[ ]:


# drop rows in inducer2 that are 'nigericin'
nELISA_plate_430420_no_inhibitors = nELISA_plate_430420_no_inhibitors[
    nELISA_plate_430420_no_inhibitors["inducer2"] != "Nigericin"
]


# In[ ]:


# save nELISA_plate_430420_no_inhibitors dataframe to csv file
nELISA_plate_430420_no_inhibitors.to_csv(
    "../../Data/filtered/nELISA_plate_430420_no_inhibitors.csv", index=False
)


# ## Graphing the data

# In[ ]:


# define cell type
cell_type = nELISA_plate_430420["cell_type"].unique()[0]
# open pdf file
with PdfPages(f"figures/inducers_and_dose_{cell_type}.pdf") as pdf:
    # plot all cytokines and selected inducer and plot them in a pdf
    # plot with x axis as inducer and dose order by list treatment
    for i in nELISA_plate_430420.filter(like="NSU", axis=1).columns:
        plt.figure()
        plt.tight_layout()
        sns.set(rc={"figure.figsize": (8, 5)})
        # plot a bar chart
        # order by list treatment
        sns.barplot(
            x="Inducer1_and_dose",
            y=i,
            data=nELISA_plate_430420,
            hue="cell_type",
            estimator=np.mean,
            # standard deviation errorbars
            errorbar=("sd"),
            order=treatments,
        )
        plt.xticks(rotation=90)
        pdf.savefig(bbox_inches="tight")
        plt.close()


# In[ ]:


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
