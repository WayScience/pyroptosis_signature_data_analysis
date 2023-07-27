#!/usr/bin/env python
# coding: utf-8

# # Correlation of single cell morhologies and nELISA Cytokine/Chemokine Panel
# Each well is median aggregated and normalized.
# The correlation of the sc morphology features and nELISA features is calculated per:
# * well
# * per treatment
# * per selected treatment

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import seaborn as sns

# In[2]:


cell_type = "SHSY5Y"
data_type = "norm"


# In[3]:


# set path for sc data
sc_data_path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_{data_type}.parquet")

# set path for nomic data
nomic_data_path = pathlib.Path(
    f"../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}.csv"
)


# In[4]:


# read in sc data
sc_data = pd.read_parquet(sc_data_path)

# read in nomic data
nomic_data = pd.read_csv(nomic_data_path)


# In[5]:


# remove columns in nomic data that contain 'pgml'
nomic_data = nomic_data.loc[:, ~nomic_data.columns.str.contains("pgML")]

# drop the first 25 columns
nomic_data = nomic_data.drop(nomic_data.columns[3:25], axis=1)

# drop the first 3 columns
nomic_data = nomic_data.drop(nomic_data.columns[[0, 1]], axis=1)


# In[6]:


# subset sc data to only have columns that contain 'Metadata'
sc_data_metadata = sc_data.loc[:, sc_data.columns.str.contains("Metadata")]
sc_data_metadata

# get the non metadata columns
sc_data_features = sc_data.loc[:, ~sc_data.columns.str.contains("Metadata")]
sc_data_features

sc_data_features = pd.concat(
    [sc_data_features, sc_data_metadata["Metadata_Well"]], axis=1
)
sc_data_features

# get aggregate mean of each well in sc_data_features
sc_data_features_agg = sc_data_features.groupby("Metadata_Well").mean()
# reset index
sc_data_features_agg = sc_data_features_agg.reset_index()

# merge the two dataframes on Metadata_Well
sc_data_features_agg = pd.merge(
    sc_data_features_agg,
    sc_data_metadata[
        ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose", "Metadata_Well"]
    ].drop_duplicates(),
    on="Metadata_Well",
)
sc_data_features_agg


# In[ ]:


# In[ ]:


# In[7]:


# merge sc_data_features_agg and nomics_data on Metadata_Well and position_x
merged_data = pd.merge(
    sc_data_features_agg, nomic_data, left_on="Metadata_Well", right_on="position_x"
)
merged_data


# In[8]:


# drop position_x column and metadata_well
merged_data = merged_data.drop(["position_x", "Metadata_Well"], axis=1)
merged_data


# In[9]:


# aggregate the data by treatment dose
merged_data_agg = merged_data.groupby(
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
).mean()
print(merged_data_agg.shape)


# In[10]:


# reset index
merged_data_agg = merged_data_agg.reset_index()
merged_data_agg.head()


# In[11]:


# make a blank df
corr_matrix_df = pd.DataFrame(
    columns=["feature", "feature2", "correlation", "treatment", "well"]
)


# In[12]:


corr_matrix_df


# In[55]:


treatment = "DMSO_0.100_DMSO_0.025"
# treatment = 'DMSO_0.100_Z-VAD-FMK_100.0'
print(treatment)
# get the row
tmpdf = merged_data_agg.loc[
    merged_data_agg["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatment
]
tmpdf
# get the correlation between the sc data and nomic data for well B13
# tmpdf = merged_data_agg.drop(['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'], axis=1)
# tmpdf = tmpdf.corr()
# # drop sc data columns correlation to sc data columns
# tmpdf = tmpdf.iloc[0:sc_data_features_agg.shape[1]-2, sc_data_features_agg.shape[1]-2:]
# tmpdf = tmpdf.reset_index()
# tmpdf = tmpdf.rename(columns={'index': 'feature'})
# # save the correlation matrix to a csv

# # melt the correlation matrix
# corr_matrix_melted = tmpdf.reset_index().melt(id_vars='feature', value_vars=tmpdf.columns[1:])
# corr_matrix_melted = corr_matrix_melted.rename(columns={'index': 'feature', 'variable': 'feature2', 'value': 'correlation'})


# # find the treatment dose for well B13 in sc_data_metadata

# corr_matrix_melted['treatment'] = treatment
# corr_matrix_melted.to_csv(f'{cell_type}_{treatment}_correlation_matrix.csv')
# # change corralted values that are abs less than 0.5 to 0
# corr_matrix_melted.loc[abs(corr_matrix_melted['correlation']) < 0.5, 'correlation'] = 0
# # remove rows that contain all 0s
# corr_matrix_melted = corr_matrix_melted.loc[~(corr_matrix_melted==0).all(axis=1)]
# # remove all columns that contain all 0s
# corr_matrix_melted = corr_matrix_melted.loc[:, (corr_matrix_melted!=0).any(axis=0)]
# # set pivot table for corr_matrix_melted
# corr_matrix_melted
# # sc_data_features_agg
# corr_matrix_pivot = corr_matrix_melted.pivot_table(index=['feature'], columns='feature2', values='correlation')
# # remove rows that contain all 0s
# corr_matrix_pivot = corr_matrix_pivot.loc[~(corr_matrix_pivot==0).all(axis=1)]
# # remove all columns that contain all 0s
# corr_matrix_pivot = corr_matrix_pivot.loc[:, (corr_matrix_pivot!=0).any(axis=0)]
# # plot the pivot table
# # plt.figure(figsize=(40, 40))
# # sns.heatmap(corr_matrix_pivot, cmap='coolwarm', vmin=-1, vmax=1)
# # plt.show()
# corr_matrix_pivot['Activin A [NSU]'].unique()


# In[ ]:


# In[15]:


for treatment in merged_data_agg[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].unique():
    print(treatment)
    # get the row
    tmpdf = merged_data_agg.loc[
        merged_data_agg["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] == treatment
    ]
    # get the correlation between the sc data and nomic data for well B13
    tmpdf = merged_data_agg.drop(
        ["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"], axis=1
    )
    tmpdf = tmpdf.corr()
    # drop sc data columns correlation to sc data columns
    tmpdf = tmpdf.iloc[
        0 : sc_data_features_agg.shape[1], sc_data_features_agg.shape[1] :
    ]
    tmpdf = tmpdf.reset_index()
    tmpdf = tmpdf.rename(columns={"index": "feature"})
    # save the correlation matrix to a csv

    # # melt the correlation matrix
    corr_matrix_melted = tmpdf.reset_index().melt(
        id_vars="feature", value_vars=tmpdf.columns[1:]
    )
    corr_matrix_melted = corr_matrix_melted.rename(
        columns={"index": "feature", "variable": "feature2", "value": "correlation"}
    )

    # # find the treatment dose for well B13 in sc_data_metadata

    corr_matrix_melted["treatment"] = treatment
    corr_matrix_melted.to_csv(f"{cell_type}_{treatment}_correlation_matrix.csv")
    # filter out the correlation values that are less than abs 0.5
    corr_matrix_melted = corr_matrix_melted.loc[
        abs(corr_matrix_melted["correlation"]) > 0.5
    ]
    # remove rows and columns that contain all NaNs
    corr_matrix_melted = corr_matrix_melted.dropna(axis=0, how="all")
    # set pivot table for corr_matrix_melted
    corr_matrix_pivot = corr_matrix_melted.pivot_table(
        index=["feature"], columns="feature2", values="correlation"
    )
    # plot the pivot table
    plt.figure(figsize=(40, 40))
    sns.heatmap(corr_matrix_pivot, cmap="coolwarm", vmin=-1, vmax=1)
    plt.show()


# In[ ]:


corr_matrix_df
# aggreate corr_matrix_df by treatment
corr_matrix_df_agg = corr_matrix_df.groupby("treatment").mean()


# In[17]:


corr_matrix_df


# In[20]:


# In[ ]:


# get the row that is well B13
tmpdf = merged_data.loc[merged_data["position_x"] == "B13"]
# get the correlation between the sc data and nomic data for well B13
tmpdf = merged_data.drop(["position_x"], axis=1)
tmpdf = tmpdf.corr()
# drop sc data columns correlation to sc data columns
tmpdf = tmpdf.iloc[0 : sc_data_features_agg.shape[1], sc_data_features_agg.shape[1] :]
tmpdf = tmpdf.reset_index()
tmpdf = tmpdf.rename(columns={"index": "feature"})
# save the correlation matrix to a csv
# tmpdf.to_csv(f'{cell_type}_B13_correlation_matrix.csv')


# In[ ]:


# melt the correlation matrix
corr_matrix_melted = tmpdf.reset_index().melt(id_vars="index")
corr_matrix_melted = corr_matrix_melted.rename(
    columns={"index": "feature", "variable": "feature2", "value": "correlation"}
)
corr_matrix_melted


# find the treatment dose for well B13 in sc_data_metadata
corr_matrix_melted["treatment"] = sc_data_metadata.loc[
    sc_data_metadata["Metadata_Well"] == "B13"
]["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].values[0]
corr_matrix_melted["well"] = "B13"
# add the treatment dose to merged_data
# tmpdf['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'] = sc_data_metadata['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].values
corr_matrix_melted.head()


# In[ ]:


# In[ ]:


# In[ ]:


# drop sc data columns correlation to sc data columns
corr_matrix = corr_matrix.iloc[
    0 : sc_data_features_agg.shape[1], sc_data_features_agg.shape[1] :
]


# In[ ]:


# fig save path
fig_save_path = pathlib.Path("./Figures/feature_nomic_correlation")
# create the directory if it doesn't exist
fig_save_path.mkdir(parents=True, exist_ok=True)


# In[ ]:


# make correlation values that are less than abs 0.5 to 0
corr_matrix[abs(corr_matrix) < 0.5] = 0
# drop rows that are all 0
corr_matrix = corr_matrix.loc[(corr_matrix != 0).any(axis=1)]
# drop columns that are all 0
corr_matrix = corr_matrix.loc[:, (corr_matrix != 0).any(axis=0)]
# sort rows by highest correlation to lowest
corr_matrix = corr_matrix.reindex(
    corr_matrix.abs().sort_values(by="TNF alpha [NSU]", ascending=False).index
)
corr_matrix


# In[ ]:


# heatmap of corr_matrix
plt.figure(figsize=(20, 20))
sns.heatmap(corr_matrix, annot=False, cmap="coolwarm")


# In[ ]:


# heatmap of corr_matrix with seaborn clustermap function
# hierarchical clustering is performed on rows and columns
# sorted by similarity of their profiles (rows and columns)
# set fig size 8.5 x 11
# use fastcluster to speed up clustering
# use euclidean distance metric
sns.clustermap(
    corr_matrix,
    figsize=(50, 50),
    method="average",
    metric="euclidean",
    cmap="RdYlBu_r",
    linewidths=0,
    linecolor="black",
)
plt.savefig(pathlib.Path(f"{fig_save_path}/{cell_type}_clustermap.png"))
plt.show()


# In[ ]:
