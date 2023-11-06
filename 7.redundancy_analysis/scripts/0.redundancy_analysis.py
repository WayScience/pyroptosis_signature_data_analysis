#!/usr/bin/env python
# coding: utf-8

# # Canonical Correlation Analysis (CCA)
# Here I calculate the Canonical Correlation Coefficients and the canonical variables for the two datasets.
# I also plot the correlation coefficients and the canonical variables.

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cross_decomposition import CCA
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# In[2]:


# Parameters
cell_type = "SHSY5Y"
Shuffle = False


# In[3]:


# set paths to data
morphology_data_path = pathlib.Path(
    f"../../data/{cell_type}_preprocessed_sc_norm_aggregated.parquet"
).resolve(strict=True)
nomic_data_path = pathlib.Path(
    f"../../2.Nomic_nELISA_Analysis/Data/clean/Plate2/nELISA_plate_430420_{cell_type}_clean.parquet"
).resolve(strict=True)

# output path
results_file_path = pathlib.Path(f"../results/{cell_type}_redundancy_analysis.csv")
results_file_path.parent.mkdir(parents=True, exist_ok=True)

# read data
morphology_data = pd.read_parquet(morphology_data_path)
nomic_data = pd.read_parquet(nomic_data_path)


# In[4]:


# get the columns that contain metadata
morphology_metadata = morphology_data[
    morphology_data.columns[morphology_data.columns.str.contains("Metadata")]
]
morphology_data = morphology_data.drop(morphology_metadata.columns, axis=1)

nomic_data_values = nomic_data[
    nomic_data.columns[nomic_data.columns.str.contains("[NSU]", regex=True)]
]
nomic_metadata = nomic_data.drop(nomic_data_values.columns, axis=1)


# In[5]:


# standardize the data for nomic standard scalar
scaler = StandardScaler()
nomic_data_values = scaler.fit_transform(nomic_data_values)
nomic_data_values = pd.DataFrame(
    nomic_data_values,
    columns=nomic_data.columns[nomic_data.columns.str.contains("[NSU]", regex=True)],
)


# In[6]:


# check the scale of the data
nomic_data_values.describe()


# In[7]:


# shuffle the data both rows and columns
if Shuffle:
    for column in nomic_data_values:
        np.random.shuffle(nomic_data_values[column].values)
    for column in morphology_data:
        np.random.shuffle(morphology_data[column].values)


# ### Variables
# $Y_{M \times P} = MorphologyData$
# $X_{N \times Q} = NomicData$
# Where
# $M = Rows of MorphologyData$
# $P = Columns of MorphologyData$
# $N = Rows of NomicData$
# $Q = Columns of NomicData$

# In[8]:


# define the variables
N = morphology_data.shape[0]
P = morphology_data.shape[1]

N = nomic_data_values.shape[0]
Q = nomic_data_values.shape[1]
print("N:", N, "P:", P, "Q:", Q)
K = min(N, P, Q)
print("K:", K)


# In[9]:


cca = CCA(n_components=K)
cca.fit(morphology_data, nomic_data_values)
X_c, Y_c = cca.transform(morphology_data, nomic_data_values)
r2 = [cca.score(morphology_data, nomic_data_values), X_c, Y_c][0]
print("The R2 score for the Canonical Correlation is:", r2)


# In[10]:


A_tilde = cca.x_loadings_.T
B_tilde = cca.y_loadings_.T


# From the canonical coefficients we can calculate the variance extracted by each canonical variable.
# $u_k = \frac{1}{P} \sum^P_{p=1} \tilde a^2_{pk}$
#
# $v_k = \frac{1}{Q} \sum^Q_{q=1} \tilde a^2_{qk}$

# In[11]:


u_k = []
v_k = []
for i in A_tilde:
    u_k.append(np.mean(i**2))
for i in B_tilde:
    v_k.append(np.mean(i**2))


# In[12]:


# coefficients of determination for each canonical variable
r2 = r2_score(u_k, v_k)


# In[13]:


# calculate the redundancy index for each canonical variable
RI_u = []
RI_v = []

for i in u_k:
    RI_u.append(i * r2)
for i in v_k:
    RI_v.append(i * r2)


# In[14]:


RI_u_min = np.min(RI_u)
RI_v_min = np.min(RI_v)
RI_u_max = np.max(RI_u)
RI_v_max = np.max(RI_v)
global_min = np.min([RI_u_min, RI_v_min])
global_max = np.max([RI_u_max, RI_v_max])

# Calulate the global redundancy index
global_RI_u_v = np.sum(RI_u) + np.sum(RI_v)
global_RI_u = np.sum(RI_u) / global_RI_u_v * 100
global_RI_v = np.sum(RI_v) / global_RI_u_v * 100

# plot RI_u and RI_v
sns.set_theme(style="whitegrid")
plt.figure(figsize=(10, 10))
plt.scatter(RI_u, RI_v)
plt.xlabel(f"RI of u ( {round(global_RI_u,2)}%)")
plt.ylabel(f"RI of v ( {round(global_RI_v,2)}%)")
plt.xlim(global_min, global_max)
plt.ylim(global_min, global_max)
# add a line of y=x to the plot
plt.plot(
    [global_min, global_max],
    [global_min, global_max],
    color="red",
    linestyle="-",
    linewidth=2,
)


# In[15]:


# make a dataframe of the results
results_df = pd.DataFrame(
    {
        "RI_u": RI_u,
        "RI_v": RI_v,
        "u_k": u_k,
        "v_k": v_k,
        "r2": r2,
        "Shuffle": Shuffle,
    }
)

results_df.head(5)


# In[16]:


# check for file existence
if results_file_path.is_file():
    print("The results file exists.")
    #  read the results file
    existing_file_df = pd.read_csv(results_file_path)
    # check for if it is full for shuffle type
    if len(existing_file_df["Shuffle"].unique()) > 1:
        # delete the existing file
        results_file_path.unlink()
    elif not existing_file_df["Shuffle"].unique() == Shuffle:
        pd.concat([existing_file_df, results_df]).to_csv(results_file_path, index=False)
else:
    results_df.to_csv(results_file_path, index=False)
    print("The results file is created.")
