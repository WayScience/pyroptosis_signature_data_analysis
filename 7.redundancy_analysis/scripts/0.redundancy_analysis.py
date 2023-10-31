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
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# In[2]:


# Parameters
cell_type = "PBMC"
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


# shuffle the data both rows and columns
if Shuffle:
    morphology_data = morphology_data.sample(frac=1).reset_index(drop=True)
    morphology_data = morphology_data.sample(frac=1, axis=1).reset_index(drop=True)
    nomic_data_values = nomic_data_values.sample(frac=1).reset_index(drop=True)
    nomic_data_values = nomic_data_values.sample(frac=1, axis=1).reset_index(drop=True)


# ### Variables
# $Y_{M \times P} = MorphologyData$
# $X_{N \times Q} = NomicData$
# Where
# $M = Rows of MorphologyData$
# $P = Columns of MorphologyData$
# $N = Rows of NomicData$
# $Q = Columns of NomicData$

# In[7]:


# define the variables
N = morphology_data.shape[0]
P = morphology_data.shape[1]

N = nomic_data_values.shape[0]
Q = nomic_data_values.shape[1]
print("N:", N, "P:", P, "Q:", Q)
K = min(N, P, Q)
print("K:", K)


# In[8]:


cca = CCA(n_components=K)
cca.fit(morphology_data, nomic_data_values)
X_c, Y_c = cca.transform(morphology_data, nomic_data_values)
ccascore = [cca.score(morphology_data, nomic_data_values), X_c, Y_c]


# In[9]:


# make a dataframe of the coefficients
coef_df = pd.DataFrame(
    cca.coef_, columns=morphology_data.columns, index=nomic_data_values.columns
)
# get the X and Y coefficients as np arrays
X_coef = cca.x_weights_
Y_coef = cca.y_weights_


# From the canonical coefficients we can calculate the variance extracted by each canonical variable.
# $u_k = \frac{1}{P} \sum^P_{p=1} \tilde a^2_{pk}$
#
# $v_k = \frac{1}{Q} \sum^Q_{q=1} \tilde a^2_{qk}$

# In[10]:


# get the variance explained by each canonical variate
u_k = sum(cca.x_loadings_**2) / P
v_k = sum(cca.y_loadings_**2) / Q


# In[11]:


# scatter plot of the canonical variates
plt.scatter(u_k, v_k)
# fit a line to the scatter plot
m, b = np.polyfit(u_k, v_k, 1)
plt.plot(u_k, m * u_k + b)
plt.xlabel("u_k")
plt.ylabel("v_k")
plt.title("Canonical Variates")
# add r2 of regression line to plot
plt.show()

# calculate r2 of k
from sklearn.metrics import r2_score

k_r2 = r2_score(v_k, m * u_k + b)
print(k_r2)


# In[12]:


# calculate the redundancy index of each variable
u_k_RI = []
v_k_RI = []
for i in enumerate(u_k):
    # add to list
    u_k_RI.append(i[1] / k_r2)
for i in enumerate(v_k):
    # add to list
    v_k_RI.append(i[1] / k_r2)

sum_u_k_RI = sum(u_k_RI)
sum_v_k_RI = sum(v_k_RI)
print(
    f"The variance of the canonical variates explained by the variables is {sum_u_k_RI} for X and {sum_v_k_RI} for Y"
)

# plot the redundancy index
plt.scatter(u_k_RI, v_k_RI)
plt.xlabel("u_k_RI")
plt.ylabel("v_k_RI")
plt.title("Redundancy Index")
plt.show()


# In[13]:


# df from the canonical variates
RI_df = pd.DataFrame([u_k_RI, v_k_RI], index=["u_k_RI", "v_k_RI"]).T
RI_df["Shuffle"] = Shuffle
# set output path
results_file_path_RI = pathlib.Path(f"../results/{cell_type}_redundancy_index.csv")
results_file_path_RI.parent.mkdir(parents=True, exist_ok=True)

if not results_file_path_RI.exists():
    RI_df.to_csv(results_file_path_RI, index=False)
elif results_file_path_RI.exists():
    # read in the existing file
    existing_df = pd.read_csv(results_file_path_RI)
    if len(existing_df["Shuffle"].unique()) > 1:
        # overwrite the file
        RI_df.to_csv(results_file_path_RI, index=False)
    elif existing_df["Shuffle"].unique() != Shuffle:
        # append to the file
        RI_df.to_csv(results_file_path_RI, mode="a", header=False, index=False)
    else:
        print("The file already exists and the shuffle value is the same")
        print("No write occured")
else:
    print("Something went wrong: check path for the redundancy index file")


# In[14]:


out_dict = {}
for i in tqdm(range(2, K)):
    cca = CCA(n_components=i)
    cca.fit(morphology_data, nomic_data_values)
    X_c, Y_c = cca.transform(morphology_data, nomic_data_values)
    cca.score(morphology_data, nomic_data_values), X_c, Y_c
    coef_df = pd.DataFrame(
        cca.coef_, columns=morphology_data.columns, index=nomic_data_values.columns
    )
    # get the X and Y coefficients as np arrays
    X_coef = cca.x_weights_
    Y_coef = cca.y_weights_
    # get the variance explained by each canonical variate
    u_k = sum(cca.x_loadings_**2) / P
    v_k = sum(cca.y_loadings_**2) / Q
    k_r2 = r2_score(v_k, m * u_k + b)
    # calculate the redundancy index of each variable
    u_k_RI = []
    v_k_RI = []
    for i in enumerate(u_k):
        # add to list
        u_k_RI.append(i[1] / k_r2)
    for i in enumerate(v_k):
        # add to list
        v_k_RI.append(i[1] / k_r2)

    sum_u_k_RI = sum(u_k_RI)
    sum_v_k_RI = sum(v_k_RI)
    out_dict[i[0]] = [sum_u_k_RI, sum_v_k_RI, k_r2]


# In[15]:


# dict to df
out_df = pd.DataFrame.from_dict(
    out_dict, orient="index", columns=["X_RI", "Y_RI", "r2"]
)
# reset index
out_df = out_df.reset_index()
# rename index column
out_df = out_df.rename(columns={"index": "K"})
out_df["Shuffle"] = Shuffle
# plot the redundancy index
plt.plot(out_df["K"], out_df["r2"], label="X")
plt.xlabel("K")
plt.ylabel("r2")
plt.title("Skree Plot of r2 against K")
plt.show()

# plot the redundancy index
plt.plot(out_df["K"], out_df["X_RI"], label="X")
plt.plot(out_df["K"], out_df["Y_RI"], label="Y")
plt.xlabel("K")
plt.ylabel("Redundancy Index")
plt.title("Redundancy Index against K")
plt.legend()
plt.show()


# In[16]:


# check if file exists
if not results_file_path.exists():
    # write to file
    out_df.to_csv(results_file_path, index=False)
    pass
elif results_file_path.exists():
    # read in the file
    old_df = pd.read_csv(results_file_path)
    if len(old_df["Shuffle"] > 1):
        # overwrite the file
        out_df.to_csv(results_file_path, index=False)
elif results_file_path.exists():
    # read in the file
    old_df = pd.read_csv(results_file_path)
    if old_df["Shuffle"].unique() == Shuffle:
        pass
    else:
        # concat the dfs and write to file
        out_df = pd.concat([old_df, out_df])
        out_df.to_csv(results_file_path, index=False)
