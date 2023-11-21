#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import toml
from matplotlib import rcParams
from tqdm import tqdm

rcParams.update({"figure.autolayout": True})

# create a venn diagram of the features that are significant in all conditions

warnings.filterwarnings("ignore")
from pycytominer.cyto_utils import infer_cp_features
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd


# In[2]:


cell_type = "PBMC"


# In[3]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_norm.parquet")
df = pd.read_parquet(file_path)


# In[4]:


# toml file path
ground_truth_file = pathlib.Path(
    "../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
).resolve(strict=True)
# read toml file
ground_truth = toml.load(ground_truth_file)
apopotosis_trts = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_trts = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
healthy_trts = ground_truth["Healthy"]["healthy_groups_list"]

# make a column that has the class of each treatment


df["apoptosis"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] in apopotosis_trts,
    axis=1,
)
df["pyroptosis"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] in pyroptosis_trts,
    axis=1,
)
df["healthy"] = df.apply(
    lambda row: row["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] in healthy_trts,
    axis=1,
)

# merge apoptosis, pyroptosis, and healthy columns into one column

df["labels"] = df.apply(
    lambda row: "apoptosis"
    if row["apoptosis"]
    else "pyroptosis"
    if row["pyroptosis"]
    else "healthy",
    axis=1,
)
# drop apoptosis, pyroptosis, and healthy columns
df.drop(columns=["apoptosis", "pyroptosis", "healthy"], inplace=True)


# In[5]:


df_metadata = df.filter(regex="Metadata")
df_data = df.drop(df_metadata.columns, axis=1)
df_data["Metadata_number_of_singlecells"] = df_metadata[
    "Metadata_number_of_singlecells"
]
cp_features = infer_cp_features(df)


# In[6]:


# anova for each feature in the dataframe with posthoc tukey test to determine which groups are different from each other
lst = []
# for i in cp_features:
for i in tqdm(cp_features):
    formula = f"{i} ~ C(labels) + C(Metadata_number_of_singlecells)"
    model = ols(formula, df_data).fit()
    aov_table = sm.stats.anova_lm(model, typ=2)
    posthoc = pairwise_tukeyhsd(
        df_data[i],
        df_data["labels"],
        alpha=0.001,
    )
    # print(posthoc)
    lst.append([posthoc, i])


# In[ ]:


tukey_df = pd.DataFrame()
for i in lst:
    j = pd.DataFrame(i[0]._results_table.data[1:])
    j["features"] = np.repeat(i[1], len(j))
    tukey_df = pd.concat([tukey_df, j], axis=0)

    np.repeat(i[1], len(j))

tukey_df.columns = [
    "group1",
    "group2",
    "meandiff",
    "lower",
    "upper",
    "p-adj",
    "reject",
    "features",
]
# drop the other organelle
# make new column with the absolute value of the p-adj
tukey_df["p-adj_abs"] = abs(tukey_df["p-adj"])
# make new column that states if the relationship is positive or negative
tukey_df["pos_neg"] = np.where(tukey_df["p-adj"] > 0, "positive", "negative")
# order the features by p-adj value


# In[ ]:


tukey_df.head()


# In[ ]:


# save the dataframe as a parquet file
anova_results_path = pathlib.Path(
    f"./results/{cell_type}_anova_results_all_treatments.parquet"
)
tukey_df.to_parquet(anova_results_path)

