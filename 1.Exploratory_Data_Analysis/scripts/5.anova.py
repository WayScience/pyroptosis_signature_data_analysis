#!/usr/bin/env python
# coding: utf-8

# ## Imports

# In[1]:


import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import plotly
import pyarrow.parquet as pq
import seaborn as sns
import statsmodels.api as sm
import toml

# create a venn diagram of the features that are significant in all conditions
from matplotlib_venn import venn2, venn3

warnings.filterwarnings("ignore")
from pycytominer.cyto_utils import infer_cp_features
from scipy.stats import f_oneway
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import MultiComparison, pairwise_tukeyhsd

# In[2]:


# Parameters
cell_type = "PBMC"
treatment1 = "DMSO_0.100_DMSO_0.025"
treatment2 = "LPS_0.010_DMSO_0.025"
treatment3 = "H202_100_DMSO_0.025"


# In[3]:


# Import Data
# set data file path under pathlib path for multi-system use
file_path = pathlib.Path(f"../data/{cell_type}_preprocessed_sc_norm.parquet")
df = pq.read_table(file_path).to_pandas()
df_metadata = df.filter(regex="Metadata")
df_data = df.drop(df_metadata.columns, axis=1)
cp_features = infer_cp_features(df)


# In[4]:


trt1 = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.contains(f"{treatment1}")
]
trt2 = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.contains(f"{treatment2}")
]
trt3 = df[
    df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].str.contains(f"{treatment3}")
]


# In[5]:


# df_LPS_100 = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('LPS_100.000_DMSO_0.025')]
# df_LPS_100['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_LPS_10 = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('LPS_10.000_DMSO_0.025')]
# df_LPS_10['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_flag = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('Flagellin_1.000_DMSO_0.025')]
# df_flag['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_DMSO = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('DMSO_0.100_DMSO_0.025')]
# df_DMSO['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_thapsi10 = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('Thapsigargin_10.000_DMSO_0.025')]
# df_thapsi10['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_thapsi_DMSO = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('Thapsigargin_10.000_DMSO_0.025','DMSO_0.100_DMSO_0.025')]
# df_thapsi_DMSO['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()
# df_thapsi1 = df[df['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].str.contains('Thapsigargin_1.000_DMSO_0.025')]
# df_thapsi1['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'].unique()


# ### Set up DF

# In[6]:


combined_df = pd.concat([trt1, trt2, trt3], axis=0)
print(len(trt1), len(trt2), len(trt3), len(combined_df))
combined_df.head(3)


# ## Anova + Post Hoc testing (tukeyHSD)

# In[7]:


# anova for each feature in the dataframe with posthoc tukey test to determine which groups are different from each other
lst = []
for i in cp_features:
    formula = f"{i} ~ C(oneb_Metadata_Treatment_Dose_Inhibitor_Dose) + C(Metadata_number_of_singlecells)"
    model = ols(formula, combined_df).fit()
    aov_table = sm.stats.anova_lm(model, typ=2)
    # posthoc = MultiComparison(df_DMSO_thapsi[i], df_DMSO_thapsi['oneb_Metadata_Treatment_Dose_Inhibitor_Dose'])
    posthoc = pairwise_tukeyhsd(
        combined_df[i],
        combined_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
        alpha=0.001,
    )
    # print(posthoc)
    lst.append([posthoc, i])


# In[8]:


# add all tukey test results to a dataframe from the list of tukey test results
tukey_df = pd.DataFrame()
for i in lst:
    # create data frame from tukey test results with feature name as column name and the results as the rows
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
tukey_df.head(3)


# In[9]:


# make new column with the absolute value of the p-adj
tukey_df["p-adj_abs"] = abs(tukey_df["p-adj"])
# get all p-adj values that are less than 0.05
tukey_df_sig = tukey_df[tukey_df["p-adj_abs"] < 0.01]
# make new column that states if the relationship is positive or negative
tukey_df_sig["pos_neg"] = np.where(tukey_df_sig["p-adj"] > 0, "positive", "negative")
# order the features by p-adj value
tukey_df_sig.head(2)


# ## Venn Diagrams

# #### Venn diagram prep

# In[10]:


# get all group1 rows that are DMSO and group2 rows that are LPS treatment
tukey_df_sig_trt_1v2 = tukey_df_sig[
    (tukey_df_sig["group1"] == f"{treatment1}")
    & (tukey_df_sig["group2"] == f"{treatment2}")
]
# get all group1 rows that are DMSO and group2 rows that are Thapsigargin treatment
tukey_df_sig_trt_1v3 = tukey_df_sig[
    (tukey_df_sig["group1"] == f"{treatment1}")
    & (tukey_df_sig["group2"] == f"{treatment3}")
]
# get all group1 rows that are LPS treatment and group2 rows that are Thapsigargin treatment
tukey_df_sig_trt_2v3 = tukey_df_sig[
    (tukey_df_sig["group1"] == f"{treatment2}")
    & (tukey_df_sig["group2"] == f"{treatment3}")
]


# In[11]:


from matplotlib import rcParams

rcParams.update({"figure.autolayout": True})


# ### Venn Diagram 2 groups

# In[12]:


venn2(
    [set(tukey_df_sig_trt_1v3["features"]), set(tukey_df_sig_trt_2v3["features"])],
    set_labels=(f"{treatment1} vs {treatment2}", f"{treatment1} vs {treatment3}"),
)
plt.title("Number of Significant (p-adj < 0.01) Features \n per Organelle", size=24)
save_path = pathlib.Path(
    f"./Figures/anova_of_features/{treatment1}_vs_{treatment2}_vs_{treatment3}"
)
save_path.mkdir(parents=True, exist_ok=True)
# plt.tight_layout()
# plt.autoscale()
plt.subplots_adjust(left=0.15)
plt.savefig(
    f"{save_path}/{treatment1}_vs_{treatment2}_number_sig_overlaping_features_per_organelle.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()


# In[ ]:


# In[13]:


# create a venn diagram of the features that are significant in all conditions
from matplotlib_venn import venn2, venn3_unweighted

venn3_unweighted(
    [
        set(tukey_df_sig_trt_1v2["features"]),
        set(tukey_df_sig_trt_1v3["features"]),
        set(tukey_df_sig_trt_2v3["features"]),
    ],
    set_labels=(
        f"{treatment1} vs {treatment2}",
        f"{treatment1} vs {treatment3}",
        f"{treatment2} vs {treatment3}",
    ),
)
plt.title("Number of Significant (p-adj < 0.01) Features \n per Organelle", size=24)
save_path = pathlib.Path(
    f"./Figures/anova_of_features/{treatment1}_vs_{treatment2}_vs_{treatment3}"
)
plt.tight_layout()
save_path.mkdir(parents=True, exist_ok=True)
plt.savefig(
    f"{save_path}/{treatment1}_vs_{treatment2}_{treatment3}_number_sig_overlaping_features_per_organelle.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()


# In[14]:


# get the features that are only in the DMSO vs LPS 100 condition and not in the other two conditions
tukey_df_sig_trt_1v2_unique = tukey_df_sig_trt_1v2[
    ~tukey_df_sig_trt_1v2["features"].isin(tukey_df_sig_trt_1v3["features"])
]
tukey_df_sig_trt_1v2_unique = tukey_df_sig_trt_1v2_unique[
    ~tukey_df_sig_trt_1v2_unique["features"].isin(tukey_df_sig_trt_2v3["features"])
]
# get the features that are only in the DMSO vs Thapsigargin 10 condition and not in the other two conditions
tukey_df_sig_trt_1v3_unique = tukey_df_sig_trt_1v3[
    ~tukey_df_sig_trt_1v3["features"].isin(tukey_df_sig_trt_1v2["features"])
]
tukey_df_sig_trt_1v3_unique = tukey_df_sig_trt_1v3_unique[
    ~tukey_df_sig_trt_1v3_unique["features"].isin(tukey_df_sig_trt_2v3["features"])
]
# get the features that are only in the LPS 100 vs Thapsigargin 10 condition and not in the other two conditions
tukey_df_sig_trt_2v3_unique = tukey_df_sig_trt_2v3[
    ~tukey_df_sig_trt_2v3["features"].isin(tukey_df_sig_trt_1v2["features"])
]
tukey_df_sig_trt_2v3_unique = tukey_df_sig_trt_2v3_unique[
    ~tukey_df_sig_trt_2v3_unique["features"].isin(
        tukey_df_sig_trt_1v3_unique["features"]
    )
]
# print the number of features that are only in each condition
print(
    len(tukey_df_sig_trt_1v2_unique),
    len(tukey_df_sig_trt_1v3_unique),
    len(tukey_df_sig_trt_2v3_unique),
)


# ## Get organelle names

# In[15]:


# split each feature by "_" and get the organelle name
tukey_df_sig_trt_1v2_unique["organelle"] = tukey_df_sig_trt_1v2_unique[
    "features"
].str.split("_", expand=True)[3]
tukey_df_sig_trt_1v2_unique["organelle"].unique()


# In[16]:


tukey_df_sig_trt_1v2_unique["organelle"] = tukey_df_sig_trt_1v2_unique[
    "organelle"
].replace(
    [
        None,
        "1" "0",
        "CorrGasdermin",
        "CorrER",
        "CorrMito",
        "CorrDNA",
        "CorrPM",
        "2",
        "4",
        "9",
        "Adjacent",
        "8",
        "6",
        "5",
        "3",
        "7",
    ],
    [
        "Other",
        "Other" "Other",
        "GasderminD",
        "ER",
        "Mito",
        "DNA",
        "PM",
        "Other",
        "Other",
        "Other",
        "Other",
        "Other",
        "Other",
        "Other",
        "Other",
        "Other",
    ],
)
print(tukey_df_sig_trt_1v2_unique["organelle"].unique())


# #### Plot number of significant features by organelle

# In[17]:


# drop the other organelle
tukey_df_sig_trt_1v2_unique = tukey_df_sig_trt_1v2_unique[
    tukey_df_sig_trt_1v2_unique["organelle"] != "Other"
]
# seaborn plot of each organelle and the counts of features that are significant
# set the order of the organelles
organelle_order = ["Mito", "PM", "DNA", "ER", "GasderminD"]
plt.figure(figsize=(10, 5))
sns.countplot(x="organelle", data=tukey_df_sig_trt_1v2_unique, order=organelle_order)
plt.title(
    f"Number of Significant Features per Organelle \n LPS 100 ug/mL vs DMSO", size=24
)
# add bar labels
for p in plt.gca().patches:
    plt.gca().text(
        p.get_x() + p.get_width() / 2,
        p.get_height(),
        "%d" % int(p.get_height()),
        fontsize=14,
        color="black",
        ha="center",
        va="bottom",
    )
plt.xlabel("Organelle", size=20)
plt.ylabel("Number of Significant Features \n (p-adj. < 0.01)", size=20)
plt.ylim(0, 175)
# rotate the x axis labels
# plt.xticks(rotation=45, ha='right')
# plt.tight_layout()
save_path = pathlib.Path(
    f"./Figures/anova_of_features/{treatment1}_vs_{treatment2}_vs_{treatment3}"
)
save_path.mkdir(parents=True, exist_ok=True)
plt.savefig(
    f"{save_path}/{treatment1}_vs_{treatment2}_number_sig_features_per_organelle.png",
    dpi=300,
)
plt.show()


# In[ ]:
