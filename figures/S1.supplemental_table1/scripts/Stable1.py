#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import pandas as pd

# In[2]:


table_path = pathlib.Path(
    f"../../../2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/tukey_unfiltered_nomic_results.csv"
)
df = pd.read_csv(table_path)


# In[3]:


# sort the df by cytokine
df = df.sort_values(by=["cytokine", "group"])
df.head()


# In[4]:


# make a new column that is boolean for if the cytokine is significant or not
df["significant"] = df["p-adj"] < 0.05

# make column for if the cytokine is significant and upregulated
df["apoptotic"] = df["significant"] & (df["group"] == "apoptosis_healthy")
df["pyroptotic"] = df["significant"] & (df["group"] == "healthy_pyroptosis")
df["apoptotic_vs_pyroptotic"] = df["significant"] & (
    df["group"] == "apoptosis_pyroptosis"
)
df

# make a consensus column that says if the cytokine is all three, two, or one or none of the above
df["consensus"] = "none"
df.loc[df["apoptotic_vs_pyroptotic"], "consensus"] = "apoptotic_vs_pyroptotic"
df.loc[df["apoptotic"], "consensus"] = "apoptotic"
df.loc[df["pyroptotic"], "consensus"] = "pyroptotic"
df.loc[df["apoptotic"] & df["pyroptotic"], "consensus"] = "apoptotic_and_pyroptotic"
df.value_counts("consensus")
# output the raw data
# make a dir if it doesn't exist
pathlib.Path("../results").mkdir(exist_ok=True, parents=True)
df.to_csv("../results/consensus_cytokine_results.csv", index=False)


# In[5]:


# unmelt the df
df1 = df.pivot(index="cytokine", columns="group", values="p-adj")
df1
# if apoptosis vs healthy is significant
df1["apoptotic_vs_healthy"] = df1["apoptosis_healthy"] < 0.05
df1["apoptosic_pyroptosic"] = df1["apoptosis_pyroptosis"] < 0.05
df1["pyroptotic_vs_healthy"] = df1["healthy_pyroptosis"] < 0.05
df1["all_significant"] = (
    df1["apoptotic_vs_healthy"]
    & df1["apoptosic_pyroptosic"]
    & df1["pyroptotic_vs_healthy"]
)
df1["none_significant"] = (
    ~df1["apoptotic_vs_healthy"]
    & ~df1["apoptosic_pyroptosic"]
    & ~df1["pyroptotic_vs_healthy"]
)
df1["healthy"] = df1["apoptotic_vs_healthy"] & df1["pyroptotic_vs_healthy"]
df1["pyroptotic"] = df1["apoptosic_pyroptosic"] & df1["pyroptotic_vs_healthy"]
df1["apoptotic"] = df1["apoptotic_vs_healthy"] & df1["apoptosic_pyroptosic"]
# change the apoptosis_healthy column name
df1 = df1.rename(columns={"apoptosis_healthy": "apoptosis_healthy_p-adj"})
df1 = df1.rename(columns={"apoptosis_pyroptosis": "apoptosis_pyroptosis_p-adj"})
df1 = df1.rename(columns={"healthy_pyroptosis": "healthy_pyroptosis_p-adj"})


# In[6]:


df1.drop(
    columns=["apoptotic_vs_healthy", "apoptosic_pyroptosic", "pyroptotic_vs_healthy"],
    inplace=True,
)
df1


# In[7]:


# flatten the multiindex
df1.columns = ["".join(col) for col in df1.columns]
df1.reset_index(inplace=True)
df1


# In[8]:


# replace values in the apoptotic column
df1["apoptotic"] = df1["apoptotic"].replace({True: "apoptotic", False: " "})
df1["pyroptotic"] = df1["pyroptotic"].replace({True: "pyroptotic", False: " "})
df1["healthy"] = df1["healthy"].replace({True: "healthy", False: ""})
df1["all_significant"] = df1["all_significant"].replace(
    {True: "all_significant", False: " "}
)
df1["none_significant"] = df1["none_significant"].replace(
    {True: "none_significant", False: " "}
)
# change the order of the columns
df1 = df1[
    [
        "cytokine",
        "apoptosis_healthy_p-adj",
        "apoptosis_pyroptosis_p-adj",
        "healthy_pyroptosis_p-adj",
        "none_significant",
        "apoptotic",
        "pyroptotic",
        "healthy",
        "all_significant",
    ]
]


# In[9]:


# remove NSU from the cytokine names
df1["cytokine"] = df1["cytokine"].str.replace(" \\[NSU\\]", "")
df1.head()


# In[10]:


# list of inflammatory cytokines
inflammatory_cytokines = ["IL-1 beta", "IL-6", "IL-18", "TNF-alpha"]


# In[11]:


df1["putative_function"] = "Not Annotated"
df1.loc[
    df1["cytokine"].isin(inflammatory_cytokines), "putative_function"
] = "Inflammatory"


# In[12]:


# print the table to a csv
df1.to_csv("../results/2023_Interstellar_Table_S1.csv", index=False)


# In[13]:


# csv to markdown
df1.to_markdown("../results/consensus_cytokine_results_cleaned.md", index=False)


# In[14]:


df1.head()


# In[ ]:
