#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import toml
from plotnine import (
    aes,
    element_text,
    facet_grid,
    geom_bar,
    geom_point,
    ggplot,
    ggsave,
    theme,
    theme_bw,
    xlim,
    ylim,
)

# In[2]:


# set paths and load data
path = pathlib.Path(
    "../../Data/clean/Plate2/nELISA_plate_430420_SHSY5Y_cleanup4correlation.csv"
)
toml_path = pathlib.Path("../../../1.Exploratory_Data_Analysis/utils/params.toml")

df = pd.read_csv(path)
params = toml.load(toml_path)
list_of_treatments = params["list_of_treatments"]["treatments"]


# In[3]:


print(df.columns.to_list())


# In[4]:


# output path for the treatment df
output_path = pathlib.Path(
    f"./results/SHSY5Y_all_cytokine_values_per_treatment_per_well.csv"
)
df.to_csv(output_path, index=False)


# In[12]:


# plot scatter plot of all the treatment groups for IL-1 beta

p = (
    ggplot(
        df,
        aes(
            x="TNF alpha [NSU]",
            y="IL-1 beta [NSU]",
            color="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
        ),
    )
    + geom_point(size=3)
    + theme_bw()
    + ylim(0, 1)
    + xlim(0, 1)
)

ggplot.save(
    p,
    filename="./figures/TNF_alpha_IL-1_beta_scatter_plot.png",
    width=6,
    height=4,
    units="in",
    dpi=300,
)
p = p + theme(figure_size=(16, 8))
p


# In[6]:


df_treatment = df.drop(
    columns=["Metadata_position_x", "fourb_Metadata_Treatment_Dose_Inhibitor_Dose"]
)
df_treatment = df_treatment.melt(
    id_vars=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"],
    value_vars=df_treatment.columns.to_list()[1:],
    var_name="Cytokine",
    value_name="Cytokine_Value",
)


# In[7]:


# outpath for the melted df
output_path = pathlib.Path(
    f"./results/SHSY5Y_all_cytokine_values_per_treatment_per_well_melted.csv"
)
df_treatment.to_csv(output_path, index=False)


# In[8]:


# facet grid of treatment and dose with TNF alpha
p = (
    ggplot(
        df_treatment,
        aes(
            x="oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
            y="Cytokine_Value",
            fill="Cytokine",
        ),
    )
    + geom_bar(stat="identity", position="dodge")
    + facet_grid(". ~ Cytokine")
    + theme_bw()
    + theme(axis_text_x=element_text(angle=90, hjust=1))
)

p = p + theme(figure_size=(16, 8))
p
