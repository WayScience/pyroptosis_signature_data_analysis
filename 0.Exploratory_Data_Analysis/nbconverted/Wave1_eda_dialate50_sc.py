#!/usr/bin/env python
# coding: utf-8

# # Wave1 Exploratory Data Analysis on dilated 50 data

# Here we understand more about the data and how it is distrubuted.
Wave1 data are extracted features from raw images.
These images were processed via Cellprofiler pipelines

Specifically wave1 is looking at Gasdermin-D and Nuclei Staining from a cell painting experiment.

Further, nuclei were dilated using multiple values of pixel dilation. Here we use data for the 50 pixel dialation

# In[16]:


# Import Packages
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# In[17]:


# Import data with low memory arg as the data are large 
df = pd.read_csv("../../Extracted Features (CSV files)/interstellar_wave1_dilate50_sc.csv.gz",low_memory=False)


# In[18]:


# Function to display df shape and # of replicates
def df_stats(df):
    # Print the dimensions of the data
    print('The dimensions of the data are:', df.shape)

    # Print the number of missing values in each column
    print('Number of total missing values across all columns:', (df.isnull().sum()).sum())
    pd.options.display.max_columns = None
    return df.head()


# In[19]:


df_stats(df)


# In[20]:


# Drop na and reindex accordingly
df = df.dropna()
df.reindex()
# Check for Nans again
df_stats(df)


# In[21]:


# Understand categorical data such as treatment and dosing 
df[['Metadata_treatment','Metadata_dose']].drop_duplicates()


# In[22]:


# Subset df where n = number of subset datam points 
df_subset = df.sample(n=1500)

# Define which columns are data and which are descriptive...
df_descriptive = df_subset[["Metadata_wellName",
                     "Metadata_row",
                     "Metadata_col",
                     "Metadata_alias",
                     "Metadata_treatment",
                     "Metadata_dose",
                     "Metadata_ImageNumber",
                     "Metadata_Plate",
                     "Metadata_Well",
                     "Metadata_TranslocatedNuclei_Parent_DilatedNuclei",
                     "Metadata_TranslocatedNuclei_Parent_Nuclei",
                     "Metadata_DilatedNuclei_Number_Object_Number",
                     "Metadata_Nuclei_Number_Object_Number"]]
df_values = df_subset.drop(columns=["Metadata_wellName",
                     "Metadata_row",
                     "Metadata_col",
                     "Metadata_alias",
                     "Metadata_treatment",
                     "Metadata_dose",
                     "Metadata_ImageNumber",
                     "Metadata_Plate",
                     "Metadata_Well",
                     "Metadata_TranslocatedNuclei_Parent_DilatedNuclei",
                     "Metadata_TranslocatedNuclei_Parent_Nuclei",
                     "Metadata_DilatedNuclei_Number_Object_Number",
                     "Metadata_Nuclei_Number_Object_Number"], 
                     axis=1)



treatment_ids = df_descriptive['Metadata_treatment']


# In[23]:


# clustering code adapted from https://www.kaggle.com/code/aussie84/clustering-with-kmeans-pca-tsne
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import plotly_express as px
import plotly.graph_objs as go
import chart_studio.plotly as py
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans


# In[24]:


# Cluster data 
kmeans = KMeans(n_clusters=9)
clustering_ori = kmeans.fit_predict(df_values)

X = df_values
Xtsne = TSNE(n_components=2).fit_transform(X)
dftsneFull = pd.DataFrame(Xtsne)

dftsneFull['cluster'] = clustering_ori
dftsneFull.columns = ['x1','x2','cluster']
dftsneFull['Treatment'] = df_descriptive['Metadata_treatment'].reset_index().drop('index',axis=1)


# In[25]:


# Figure Showing tSNE of Clusters vs Treatment
fig, ax = plt.subplots(1, 2, figsize=(12,6))
plot = sns.scatterplot(data=dftsneFull,x='x1',y='x2',hue='cluster',legend="full",alpha=0.7, ax=ax[0])
ax[0].set_title('Visualized on TSNE')
plot = sns.scatterplot(data=dftsneFull,x='x1',y='x2',hue='Treatment',legend="full",alpha=0.7,ax=ax[1])
ax[1].set_title('Visualized on TSNE')
fig.suptitle('Comparing Clusters vs Treatment tSNE')
display(fig)
df_values['cluster'] = clustering_ori


# Above tSNE shows that based on dimensionality reduction, there is no observable difference in treated cells. More sensitive methods such as machine learning models will need to be employed to achieve such.

# In[26]:


# Plot contributing features

# helper
def outside_limit(df, label_col, label, sensitivity):
  feature_list = df.columns[:-1]
  
  plot_list = []
  mean_overall_list = []
  mean_cluster_list = []
  
  for i,varname in enumerate(feature_list):
    
    #     get overall mean for a variable, set lower and upper limit
    mean_overall = df[varname].mean()
    lower_limit = mean_overall - (mean_overall*sensitivity)
    upper_limit = mean_overall + (mean_overall*sensitivity)

    #     get cluster mean for a variable
    cluster_filter = df[label_col]==label
    pd_cluster = df[cluster_filter]
    mean_cluster = pd_cluster[varname].mean()
    
    #     create filter to display graph with 0.5 deviation from the mean
    if mean_cluster <= lower_limit or mean_cluster >= upper_limit:
      plot_list.append(varname)
      mean_overall_std = mean_overall/mean_overall
      mean_cluster_std = mean_cluster/mean_overall
      mean_overall_list.append(mean_overall_std)
      mean_cluster_list.append(mean_cluster_std)
   
  mean_df = pd.DataFrame({'feature_list':plot_list,
                         'mean_overall_list':mean_overall_list,
                         'mean_cluster_list':mean_cluster_list})
  mean_df = mean_df.sort_values(by=['mean_cluster_list'], ascending=False)
  
  return mean_df
# helped
def plot_barchart_all_unique_features(df, label_col, label, ax, sensitivity):
  
  mean_df = outside_limit(df, label_col, label, sensitivity)
  mean_df_to_plot = mean_df.drop(['mean_overall_list'], axis=1)
  
  if len(mean_df.index) != 0:
    sns.barplot(y='feature_list', x='mean_cluster_list', data=mean_df_to_plot, palette=sns.cubehelix_palette(20, start=.5, rot=-.75, reverse=True),                 alpha=0.75, dodge=True, ax=ax)

    for i,p in enumerate(ax.patches):
      ax.annotate("{:.02f}".format((p.get_width())), 
                  (1, p.get_y() + p.get_height() / 2.), xycoords=('axes fraction', 'data'),
                  ha='right', va='top', fontsize=10, color='black', rotation=0, 
                  xytext=(0, 0),
                  textcoords='offset pixels')
  
  ax.set_title('Unique Characteristics of Cluster ' + str(label))
  ax.set_xlabel('Standardized Mean')
  ax.axvline(x=1, color='k')
# callable function for graphing features that contribute most to each cluster's grouping
# Though the clusters arent grouped via treatment
def plot_features_all_cluster(df, label_col, n_clusters, sensitivity):
  n_plot = n_clusters
  
  fig, ax = plt.subplots(n_plot, 1, figsize=(15, n_plot*6), sharex='col')
  plt.rc('xtick', labelsize=6)
  plt.rc('axes', labelsize=6)
  plt.tick_params(labelsize=4)
  
  ax= ax.ravel()


  label = np.arange(n_clusters)
  for i in label:
    
    plot_barchart_all_unique_features(df, label_col, label=i, ax=ax[i], sensitivity=sensitivity)
    ax[i].xaxis.set_tick_params(labelbottom=True)
    ax[i].yaxis.set_tick_params(labelsize=4)
   

    
  plt.tight_layout()
  
  display(fig)


# In[15]:


plot_features_all_cluster(df=df_values, label_col='cluster', n_clusters=6, sensitivity=0.2)


# In[ ]:




