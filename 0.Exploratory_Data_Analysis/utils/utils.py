#!/home/lippincm/miniconda3/envs/way/bin

import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Function to display df shape and # of replicates
def df_stats(df=pd.DataFrame):
    """Funtion displays shape of df and retruns # of replicates present

    Parameters:
    df (Pandas df): DataFrame of choice

    return:
    string stdout: Shape dimensions of df
    string stdout: # of duplicates present in the df
    df.head: the first 5 rows of the df
    """
    # Print the dimensions of the data
    print("The dimensions of the data are:", df.shape)

    # Print the number of missing values in each column
    print(
        "Number of total missing values across all columns:", (df.isnull().sum()).sum()
    )
    pd.options.display.max_columns = None
    return df.head()


# Plot contributing features
# helper function for plot_features_all_cluster
def outside_limit(df, label_col, label, sensitivity):
    feature_list = df.columns[:-1]

    plot_list = []
    mean_overall_list = []
    mean_cluster_list = []

    for i, varname in enumerate(feature_list):

        # get overall mean for a variable, set lower and upper limit
        mean_overall = df[varname].mean()
        lower_limit = mean_overall - (mean_overall * sensitivity)
        upper_limit = mean_overall + (mean_overall * sensitivity)

        # get cluster mean for a variable
        cluster_filter = df[label_col] == label
        pd_cluster = df[cluster_filter]
        mean_cluster = pd_cluster[varname].mean()

        # create filter to display graph with 0.5 deviation from the mean
        if mean_cluster <= lower_limit or mean_cluster >= upper_limit:
            plot_list.append(varname)
            mean_overall_std = mean_overall / mean_overall
            mean_cluster_std = mean_cluster / mean_overall
            mean_overall_list.append(mean_overall_std)
            mean_cluster_list.append(mean_cluster_std)

    mean_df = pd.DataFrame(
        {
            "feature_list": plot_list,
            "mean_overall_list": mean_overall_list,
            "mean_cluster_list": mean_cluster_list,
        }
    )
    mean_df = mean_df.sort_values(by=["mean_cluster_list"], ascending=False)

    return mean_df


# helper function for plot_features_all_cluster
def plot_barchart_all_unique_features(df, label_col, label, ax, sensitivity):

    mean_df = outside_limit(df, label_col, label, sensitivity)
    mean_df_to_plot = mean_df.drop(["mean_overall_list"], axis=1)

    if len(mean_df.index) != 0:
        sns.barplot(
            y="feature_list",
            x="mean_cluster_list",
            data=mean_df_to_plot,
            palette=sns.cubehelix_palette(20, start=0.5, rot=-0.75, reverse=True),
            alpha=0.75,
            dodge=True,
            ax=ax,
        )

        for i, p in enumerate(ax.patches):
            ax.annotate(
                "{:.02f}".format((p.get_width())),
                (1, p.get_y() + p.get_height() / 2.0),
                xycoords=("axes fraction", "data"),
                ha="right",
                va="top",
                fontsize=10,
                color="black",
                rotation=0,
                xytext=(0, 0),
                textcoords="offset pixels",
            )

    ax.set_title("Unique Characteristics of Cluster " + str(label))
    ax.set_xlabel("Standardized Mean")
    ax.axvline(x=1, color="k")


# callable function for graphing features that contribute most to each cluster's grouping
# Though the clusters arent grouped via treatment
def plot_features_all_cluster(
    df=pd.DataFrame, label_col="cluster", n_clusters=int, sensitivity=float
):
    """Function plots features that most influence cluster

    Parameters:
    df: pandas df with column indicating cluster
    label_col: color for bars
    n_clusters: number of clusters used
    sensitivity: optimizable parameter

    Return:
    type: graph image
    """
    n_plot = n_clusters

    fig, ax = plt.subplots(n_plot, 1, figsize=(15, n_plot * 6), sharex="col")
    plt.rc("xtick", labelsize=6)
    plt.rc("axes", labelsize=6)
    plt.tick_params(labelsize=4)

    ax = ax.ravel()

    label = np.arange(n_clusters)
    for i in label:
        plot_barchart_all_unique_features(
            df, label_col, label=i, ax=ax[i], sensitivity=sensitivity
        )
        ax[i].xaxis.set_tick_params(labelbottom=True)
        ax[i].yaxis.set_tick_params(labelsize=4)

    plt.tight_layout()
