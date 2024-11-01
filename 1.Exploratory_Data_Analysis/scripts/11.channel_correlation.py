#!/usr/bin/env python
# coding: utf-8

# # Correlate Channels

# This notebook calculates the correlations between channels for this Cell Painting data.

# In[1]:


import itertools
import pathlib
import warnings
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import tqdm

# ## Import data and set paths

# In[2]:


output_dir = pathlib.Path("../output").resolve()
output_dir.mkdir(exist_ok=True, parents=True)

bulk_profile_path = pathlib.Path(
    "../../data/PBMC_preprocess_sc_norm_no_fs_aggregated_nomic.parquet"
).resolve(strict=True)
# import jump data
jump_data_path = pathlib.Path("/home/lippincm/Desktop/18TB/normalized_sc_data").resolve(
    strict=True
)
# get a list of all the jump data files
plates_to_get_channel_correlations_from = list(jump_data_path.glob("*agg*.parquet"))
# set up dictionary for each plate's channels
plates_dict = {}
for plate in plates_to_get_channel_correlations_from:
    plates_dict[plate] = {"channels": ["ER", "Mito", "DNA", "AGP", "RNA"]}
plates_dict[bulk_profile_path] = {
    "channels": ["CorrER", "CorrMito", "CorrDNA", "CorrPM", "CorrGasdermin"]
}
print(len(plates_dict))


# ## Define class to calculate correlations

# In[3]:


class GetChannelCorrelations:
    """
    Class to get the correlations between channels in a profile dataframe

    Parameters
    ----------
    profile_df : pd.DataFrame
        profile dataframe of features
    channels : list, optional
        list of channels to calculate correlations for, by default ["CorrER", "CorrMito", "CorrDNA", "CorrPM", "CorrGasdermin"]
    compartments : list, optional
        list of compartments to calculate correlations for, by default ["Cytoplasm", "Nuclei", "Cell"]

    Methods
    -------
    extract_feature_names()
        Extract the feature names from the profile dataframe
        Returns
        -------
        Tuple[AssertionError, None]
            AssertionError if the length of the feature names is not equal to the length of the compartments, feature types, feature info, and channels
    calculate_correlation()
        Calculate the correlation between all the channels in the profile dataframe
        Returns
        -------
        None
    call_all_class_methods()
        Call all class methods and return the correlation dataframe
        Returns
        -------
        pd.DataFrame
            correlation dataframe
    """

    def __init__(
        self,
        profile_df: pd.DataFrame,
        channels: list = ["CorrER", "CorrMito", "CorrDNA", "CorrPM", "CorrGasdermin"],
        compartments: list = ["Cytoplasm", "Nuclei", "Cell"],
    ) -> None:
        """
        Initialize the class

        Parameters
        ----------
        profile_df : pd.DataFrame
            profile dataframe of features
        channels : list, optional
            list of channels to calculate correlations for, by default ["CorrER", "CorrMito", "CorrDNA", "CorrPM", "CorrGasdermin"]
        compartments : list, optional
            list of compartments to calculate correlations for, by default ["Cytoplasm", "Nuclei", "Cell"]

        Returns
        -------
        None
        """
        self.profile_df = profile_df
        self.channels = channels
        self.compartments = compartments

    def extract_feature_names(self) -> Tuple[AssertionError, None]:
        """
        Extract the feature names from the profile dataframe

        Returns
        -------
        Tuple[AssertionError, None]
            AssertionError if the length of the feature names is not equal to the length of the compartments, feature types, feature info, and channels
        """
        cols = [col for col in self.profile_df.columns]
        dict_feature_names = {
            "feature": [],
            "compartment": [],
            "feature_type": [],
            "feature_info": [],
            "channel": [],
        }
        # iterate through the columns and extract the feature names
        for col in cols:
            if "Metadata" in col:
                continue
            dict_feature_names["feature"].append(col)
            dict_feature_names["compartment"].append(col.split("_")[0])
            feature_type = col.split("_")[1]
            # do not include correlation features in channel comparisions
            if feature_type != "Correlation":
                dict_feature_names["feature_type"].append(feature_type)
                feature_info = col
                dict_feature_names["feature_info"].append(feature_info)
                for channel in self.channels:
                    if channel in feature_info:
                        dict_feature_names["channel"].append(channel)
                        break
                else:
                    dict_feature_names["channel"].append(feature_info)

            else:
                dict_feature_names["feature_type"].append("Correlation")
                dict_feature_names["channel"].append("Correlation")
                feature_info = "_".join(col.split("_")[2:])
                dict_feature_names["feature_info"].append(feature_info)

        assert (
            len(dict_feature_names["feature"])
            == len(dict_feature_names["compartment"])
            == len(dict_feature_names["feature_type"])
            == len(dict_feature_names["feature_info"])
            == len(dict_feature_names["channel"])
        )
        self.extracted_feature_names = pd.DataFrame(dict_feature_names)

    def calculate_correlation(self) -> None:
        """
        Calculate the correlation between all the channels in the profile dataframe

        Returns
        -------
        None
        """
        bulk_profile_long = pd.melt(
            self.profile_df,
            id_vars=["Metadata_Well"],
            value_vars=self.profile_df,
            var_name="feature",
            value_name="feature_value",
        )
        bulk_profile_long = bulk_profile_long.merge(
            self.extracted_feature_names, on="feature"
        )
        # keep well, feature_value, compartment, feature_type, and channel
        # bulk_profile_long = bulk_profile_long[["Metadata_Well", "feature_value", "compartment", "channel"]]
        # drop the channels that are not in the channel list
        bulk_profile_long = bulk_profile_long[
            bulk_profile_long["channel"].isin(self.channels)
        ]
        dict_of_correlations = {}

        index = 0
        for channel1, channel2 in itertools.combinations(self.channels, 2):
            c1_array = np.array(
                bulk_profile_long[bulk_profile_long["channel"] == channel1][
                    "feature_value"
                ]
            )
            c2_array = np.array(
                bulk_profile_long[bulk_profile_long["channel"] == channel2][
                    "feature_value"
                ]
            )
            corr, p = scipy.stats.pearsonr(c1_array, c2_array)  # corr and p value

            # write to dictionary
            dict_of_correlations[index] = {
                "channel1": channel1,
                "channel2": channel2,
                "correlation": corr,
            }
            index += 1
        # add self correlations
        for channel in self.channels:
            dict_of_correlations[index] = {
                "channel1": channel,
                "channel2": channel,
                "correlation": 1.0,
            }
            index += 1
        self.correlation_df = pd.DataFrame(dict_of_correlations).T

    def call_all_class_methods(self) -> pd.DataFrame:
        """
        Call all class methods and return the correlation dataframe

        Returns
        -------
        pd.DataFrame
            correlation dataframe
        """
        self.extract_feature_names()
        self.calculate_correlation()
        return self.correlation_df


# ## Calculate correlations for each plate and save to file

# In[4]:


# supress future warnings


warnings.simplefilter(action="ignore", category=FutureWarning)

df_list = []
for file in tqdm.tqdm(plates_dict):
    get_corr = GetChannelCorrelations(
        profile_df=pd.read_parquet(file),
        channels=plates_dict[file]["channels"],
        compartments=["Cytoplasm", "Nuclei", "Cell"],
    )
    corr_df = get_corr.call_all_class_methods()
    corr_df["plate"] = file.stem
    df_list.append(corr_df)


# In[5]:


all_plate_channel_correlations = pd.concat(df_list)
all_plate_channel_correlations.to_parquet(
    pathlib.Path(output_dir / "all_plate_channel_correlations.parquet"), index=False
)
print(all_plate_channel_correlations.shape)
all_plate_channel_correlations.head()
