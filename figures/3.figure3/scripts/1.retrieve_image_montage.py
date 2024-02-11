#!/usr/bin/env python
# coding: utf-8

# ## Imports

# In[1]:


import pathlib

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tifffile as tf  # write tiff files
from cytocherrypick.calculations import find_median
from PIL import Image  # read tiff files
from toml import load
from tqdm import tqdm  # progress bar

# In[2]:


CELL_TYPE = "PBMC"


# In[3]:


sc_cell_path = pathlib.Path(f"../../../data/{CELL_TYPE}_preprocessed_sc_norm.parquet")
sc_cell_df = pd.read_parquet(
    sc_cell_path, columns=["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
)

columns_to_load = [
    "Nuclei_Location_Center_Y",
    "Nuclei_Location_Center_X",
]
# get the unfeature selected data
unselected_df_path = pathlib.Path(
    f"../../../data/{CELL_TYPE}_sc.parquet",
)
unselected_df = pd.read_parquet(unselected_df_path, columns=columns_to_load)
# reanme the columns to start with "Metadata_"
unselected_df.columns = [f"Metadata_{x}" for x in unselected_df.columns]
unselected_df.head()


# In[4]:


# add the cell df to the unselected df
sc_cell_df = pd.concat([sc_cell_df, unselected_df], axis="columns")
sc_cell_df


# In[5]:


# Get the current working directory of the repository
cwd = pathlib.Path.cwd()

if (cwd / ".git").is_dir():
    root_dir = cwd

else:
    root_dir = None
    for parent in cwd.parents:
        if (parent / ".git").is_dir():
            root_dir = parent
            break

# Check if a Git root directory was found
if root_dir is None:
    raise FileNotFoundError("No Git root directory found.")
root_dir


# In[6]:


image_out_dir_path = pathlib.Path(
    f"{root_dir}/figures/3.figure3/figures/images/{CELL_TYPE}/"
)
composite_image_out_dir_path = pathlib.Path(
    f"{root_dir}/figures/3.figure3/figures/composite_images/{CELL_TYPE}/"
)

image_out_dir_path.mkdir(parents=True, exist_ok=True)
composite_image_out_dir_path.mkdir(parents=True, exist_ok=True)


# In[7]:


# define directories
# where the images are on a local machine
# this is a hard coded path to the 1TB image directory

#####
# THIS PATH NEEDS TO BE CHANGED TO THE LOCAL IMAGE DIRECTORY ON YOUR MACHINE
#####

image_dir_path = pathlib.Path(
    "/media/lippincm/18T/interstellar_data/70117_20230210MM1_Gasdermin514_CP_BC430856__2023-03-22T15_42_38-Measurement1/2.IC/"
).resolve(strict=True)


# In[8]:


# path
anova_path = pathlib.Path(
    f"../../../1.Exploratory_Data_Analysis/results/{CELL_TYPE}_combined.parquet"
)
# read in the anova results
anova_results = pd.read_parquet(anova_path)


# ## define the groups

# In[9]:


# read in the ground truth data
data_path_ground_truth = (
    "../../../4.sc_Morphology_Neural_Network_MLP_Model/MLP_utils/ground_truth.toml"
)
ground_truth = load(data_path_ground_truth)

# make a a list of the treatments that are in the ground truth data
apoptosis_ground_truth_list = ground_truth["Apoptosis"]["apoptosis_groups_list"]
pyroptosis_ground_truth_list = ground_truth["Pyroptosis"]["pyroptosis_groups_list"]
control_ground_truth_list = ground_truth["Healthy"]["healthy_groups_list"]


# replace Flagellin_1.000_0_DMSO_0.0_% with Flagellin_1.000_ug_per_ml_DMSO_0.025_%
sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_0.100_ug_per_ml_DMSO_0.000_%", "Flagellin_0.100_ug_per_ml_DMSO_0.025_%"
)
# replace Flagellin_1.000_0_DMSO_0.0_% with Flagellin_1.000_ug_per_ml_DMSO_0.025_%
sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_1.000_ug_per_ml_DMSO_0.000_%", "Flagellin_1.000_ug_per_ml_DMSO_0.025_%"
)
sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("Flagellin_1.000_0_DMSO_0.025_%", "Flagellin_1.000_ug_per_ml_DMSO_0.025_%")
# convert media_ctr_0.0_ug_per_ml_Media_ctr_0_0 to media_ctr_0.0_ug_per_ml_Media_ctr_0_025
sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("media_ctr_0.0_ug_per_ml_Media_ctr_0_0", "media_ctr_0.0_0_Media_ctr_0.0_0")


sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace("media_ctr_0.0_0_Media_0_0", "media_ctr_0.0_0_Media_ctr_0.0_0")

sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"] = sc_cell_df[
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose"
].replace(
    "Flagellin_1.000_0_Disulfiram_1.000_uM",
    "Flagellin_1.000_ug_per_ml_Disulfiram_1.000_uM",
)

# make a new column that is the treatment group based on the ground truth data
sc_cell_df["group"] = "NA"
sc_cell_df.loc[
    sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        apoptosis_ground_truth_list
    ),
    "group",
] = "Apoptosis"
sc_cell_df.loc[
    sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        pyroptosis_ground_truth_list
    ),
    "group",
] = "Pyroptosis"
sc_cell_df.loc[
    sc_cell_df["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"].isin(
        control_ground_truth_list
    ),
    "group",
] = "Control"

# make the group column a category
sc_cell_df["group"] = pd.Categorical(
    sc_cell_df["group"],
    categories=["Control", "Apoptosis", "Pyroptosis"],
    ordered=True,
)

print(sc_cell_df["group"].unique())


# In[10]:


# create a column that adds group1 and group2 together
anova_results["group"] = anova_results["group1"] + "_" + anova_results["group2"]
print(anova_results.shape)

# filter out rows that have p-adj_abs > 0.05
anova_results = anova_results[anova_results["p-adj_abs"] < 0.05]
print(anova_results.shape)

# change the group names to replace healthy with control
anova_results["group"] = anova_results["group"].str.replace("healthy", "control")
# make a -log10(p-adj) column
anova_results["neg-log10(p-adj_abs)"] = -np.log10(anova_results["p-adj_abs"])
# sort by neg-log10(p-adj_abs)
anova_results = anova_results.sort_values(by="neg-log10(p-adj_abs)", ascending=False)
# split the dfs into comparisons
c_p_df = anova_results[anova_results["group"] == "control_pyroptosis"]
a_c_df = anova_results[anova_results["group"] == "apoptosis_control"]
a_p_df = anova_results[anova_results["group"] == "apoptosis_pyroptosis"]
# sort by neg-log10(p-adj_abs)
c_p_df = c_p_df.sort_values(by="neg-log10(p-adj_abs)", ascending=False)
a_c_df = a_c_df.sort_values(by="neg-log10(p-adj_abs)", ascending=False)
a_p_df = a_p_df.sort_values(by="neg-log10(p-adj_abs)", ascending=False)


# In[11]:


# get the top 1 features for each comparison
c_p_top1 = c_p_df.iloc[:1, :]
a_c_top1 = a_c_df.iloc[:1, :]
a_p_top1 = a_p_df.iloc[:1, :]

c_p_top1["features"].to_list()
a_c_top1["features"].to_list()
a_p_top1["features"].to_list()
dict_of_top_all = {}
dict_of_top_all["control_pyroptosis"] = c_p_top1["features"].to_list()
dict_of_top_all["apoptosis_control"] = a_c_top1["features"].to_list()
dict_of_top_all["apoptosis_pyroptosis"] = a_p_top1["features"].to_list()

# get list of all the top features
top_features = []
for key in dict_of_top_all:
    top_features.extend(dict_of_top_all[key])
print(len(top_features))
# remove duplicates from the list
top_features = list(set(top_features))
print(len(top_features))
top_features


# In[12]:


# add columns
top_features = top_features + [
    "Metadata_Well",
    "Metadata_Site",
    "Metadata_ImageNumber",
    "Metadata_Cells_Number_Object_Number",
]


# In[13]:


# get features from df
top_features_df = pd.read_parquet(
    sc_cell_path,
    columns=top_features,
)
top_features_df
# merge the top features df with the sc_cell_df
sc_cell_df = pd.concat([sc_cell_df, top_features_df], axis="columns")


# In[14]:


# seperate the data into the different groups
control_df = sc_cell_df[sc_cell_df["group"] == "Control"]
apoptosis_df = sc_cell_df[sc_cell_df["group"] == "Apoptosis"]
pyroptosis_df = sc_cell_df[sc_cell_df["group"] == "Pyroptosis"]


# In[15]:


# define empty dictionary
final_dict = {}


# In[16]:


control_df.head()
# sort the control df by Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1
control_df = control_df.sort_values(
    by="Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1", ascending=False
)
apoptosis_df = apoptosis_df.sort_values(
    by="Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1", ascending=False
)
pyroptosis_df = pyroptosis_df.sort_values(
    by="Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1", ascending=False
)

control_df.reset_index(drop=True, inplace=True)
apoptosis_df.reset_index(drop=True, inplace=True)
pyroptosis_df.reset_index(drop=True, inplace=True)

print(
    control_df["Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1"][
        control_df.last_valid_index()
    ],
    apoptosis_df["Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1"][
        apoptosis_df.last_valid_index()
    ],
    pyroptosis_df["Cytoplasm_RadialDistribution_ZernikePhase_CorrGasdermin_9_1"][0],
)
# get the last item in the control df

dict_of_dfs = {}
dict_of_dfs["control"] = control_df
dict_of_dfs["apoptosis"] = apoptosis_df
dict_of_dfs["pyroptosis"] = pyroptosis_df


# In[17]:


for group in tqdm(dict_of_top_all):
    print(group)
    for dataset in dict_of_dfs:
        feature = dict_of_top_all[group][0]
        key = f"{dataset}__{group}__{feature}"
        df = dict_of_dfs[dataset]
        df = df.sort_values(by=feature, ascending=False, inplace=False)
        df.reset_index(inplace=True, drop=True)
        # get the first and last 3 items in the df
        first_3 = df.head(3)
        last_3 = df.tail(3)
        # add the first and last 3 items to the final dict
        df = pd.concat([first_3, last_3], axis=0)
        print(len(df))
        final_dict[key] = df


# ## Get the images

# In[18]:


# define a dictionary for coding the wells and FOVs correctly
well_dict = {
    "A": "01",
    "B": "02",
    "C": "03",
    "D": "04",
    "E": "05",
    "F": "06",
    "G": "07",
    "H": "08",
    "I": "09",
    "J": "10",
    "K": "11",
    "L": "12",
    "M": "13",
    "N": "14",
    "O": "15",
    "P": "16",
}
column_dict = {
    "1": "01",
    "2": "02",
    "3": "03",
    "4": "04",
    "5": "05",
    "6": "06",
    "7": "07",
    "8": "08",
    "9": "09",
    "10": "10",
    "11": "11",
    "12": "12",
    "13": "13",
    "14": "14",
    "15": "15",
    "16": "16",
    "17": "17",
    "18": "18",
    "19": "19",
    "20": "20",
    "21": "21",
    "22": "22",
    "23": "23",
    "24": "24",
}
fov_dict = {
    "1": "01",
    "2": "02",
    "3": "03",
    "4": "04",
    "5": "05",
    "6": "06",
    "7": "07",
    "8": "08",
    "9": "09",
    "10": "10",
    "11": "11",
    "12": "12",
    "13": "13",
    "14": "14",
    "15": "15",
    "16": "16",
}


# In[19]:


image_basename_1 = "p04-ch1sk1fk1fl1_IC.tiff"
image_basename_2 = "p04-ch2sk1fk1fl1_IC.tiff"
image_basename_3 = "p04-ch3sk1fk1fl1_IC.tiff"
image_basename_4 = "p04-ch4sk1fk1fl1_IC.tiff"
image_basename_5 = "p04-ch5sk1fk1fl1_IC.tiff"


# In[20]:


# set constants for the loop
radius = 50
# define the number of cells to select
n = 5


# In[21]:


# define an empty df
main_df = apoptosis_df.drop(apoptosis_df.index)


# In[22]:


for i in tqdm(final_dict):
    for j in range(len(final_dict[i])):
        tmp_df = pd.DataFrame(final_dict[i].iloc[j]).T
        image_id = tmp_df["Metadata_ImageNumber"].values[0]
        fov_id = tmp_df["Metadata_Site"].values[0]
        cell_id = tmp_df["Metadata_Cells_Number_Object_Number"].values[0]
        well_id = tmp_df["Metadata_Well"].values[0]
        row_id = well_id[0]
        column_id = well_id[1:]
        center_x = tmp_df["Metadata_Nuclei_Location_Center_X"].values[0]
        center_y = tmp_df["Metadata_Nuclei_Location_Center_Y"].values[0]
        # make each of the ids a string
        fov_id = str(fov_id)
        cell_id = str(cell_id)
        well_id = str(well_id)
        row_id = str(row_id)
        column_id = str(column_id)
        center_x = int(center_x)
        center_y = int(center_y)
        treatment = i.split("__")[0]
        comparison = i.split("__")[1]
        feature = i.split("__")[2]
        print(well_id, treatment, comparison, feature)
        # create a custom and contstant bounding box for the images
        # this is made from the extracted center_x and center_y of the cell (nucleus)
        min_x_box = center_x - radius
        max_x_box = center_x + radius
        min_y_box = center_y - radius
        max_y_box = center_y + radius
        print(group, fov_id, cell_id, row_id, column_id, center_x, center_y)
        # create the image paths for each channel of the image
        image_name1 = (
            f"r{well_dict[row_id]}c{column_id}f{fov_dict[fov_id]}{image_basename_1}"
        )
        image_path1 = image_dir_path.joinpath(image_name1)

        image_name2 = (
            f"r{well_dict[row_id]}c{column_id}f{fov_dict[fov_id]}{image_basename_2}"
        )
        image_path2 = image_dir_path.joinpath(image_name2)

        image_name3 = (
            f"r{well_dict[row_id]}c{column_id}f{fov_dict[fov_id]}{image_basename_3}"
        )
        image_path3 = image_dir_path.joinpath(image_name3)

        image_name4 = (
            f"r{well_dict[row_id]}c{column_id}f{fov_dict[fov_id]}{image_basename_4}"
        )
        image_path4 = image_dir_path.joinpath(image_name4)

        image_name5 = (
            f"r{well_dict[row_id]}c{column_id}f{fov_dict[fov_id]}{image_basename_5}"
        )
        image_path5 = image_dir_path.joinpath(image_name5)

        # crop all 5 channels of the image
        im1 = cv2.imread(image_path1.as_posix(), cv2.IMREAD_UNCHANGED)

        im2 = cv2.imread(image_path2.as_posix(), cv2.IMREAD_UNCHANGED)

        im3 = cv2.imread(image_path3.as_posix(), cv2.IMREAD_UNCHANGED)

        im4 = cv2.imread(image_path4.as_posix(), cv2.IMREAD_UNCHANGED)

        im5 = cv2.imread(image_path5.as_posix(), cv2.IMREAD_UNCHANGED)

        # check for non-edge cells

        ### channels ###
        # * Channel 1: DAPI
        # * Channel 2: ER
        # * Channel 3: GasderminD
        # * Channel 4: AGP (Actin, Golgi, and Plasma membrane)
        # * Channel 5: Mitochondria

        # prior to merging adjust the brightness of the image to make it easier to see
        # adjust the brightness of the image to make it easier to see
        alpha = 0.05  # Contrast control (1.0-3.0)
        beta = 0  # Brightness control (0-100)
        im3 = cv2.convertScaleAbs(im3, alpha=alpha, beta=beta)
        im4 = cv2.convertScaleAbs(im4, alpha=alpha, beta=beta)
        # blue channel does not need to be adjusted as it is the DAPI channel and is already bright

        blue_channel_stack = np.stack(im1, axis=-1)
        yellow_channel_stack = np.stack(im2, axis=-1)
        green_channel_stack = np.stack(im3, axis=-1)
        red_channel_stack = np.stack(im4, axis=-1)
        magenta_channel_stack = np.stack(im5, axis=-1)

        channel1 = "im1"
        channel2 = "im3"
        channel3 = "im4"
        channel4 = "im5"
        channel5 = "im2"

        # Scale the pixel values to fit within the 16-bit range (0-65535)
        blue_channel = (blue_channel_stack / np.max(blue_channel_stack) * 65535).astype(
            np.uint16
        )
        yellow_channel = (
            yellow_channel_stack / np.max(yellow_channel_stack) * 65535
        ).astype(np.uint16)
        green_channel = (
            green_channel_stack / np.max(green_channel_stack) * 65535
        ).astype(np.uint16)
        red_channel = (red_channel_stack / np.max(red_channel_stack) * 65535).astype(
            np.uint16
        )
        magenta_channel = (
            magenta_channel_stack / np.max(magenta_channel_stack) * 65535
        ).astype(np.uint16)

        # merge the channels together

        composite_image = cv2.merge((red_channel, green_channel, blue_channel)).astype(
            np.uint16
        )

        # The images end up being `wonky` so we need to do some post processing prior to saving
        # where wonky means that the image is not oriented correctly
        # the image is rotated 90 degrees clockwise and flipped vertically

        # this will ensure that the images are oriented correctly with X and Y centers prior to cropping
        # transformations of the image to fix the orientation post pixel scaling
        # flip the image vertically
        composite_image = cv2.flip(composite_image, 0)
        # rotate the image 90 degrees clockwise
        composite_image = cv2.rotate(composite_image, cv2.ROTATE_90_CLOCKWISE)

        # flip the channels vertically
        blue_channel = cv2.flip(blue_channel, 0)
        yellow_channel = cv2.flip(yellow_channel, 0)
        green_channel = cv2.flip(green_channel, 0)
        red_channel = cv2.flip(red_channel, 0)
        magenta_channel = cv2.flip(magenta_channel, 0)
        # rotate the channels 90 degrees clockwise
        blue_channel = cv2.rotate(blue_channel, cv2.ROTATE_90_CLOCKWISE)
        yellow_channel = cv2.rotate(yellow_channel, cv2.ROTATE_90_CLOCKWISE)
        green_channel = cv2.rotate(green_channel, cv2.ROTATE_90_CLOCKWISE)
        red_channel = cv2.rotate(red_channel, cv2.ROTATE_90_CLOCKWISE)
        magenta_channel = cv2.rotate(magenta_channel, cv2.ROTATE_90_CLOCKWISE)

        composite_image_crop = composite_image[min_y_box:max_y_box, min_x_box:max_x_box]
        # crop the individual channels
        blue_channel_crop = blue_channel[min_y_box:max_y_box, min_x_box:max_x_box]
        yellow_channel_crop = yellow_channel[min_y_box:max_y_box, min_x_box:max_x_box]
        green_channel_crop = green_channel[min_y_box:max_y_box, min_x_box:max_x_box]
        red_channel_crop = red_channel[min_y_box:max_y_box, min_x_box:max_x_box]
        magenta_channel_crop = magenta_channel[min_y_box:max_y_box, min_x_box:max_x_box]

        if composite_image_crop.shape[0] == 0 or composite_image_crop.shape[1] == 0:
            print("Cell is on the edge of the image, skipping")
            continue

            # image_out_dir_path updated to include the feature name
        # write images
        tf.imwrite(
            pathlib.Path(
                f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_cell_{j}.tiff"
            ),
            composite_image,
            compression=None,
        )
        # write each channel as a tiff file
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_blue_channel_cell_{j}.tiff"),
            blue_channel,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_yellow_channel_cell_{j}.tiff"),
            yellow_channel,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_green_channel_cell_{j}.tiff"),
            green_channel,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_red_channel_cell_{j}.tiff"),
            red_channel,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_magenta_channel_cell_{j}.tiff"),
            magenta_channel,
            compression=None,
        )

        # write crops
        tf.imwrite(
            pathlib.Path(
                f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{j}.tiff"
            ),
            composite_image_crop,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_blue_channel_crop_cell_{j}.tiff"),
            blue_channel_crop,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_yellow_channel_crop_cell_{j}.tiff"),
            yellow_channel_crop,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_green_channel_crop_cell_{j}.tiff"),
            green_channel_crop,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(f"{image_out_dir_path}/{i}_red_channel_crop_cell_{j}.tiff"),
            red_channel_crop,
            compression=None,
        )
        tf.imwrite(
            pathlib.Path(
                f"{image_out_dir_path}/{i}_magenta_channel_crop_cell_{j}.tiff"
            ),
            magenta_channel_crop,
            compression=None,
        )

        composite_image = cv2.cvtColor(composite_image, cv2.COLOR_BGR2RGB)
        composite_image_crop = cv2.cvtColor(composite_image_crop, cv2.COLOR_BGR2RGB)

        # save the image as a png file
        cv2.imwrite(
            f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_cell_{j}.png",
            composite_image,
        )
        cv2.imwrite(
            f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{j}.png",
            composite_image_crop,
        )

        tmp_df["comparison"] = comparison
        tmp_df["treatment"] = treatment
        tmp_df["feature"] = feature

        # tmp_df = tmp_df.to_frame().T
        tmp_df[
            "image_compsite_path"
        ] = f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_cell_{j}.png"
        tmp_df[
            "image_composite_crop_path"
        ] = f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{j}.png"

        tmp_df[
            "image_DAPI_path"
        ] = f"{image_out_dir_path}/{i}_blue_channel_cell_{j}.png"
        tmp_df[
            "image_ER_path"
        ] = f"{image_out_dir_path}/{i}_yellow_channel_cell_{j}.png"
        tmp_df[
            "image_GasderminD_path"
        ] = f"{image_out_dir_path}/{i}_green_channel_cell_{j}.png"
        tmp_df["image_AGP_path"] = f"{image_out_dir_path}/{i}_red_channel_cell_{j}.png"
        tmp_df[
            "image_Mitochondria_path"
        ] = f"{image_out_dir_path}/{i}_magenta_channel_cell_{j}.png"

        # crops path
        tmp_df[
            "image_compsite_crop_path"
        ] = f"{composite_image_out_dir_path}/{i}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{j}.png"
        tmp_df[
            "image_DAPI_crop_path"
        ] = f"{image_out_dir_path}/{i}_blue_channel_crop_cell_{j}.png"
        tmp_df[
            "image_ER_crop_path"
        ] = f"{image_out_dir_path}/{i}_yellow_channel_crop_cell_{j}.png"
        tmp_df[
            "image_GasderminD_crop_path"
        ] = f"{image_out_dir_path}/{i}_green_channel_crop_cell_{j}.png"
        tmp_df[
            "image_AGP_crop_path"
        ] = f"{image_out_dir_path}/{i}_red_channel_crop_cell_{j}.png"
        tmp_df[
            "image_Mitochondria_crop_path"
        ] = f"{image_out_dir_path}/{i}_magenta_channel_crop_cell_{j}.png"

        main_df = pd.concat([main_df, tmp_df], ignore_index=True)


# In[23]:


# define main_df_path
main_df_path = pathlib.Path(f"../results/{CELL_TYPE}/")
# if path does not exist, create it
main_df_path.mkdir(parents=True, exist_ok=True)
# save the dataframe
main_df.to_parquet(f"{main_df_path}/single_cell_predictions.parquet")


# In[24]:


# print the number of rows in the df
print(main_df.shape)
main_df
