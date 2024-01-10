#!/usr/bin/env python
# coding: utf-8

# This notebook finds random cells from each prediction category and displays them. The purpose is to get representative images examples of each category.

# In[1]:


import pathlib

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile as tf  # write tiff files
from PIL import Image  # read tiff files
from tqdm import tqdm  # progress bar

# In[2]:


# function that selects a random image from the dataframe


def random_cell_select(
    df: pd.DataFrame,
    n: int = 1,
) -> pd.DataFrame:
    """
    Selects a random cell from the dataframe

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the cell features
    n : int, optional
        Number of random cells to select, by default 1

    Returns
    -------
    pd.DataFrame
        The return dataframe with the random cell selected
    """

    # select a random cell
    random_cell = df.sample(n=n, random_state=0)
    return random_cell


# In[3]:


# parameters
CELL_TYPE = "PBMC"


# In[4]:


# Get the current working directory
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


# In[5]:


image_out_dir_path = pathlib.Path(f"{root_dir}/8.cytopick_analysis/figures/PBMC/")


# In[6]:


# define directories
# where the images are
image_dir_path = pathlib.Path(
    "/media/lippincm/18T/interstellar_data/70117_20230210MM1_Gasdermin514_CP_BC430856__2023-03-22T15_42_38-Measurement1/2.IC/"
).resolve(strict=True)


# if path does not exist, create it
image_out_dir_path.mkdir(parents=True, exist_ok=True)


# In[7]:


df_path = pathlib.Path(
    f"../../4.sc_Morphology_Neural_Network_MLP_Model/results/Multi_Class/MultiClass_MLP/{CELL_TYPE}/single_cell_predictions.parquet"
)
# read in the data
df = pd.read_parquet(df_path)
df.head()


# In[8]:


# add column for if the prediction was correct
df["correct"] = df.apply(lambda x: x["true_label"] == x["predicted_label"], axis=1)
# split the data into correct and incorrect
df_correct = df[df["correct"] == True]
df_incorrect = df[df["correct"] == False]
assert len(df_correct) + len(df_incorrect) == len(df)


# In[9]:


# split the data into the different classes
pyroptosis_df = df_correct[df_correct["labels"] == "pyroptosis"]
apoptosis_df = df_correct[df_correct["labels"] == "apoptosis"]
control_df = df_correct[df_correct["labels"] == "healthy"]

# split the data classes by shuffled and unshuffled
pyroptosis_shuffled_df = pyroptosis_df[pyroptosis_df["shuffle"] == True]
pyroptosis_unshuffled_df = pyroptosis_df[pyroptosis_df["shuffle"] == False]
apoptosis_shuffled_df = apoptosis_df[apoptosis_df["shuffle"] == True]
apoptosis_unshuffled_df = apoptosis_df[apoptosis_df["shuffle"] == False]
control_shuffled_df = control_df[control_df["shuffle"] == True]
control_unshuffled_df = control_df[control_df["shuffle"] == False]

# split the shuffled/unshuffled data by the data splits
pyroptosis_shuffled_train_df = pyroptosis_shuffled_df[
    pyroptosis_shuffled_df["data_split"] == "train"
]
pyroptosis_shuffled_test_df = pyroptosis_shuffled_df[
    pyroptosis_shuffled_df["data_split"] == "test"
]
pyroptosis_shuffled_validation_df = pyroptosis_shuffled_df[
    pyroptosis_shuffled_df["data_split"] == "validation"
]
pyroptosis_shuffled_treatment_holdout_df = pyroptosis_shuffled_df[
    pyroptosis_shuffled_df["data_split"] == "treatment_holdout"
]
pyroptosis_shuffled_holdout_df = pyroptosis_shuffled_df[
    pyroptosis_shuffled_df["data_split"] == "holdout"
]

pyroptosis_unshuffled_train_df = pyroptosis_unshuffled_df[
    pyroptosis_unshuffled_df["data_split"] == "train"
]
pyroptosis_unshuffled_test_df = pyroptosis_unshuffled_df[
    pyroptosis_unshuffled_df["data_split"] == "test"
]
pyroptosis_unshuffled_validation_df = pyroptosis_unshuffled_df[
    pyroptosis_unshuffled_df["data_split"] == "validation"
]
pyroptosis_unshuffled_treatment_holdout_df = pyroptosis_unshuffled_df[
    pyroptosis_unshuffled_df["data_split"] == "treatment_holdout"
]
pyroptosis_unshuffled_holdout_df = pyroptosis_unshuffled_df[
    pyroptosis_unshuffled_df["data_split"] == "holdout"
]

apoptosis_shuffled_train_df = apoptosis_shuffled_df[
    apoptosis_shuffled_df["data_split"] == "train"
]
apoptosis_shuffled_test_df = apoptosis_shuffled_df[
    apoptosis_shuffled_df["data_split"] == "test"
]
apoptosis_shuffled_validation_df = apoptosis_shuffled_df[
    apoptosis_shuffled_df["data_split"] == "validation"
]
apoptosis_shuffled_treatment_holdout_df = apoptosis_shuffled_df[
    apoptosis_shuffled_df["data_split"] == "treatment_holdout"
]
apoptosis_shuffled_holdout_df = apoptosis_shuffled_df[
    apoptosis_shuffled_df["data_split"] == "holdout"
]

apoptosis_unshuffled_train_df = apoptosis_unshuffled_df[
    apoptosis_unshuffled_df["data_split"] == "train"
]
apoptosis_unshuffled_test_df = apoptosis_unshuffled_df[
    apoptosis_unshuffled_df["data_split"] == "test"
]
apoptosis_unshuffled_validation_df = apoptosis_unshuffled_df[
    apoptosis_unshuffled_df["data_split"] == "validation"
]
apoptosis_unshuffled_treatment_holdout_df = apoptosis_unshuffled_df[
    apoptosis_unshuffled_df["data_split"] == "treatment_holdout"
]
apoptosis_unshuffled_holdout_df = apoptosis_unshuffled_df[
    apoptosis_unshuffled_df["data_split"] == "holdout"
]

control_shuffled_train_df = control_shuffled_df[
    control_shuffled_df["data_split"] == "train"
]
control_shuffled_test_df = control_shuffled_df[
    control_shuffled_df["data_split"] == "test"
]
control_shuffled_validation_df = control_shuffled_df[
    control_shuffled_df["data_split"] == "validation"
]
control_shuffled_treatment_holdout_df = control_shuffled_df[
    control_shuffled_df["data_split"] == "treatment_holdout"
]
control_shuffled_holdout_df = control_shuffled_df[
    control_shuffled_df["data_split"] == "holdout"
]

control_unshuffled_train_df = control_unshuffled_df[
    control_unshuffled_df["data_split"] == "train"
]
control_unshuffled_test_df = control_unshuffled_df[
    control_unshuffled_df["data_split"] == "test"
]
control_unshuffled_validation_df = control_unshuffled_df[
    control_unshuffled_df["data_split"] == "validation"
]
control_unshuffled_treatment_holdout_df = control_unshuffled_df[
    control_unshuffled_df["data_split"] == "treatment_holdout"
]
control_unshuffled_holdout_df = control_unshuffled_df[
    control_unshuffled_df["data_split"] == "holdout"
]

# add each df to a dictionary
dict_of_dfs = {}
dict_of_dfs["pyroptosis_shuffled_train_df"] = pyroptosis_shuffled_train_df
dict_of_dfs["pyroptosis_shuffled_test_df"] = pyroptosis_shuffled_test_df
dict_of_dfs["pyroptosis_shuffled_validation_df"] = pyroptosis_shuffled_validation_df
dict_of_dfs[
    "pyroptosis_shuffled_treatment_holdout_df"
] = pyroptosis_shuffled_treatment_holdout_df
dict_of_dfs["pyroptosis_shuffled_holdout_df"] = pyroptosis_shuffled_holdout_df

dict_of_dfs["pyroptosis_unshuffled_train_df"] = pyroptosis_unshuffled_train_df
dict_of_dfs["pyroptosis_unshuffled_test_df"] = pyroptosis_unshuffled_test_df
dict_of_dfs["pyroptosis_unshuffled_validation_df"] = pyroptosis_unshuffled_validation_df
dict_of_dfs[
    "pyroptosis_unshuffled_treatment_holdout_df"
] = pyroptosis_unshuffled_treatment_holdout_df
dict_of_dfs["pyroptosis_unshuffled_holdout_df"] = pyroptosis_unshuffled_holdout_df

dict_of_dfs["apoptosis_shuffled_train_df"] = apoptosis_shuffled_train_df
dict_of_dfs["apoptosis_shuffled_test_df"] = apoptosis_shuffled_test_df
dict_of_dfs["apoptosis_shuffled_validation_df"] = apoptosis_shuffled_validation_df
dict_of_dfs[
    "apoptosis_shuffled_treatment_holdout_df"
] = apoptosis_shuffled_treatment_holdout_df
dict_of_dfs["apoptosis_shuffled_holdout_df"] = apoptosis_shuffled_holdout_df

dict_of_dfs["apoptosis_unshuffled_train_df"] = apoptosis_unshuffled_train_df
dict_of_dfs["apoptosis_unshuffled_test_df"] = apoptosis_unshuffled_test_df
dict_of_dfs["apoptosis_unshuffled_validation_df"] = apoptosis_unshuffled_validation_df
dict_of_dfs[
    "apoptosis_unshuffled_treatment_holdout_df"
] = apoptosis_unshuffled_treatment_holdout_df
dict_of_dfs["apoptosis_unshuffled_holdout_df"] = apoptosis_unshuffled_holdout_df

dict_of_dfs["control_shuffled_train_df"] = control_shuffled_train_df
dict_of_dfs["control_shuffled_test_df"] = control_shuffled_test_df
dict_of_dfs["control_shuffled_validation_df"] = control_shuffled_validation_df
dict_of_dfs[
    "control_shuffled_treatment_holdout_df"
] = control_shuffled_treatment_holdout_df
dict_of_dfs["control_shuffled_holdout_df"] = control_shuffled_holdout_df

dict_of_dfs["control_unshuffled_train_df"] = control_unshuffled_train_df
dict_of_dfs["control_unshuffled_test_df"] = control_unshuffled_test_df
dict_of_dfs["control_unshuffled_validation_df"] = control_unshuffled_validation_df
dict_of_dfs[
    "control_unshuffled_treatment_holdout_df"
] = control_unshuffled_treatment_holdout_df
dict_of_dfs["control_unshuffled_holdout_df"] = control_unshuffled_holdout_df

# check the length of each df
for key, value in dict_of_dfs.items():
    if not len(dict_of_dfs[key]) == 0:
        pass
    else:
        print(key)


# In[10]:


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


# In[11]:


image_basename_1 = "p04-ch1sk1fk1fl1_IC.tiff"
image_basename_2 = "p04-ch2sk1fk1fl1_IC.tiff"
image_basename_3 = "p04-ch3sk1fk1fl1_IC.tiff"
image_basename_4 = "p04-ch4sk1fk1fl1_IC.tiff"
image_basename_5 = "p04-ch5sk1fk1fl1_IC.tiff"


# In[12]:


# set constants for the loop
radius = 50
# define the number of cells to select
n = 5


# In[13]:


dict_of_subset_dfs = {}
for key in tqdm(dict_of_dfs):
    df = dict_of_dfs[key]
    if len(df) == 0:
        pass
    else:
        # select n random cells from the dataframe
        df = random_cell_select(df, n)
        # add the df to the dictionary
        dict_of_subset_dfs[key] = df


# In[14]:


# create a blank df to append the data to
main_df = dict_of_subset_dfs["pyroptosis_shuffled_train_df"]
# drop all rows from the df
main_df = main_df.drop(main_df.index)


# In[15]:


for key in tqdm(dict_of_subset_dfs):
    if len(dict_of_subset_dfs[key]) >= 1:
        # loop through the dataframe
        for cell in range(len(dict_of_subset_dfs[key])):
            # get the first row of the dataframe
            df = dict_of_subset_dfs[key].iloc[cell]
            image_id = df["Metadata_ImageNumber"]
            fov_id = df["Metadata_Site"].astype(str)
            cell_id = df["Metadata_Cells_Number_Object_Number"]
            well_id = df["Metadata_Well"]
            row_id = well_id[0]
            column_id = well_id[1:]
            center_x = df["Metadata_Nuclei_Location_Center_X"].astype(int)
            center_y = df["Metadata_Nuclei_Location_Center_Y"].astype(int)
            # create a custom and contstant bounding box for the images
            # this is made from the extracted center_x and center_y of the cell (nucleus)
            min_x_box = center_x - radius
            max_x_box = center_x + radius
            min_y_box = center_y - radius
            max_y_box = center_y + radius
            print(cell + 1, key, row_id, column_id, fov_id, cell_id, center_x, center_y)

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
            # im1_crop = im1[min_y_box:max_y_box, min_x_box:max_x_box]

            im2 = cv2.imread(image_path2.as_posix(), cv2.IMREAD_UNCHANGED)
            # im2_crop = im2[min_y_box:max_y_box, min_x_box:max_x_box]

            im3 = cv2.imread(image_path3.as_posix(), cv2.IMREAD_UNCHANGED)
            # im3_crop = im3[min_y_box:max_y_box, min_x_box:max_x_box]

            im4 = cv2.imread(image_path4.as_posix(), cv2.IMREAD_UNCHANGED)
            # im4_crop = im4[min_y_box:max_y_box, min_x_box:max_x_box]

            im5 = cv2.imread(image_path5.as_posix(), cv2.IMREAD_UNCHANGED)
            # im5_crop = im5[min_y_box:max_y_box, min_x_box:max_x_box]

            # check for non-edge cells

            ### channels ###
            # * Channel 1: DAPI
            # * Channel 2: ER
            # * Channel 3: GasderminD
            # * Channel 4: AGP (Actin, Golgi, and Plasma membrane)
            # * Channel 5: Mitochondria

            blue_channel_stack = np.stack(im1, axis=-1)
            green_channel_stack = np.stack(im3, axis=-1)
            red_channel_stack = np.stack(im4, axis=-1)

            # blue_channel_stack_crop = np.stack(im1_crop, axis=-1)
            # green_channel_stack_crop = np.stack(im3_crop, axis=-1)
            # red_channel_stack_crop = np.stack(im4_crop, axis=-1)

            channel1 = "im1"
            channel2 = "im3"
            channel3 = "im4"

            # Scale the pixel values to fit within the 16-bit range (0-65535)
            blue_channel = (
                blue_channel_stack / np.max(blue_channel_stack) * 65535
            ).astype(np.uint16)
            green_channel = (
                green_channel_stack / np.max(green_channel_stack) * 65535
            ).astype(np.uint16)
            red_channel = (
                red_channel_stack / np.max(red_channel_stack) * 65535
            ).astype(np.uint16)

            # blue_channel_crop = (
            #     blue_channel_stack_crop / np.max(blue_channel_stack_crop) * 65535
            # ).astype(np.uint16)
            # green_channel_crop = (
            #     green_channel_stack_crop / np.max(green_channel_stack_crop) * 65535
            # ).astype(np.uint16)
            # red_channel_crop = (
            #     red_channel_stack_crop / np.max(red_channel_stack_crop) * 65535
            # ).astype(np.uint16)

            # merge the channels together

            composite_image = cv2.merge(
                (red_channel, green_channel, blue_channel)
            ).astype(np.uint16)
            # composite_image = cv2.cvtColor(composite_image, cv2.COLOR_BGR2RGB)

            composite_image_crop = composite_image[
                min_y_box:max_y_box, min_x_box:max_x_box
            ]

            # im_crop = composite_image[min_y_box:max_y_box, min_x_box:max_x_box]
            if composite_image_crop.shape[0] == 0 or composite_image_crop.shape[1] == 0:
                print("Cell is on the edge of the image, skipping")
                continue

            # composite_image_crop = cv2.merge(
            #     (blue_channel_crop, green_channel_crop, red_channel_crop)
            # ).astype(np.uint16)
            composite_image = cv2.cvtColor(composite_image, cv2.COLOR_BGR2RGB)
            composite_image_crop = cv2.cvtColor(composite_image_crop, cv2.COLOR_BGR2RGB)

            # The images end up being `wonky` so we need to do some post processing prior to saving
            # this will ensure that the images are oriented correctly with X and Y centers prior to cropping
            # transformations of the image to fix the orientation post pixel scaling
            # flip the image vertically
            composite_image = cv2.flip(composite_image, 0)
            composite_image_crop = cv2.flip(composite_image_crop, 0)
            # rotate the image 90 degrees clockwise
            composite_image = cv2.rotate(composite_image, cv2.ROTATE_90_CLOCKWISE)
            composite_image_crop = cv2.rotate(
                composite_image_crop, cv2.ROTATE_90_CLOCKWISE
            )

            print(composite_image.shape)

            # save the image as a png file
            cv2.imwrite(
                f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_cell_{cell}.png",
                composite_image,
            )
            cv2.imwrite(
                f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{cell}.png",
                composite_image_crop,
            )

            # image_out_dir_path updated to include the feature name
            # write images
            tf.imwrite(
                pathlib.Path(
                    f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_cell_{cell}.tiff"
                ),
                composite_image,
                compression=None,
            )
            tf.imwrite(
                pathlib.Path(
                    f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{cell}.tiff"
                ),
                composite_image_crop,
                compression=None,
            )
            df = df.to_frame().T
            df[
                "image_path"
            ] = f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_cell_{cell}.png"
            df[
                "image_crop_path"
            ] = f"{image_out_dir_path}/{key}_{channel1}_{channel2}_{channel3}_composite_image_crop_cell_{cell}.png"
            main_df = pd.concat([main_df, df], ignore_index=True)


# In[16]:


# define main_df_path
main_df_path = pathlib.Path(f"../results/{CELL_TYPE}/")
# if path does not exist, create it
main_df_path.mkdir(parents=True, exist_ok=True)
# save the dataframe
main_df.to_parquet(f"{main_df_path}/single_cell_predictions.parquet")


# In[17]:


main_df


# In[ ]:
