#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pathlib

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# import pillow and open cv
import PIL
import seaborn as sns
import tifffile as tf
from cytocherrypick import cherrypick
from PIL import Image, ImageEnhance

# In[2]:


# parameters
cell_type = "SHSY5Y"
feature = "Nuclei_Texture_SumVariance_CorrGasdermin_3_01_256"


# In[4]:


# define directories
# where the images are
image_dir_path = pathlib.Path(
    "/media/lippincm/c58d4f19-ae4d-4b78-8370-2c2639886da0/interstellar_data/70117_20230210MM1_Gasdermin514_CP_BC430856__2023-03-22T15_42_38-Measurement1/2.IC/"
)
# if path does not exist, create it
image_dir_path.mkdir(parents=True, exist_ok=True)

image_out_dir_path = pathlib.Path("../figures/")
# if path does not exist, create it
image_out_dir_path.mkdir(parents=True, exist_ok=True)


# In[6]:


df_path = pathlib.Path(f"../../data/{cell_type}_sc_norm_fs.parquet")
# read in the data
df = pd.read_parquet(df_path)

df_no_fs_path = pathlib.Path(f"../../data/{cell_type}_sc.parquet")
# read in the data
df_no_fs = pd.read_parquet(df_no_fs_path)


# In[7]:


df["Nuclei_Location_Center_X"] = df_no_fs["Nuclei_Location_Center_X"]
df["Nuclei_Location_Center_Y"] = df_no_fs["Nuclei_Location_Center_Y"]
df["Cytoplasm_AreaShape_BoundingBoxMaximum_X"] = df_no_fs[
    "Cytoplasm_AreaShape_BoundingBoxMaximum_X"
]
df["Cytoplasm_AreaShape_BoundingBoxMaximum_Y"] = df_no_fs[
    "Cytoplasm_AreaShape_BoundingBoxMaximum_Y"
]
df["Cytoplasm_AreaShape_BoundingBoxMinimum_X"] = df_no_fs[
    "Cytoplasm_AreaShape_BoundingBoxMinimum_X"
]
df["Cytoplasm_AreaShape_BoundingBoxMinimum_Y"] = df_no_fs[
    "Cytoplasm_AreaShape_BoundingBoxMinimum_Y"
]


# In[8]:


median, median_index = cherrypick.find_median(df=df, feature_name=feature)
# get the row of the median
median_row = df.iloc[median_index]
image_id = median_row["Metadata_ImageNumber"]
fov_id = median_row["Metadata_Site"].astype(str)
cell_id = median_row["Metadata_Cells_Number_Object_Number"]
well_id = median_row["Metadata_Well"]
row_id = well_id[0]
column_id = well_id[1:]
center_x = median_row["Nuclei_Location_Center_X"]
center_y = median_row["Nuclei_Location_Center_Y"]
# median_row['Cytoplasm_']
# 'Cytoplasm_AreaShape_BoundingBoxMaximum_X', 'Cytoplasm_AreaShape_BoundingBoxMaximum_Y', 'Cytoplasm_AreaShape_BoundingBoxMinimum_X', 'Cytoplasm_AreaShape_BoundingBoxMinimum_Y'
max_x_box = median_row["Cytoplasm_AreaShape_BoundingBoxMaximum_X"]
max_y_box = median_row["Cytoplasm_AreaShape_BoundingBoxMaximum_Y"]
min_x_box = median_row["Cytoplasm_AreaShape_BoundingBoxMinimum_X"]
min_y_box = median_row["Cytoplasm_AreaShape_BoundingBoxMinimum_Y"]
print(median_index)
print(
    f"Median: {median}",
    f"Image ID: {image_id}",
    f"Cell ID: {cell_id}",
    f"Well ID: {well_id}",
    f"row_id: {row_id}",
    f"fov_id: {fov_id}",
    f"column id: {column_id}",
    f"Center X: {center_x}",
    f"Center Y: {center_y}",
    f"Max X: {max_x_box}",
    f"Max Y: {max_y_box}",
    f"Min X: {min_x_box}",
    f"Min Y: {min_y_box}",
    sep="\n",
)


# In[9]:


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


# In[10]:


image_basename_1 = "p04-ch1sk1fk1fl1_IC.tiff"
image_basename_2 = "p04-ch2sk1fk1fl1_IC.tiff"
image_basename_3 = "p04-ch3sk1fk1fl1_IC.tiff"
image_basename_4 = "p04-ch4sk1fk1fl1_IC.tiff"
image_basename_5 = "p04-ch5sk1fk1fl1_IC.tiff"


# In[11]:


print(row_id, column_id, image_id, fov_id)
print(f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}")


# In[12]:


image_name1 = f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}{image_basename_1}"
image_path1 = image_dir_path.joinpath(image_name1)
print(image_name1, "\n", image_path1)

image_name2 = f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}{image_basename_2}"
image_path2 = image_dir_path.joinpath(image_name2)
print(image_name2, "\n", image_path2)

image_name3 = f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}{image_basename_3}"
image_path3 = image_dir_path.joinpath(image_name3)
print(image_name3, "\n", image_path3)

image_name4 = f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}{image_basename_4}"
image_path4 = image_dir_path.joinpath(image_name4)
print(image_name4, "\n", image_path4)

image_name5 = f"r{well_dict[row_id]}c{column_dict[column_id]}f{fov_dict[fov_id]}{image_basename_5}"
image_path5 = image_dir_path.joinpath(image_name5)
print(image_name5, "\n", image_path5)


# In[13]:


# define the radius of the circle to be drawn (actually a square) but a circle fits inside a square so it's all good
# side of square - diameter of circle
# radius = side / 2
radius = 75


# In[14]:


# make each coordinate an integer
center_x = int(center_x)
center_y = int(center_y)

# calculate the bounding box coordinates
min_x_box = center_x - radius
min_y_box = center_y - radius
max_x_box = center_x + radius
max_y_box = center_y + radius

print(center_x, center_y, min_x_box, min_y_box, max_x_box, max_y_box)


# In[15]:


# crop all 5 channels of the image
im1 = cv2.imread(image_path1.as_posix(), cv2.IMREAD_GRAYSCALE)
# im1 = cv2.convertScaleAbs(im1, alpha=alpha, beta=beta)
im_crop1 = im1[min_y_box:max_y_box, min_x_box:max_x_box]

im2 = cv2.imread(image_path2.as_posix(), cv2.IMREAD_GRAYSCALE)
# im2 = cv2.convertScaleAbs(im2, alpha=alpha, beta=beta)
im_crop2 = im2[min_y_box:max_y_box, min_x_box:max_x_box]

im3 = cv2.imread(image_path3.as_posix(), cv2.IMREAD_GRAYSCALE)
# im3 = cv2.convertScaleAbs(im3, alpha=alpha, beta=beta)
im_crop3 = im3[min_y_box:max_y_box, min_x_box:max_x_box]

im4 = cv2.imread(image_path4.as_posix(), cv2.IMREAD_GRAYSCALE)
# im4 = cv2.convertScaleAbs(im4, alpha=alpha, beta=beta)
im_crop4 = im4[min_y_box:max_y_box, min_x_box:max_x_box]

im5 = cv2.imread(image_path5.as_posix(), cv2.IMREAD_GRAYSCALE)
# im5 = cv2.convertScaleAbs(im5, alpha=alpha, beta=beta)
im_crop5 = im5[min_y_box:max_y_box, min_x_box:max_x_box]


# In[16]:


# pick three channels to stack
# nuclei = blue
# Gasdermin = green
# Actin = red

blue_channel_stack = np.stack(im1, axis=-1)
green_channel_stack = np.stack(im3, axis=-1)
red_channel_stack = np.stack(im5, axis=-1)

channel1 = "im1"
channel2 = "im3"
channel3 = "im5"

# Scale the pixel values to fit within the 16-bit range (0-65535)
blue_channel = (blue_channel_stack / np.max(blue_channel_stack) * 65535).astype(
    np.uint16
)
green_channel = (green_channel_stack / np.max(green_channel_stack) * 65535).astype(
    np.uint16
)
red_channel = (red_channel_stack / np.max(red_channel_stack) * 65535).astype(np.uint16)


# In[17]:


composite_image = cv2.merge((blue_channel, green_channel, red_channel)).astype(
    np.uint16
)
composite_image.shape
composite_image = cv2.cvtColor(composite_image, cv2.COLOR_BGR2RGB)


# In[18]:


# transformations of the image to fix the orientation post pixel scaling
# flip the image vertically
composite_image = cv2.flip(composite_image, 0)
# rotate the image 90 degrees clockwise
composite_image = cv2.rotate(composite_image, cv2.ROTATE_90_CLOCKWISE)


# In[19]:


# cv2.imshow("Composite", composite_image)
# cv2.waitKey(0)
# cv2.destroyAllWindows()


# In[20]:


# crop the composite image
im_crop = composite_image[min_y_box:max_y_box, min_x_box:max_x_box]
# cv2.imshow("Composite", im_crop)
# cv2.waitKey(0)
# cv2.destroyAllWindows()


# In[21]:


# image_out_dir_path updated to include the feature name
image_out_dir_path = pathlib.Path(f"{image_out_dir_path}/{feature}")
image_out_dir_path.mkdir(parents=True, exist_ok=True)
# write images
tf.imwrite(
    pathlib.Path(
        f"{image_out_dir_path}/{channel1}_{channel2}_{channel3}_composite_image.tiff"
    ),
    composite_image,
    compression=None,
)
tf.imwrite(
    pathlib.Path(
        f"{image_out_dir_path}/{channel1}_{channel2}_{channel3}_composite_image_crop.tiff"
    ),
    im_crop,
    compression=None,
)


# In[ ]:
