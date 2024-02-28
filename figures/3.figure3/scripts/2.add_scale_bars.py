#!/usr/bin/env python
# coding: utf-8

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


path = pathlib.Path(f"../results/{CELL_TYPE}/single_cell_predictions.parquet").resolve()

# read in the data
df = pd.read_parquet(path)
# show all columns
pd.set_option("display.max_columns", None)
df.head()


# In[4]:


df.columns
image_path_columns = [
    "image_compsite_path",
    "image_DAPI_path",
    "image_ER_path",
    "image_GasderminD_path",
    "image_AGP_path",
    "image_Mitochondria_path",
    "image_composite_crop_path",
    "image_DAPI_crop_path",
    "image_ER_crop_path",
    "image_GasderminD_crop_path",
    "image_AGP_crop_path",
    "image_Mitochondria_crop_path",
]
for path_type in image_path_columns:
    print(path_type)
    for image_path in df[path_type]:
        print(image_path)
        # open the image
        image = cv2.imread(image_path, cv2.IMREAD_UNCHANGED)
        # image to array
        image_array = np.array(image)
        resolution = 3.3446  # pixels per micrometer
        # get the image size
        image_size = image_array.shape
        if image_array.shape[1] <= 100:
            scale_bar_length = 5  # um
            scale_bar_height = 1  # pixels
            padding = 3  # pixels
        elif image_array.shape[1] > 100:
            scale_bar_length = 100
            scale_bar_height = 10
            padding = 10
        # get the bottom right most corner based on scale 1 % of pixels
        scale_bar_x = (
            image_size[1] - (scale_bar_length * resolution) - padding - padding
        )
        scale_bar_y = image_size[0] - (scale_bar_height) - padding
        print(scale_bar_x, scale_bar_y)
        # draw the scale bar
        new = cv2.rectangle(
            image_array,
            (int(scale_bar_x), int(scale_bar_y)),
            (image_size[1] - padding, image_size[0] - padding),
            (255, 255, 255),
            -1,
        )
        # show the image
        # plt.imshow(new)
        # plt.show()
        # save the image with the scale bar via cv2
        cv2.imwrite(image_path, new)
#
