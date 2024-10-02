#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import pathlib
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import toml
import tqdm
from joblib import dump
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV, LogisticRegression, MultiTaskElasticNetCV

# import RepeatedKFold
from sklearn.model_selection import (
    GridSearchCV,
    LeaveOneOut,
    RepeatedKFold,
    StratifiedKFold,
    cross_val_score,
    train_test_split,
)
from sklearn.utils import parallel_backend

# In[2]:


# argparser = argparse.ArgumentParser()
# argparser.add_argument("--cell_type", type=str, default="all")
# argparser.add_argument("--shuffle", type=str, default=False)
# argparser.add_argument("--cytokine", type=str, default="cytokine")

# args = argparser.parse_args()

# cell_type = args.cell_type
# cytokine = args.cytokine
# shuffle = args.shuffle

cell_type = "PBMC"
cytokine = "IL-1 beta [NSU]"
shuffle = "False"
columns = [
    "Cytoplasm_AreaShape_Compactness",
    "Cytoplasm_AreaShape_FormFactor",
    "Cytoplasm_AreaShape_MajorAxisLength",
    "Cytoplasm_AreaShape_MinorAxisLength",
    "Cytoplasm_AreaShape_Orientation",
    "Cytoplasm_AreaShape_Zernike_0_0",
    "Cytoplasm_AreaShape_Zernike_1_1",
    "Cytoplasm_AreaShape_Zernike_2_0",
    "Cytoplasm_AreaShape_Zernike_2_2",
    "Cytoplasm_AreaShape_Zernike_3_1",
    "Cytoplasm_AreaShape_Zernike_3_3",
    "Cytoplasm_AreaShape_Zernike_4_0",
    "Cytoplasm_AreaShape_Zernike_4_2",
    "Cytoplasm_AreaShape_Zernike_4_4",
    "Cytoplasm_AreaShape_Zernike_5_1",
    "Cytoplasm_AreaShape_Zernike_5_3",
    "Cytoplasm_AreaShape_Zernike_5_5",
    "Cytoplasm_AreaShape_Zernike_6_0",
    "Cytoplasm_AreaShape_Zernike_6_2",
    "Cytoplasm_AreaShape_Zernike_6_4",
    "Cytoplasm_AreaShape_Zernike_6_6",
    "Cytoplasm_AreaShape_Zernike_7_1",
    "Cytoplasm_AreaShape_Zernike_7_3",
    "Cytoplasm_AreaShape_Zernike_7_5",
    "Cytoplasm_AreaShape_Zernike_7_7",
    "Cytoplasm_AreaShape_Zernike_8_0",
    "Cytoplasm_AreaShape_Zernike_8_2",
    "Cytoplasm_AreaShape_Zernike_8_4",
    "Cytoplasm_AreaShape_Zernike_8_6",
    "Cytoplasm_AreaShape_Zernike_9_1",
    "Cytoplasm_AreaShape_Zernike_9_3",
    "Cytoplasm_AreaShape_Zernike_9_5",
    "Cytoplasm_AreaShape_Zernike_9_7",
    "Cells_AreaShape_Compactness",
    "Cells_AreaShape_Eccentricity",
    "Cells_AreaShape_FormFactor",
    "Cells_AreaShape_MeanRadius",
    "Cells_AreaShape_Orientation",
    "Cells_AreaShape_Zernike_0_0",
    "Cells_AreaShape_Zernike_1_1",
    "Cells_AreaShape_Zernike_2_0",
    "Cells_AreaShape_Zernike_2_2",
    "Cells_AreaShape_Zernike_3_1",
    "Cells_AreaShape_Zernike_3_3",
    "Cells_AreaShape_Zernike_4_0",
    "Cells_AreaShape_Zernike_4_2",
    "Cells_AreaShape_Zernike_4_4",
    "Cells_AreaShape_Zernike_5_1",
    "Cells_AreaShape_Zernike_5_3",
    "Cells_AreaShape_Zernike_5_5",
    "Cells_AreaShape_Zernike_6_0",
    "Cells_AreaShape_Zernike_6_2",
    "Cells_AreaShape_Zernike_6_4",
    "Cells_AreaShape_Zernike_6_6",
    "Cells_AreaShape_Zernike_7_1",
    "Cells_AreaShape_Zernike_7_3",
    "Cells_AreaShape_Zernike_7_5",
    "Cells_AreaShape_Zernike_7_7",
    "Cells_AreaShape_Zernike_8_0",
    "Cells_AreaShape_Zernike_8_2",
    "Cells_AreaShape_Zernike_8_4",
    "Cells_AreaShape_Zernike_8_6",
    "Cells_AreaShape_Zernike_8_8",
    "Cells_AreaShape_Zernike_9_1",
    "Cells_AreaShape_Zernike_9_3",
    "Cells_AreaShape_Zernike_9_5",
    "Cells_AreaShape_Zernike_9_7",
    "Cells_AreaShape_Zernike_9_9",
    "Cells_Neighbors_AngleBetweenNeighbors_Adjacent",
    "Cells_Neighbors_FirstClosestDistance_Adjacent",
    "Cells_Neighbors_SecondClosestDistance_Adjacent",
    "Nuclei_AreaShape_MajorAxisLength",
    "Nuclei_AreaShape_MinorAxisLength",
    "Nuclei_AreaShape_Orientation",
    "Nuclei_AreaShape_Zernike_0_0",
    "Nuclei_AreaShape_Zernike_1_1",
    "Nuclei_AreaShape_Zernike_2_0",
    "Nuclei_AreaShape_Zernike_2_2",
    "Nuclei_AreaShape_Zernike_3_1",
    "Nuclei_AreaShape_Zernike_3_3",
    "Nuclei_AreaShape_Zernike_4_0",
    "Nuclei_AreaShape_Zernike_4_2",
    "Nuclei_AreaShape_Zernike_4_4",
    "Nuclei_AreaShape_Zernike_5_1",
    "Nuclei_AreaShape_Zernike_5_3",
    "Nuclei_AreaShape_Zernike_5_5",
    "Nuclei_AreaShape_Zernike_6_0",
    "Nuclei_AreaShape_Zernike_6_2",
    "Nuclei_AreaShape_Zernike_6_4",
    "Nuclei_AreaShape_Zernike_6_6",
    "Nuclei_AreaShape_Zernike_7_1",
    "Nuclei_AreaShape_Zernike_7_3",
    "Nuclei_AreaShape_Zernike_7_5",
    "Nuclei_AreaShape_Zernike_7_7",
    "Nuclei_AreaShape_Zernike_8_0",
    "Nuclei_AreaShape_Zernike_8_2",
    "Nuclei_AreaShape_Zernike_8_4",
    "Nuclei_AreaShape_Zernike_8_6",
    "Nuclei_AreaShape_Zernike_8_8",
    "Nuclei_AreaShape_Zernike_9_1",
    "Nuclei_AreaShape_Zernike_9_3",
    "Nuclei_AreaShape_Zernike_9_5",
    "Nuclei_AreaShape_Zernike_9_7",
    "Nuclei_AreaShape_Zernike_9_9",
    "Nuclei_Neighbors_AngleBetweenNeighbors_Adjacent",
    "Nuclei_Neighbors_FirstClosestDistance_Adjacent",
    "Nuclei_Neighbors_SecondClosestDistance_Adjacent",
    "oneb_Treatment_Dose_Inhibitor_Dose",
    "Cytoplasm_AreaShape_Compactness",
    "Cytoplasm_AreaShape_FormFactor",
    "Cytoplasm_AreaShape_MajorAxisLength",
    "Cytoplasm_AreaShape_MinorAxisLength",
    "Cytoplasm_AreaShape_Orientation",
    "Cytoplasm_AreaShape_Zernike_0_0",
    "Cytoplasm_AreaShape_Zernike_1_1",
    "Cytoplasm_AreaShape_Zernike_2_0",
    "Cytoplasm_AreaShape_Zernike_2_2",
    "Cytoplasm_AreaShape_Zernike_3_1",
    "Cytoplasm_AreaShape_Zernike_3_3",
    "Cytoplasm_AreaShape_Zernike_4_0",
    "Cytoplasm_AreaShape_Zernike_4_2",
    "Cytoplasm_AreaShape_Zernike_4_4",
    "Cytoplasm_AreaShape_Zernike_5_1",
    "Cytoplasm_AreaShape_Zernike_5_3",
    "Cytoplasm_AreaShape_Zernike_5_5",
    "Cytoplasm_AreaShape_Zernike_6_0",
    "Cytoplasm_AreaShape_Zernike_6_2",
    "Cytoplasm_AreaShape_Zernike_6_4",
    "Cytoplasm_AreaShape_Zernike_6_6",
    "Cytoplasm_AreaShape_Zernike_7_1",
    "Cytoplasm_AreaShape_Zernike_7_3",
    "Cytoplasm_AreaShape_Zernike_7_5",
    "Cytoplasm_AreaShape_Zernike_7_7",
    "Cytoplasm_AreaShape_Zernike_8_0",
    "Cytoplasm_AreaShape_Zernike_8_2",
    "Cytoplasm_AreaShape_Zernike_8_4",
    "Cytoplasm_AreaShape_Zernike_8_6",
    "Cytoplasm_AreaShape_Zernike_9_1",
    "Cytoplasm_AreaShape_Zernike_9_3",
    "Cytoplasm_AreaShape_Zernike_9_5",
    "Cytoplasm_AreaShape_Zernike_9_7",
    "Cells_AreaShape_Compactness",
    "Cells_AreaShape_Eccentricity",
    "Cells_AreaShape_FormFactor",
    "Cells_AreaShape_MeanRadius",
    "Cells_AreaShape_Orientation",
    "Cells_AreaShape_Zernike_0_0",
    "Cells_AreaShape_Zernike_1_1",
    "Cells_AreaShape_Zernike_2_0",
    "Cells_AreaShape_Zernike_2_2",
    "Cells_AreaShape_Zernike_3_1",
    "Cells_AreaShape_Zernike_3_3",
    "Cells_AreaShape_Zernike_4_0",
    "Cells_AreaShape_Zernike_4_2",
    "Cells_AreaShape_Zernike_4_4",
    "Cells_AreaShape_Zernike_5_1",
    "Cells_AreaShape_Zernike_5_3",
    "Cells_AreaShape_Zernike_5_5",
    "Cells_AreaShape_Zernike_6_0",
    "Cells_AreaShape_Zernike_6_2",
    "Cells_AreaShape_Zernike_6_4",
    "Cells_AreaShape_Zernike_6_6",
    "Cells_AreaShape_Zernike_7_1",
    "Cells_AreaShape_Zernike_7_3",
    "Cells_AreaShape_Zernike_7_5",
    "Cells_AreaShape_Zernike_7_7",
    "Cells_AreaShape_Zernike_8_0",
    "Cells_AreaShape_Zernike_8_2",
    "Cells_AreaShape_Zernike_8_4",
    "Cells_AreaShape_Zernike_8_6",
    "Cells_AreaShape_Zernike_8_8",
    "Cells_AreaShape_Zernike_9_1",
    "Cells_AreaShape_Zernike_9_3",
    "Cells_AreaShape_Zernike_9_5",
    "Cells_AreaShape_Zernike_9_7",
    "Cells_AreaShape_Zernike_9_9",
    "Cells_Neighbors_AngleBetweenNeighbors_Adjacent",
    "Cells_Neighbors_FirstClosestDistance_Adjacent",
    "Cells_Neighbors_SecondClosestDistance_Adjacent",
    "Nuclei_AreaShape_MajorAxisLength",
    "Nuclei_AreaShape_MinorAxisLength",
    "Nuclei_AreaShape_Orientation",
    "Nuclei_AreaShape_Zernike_0_0",
    "Nuclei_AreaShape_Zernike_1_1",
    "Nuclei_AreaShape_Zernike_2_0",
    "Nuclei_AreaShape_Zernike_2_2",
    "Nuclei_AreaShape_Zernike_3_1",
    "Nuclei_AreaShape_Zernike_3_3",
    "Nuclei_AreaShape_Zernike_4_0",
    "Nuclei_AreaShape_Zernike_4_2",
    "Nuclei_AreaShape_Zernike_4_4",
    "Nuclei_AreaShape_Zernike_5_1",
    "Nuclei_AreaShape_Zernike_5_3",
    "Nuclei_AreaShape_Zernike_5_5",
    "Nuclei_AreaShape_Zernike_6_0",
    "Nuclei_AreaShape_Zernike_6_2",
    "Nuclei_AreaShape_Zernike_6_4",
    "Nuclei_AreaShape_Zernike_6_6",
    "Nuclei_AreaShape_Zernike_7_1",
    "Nuclei_AreaShape_Zernike_7_3",
    "Nuclei_AreaShape_Zernike_7_5",
    "Nuclei_AreaShape_Zernike_7_7",
    "Nuclei_AreaShape_Zernike_8_0",
    "Nuclei_AreaShape_Zernike_8_2",
    "Nuclei_AreaShape_Zernike_8_4",
    "Nuclei_AreaShape_Zernike_8_6",
    "Nuclei_AreaShape_Zernike_8_8",
    "Nuclei_AreaShape_Zernike_9_1",
    "Nuclei_AreaShape_Zernike_9_3",
    "Nuclei_AreaShape_Zernike_9_5",
    "Nuclei_AreaShape_Zernike_9_7",
    "Nuclei_AreaShape_Zernike_9_9",
    "Nuclei_Neighbors_AngleBetweenNeighbors_Adjacent",
    "Nuclei_Neighbors_FirstClosestDistance_Adjacent",
    "Nuclei_Neighbors_SecondClosestDistance_Adjacent",
    "oneb_Treatment_Dose_Inhibitor_Dose",
    "Activin A [NSU]",
    "AITRL (GITR Ligand) [NSU]",
    "Amphiregulin [NSU]",
    "Amyloid beta [NSU]",
    "APRIL [NSU]",
    "BAFF [NSU]",
    "BCMA (TNFRSF17) [NSU]",
    "BDNF [NSU]",
    "BMP2 [NSU]",
    "BMP3 [NSU]",
    "BMP4 [NSU]",
    "BMP6 [NSU]",
    "BMP7 [NSU]",
    "BMP9 [NSU]",
    "C5_C5a [NSU]",
    "Calbindin [NSU]",
    "CCL1 [NSU]",
    "CCL11 [NSU]",
    "CCL13 [NSU]",
    "CCL15 [NSU]",
    "CCL16 [NSU]",
    "CCL17 [NSU]",
    "CCL18 [NSU]",
    "CCL19 [NSU]",
    "CCL2 [NSU]",
    "CCL20 [NSU]",
    "CCL21 [NSU]",
    "CCL22 [NSU]",
    "CCL23 [NSU]",
    "CCL24 [NSU]",
    "CCL25 [NSU]",
    "CCL27 [NSU]",
    "CCL28 [NSU]",
    "CCL3 [NSU]",
    "CCL4 [NSU]",
    "CCL5 [NSU]",
    "CCL7 [NSU]",
    "CCL8 [NSU]",
    "CD14 [NSU]",
    "CD163 [NSU]",
    "CD276 (B7-H3) [NSU]",
    "CD27L [NSU]",
    "CD30 [NSU]",
    "CD40L [NSU]",
    "CNTF [NSU]",
    "CRP [NSU]",
    "CX3CL1 [NSU]",
    "CXCL1 [NSU]",
    "CXCL10 [NSU]",
    "CXCL11 [NSU]",
    "CXCL12 (alpha) [NSU]",
    "CXCL12 (beta) [NSU]",
    "CXCL13 [NSU]",
    "CXCL14 [NSU]",
    "CXCL16 [NSU]",
    "CXCL17 [NSU]",
    "CXCL3 [NSU]",
    "CXCL4 [NSU]",
    "CXCL5 [NSU]",
    "CXCL6 [NSU]",
    "CXCL7 [NSU]",
    "CXCL9 [NSU]",
    "Cytochrome C [NSU]",
    "EGF [NSU]",
    "EGFR [NSU]",
    "EMMPRIN [NSU]",
    "FAS-L [NSU]",
    "FGF-1 [NSU]",
    "FGF-19 [NSU]",
    "FGF-2 [NSU]",
    "FGF-21 [NSU]",
    "FGF-4 [NSU]",
    "FGF-6 [NSU]",
    "FGF-7 (KGF) [NSU]",
    "FGF-9 [NSU]",
    "FGFR3 (IIIc) [NSU]",
    "FLRG (FSTL3) [NSU]",
    "Flt-3 Ligand [NSU]",
    "G-CSF [NSU]",
    "GDF-11 (BMP-11) [NSU]",
    "GDF-15 (MIC-1) [NSU]",
    "GDNF [NSU]",
    "GM-CSF [NSU]",
    "Granzyme B [NSU]",
    "Growth Hormone (Somatotropin) [NSU]",
    "HGF [NSU]",
    "HVEM [NSU]",
    "ICAM-1 [NSU]",
    "ICAM-2 [NSU]",
    "IFN alpha 2 (alpha 2b) [NSU]",
    "IFN beta [NSU]",
    "IFN gamma [NSU]",
    "IFN-epsilon [NSU]",
    "IGF-1 [NSU]",
    "IL-1 alpha [NSU]",
    "IL-1 beta [NSU]",
    "IL-1 R1 [NSU]",
    "IL-1 RA_RN [NSU]",
    "IL-10 [NSU]",
    "IL-11 [NSU]",
    "IL-12 p35 [NSU]",
    "IL-12 p40 [NSU]",
    "IL-12 p70 [NSU]",
    "IL-15 [NSU]",
    "IL-15_IL-15R alpha complex [NSU]",
    "IL-16 [NSU]",
    "IL-17A [NSU]",
    "IL-17B [NSU]",
    "IL-17C [NSU]",
    "IL-17D [NSU]",
    "IL-17E (IL-25) [NSU]",
    "IL-17F [NSU]",
    "IL-18 [NSU]",
    "IL-2 [NSU]",
    "IL-2 RA [NSU]",
    "IL-21 [NSU]",
    "IL-22 [NSU]",
    "IL-22 BP [NSU]",
    "IL-23 [NSU]",
    "IL-24 [NSU]",
    "IL-27 [NSU]",
    "IL-28A [NSU]",
    "IL-29 [NSU]",
    "IL-3 [NSU]",
    "IL-31 [NSU]",
    "IL-32 (alpha) [NSU]",
    "IL-33 [NSU]",
    "IL-35 [NSU]",
    "IL-4 [NSU]",
    "IL-5 [NSU]",
    "IL-6 [NSU]",
    "IL-6 R alpha [NSU]",
    "IL-7 [NSU]",
    "IL-8 [NSU]",
    "IL-9 [NSU]",
    "Leptin [NSU]",
    "LIF [NSU]",
    "LOX1 (OLR1) [NSU]",
    "M-CSF [NSU]",
    "M-CSF R (CD115) [NSU]",
    "Mesothelin [NSU]",
    "MIF [NSU]",
    "MMP-1 [NSU]",
    "MMP-10 [NSU]",
    "MMP-12 [NSU]",
    "MMP-2 [NSU]",
    "MMP-3 [NSU]",
    "MMP-7 [NSU]",
    "MMP-9 [NSU]",
    "NF-L [NSU]",
    "NGF beta [NSU]",
    "NRG1 beta 1 [NSU]",
    "Oncostatin M (OSM) [NSU]",
    "Osteopontin (OPN) [NSU]",
    "PCSK9 [NSU]",
    "PDGF-BB [NSU]",
    "PLGF [NSU]",
    "PTX3 (Pentraxin 3) [NSU]",
    "Resistin [NSU]",
    "SAA [NSU]",
    "SCF [NSU]",
    "ST2 (IL-33R) [NSU]",
    "TGF-beta 1 (LAP domain in precursor) [NSU]",
    "TGF-beta 1 (total) [NSU]",
    "TGF-beta 2 [NSU]",
    "TGF-beta 3 [NSU]",
    "Tie-2 [NSU]",
    "TIMP1 [NSU]",
    "Tissue Factor (TF) [NSU]",
    "TNF alpha [NSU]",
    "TNF RI [NSU]",
    "TNF RII [NSU]",
    "TNF RIII (Lymphotoxin Beta R) [NSU]",
    "TPO (Thrombopoietin) [NSU]",
    "TRAIL [NSU]",
    "TREM2 [NSU]",
    "TSLP [NSU]",
    "TWEAK [NSU]",
    "uPA [NSU]",
    "VCAM-1 [NSU]",
    "VEGF Receptor 2 (Flk-1) [NSU]",
    "VEGF-A (165) [NSU]",
    "VEGF-C [NSU]",
    "VEGF-D [NSU]",
    "VEGFR-1 [NSU]",
    "WISP-1 (CCN4) [NSU]",
    "XCL1 (Lymphotactin) [NSU]",
    "Metadata_Well",
    "oneb_Metadata_Treatment_Dose_Inhibitor_Dose",
]
print(cell_type, shuffle, cytokine)
if shuffle == "True":
    shuffle = True
elif shuffle == "False":
    shuffle = False
else:
    raise ValueError("shuffle must be True or False")
print(f"shuffle: {shuffle}")


# In[3]:


aggregation = True
nomic = True


# In[4]:


# set shuffle value
if shuffle:
    shuffle = "shuffled_baseline"
else:
    shuffle = "final"


# In[5]:


MODEL_TYPE = "regression"


# In[6]:


# load training data from indexes and features dataframe
data_split_path = pathlib.Path(
    f"../../0.split_data/indexes/{cell_type}/regression/aggregated_sc_and_nomic_data_split_indexes.tsv"
)
feature_combinations_path = pathlib.Path(
    f"../../0.split_data/results/feature_combinations_{cell_type}.toml"
).resolve()
data_path = pathlib.Path(
    f"../../../data/{cell_type}_preprocessed_sc_norm_aggregated_nomic.parquet"
)

# dataframe with only the labeled data we want (exclude certain phenotypic classes)
data_df = pd.read_parquet(data_path)
data_df = data_df[columns]

data_split_indexes = pd.read_csv(data_split_path, sep="\t")


# In[7]:


# select tht indexes for the training and test set
train_indexes = data_split_indexes.loc[data_split_indexes["label"] == "train"]
# subset data_df by indexes in data_split_indexes
training_data = data_df.loc[train_indexes["labeled_data_index"]]
# define metadata columns
# subset each column that contains metadata
metadata = training_data.filter(regex="Metadata")
# drop all metadata columns
data_x = training_data.drop(metadata.columns, axis=1)
labeled_data = training_data["oneb_Metadata_Treatment_Dose_Inhibitor_Dose"]
# get all columns that contain "NSU" in the column name
data_y_cols = data_x.filter(regex="NSU").columns
train_y = training_data[data_y_cols]
train_x = data_x.drop(data_y_cols, axis=1)
train_x = train_x.drop(columns="oneb_Treatment_Dose_Inhibitor_Dose")
loo = LeaveOneOut()
loo.get_n_splits(train_y)

train_data_y = train_y[cytokine]
model = ElasticNetCV(
    random_state=0,
    max_iter=10000,
    cv=loo,
    l1_ratio=[0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.99],
    alphas=[0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],
    fit_intercept=True,
    selection="random",
)
# train model on training data on all combinations of model types, feature types, and phenotypic classes


# In[8]:


if shuffle == "shuffled_baseline":
    print("Shuffling data")
    for column in train_x:
        np.random.shuffle(train_x[column].values)
else:
    print("Not shuffling data")
# define parameters to search over
with parallel_backend("multiprocessing"):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning, module="sklearn")
        # create a logistic regression model
        model.fit(train_x, train_data_y)
        scores = cross_val_score(
            model,
            train_x,
            train_data_y,
            scoring="neg_mean_absolute_error",
            cv=loo,
            n_jobs=-1,
        )
        print(scores)
        print(f"Mean MAE: {scores.mean()}")
        print(f"Std MAE: {scores.std()}")
        print(f"R2: {model.score(train_x, train_data_y)}")

if (aggregation == True) and (nomic == True):
    results_dir = f"../models/regression/{cell_type}/aggregated_with_nomic/"
elif (aggregation == True) and (nomic == False):
    results_dir = f"../models/regression/{cell_type}/aggregated/"
elif (aggregation == False) and (nomic == True):
    results_dir = f"../models/regression/{cell_type}/sc_with_nomic/"
elif (aggregation == False) and (nomic == False):
    results_dir = f"../models/regression/{cell_type}/sc/"
else:
    print("Error")

# create results directory if it doesn't exist
pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

# save final estimator
if shuffle == "shuffled_baseline":
    dump(
        model,
        f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.joblib",
    )
elif shuffle == "final":
    dump(
        model,
        f"{results_dir}/{cytokine}_final__all_nomic.joblib",
    )
else:
    print("Error")

# save condfig copy specific to this model to the folder with the results
# use pathlib
if shuffle == "shuffled_baseline":
    config_copy_path = pathlib.Path(
        f"{results_dir}/{cytokine}_shuffled_baseline__all_nomic.toml"
    )
elif shuffle == "final":
    config_copy_path = pathlib.Path(f"{results_dir}/{cytokine}_final__all_nomic.toml")
else:
    print("Error")

# write toml file with parameters used from injected parameters

with open(config_copy_path, "w") as f:
    f.write(f"model_type='{shuffle}'\n")
    f.write(f"aggregation={aggregation}\n")
    f.write(f"nomic={nomic}\n")
    f.write(f"cell_type='{cell_type}'\n")
    f.write(f"feature=all\n")
