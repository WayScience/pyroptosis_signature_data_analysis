#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# Merge first and second run of the SH-SY5Y cell line Cell profiler analysis
# _sc_ = single cell resolved features
# _fs_ = feature selection performed on the singled cell resolved feature space
#
# Each run need to be merged to return one continuous run

# In[2]:


# define input paths
first_run_sc_fs_path = Path(
    "/media/lippincm/c58d4f19-ae4d-4b78-8370-2c2639886da0/Interstellar_Plate2/SHSY5Y_first_run_sc_norm_fs.parquet"
)
first_run_sc_path = Path(
    "/media/lippincm/c58d4f19-ae4d-4b78-8370-2c2639886da0/Interstellar_Plate2/SHSY5Y_first_run_sc_norm.parquet"
)

second_run_sc_fs_path = Path(
    "/media/lippincm/c58d4f19-ae4d-4b78-8370-2c2639886da0/Interstellar_Plate2/SHSY5Y_second_run_sc_norm_fs.parquet"
)
second_run_sc_path = Path(
    "/media/lippincm/c58d4f19-ae4d-4b78-8370-2c2639886da0/Interstellar_Plate2/SHSY5Y_second_run_sc_norm.parquet"
)

# define output paths
SHSY5Y_sc_norm_fs_path = Path(
    "../../Extracted_Features_(CSV_files)/SHSY5Y_run_sc_norm_fs.parquet"
)
SHSY5Y_sc_norm_path = Path(
    "../../Extracted_Features_(CSV_files)/SHSY5Y_run_sc_norm.parquet"
)


# In[3]:


# read parquet files into pandas dataframes
first_run_sc_fs = pq.read_table(first_run_sc_fs_path).to_pandas()
first_run_sc = pq.read_table(first_run_sc_path).to_pandas()
second_run_sc_fs = pq.read_table(second_run_sc_fs_path).to_pandas()
second_run_sc = pq.read_table(second_run_sc_fs_path).to_pandas()


# In[ ]:


# concatenate dataframes (append second run to first run)
SHSY5Y_sc_norm_fs = pd.concat([first_run_sc_fs, second_run_sc_fs], ignore_index=True)
SHSY5Y_sc_norm = pd.concat([first_run_sc, second_run_sc], ignore_index=True)


# In[ ]:


# convert dataframes to pyarrow tables
SHSY5Y_sc_norm_fs = pa.Table.from_pandas(SHSY5Y_sc_norm_fs)
SHSY5Y_sc_norm = pa.Table.from_pandas(SHSY5Y_sc_norm)


# In[ ]:


# write pyarrow tables to parquet files
pq.write_table(SHSY5Y_sc_norm_fs, SHSY5Y_sc_norm_fs_path)
pq.write_table(SHSY5Y_sc_norm, SHSY5Y_sc_norm_path)
