# Exploratory data analysis

This folder contains notebooks that demonstrate how to perform exploratory data analysis (EDA) on a dataset.
EDA is an essential step in the data analysis process.
It allows you to summarize the main characteristics of the data, often with visual methods. This helps you understand the data and the relationships between the variables.

Due to the exploratory nature of the analysis, some analysis are secondary and were not included in the manuscript that this repository is for.
However, they are still useful for understanding the data and the relationships between the variables.

The main analysis that are needed for the manuscript are:

1. [UMAP](1.umap_analysis_plate2.ipynb)
    This notebook extracts the UMAP embeddings of the plate 2 data and visualizes the embeddings.
2. [Cell Counts](4.cell_count_analysis.ipynb)
    This notebook extracts the cell counts of the plate 2 data and visualizes the cell counts.
3. [Feature index](8.0_create_feature_index.ipynb)
    This notebook creates a feature index for the plate 2 data.
4. [ANOVA](8.1_anova_all_groupings.ipynb)
    This notebook performs ANOVA on the plate 2 data.
5. [Combine ANOVA](8.2_combine_files.ipynb)
    This notebook combines the ANOVA results from the previous notebook.
6. [Subset ANOVA](9.subset_umap.ipynb)
    This notebook subsets the UMAP embeddings based on the ANOVA results.
