suppressWarnings(suppressPackageStartupMessages(library(ggplot2))) # plotting
suppressWarnings(suppressPackageStartupMessages(library(platetools))) # plotting
suppressWarnings(suppressPackageStartupMessages(library(gridExtra))) # plot assembly
suppressWarnings(suppressPackageStartupMessages(library(cowplot))) # plot assembly
suppressWarnings(suppressPackageStartupMessages(library(viridis))) # color palettes
suppressWarnings(suppressPackageStartupMessages(library(patchwork))) # plot assembly
suppressWarnings(suppressPackageStartupMessages(library(arrow))) # reading parquet files
suppressWarnings(suppressPackageStartupMessages(library(dplyr))) # data manipulation

# define the base directory
base_dir <- file.path(
    "..",
    "..",
    ".."
)

# Set the path to the index file
index_file_path_PBMC <- file.path(
    base_dir,
    "4.sc_Morphology_Neural_Network_MLP_Model/0.hyperparameter_optimization/indexes/PBMC/multi_class/MultiClass_MLP_data_split_indexes.tsv"
)
index_file_path_SHSY5Y <- file.path(
    base_dir,
    "4.sc_Morphology_Neural_Network_MLP_Model/0.hyperparameter_optimization/indexes/SHSY5Y/multi_class/MultiClass_MLP_data_split_indexes.tsv"
)
# load the index file
index_file_PBMC <- read.table(
    file = index_file_path_PBMC,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)
index_file_SHSY5Y <- read.table(
    file = index_file_path_SHSY5Y,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)


# create the figure directory if it does not exist
figure_dir <- file.path(
    "../",
    "figures"
)
if (!dir.exists(figure_dir)) {
    dir.create(figure_dir)
}

# sort the index file by labeled_data_index
index_file_PBMC <- index_file_PBMC[order(index_file_PBMC$labeled_data_index),]
head(index_file_PBMC,2)
index_file_SHSY5Y <- index_file_SHSY5Y[order(index_file_SHSY5Y$labeled_data_index),]
head(index_file_SHSY5Y,2)


# set path to the PBMC metadata file
metadata_file_path_PBMC <- file.path(
    "..",
    "..",
    "..",
    "data",
    "PBMC_preprocessed_sc_norm.parquet"
)

# load via arrow with certain columns
metadata_file_PBMC <- arrow::read_parquet(
    metadata_file_path_PBMC,
    col_select = c(
        'Metadata_cell_type',
        'Metadata_Well'	,
        'Metadata_number_of_singlecells',
    )
)
head(metadata_file_PBMC)

# SHSY5Y
# set path to the SHSY5Y metadata file
metadata_file_path_SHSY5Y <- file.path(
    "..",
    "..",
    "..",
    "data",
    "SHSY5Y_preprocessed_sc_norm.parquet"
)
# load the metadata file via arrow
metadata_file_SHSY5Y <- arrow::read_parquet(
    metadata_file_path_SHSY5Y,
    col_select = c(
        'Metadata_cell_type',
        'Metadata_Well'	,
        'Metadata_number_of_singlecells',
    )
)
head(metadata_file_SHSY5Y)

## PBMC

# add column to metadata file with the labeled_data_index
metadata_file_PBMC$label <- index_file_PBMC$label
# replace test, train, or validation with "pool" in label column
metadata_file_PBMC$label <- gsub(
    pattern = "test",
    replacement = "pool",
    x = metadata_file_PBMC$label
)
metadata_file_PBMC$label <- gsub(
    pattern = "train",
    replacement = "pool",
    x = metadata_file_PBMC$label
)
metadata_file_PBMC$label <- gsub(
    pattern = "val",
    replacement = "pool",
    x = metadata_file_PBMC$label
)
# aggregate the metadata via Metadata_Well
metadata_file_PBMC <- metadata_file_PBMC %>%
    group_by(Metadata_Well, Metadata_cell_type, label) %>%
    summarize(
        n = n()
    )
head(metadata_file_PBMC)

# SHSY5Y
# add column to metadata file with the labeled_data_index
metadata_file_SHSY5Y$label <- index_file_SHSY5Y$label
# replace test, train, or validation with "pool" in label column
metadata_file_SHSY5Y$label <- gsub(
    pattern = "test",
    replacement = "pool",
    x = metadata_file_SHSY5Y$label
)
metadata_file_SHSY5Y$label <- gsub(
    pattern = "train",
    replacement = "pool",
    x = metadata_file_SHSY5Y$label
)
metadata_file_SHSY5Y$label <- gsub(
    pattern = "val",
    replacement = "pool",
    x = metadata_file_SHSY5Y$label
)
# aggregate the metadata via Metadata_Well
metadata_file_SHSY5Y <- metadata_file_SHSY5Y %>%
    group_by(Metadata_Well, Metadata_cell_type, label) %>%
    summarize(
        n = n()
    )
head(metadata_file_SHSY5Y)


# combine the SHSY5Y and PBMC metadata files
combined_metadata <- rbind(
    metadata_file_PBMC,
    metadata_file_SHSY5Y
)
head(combined_metadata)
unique(combined_metadata$Metadata_cell_type)
unique(combined_metadata$label)

# load in the platemap
platemap_file_path <- file.path(
    base_dir,
    "data/",
    "Interstellar_plate2_platemap.csv"
)
# read in the platemap
platemap_file <- read.csv(
    file = platemap_file_path,
    header = TRUE,
    stringsAsFactors = FALSE
)
head(platemap_file)

# merge the combined metadata file with plate map on Metadata_Well
updated_platemap <- merge(
    x = platemap_file,
    y = combined_metadata,
    by.x = "well_id",
    by.y = "Metadata_Well",
    all.x = TRUE
)
head(updated_platemap)
unique(updated_platemap$Metadata_cell_type)
# replace "" in cell_type with "Blank" in label
updated_platemap$Metadata_cell_type[is.na(updated_platemap$Metadata_cell_type)] <- "Blank"
# replace NA with "Blank" in label
updated_platemap$label[is.na(updated_platemap$label)] <- "Blank"
# replace pool with Pool in label
updated_platemap$label[updated_platemap$label == "pool"] <- "Pool"
# replace holdout with Holdout in label
updated_platemap$label[updated_platemap$label == "holdout"] <- "Holdout well"
# replace treatment_holdout with Treatment Holdout in label
updated_platemap$label[updated_platemap$label == "treatment_holdout"] <- "Treatment holdout"
unique(updated_platemap$Metadata_cell_type)
unique(updated_platemap$label)

# Make the cell type factor ordered
updated_platemap$Metadata_cell_type <- factor(
    x = updated_platemap$Metadata_cell_type,
    levels = c(
        "Blank",
        "PBMC",
        "SH-SY5Y"
    )
)
# make the label factor ordered
updated_platemap$label <- factor(
    x = updated_platemap$label,
    levels = c(
        'Blank',
        'Pool',
        'Holdout well',
        'Treatment holdout'
    )
)

width <- 14
height <- 14
options(repr.plot.width = width, repr.plot.height = height, units = "cm")
# set pallette
viridis_pal_custom <- viridis::viridis_pal(option = "C")(5)

data_split_plate_map_full <- (
    raw_map(
        data = updated_platemap$label,
        well = updated_platemap$well_id,
        plate = 384, # number of wells in plate apriori known
        size = 14 # size of the wells displayed
        )
    + theme_dark()
        + ggplot2::geom_point(
        aes(shape = updated_platemap$Metadata_cell_type),
        size = 5
        )
        + labs(fill = "Data Split", shape = "Cell Type")
    # change legend text size for fill

    + guides(shape = guide_legend(override.aes = list(size = 12), nrow = 1))
    + guides(fill = guide_legend(override.aes = list(size = 12),ncol = 2))
    + theme(
        legend.title = element_text(size = 18,hjust = 0.5),
        legend.text = element_text(size = 16),
    )
    # cell type legend
    + scale_shape_manual(
        values = c(
            'Blank' = 0,
            'PBMC' = 19,
            'SH-SY5Y' = 8
        )
    )
    # data split legend
    + scale_fill_manual(
    values = c(
        'Blank' = "grey",
        'Pool' = viridis_pal_custom[3],
        'Holdout well' = viridis_pal_custom[4],
        'Treatment holdout' = viridis_pal_custom[5]

    )
    )
    + theme(
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)
    )

)
data_split_plate_map_full
ggsave(
    filename = "../figures/data_split_plate_map_full.png",
    plot = data_split_plate_map_full,
    width = width,
    height = height,
    units = "in",
    dpi = 600
)


# remove SHSY5Y from the platemap
updated_platemap <- updated_platemap[updated_platemap$Metadata_cell_type != "SH-SY5Y",]
head(updated_platemap)

width <- 14
height <- 14
options(repr.plot.width = width, repr.plot.height = height, units = "cm")
# set pallette
viridis_pal_custom <- viridis::viridis_pal(option = "C")(5)

data_split_plate_map <- (
    raw_map(
        data = updated_platemap$label,
        well = updated_platemap$well_id,
        plate = 384, # number of wells in plate apriori known
        size = 14 # size of the wells displayed
        )
    + theme_dark()
        + ggplot2::geom_point(
        aes(shape = updated_platemap$Metadata_cell_type),
        size = 5
        )
        + labs(fill = "Data Split", shape = "Cell Type")
    # change legend text size for fill

    + guides(shape = guide_legend(override.aes = list(size = 12), nrow = 1))
    + guides(fill = guide_legend(override.aes = list(size = 12),ncol = 2))
    + theme(
        legend.title = element_text(size = 18,hjust = 0.5),
        legend.text = element_text(size = 16),
    )
    # cell type legend
    + scale_shape_manual(
        values = c(
            'Blank' = 0,
            'PBMC' = 19,
            'SH-SY5Y' = 8
        )
    )
    # data split legend
    + scale_fill_manual(
    values = c(
        'Blank' = "grey",
        'Pool' = viridis_pal_custom[3],
        'Holdout well' = viridis_pal_custom[4],
        'Treatment holdout' = viridis_pal_custom[5]

    )
    )
    + theme(
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)
    )

)
data_split_plate_map
ggsave(
    filename = "../figures/data_split_plate_map.png",
    plot = data_split_plate_map,
    width = width,
    height = height,
    units = "in",
    dpi = 600
)

