list_of_packages <- c(
    "ggplot2", # for plotting
    "dplyr", # for data manipulation
    "patchwork", # for combining plots
    "cowplot", # for combining plots
    "RcppTOML", # for reading TOML files
    "pheatmap", # for heatmaps
    "lattice", # for heatmaps
    "RColorBrewer", # for heatmaps
    "RColorBrewer", # for heatmaps
    "ComplexHeatmap", # for heatmaps
    "ggplotify", # for heatmaps
    "viridis", # for color palettes
    "platetools", # for plate visualization
    "circlize", # for plate visualization
    "reshape2", # for data manipulation
    "stringr", # for string manipulation
    "purrr" # for data manipulation
)
for (package in list_of_packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            suppressMessages(
                library(package, character.only = TRUE)
            )
        )
    )
}

#  import the theme
source("../../utils/figure_themes.r")

# set the cell type
cell_type <- "PBMC"


# set the path to the data files
df_stats_path <- file.path(
    paste0("../../../6.bulk_Morphology_Elastic_Network/2.test_model/results/regression/",cell_type,"/aggregated_with_nomic/model_stats.csv"
    )
    )
df_variance_path <- file.path(
    paste0("../../../6.bulk_Morphology_Elastic_Network/2.test_model/results/regression/",cell_type,"/aggregated_with_nomic/variance_r2_stats.csv"
    )
)

# set the path to the figure output
enet_cp_fig_path <- paste0("../figures/regression/",cell_type,"/")
# if path does not exist, create it
if (!file.exists(dirname(enet_cp_fig_path))) {
    dir.create(dirname(enet_cp_fig_path), recursive = TRUE)
}

# read the data
df_stats <- read.csv(df_stats_path)
df_variance <- read.csv(df_variance_path)


head(df_stats)
head(df_variance)
# remove '[]' from the string in the column
df_variance$r2 <- gsub("\\[|\\]", "", df_variance$r2)
# set the column as numeric
df_variance$r2 <- as.numeric(df_variance$r2)
head(df_variance)


df_variance$shuffle <- gsub("final", "Final\n ", df_variance$shuffle)
df_variance$shuffle <- gsub("shuffled_baseline", "Shuffled\nbaseline", df_variance$shuffle)
df_variance$data_split <- gsub("test_data", "Test Data", df_variance$data_split)
df_variance$data_split <- gsub("train_data", "Train Data", df_variance$data_split)

# set plot size
options(repr.plot.width=7, repr.plot.height=5)
# set output path
global_variance_r2_path <- file.path(paste0(enet_cp_fig_path,"global_variance_r2.png"))
# if path does not exist, create it
if (!file.exists(dirname(global_variance_r2_path))) {
    print(dirname(global_variance_r2_path))
    dir.create(dirname(global_variance_r2_path), recursive = TRUE)
}
# plot df_var df
variance_r2_plot_global <- (
    ggplot(df_variance, aes(x=r2, y=actual_value,col=shuffle, shape = data_split))
    + geom_point()
    + theme_bw()
    + labs(x="R2", y="Explained variance")
    # update the legend title
    + labs(shape = "Data Split", col = "Model Shuffle")
    # alter the text size of the legend title
    + figure_theme
    + scale_shape_manual(values=c(16, 4))
    # change the color of the points
    + scale_color_manual(values=c("Darkblue", "orange"))
    # legend position
    # + theme(legend.position = "bottom", legend.box = "horizontal")
    # change legend dot size
    + guides(colour = guide_legend(override.aes = list(size=3)))
)
variance_r2_plot_global
ggsave(global_variance_r2_path, variance_r2_plot_global, width=5, height=5, dpi=500)



df_stats$shuffle_plus_data_split <- paste0(df_stats$shuffle, "_", df_stats$data_split)
# replace 'final_test_data' with 'Final + Test' and 'final_train_data' with 'Final + Train'
df_stats$shuffle_plus_data_split <- gsub("final_test_data", "Final (Test)", df_stats$shuffle_plus_data_split)
df_stats$shuffle_plus_data_split <- gsub("final_train_data", "Final (Train)", df_stats$shuffle_plus_data_split)
df_stats$shuffle_plus_data_split <- gsub("shuffled_baseline_test_data", "Shuffled (Test)", df_stats$shuffle_plus_data_split)
df_stats$shuffle_plus_data_split <- gsub("shuffled_baseline_train_data", "Shuffled (Train)", df_stats$shuffle_plus_data_split)


options(repr.plot.width=6, repr.plot.height=5)
# set output path
global_prediction_trend_path <- file.path(paste0(enet_cp_fig_path,"global_prediction_trend.png"))
# if path does not exist, create it
if (!file.exists(dirname(global_prediction_trend_path))) {
    print(dirname(global_prediction_trend_path))
    dir.create(dirname(global_prediction_trend_path), recursive = TRUE)
}
# plot the data
global_prediction_trend_scatter <- (
    ggplot(df_stats, aes(x=actual_value, y=predicted_value, col=shuffle_plus_data_split))
    + geom_point(alpha=0.5, size=0.5)
    + geom_smooth(method="lm", se=TRUE, alpha=0.5, size=0.5, color="black")

    # add geom smooth with each line being a different color
    + labs(x="Actual", y="Predicted")
    + theme_bw()
    + labs(title="Global Prediction Trends of Cytokine Concentrations")
    # add y=x line
    + geom_abline(intercept = 0, slope = 1, linetype="dashed", color="black")
    + facet_wrap(.~shuffle_plus_data_split, ncol=2)
    + labs(color="Model", hjust=0.5)
    + figure_theme
    + scale_x_continuous(breaks = seq(0, 1, 0.5))
    # legend dot size
    + guides(colour = guide_legend(override.aes = list(size=5), ncol = 2))
    # legend position
    + theme(legend.position = "bottom", legend.box = "horizontal")
    # rotate x axis text
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

# save the plot
ggsave(global_prediction_trend_path, global_prediction_trend_scatter, width=5, height=5, dpi=500)
global_prediction_trend_scatter


#  aggregate the data by the cytokine
df_stats_agg <- df_stats %>%
    group_by(cytokine) %>%
    summarise(
        mean_r2 = mean(r2),
        mean_neg_mean_squared_error = mean(neg_mean_squared_error),
        mean_log10_neg_mean_squared_error = mean(log10_neg_mean_squared_error),
        mean_actual_value = mean(actual_value),
        mean_predicted_value = mean(predicted_value)
    )
# sort the data by the log10_neg_mean_squared_error descending
df_stats_agg <- df_stats_agg[order(df_stats_agg$mean_log10_neg_mean_squared_error),]
poor_performance_cytokines <- df_stats_agg$cytokine[1:9]

# df_stats factor levels
df_stats$shuffle_plus_data_split <- factor(
    df_stats$shuffle_plus_data_split,
    levels = c(
        "Final (Train)",
        "Final (Test)",
        "Shuffled (Train)",
        "Shuffled (Test)"
    )
)

enet_cp_fig <- file.path(paste0(enet_cp_fig_path,"Predicted_vs_Actual_all_cytokines.png"))
# set plot size
width <- 10
height <- 10
options(repr.plot.width=width, repr.plot.height=height)

# subset the df_stats to only include the cytokines of interest
cytokines <- c("IL-1beta", "TFalpha", "CCL24", "IL-18", "Osteopontin(OP)", "CCL13", "IL-2", "IL-6", "CCL4")
df_stats1 <- df_stats[df_stats$cytokine %in% cytokines,]

cytokine_predictions1 <- (
    ggplot(df_stats1, aes(x=actual_value, y=predicted_value, col=shuffle_plus_data_split))
    + geom_point()
    + theme_bw()
    + geom_smooth(method=lm, se=TRUE, formula = y ~ x, alpha=0.5, size=0.5)
    + labs(x="Actual", y="Predicted ")
    + ylim(0, 1)
    + xlim(0, 1)
    + figure_theme
    + labs(color="Model", hjust=0.5)
    # change legend title
    # make kegend key background white
    + guides(color = guide_legend(override.aes = list(fill = NA, size = 5), ncol = 2),
        linetype = guide_legend(override.aes = list(fill = NA)))
    + theme(legend.key = element_rect(fill = "white"))
    + facet_wrap(.~cytokine, ncol=3)
        # x tick marks
    + scale_x_continuous(breaks = seq(0, 1, 0.5))
    # legend position
    + theme(legend.position = "bottom", legend.box = "horizontal")
    # rotate x axis text
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # make the legend background white instead of grey
    + theme(legend.key=element_blank(), legend.background = element_blank())


    )
cytokine_predictions1


df_stats2 <- df_stats[df_stats$cytokine %in% poor_performance_cytokines,]

cytokine_predictions2 <- (
    ggplot(df_stats2, aes(x=actual_value, y=predicted_value, col=shuffle_plus_data_split))
    + geom_point()
    + theme_bw()
    + geom_smooth(method=lm, se=TRUE, formula = y ~ x, alpha=0.5, size=0.5)
    + labs(x="Actual", y="Predicted ")
    + ylim(0, 1)
    + xlim(0, 1)
    + figure_theme
    + labs(color="Model", hjust=0.5)
    # change legend title
    # make kegend key background white
    + guides(color = guide_legend(override.aes = list(fill = NA, size = 5), ncol = 2),
        linetype = guide_legend(override.aes = list(fill = NA)))
    + theme(legend.key = element_rect(fill = "white"))
    + facet_wrap(.~cytokine, ncol=3)
        # x tick marks
    + scale_x_continuous(breaks = seq(0, 1, 0.5))
    # legend position
    + theme(legend.position = "bottom", legend.box = "horizontal")
    # rotate x axis text

    )
cytokine_predictions2


# path set
input_file_path <- file.path(paste0("../../../6.bulk_Morphology_Elastic_Network/3.model_coefficients/results/regression/",cell_type))
# read in the data
output_path <- file.path(paste0("../figures/","regression/",cell_type,"/"))
# create output directory if it doesn't exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


## function to process the data for visualization
process_subset_data <- function(df){
    # read in the data
    # data <- read.csv(data_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    # get the basename of the files

    data <- df %>%
        dplyr::arrange(desc(abs(coefficients))) %>%
        tidyr::separate(
            feature_names,
            into = c(
                "compartment",
                "feature_group",
                "measurement",
                "channel",
                "parameter1",
                "parameter2"
            ),
            sep = "_",
            remove = FALSE
        ) %>%
        dplyr::mutate(channel_cleaned = channel) %>%
        dplyr::arrange(desc(abs(coefficients)))

    # Clean channel for visualization
    data$channel_learned <- dplyr::recode(data$channel,
            "CorrDNA" = "nuclei",
            "CorrMito" = "Mito",
            "CorrER" = "ER",
            "CorrGasdermin" = "gasdermin",
            "CorrPM" = "PM",
            .default = "other",
            .missing="other"
    )
    data <- data %>%
        dplyr::group_by(feature_group, channel_learned, compartment) %>%
        dplyr::slice_max(order_by = coefficients, n = 1)
    return(data)
}


plot_coeffs <- function(df, cytokine, shuffle){
    # replace "[NSU]" with ""
    cytokine <- gsub("\\[NSU\\]", "", cytokine)
    # plot the data
    coef_gg <- (
        ggplot(df, aes(x = channel_learned, y = feature_group))
        + geom_point(aes(fill = abs(coefficients)), pch = 22, size = 5.75)
        + facet_wrap("~compartment", ncol = 3)
        + theme_bw()
        + scale_fill_continuous(
            name=paste0(cytokine,"\ntop Abs. val\ntreatment\nlinear model\ncoefficient"),
            low = "darkblue",
            high = "yellow",
        )
        + xlab("Channel")
        + ylab("Feature")

        + figure_theme
        + theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        )
        # rotate x axis labels
        + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        + theme(plot.title = element_text(hjust = 0.5))
        )
        return(coef_gg)
    }


# get all files in a directory
files <- list.files(path = input_file_path, pattern = "*.csv", full.names = TRUE)

# create empty list (mimics a dictionary )
nested_list <- list(
        filename = list(),
        cytokine = list(),
        shuffle = list()
    )

for (i in files){
    filename <- basename(i)
    # split the string at the first _
    filename <- strsplit(filename, "_", fixed = TRUE)[[1]]
    cytokine <- filename[1]
    shuffle <- filename[2]
    nested_list$filename <- c(nested_list$filename, i)
    nested_list$cytokine <- c(nested_list$cytokine, cytokine)
    nested_list$shuffle <- c(nested_list$shuffle, shuffle)
}


cytokine <- 'IL-1 beta [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
il1beta_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

il1beta_final <- process_subset_data(il1beta_final)
head(il1beta_final)
# factor levels
# rename the factor levels
il1beta_final$channel_learned <- factor(
    il1beta_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
il1beta_final$feature_group <- factor(
    il1beta_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
il1beta_final$feature_group <- factor(
    il1beta_final$feature_group,
    levels = rev(levels(il1beta_final$feature_group))
)

il1beta_final_plot <- plot_coeffs(il1beta_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, il1beta_final_plot, width=width, height=height, dpi=500)
il1beta_final_plot


cytokine <- 'TNF alpha [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
tnfa_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

tnfa_final <- process_subset_data(tnfa_final)
head(tnfa_final)
# factor levels
# rename the factor levels
tnfa_final$channel_learned <- factor(
    tnfa_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
tnfa_final$feature_group <- factor(
    tnfa_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
tnfa_final$feature_group <- factor(
    tnfa_final$feature_group,
    levels = rev(levels(tnfa_final$feature_group))
)

tnfa_final_plot <- plot_coeffs(tnfa_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, tnfa_final_plot, width=width, height=height, dpi=500)
tnfa_final_plot


cytokine <- 'CCL24 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
CCL24_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

CCL24_final <- process_subset_data(CCL24_final)
head(CCL24_final)
# factor levels
# rename the factor levels
CCL24_final$channel_learned <- factor(
    CCL24_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
CCL24_final$feature_group <- factor(
    CCL24_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
CCL24_final$feature_group <- factor(
    CCL24_final$feature_group,
    levels = rev(levels(CCL24_final$feature_group))
)

CCL24_final_plot <- plot_coeffs(CCL24_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, CCL24_final_plot, width=width, height=height, dpi=500)
CCL24_final_plot


cytokine <- 'IL-18 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
il18beta_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

il18beta_final <- process_subset_data(il18beta_final)
head(il18beta_final)
# factor levels
# rename the factor levels
il18beta_final$channel_learned <- factor(
    il18beta_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
il18beta_final$feature_group <- factor(
    il18beta_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
il18beta_final$feature_group <- factor(
    il18beta_final$feature_group,
    levels = rev(levels(il18beta_final$feature_group))
)

il18beta_final_plot <- plot_coeffs(il18beta_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, il18beta_final_plot, width=width, height=height, dpi=500)
il18beta_final_plot


cytokine <- 'Osteopontin (OPN) [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
op_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

op_final <- process_subset_data(op_final)
# replace (OPN with "" in the cytokine column)
# op_final$secreted_proteins <- stringr::str_replace_all(op_final$secreted_proteins, "Osteopontin (OPN) [NSU]", "Osteopontin [NSU]")
head(op_final)
# factor levels
# rename the factor levels
op_final$channel_learned <- factor(
    op_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
op_final$feature_group <- factor(
    op_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
op_final$feature_group <- factor(
    op_final$feature_group,
    levels = rev(levels(op_final$feature_group))
)

op_final_plot <- plot_coeffs(op_final, cytokine, shuffle)
op_final_plot <- (
    op_final_plot
    + scale_fill_continuous(
            name=paste0("Osteopontin","\ntop Abs. val\ntreatment\nlinear model\ncoefficient"),
            low = "darkblue",
            high = "yellow",
        )
)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, op_final_plot, width=width, height=height, dpi=500)
op_final_plot


cytokine <- 'CCL13 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
CCL13_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

CCL13_final <- process_subset_data(CCL13_final)
head(CCL13_final)
# factor levels
# rename the factor levels
CCL13_final$channel_learned <- factor(
    CCL13_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
CCL13_final$feature_group <- factor(
    CCL13_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
CCL13_final$feature_group <- factor(
    CCL13_final$feature_group,
    levels = rev(levels(CCL13_final$feature_group))
)

CCL13_final_plot <- plot_coeffs(CCL13_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, CCL13_final_plot, width=width, height=height, dpi=500)
CCL13_final_plot


cytokine <- 'IL-2 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
il2_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

il2_final <- process_subset_data(il2_final)
head(il2_final)
# factor levels
# rename the factor levels
il2_final$channel_learned <- factor(
    il2_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
il2_final$feature_group <- factor(
    il2_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
il2_final$feature_group <- factor(
    il2_final$feature_group,
    levels = rev(levels(il2_final$feature_group))
)

il2_final_plot <- plot_coeffs(il2_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, il2_final_plot, width=width, height=height, dpi=500)
il2_final_plot


cytokine <- 'IL-6 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
il6_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

il6_final <- process_subset_data(il6_final)
head(il6_final)
# factor levels
# rename the factor levels
il6_final$channel_learned <- factor(
    il6_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
il6_final$feature_group <- factor(
    il6_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
il6_final$feature_group <- factor(
    il6_final$feature_group,
    levels = rev(levels(il6_final$feature_group))
)

il6_final_plot <- plot_coeffs(il6_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, il6_final_plot, width=width, height=height, dpi=500)
il6_final_plot


cytokine <- 'CCL4 [NSU]'
shuffle <- 'final'
filename <- nested_list$filename[which(nested_list$cytokine == cytokine & nested_list$shuffle == shuffle)]
# set to string
filename <- filename[[1]]
CCL4_final <- read.csv(filename, header = TRUE, sep = ",", stringsAsFactors = FALSE)

CCL4_final <- process_subset_data(CCL4_final)
head(CCL4_final)
# factor levels
# rename the factor levels
CCL4_final$channel_learned <- factor(
    CCL4_final$channel_learned,
    levels = c(
        "nuclei",
        "Mito",
        "ER",
        "gasdermin",
        "PM",
        "other"
    ),
    labels = c(
        "Nuclei",
        "Mito",
        "ER",
        "GasderminD",
        "AGP",
        "Other"
    )
)
CCL4_final$feature_group <- factor(
    CCL4_final$feature_group,
    levels = c(
        "AreaShape",
        "Correlation",
        "Granularity",
        "Intensity",
        "Location",
        "Neighbors",
        "RadialDistribution",
        "Texture"
    )
)
# reverese the order of the factor levels
CCL4_final$feature_group <- factor(
    CCL4_final$feature_group,
    levels = rev(levels(CCL4_final$feature_group))
)

CCL4_final_plot <- plot_coeffs(CCL4_final, cytokine, shuffle)
# output
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
# set plot size
width <- 8.5
height <- 4
options(repr.plot.width=width, repr.plot.height=height)
ggsave(coef_gg_file, CCL4_final_plot, width=width, height=height, dpi=500)
CCL4_final_plot


# set cell type
cell_type <- "PBMC"


# set path for data of all models
data_path <- file.path(paste0("../../../6.bulk_Morphology_Elastic_Network/4.model_performance/results/regression/", cell_type, "/", "all_model_performance.csv"))
df <- read.csv(data_path)
# setfigure path
figure_path <- file.path(paste0("../figures/regression/", cell_type, "/"))
# make the directory if it doesn't exist
dir.create(figure_path, recursive = TRUE, showWarnings = FALSE)


# fix the col name
df <- df %>%
  mutate(secreted_proteins = case_when(
    secreted_proteins == "MMP-1 [NSU]" ~ "MMP-1",
    secreted_proteins == "VEGFR-1 [NSU]" ~ "VEGFR-1",
    secreted_proteins == "CCL4 [NSU]" ~ "CCL4",
    secreted_proteins == "MMP-12 [NSU]" ~ "MMP-12",
    secreted_proteins == "CCL18 [NSU]" ~ "CCL18",
    secreted_proteins == "IL-9 [NSU]" ~ "IL-9",
    secreted_proteins == "TWEAK [NSU]" ~ "TWEAK",
    secreted_proteins == "EGFR [NSU]" ~ "EGFR",
    secreted_proteins == "IL-21 [NSU]" ~ "IL-21",
    secreted_proteins == "FGF-1 [NSU]" ~ "FGF-1",
    secreted_proteins == "FAS-L [NSU]" ~ "FAS-L",
    secreted_proteins == "CXCL12 (beta) [NSU]" ~ "CXCL12 (beta)",
    secreted_proteins == "CXCL12 (alpha) [NSU]" ~ "CXCL12 (alpha)",
    secreted_proteins == "CXCL14 [NSU]" ~ "CXCL14",
    secreted_proteins == "HGF [NSU]" ~ "HGF",
    secreted_proteins == "IL-3 [NSU]" ~ "IL-3",
    secreted_proteins == "CXCL7 [NSU]" ~ "CXCL7",
    secreted_proteins == "CCL25 [NSU]" ~ "CCL25",
    secreted_proteins == "BMP9 [NSU]" ~ "BMP9",
    secreted_proteins == "IL-12 p35 [NSU]" ~ "IL-12 p35",
    secreted_proteins == "CCL16 [NSU]" ~ "CCL16",
    secreted_proteins == "CCL2 [NSU]" ~ "CCL2",
    secreted_proteins == "LIF [NSU]" ~ "LIF",
    secreted_proteins == "CXCL9 [NSU]" ~ "CXCL9",
    secreted_proteins == "CNTF [NSU]" ~ "CNTF",
    secreted_proteins == "TSLP [NSU]" ~ "TSLP",
    secreted_proteins == "Flt-3 Ligand [NSU]" ~ "Flt-3 Ligand",
    secreted_proteins == "CD14 [NSU]" ~ "CD14",
    secreted_proteins == "IL-16 [NSU]" ~ "IL-16",
    secreted_proteins == "FGF-21 [NSU]" ~ "FGF-21",
    secreted_proteins == "IL-29 [NSU]" ~ "IL-29",
    secreted_proteins == "IL-17C [NSU]" ~ "IL-17C",
    secreted_proteins == "IFN-epsilon [NSU]" ~ "IFN-epsilon",
    secreted_proteins == "PCSK9 [NSU]" ~ "PCSK9",
    secreted_proteins == "TPO (Thrombopoietin) [NSU]" ~ "Thrombopoietin",
    secreted_proteins == "TREM2 [NSU]" ~ "TREM2",
    secreted_proteins == "Growth Hormone (Somatotropin) [NSU]" ~ "Somatotropin",
    secreted_proteins == "CCL1 [NSU]" ~ "CCL1",
    secreted_proteins == "LOX1 (OLR1) [NSU]" ~ "LOX1 (OLR1)",
    secreted_proteins == "MMP-3 [NSU]" ~ "MMP-3",
    secreted_proteins == "IL-32 (alpha) [NSU]" ~ "IL-32 (alpha)",
    secreted_proteins == "IL-7 [NSU]" ~ "IL-7",
    secreted_proteins == "CCL21 [NSU]" ~ "CCL21",
    secreted_proteins == "CD276 (B7-H3) [NSU]" ~ "CD276 (B7-H3)",
    secreted_proteins == "IL-2 RA [NSU]" ~ "IL-2 RA",
    secreted_proteins == "Calbindin [NSU]" ~ "Calbindin",
    secreted_proteins == "CCL3 [NSU]" ~ "CCL3",
    secreted_proteins == "ICAM-1 [NSU]" ~ "ICAM-1",
    secreted_proteins == "IL-17A [NSU]" ~ "IL-17A",
    secreted_proteins == "CCL28 [NSU]" ~ "CCL28",
    secreted_proteins == "TIMP1 [NSU]" ~ "TIMP1",
    secreted_proteins == "GDF-15 (MIC-1) [NSU]" ~ "GDF-15 (MIC-1)",
    secreted_proteins == "CXCL17 [NSU]" ~ "CXCL17",
    secreted_proteins == "M-CSF R (CD115) [NSU]" ~ "M-CSF R (CD115)",
    secreted_proteins == "CCL7 [NSU]" ~ "CCL7",
    secreted_proteins == "Granzyme B [NSU]" ~ "Granzyme B",
    secreted_proteins == "CXCL4 [NSU]" ~ "CXCL4",
    secreted_proteins == "PDGF-BB [NSU]" ~ "PDGF-BB",
    secreted_proteins == "CX3CL1 [NSU]" ~ "CX3CL1",
    secreted_proteins == "FGF-6 [NSU]" ~ "FGF-6",
    secreted_proteins == "IL-35 [NSU]" ~ "IL-35",
    secreted_proteins == "MMP-7 [NSU]" ~ "MMP-7",
    secreted_proteins == "GM-CSF [NSU]" ~ "GM-CSF",
    secreted_proteins == "CCL24 [NSU]" ~ "CCL24",
    secreted_proteins == "IL-12 p40 [NSU]" ~ "IL-12 p40",
    secreted_proteins == "IL-5 [NSU]" ~ "IL-5",
    secreted_proteins == "BCMA (TNFRSF17) [NSU]" ~ "BCMA (TNFRSF17)",
    secreted_proteins == "Tissue Factor (TF) [NSU]" ~ "Tissue Factor",
    secreted_proteins == "IL-1 beta [NSU]" ~ "IL-1 beta",
    secreted_proteins == "CD30 [NSU]" ~ "CD30",
    secreted_proteins == "CCL27 [NSU]" ~ "CCL27",
    secreted_proteins == "ICAM-2 [NSU]" ~ "ICAM-2",
    secreted_proteins == "CXCL16 [NSU]" ~ "CXCL16",
    secreted_proteins == "VEGF-A (165) [NSU]" ~ "VEGF-A (165)",
    secreted_proteins == "IL-2 [NSU]" ~ "IL-2",
    secreted_proteins == "HVEM [NSU]" ~ "HVEM",
    secreted_proteins == "PTX3 (Pentraxin 3) [NSU]" ~ "PTX3",
    secreted_proteins == "IL-1 alpha [NSU]" ~ "IL-1 alpha",
    secreted_proteins == "CXCL3 [NSU]" ~ "CXCL3",
    secreted_proteins == "Oncostatin M (OSM) [NSU]" ~ "Oncostatin M",
    secreted_proteins == "CCL8 [NSU]" ~ "CCL8",
    secreted_proteins == "CCL15 [NSU]" ~ "CCL15",
    secreted_proteins == "FLRG (FSTL3) [NSU]" ~ "FLRG",
    secreted_proteins == "CXCL5 [NSU]" ~ "CXCL5",
    secreted_proteins == "CD163 [NSU]" ~ "CD163",
    secreted_proteins == "IL-17E (IL-25) [NSU]" ~ "IL-17E",
    secreted_proteins == "NF-L [NSU]" ~ "NF-L",
    secreted_proteins == "IFN alpha 2 (alpha 2b) [NSU]" ~ "IFN alpha 2",
    secreted_proteins == "TNF RI [NSU]" ~ "TNF RI",
    secreted_proteins == "CD40L [NSU]" ~ "CD40L",
    secreted_proteins == "IFN beta [NSU]" ~ "IFN beta",
    secreted_proteins == "VEGF Receptor 2 (Flk-1) [NSU]" ~ "VEGF Receptor 2",
    secreted_proteins == "BDNF [NSU]" ~ "BDNF",
    secreted_proteins == "Amyloid beta [NSU]" ~ "Amyloid beta",
    secreted_proteins == "MMP-2 [NSU]" ~ "MMP-2",
    secreted_proteins == "SAA [NSU]" ~ "SAA",
    secreted_proteins == "uPA [NSU]" ~ "uPA",
    secreted_proteins == "IL-22 BP [NSU]" ~ "IL-22 BP",
    secreted_proteins == "TRAIL [NSU]" ~ "TRAIL",
    secreted_proteins == "Mesothelin [NSU]" ~ "Mesothelin",
    secreted_proteins == "Activin A [NSU]" ~ "Activin A",
    secreted_proteins == "MMP-9 [NSU]" ~ "MMP-9",
    secreted_proteins == "CCL13 [NSU]" ~ "CCL13",
    secreted_proteins == "CXCL11 [NSU]" ~ "CXCL11",
    secreted_proteins == "IL-31 [NSU]" ~ "IL-31",
    secreted_proteins == "MIF [NSU]" ~ "MIF",
    secreted_proteins == "BMP7 [NSU]" ~ "BMP7",
    secreted_proteins == "IL-12 p70 [NSU]" ~ "IL-12 p70",
    secreted_proteins == "CCL19 [NSU]" ~ "CCL19",
    secreted_proteins == "CCL5 [NSU]" ~ "CCL5",
    secreted_proteins == "IL-33 [NSU]" ~ "IL-33",
    secreted_proteins == "IL-22 [NSU]" ~ "IL-22",
    secreted_proteins == "CCL11 [NSU]" ~ "CCL11",
    secreted_proteins == "IL-8 [NSU]" ~ "IL-8",
    secreted_proteins == "SCF [NSU]" ~ "SCF",
    secreted_proteins == "TNF RII [NSU]" ~ "TNF RII",
    secreted_proteins == "FGF-2 [NSU]" ~ "FGF-2",
    secreted_proteins == "Leptin [NSU]" ~ "Leptin",
    secreted_proteins == "CXCL13 [NSU]" ~ "CXCL13",
    secreted_proteins == "TNF alpha [NSU]" ~ "TNF alpha",
    secreted_proteins == "IL-4 [NSU]" ~ "IL-4",
    secreted_proteins == "CCL23 [NSU]" ~ "CCL23",
    secreted_proteins == "IGF-1 [NSU]" ~ "IGF-1",
    secreted_proteins == "FGF-4 [NSU]" ~ "FGF-4",
    secreted_proteins == "GDF-11 (BMP-11) [NSU]" ~ "GDF-11 (BMP-11)",
    secreted_proteins == "IL-10 [NSU]" ~ "IL-10",
    secreted_proteins == "IL-23 [NSU]" ~ "IL-23",
    secreted_proteins == "TNF RIII (Lymphotoxin Beta R) [NSU]" ~ "TNF RIII",
    secreted_proteins == "IL-17B [NSU]" ~ "IL-17B",
    secreted_proteins == "ST2 (IL-33R) [NSU]" ~ "ST2 (IL-33R)",
    secreted_proteins == "PLGF [NSU]" ~ "PLGF",
    secreted_proteins == "VEGF-D [NSU]" ~ "VEGF-D",
    secreted_proteins == "XCL1 (Lymphotactin) [NSU]" ~ "XCL1",
    secreted_proteins == "GDNF [NSU]" ~ "GDNF",
    secreted_proteins == "C5 [NSU]" ~ "C5",
    secreted_proteins == "IL-1 RA" ~ "IL-1 RA",
    secreted_proteins == "IL-17D [NSU]" ~ "IL-17D",
    secreted_proteins == "IL-27 [NSU]" ~ "IL-27",
    secreted_proteins == "Osteopontin (OPN) [NSU]" ~ "Osteopontin",
    secreted_proteins == "FGF-9 [NSU]" ~ "FGF-9",
    secreted_proteins == "BAFF [NSU]" ~ "BAFF",
    secreted_proteins == "TGF-beta 3 [NSU]" ~ "TGF-beta 3",
    secreted_proteins == "EGF [NSU]" ~ "EGF",
    secreted_proteins == "IL-5 [NSU]" ~ "IL-5",
    secreted_proteins == "FGF-7 (KGF) [NSU]" ~ "FGF-7 (KGF)",
    secreted_proteins == "APRIL [NSU]" ~ "APRIL",
    secreted_proteins == "WISP-1 (CCN4) [NSU]" ~ "WISP-1 (CCN4)",
    secreted_proteins == "CCL22 [NSU]" ~ "CCL22",
    secreted_proteins == "FGF-19 [NSU]" ~ "FGF-19",
    secreted_proteins == "M-CSF [NSU]" ~ "M-CSF",
    secreted_proteins == "CXCL10 [NSU]" ~ "CXCL10",
    secreted_proteins == "TGF-beta 1 (total) [NSU]" ~ "TGF-beta 1 ",
    secreted_proteins == "Tie-2 [NSU]" ~ "Tie-2",
    secreted_proteins == "TGF-beta 1 (LAP domain in precursor) [NSU]" ~ "TGF-beta 1",
    secreted_proteins == "FGFR3 (IIIc) [NSU]" ~ "FGFR3 (IIIc)",
    secreted_proteins == "AITRL (GITR Ligand) [NSU]" ~ "AITRL (GITR Ligand)",
    secreted_proteins == "Amphiregulin [NSU]" ~ "Amphiregulin",
    secreted_proteins == "BMP4 [NSU]" ~ "BMP4",
    secreted_proteins == "G-CSF [NSU]" ~ "G-CSF",
    secreted_proteins == "TGF-beta 2 [NSU]" ~ "TGF-beta 2",
    secreted_proteins == "IL-6 R alpha [NSU]" ~ "IL-6 R alpha",
    secreted_proteins == "BMP6 [NSU]" ~ "BMP6",
    secreted_proteins == "NGF beta [NSU]" ~ "NGF beta",
    secreted_proteins == "IL-1 R1 [NSU]" ~ "IL-1 R1",
    secreted_proteins == "MMP-10 [NSU]" ~ "MMP-10",
    secreted_proteins == "IL-17F [NSU]" ~ "IL-17F",
    secreted_proteins == "IL-18 [NSU]" ~ "IL-18",
    secreted_proteins == "CXCL6 [NSU]" ~ "CXCL6",
    secreted_proteins == "IL-6 [NSU]" ~ "IL-6",
    secreted_proteins == "CXCL1 [NSU]" ~ "CXCL1",
    secreted_proteins == "VEGF-C [NSU]" ~ "VEGF-C",
    secreted_proteins == "Resistin [NSU]" ~ "Resistin",
    secreted_proteins == "EMMPRIN [NSU]" ~ "EMMPRIN",
    secreted_proteins == "IFN gamma [NSU]" ~ "IFN gamma",
    secreted_proteins == "CCL20 [NSU]" ~ "CCL20",
    secreted_proteins == "CRP [NSU]" ~ "CRP",
    secreted_proteins == "VCAM-1 [NSU]" ~ "VCAM-1",
    secreted_proteins == "Cytochrome C [NSU]" ~ "Cytochrome C",
    secreted_proteins == "BMP3 [NSU]" ~ "BMP3",
    secreted_proteins == "IL-24 [NSU]" ~ "IL-24",
    secreted_proteins == "IL-28A [NSU]" ~ "IL-28A",
    secreted_proteins == "CCL17 [NSU]" ~ "CCL17",
    secreted_proteins == "BMP2 [NSU]" ~ "BMP2",
    secreted_proteins == "CD27L [NSU]" ~ "CD27L",
    secreted_proteins == "NRG1 beta 1 [NSU]" ~ "NRG1 beta 1",
    secreted_proteins == "IL-11 [NSU]" ~ "IL-11",
    TRUE ~ secreted_proteins
  ))


# select MMP-1 secreted protein as the target
df <- df %>% filter(shuffle == "final")


# get the feature names for color bar visualization
features <- df %>% select(feature_names)
# drop duplicate features from the feature names
features <- unique(features)
features <- features %>%
        tidyr::separate(
            feature_names,
            into = c(
                "compartment",
                "feature_group",
                "measurement",
                "channel",
                "parameter1",
                "parameter2"
            ),
            sep = "_",
            remove = FALSE
        ) %>%
        dplyr::mutate(channel_cleaned = channel)

    # Clean channel for visualization
    features$channel_learned <- dplyr::recode(features$channel,
            "CorrDNA" = "Nuclei",
            "CorrMito" = "Mito",
            "CorrER" = "ER",
            "CorrGasdermin" = "GasderminD",
            "CorrPM" = "AGP",
            .default = "Other",
            .missing="Other"
    )
# make the channel learned a factor
features$channel_learned <- factor(features$channel_learned, levels = c("Nuclei", "Mito", "ER", "Gasdermin", "AGP", "Other"))


legend_title_size <- 28
legend_text_size <- 24

r2_df <- df %>% select(r2)
r2_df <- unique(r2_df)
column_ha <- HeatmapAnnotation(
    R2 = r2_df$r2,
    show_legend = TRUE,
    annotation_name_side = "left",
    # rotate the title
    annotation_legend_param = list(
        labels_gp = gpar(fontsize = legend_text_size, angle = 0),
        title_position = "topcenter",
        title_gp = gpar(fontsize = legend_title_size, hjust = 0.5, fontface = "bold")
    ),
    annotation_name_gp = gpar(fontsize = 20),
    # set color bar for r2 continuous value with brewer palette
    col = list(R2 = colorRamp2(c(0, 0.5, 1), spectral_palette <- c(
        # white
        "#FFFFFF",
        # light blue
        "#A6CEE3",
        # dark blue
        "#1F78B4"
    )))

)


# make the df into a matrix for heatmap
mat <- dcast(df, feature_names ~ secreted_proteins, value.var = "coefficients")
row.names(mat) <- mat$feature_names
mat <- mat %>% select(-feature_names)
mat <- as.matrix(mat)
# na to 0
mat[is.na(mat)] <- 0

legend_height <- unit(4, "cm")
legend_width <- unit(10, "cm")
mat <- as.data.frame(mat)
# mat <- as.matrix(mat)
# get the feature names from the index
mat$feature_names <- row.names(mat)
# get the feature names for color bar visualization
# drop duplicate features from the feature names
features <- unique(features)
features <- features %>%
        tidyr::separate(
            feature_names,
            into = c(
                "compartment",
                "feature_group",
                "measurement",
                "channel",
                "parameter1",
                "parameter2"
            ),
            sep = "_",
            remove = FALSE
        ) %>%
        dplyr::mutate(channel_cleaned = channel)

    # Clean channel for visualization
    features$channel_learned <- dplyr::recode(features$channel,
            "CorrDNA" = "Nuclei",
            "CorrMito" = "Mito",
            "CorrER" = "ER",
            "CorrGasdermin" = "GasderminD",
            "CorrPM" = "AGP",
            .default = "Other",
            .missing="Other"
    )

# set annotations
row_ha_1 <- rowAnnotation(
    Object = features$compartment,
    show_legend = TRUE,
    # change the legend titles
    annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = legend_title_size, hjust = 0.5,fontface = "bold"),
        labels_gp = gpar(fontsize = legend_text_size),
        legend_width = unit(legend_width, "cm"),
        legend_height = unit(legend_height, "cm"),
        direction = "horizontal"

        ),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 20),
    # color
    col = list(
        Object = c(
            "Cells" = brewer.pal(12, "Accent")[7],
            "Cytoplasm" = brewer.pal(12, "Accent")[6],
            "Nuclei" = "#0000AB"
        )
    )
)

row_ha_2 <- rowAnnotation(
        FeatureType = features$feature_group,
       annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = legend_title_size, hjust = 0.5,fontface = "bold"),
        labels_gp = gpar(fontsize = legend_text_size),
        legend_width = unit(legend_width, "cm"),
        legend_height = unit(legend_height, "cm"),
        direction = "horizontal"
        ),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 20),
    col = list(
            FeatureType = c(
            "AreaShape" = brewer.pal(8, "Dark2")[1],
            "Correlation" = brewer.pal(8, "Dark2")[2],
            "Granularity" = brewer.pal(8, "Dark2")[3],
            "Intensity" = brewer.pal(8, "Dark2")[4],
            "Location" = brewer.pal(8, "Dark2")[5],
            "Neighbors" =  brewer.pal(8, "Dark2")[6],
            "RadialDistribution" = brewer.pal(8, "Dark2")[7],
            "Texture" = brewer.pal(8, "Dark2")[8]

    )
)
)

row_ha_3 <- rowAnnotation(
    Channel = features$channel_learned,
    annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = legend_title_size, hjust = 0.5,fontface = "bold"),
        labels_gp = gpar(fontsize = legend_text_size),
        # make annotation bar text bigger
        legend = gpar(fontsize = legend_text_size),
        annotation_name = gpar(fontsize = 16),
        # make legend taller
        legend_width = unit(legend_width, "cm"),
        legend_height = unit(legend_height, "cm"),
        legend_key = gpar(fontsize = legend_text_size),
        direction = "horizontal"


            ),

    annotation_name_side = "bottom",
    # make font size bigger
    annotation_name_gp = gpar(fontsize = 20),
    col = list(
    Channel = c(
            "Nuclei" = "#0000AB",
            "Mito" = "#B000B0",
            "ER" = "#00D55B",
            "GasderminD" = "#FFFF00",
            "AGP" = "#C90000",
            "Other" = "#B09FB0")
    )
)

# # drop the feature names column
mat <- mat %>% select(-feature_names)
mat <- as.matrix(mat)

# plot size
width <- 23
height <- 23
options(repr.plot.width=width, repr.plot.height=height)
# change margins
# par(mar = c(1, 1, 1, 1))

model_heatmap <- (
        Heatmap(
        (mat),
        cluster_rows = TRUE,    # Cluster rows
        cluster_columns = TRUE, # Cluster columns
        show_row_names = FALSE,  # Show row names
        show_column_names = FALSE, # Show column names
        column_names_gp = gpar(fontsize = 16), # Column name label formatting
        row_names_gp = gpar(fontsize = 14),    # Row name label formatting
        right_annotation = c(row_ha_1,row_ha_2,row_ha_3),
        bottom_annotation = column_ha,
        # rename fill legend
        heatmap_legend_param = list(
                title = "Coef",
                title_position = "topcenter",
                title_gp = gpar(fontsize = legend_title_size, hjust = 0.5,fontface = "bold"),
                labels_gp = gpar(fontsize = legend_text_size),
                legend_height = unit(legend_height, "cm"),
                legend_width = unit(legend_width, "cm")
                ),
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),


        )
)
model_heatmap <- draw(

    model_heatmap,
    merge_legend = TRUE,
    padding = unit(c(1, 1, 1, 1), "cm")
)

# save the plot
png(filename = paste0(figure_path, "pt1.png"), width = width, height = height, units = "in", res = 600)
model_heatmap
dev.off()

width <- 17
height <- 23
options(repr.plot.width=width, repr.plot.height=height)

design2 <- "
            AAACCC
            DDDEEE
            "

pt2 <- (
    wrap_elements(full = variance_r2_plot_global)
    + wrap_elements(full = global_prediction_trend_scatter)
    + cytokine_predictions1
    + cytokine_predictions2
    + plot_layout(design = design2)
    + plot_annotation(tag_levels = "A")  & theme(plot.tag = element_text(size = 20))

)
pt2

# functions to remove x or y axis labels
remove_x <- function(plot){
    plot <- plot + theme(axis.title.x = element_blank(),  # Remove x-axis label
        axis.text.x = element_blank())   # Remove x-axis text
    return(plot)
}

# remove y axis labels
remove_y <- function(plot){
    plot <- plot + theme(axis.title.y = element_blank(),  # Remove y-axis label
        axis.text.y = element_blank())   # Remove y-axis text
    return(plot)
}

# remove the x axis information from certain plots
il1beta_final_plot <- remove_x(il1beta_final_plot)
tnfa_final_plot <- remove_x(tnfa_final_plot)
CCL24_final_plot <- remove_x(CCL24_final_plot)
il18beta_final_plot <- remove_x(il18beta_final_plot)
op_final_plot <- remove_x(op_final_plot)
CCL13_final_plot <- remove_x(CCL13_final_plot)

# remove the y axis information from certain plots
tnfa_final_plot <- remove_y(tnfa_final_plot)
CCL24_final_plot <- remove_y(CCL24_final_plot)
op_final_plot <- remove_y(op_final_plot)
CCL13_final_plot <- remove_y(CCL13_final_plot)
il6_final_plot <- remove_y(il6_final_plot)
CCL4_final_plot <- remove_y(CCL4_final_plot)


height <- 15
width <- 17
options(repr.plot.width=width, repr.plot.height=height)

pt3 <- (
    il1beta_final_plot
    + tnfa_final_plot

    + CCL24_final_plot

    +il18beta_final_plot


    + op_final_plot
    + CCL13_final_plot
    + il2_final_plot
    + il6_final_plot
    + CCL4_final_plot
    + plot_annotation(tag_levels = "A")  & theme(plot.tag = element_text(size = 20))
)

ggsave(
    filename = paste0(figure_path, "pt2.png"),
    plot = pt2,
    width = width,
    height = height,
    units = "in",
    dpi = 600
)
width <- 23
height <- 15
options(repr.plot.width=width, repr.plot.height=height)
ggsave(
    filename = paste0(figure_path, "pt3.png"),
    plot = pt3,
    width = width,
    height = height,
    units = "in",
    dpi = 600
)
pt3
