# import the libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(cowplot))

figure_theme_path <- file.path(
    "..","visulaization_utils", "figure_themes.R")
source(figure_theme_path)

# define command line arguments
parser <- ArgumentParser(description = "Visualize MLP results")
# add arguments
parser$add_argument('--cell_type', type='character', help='Cell type to visualize')
parser$add_argument('--model_name', type='character', help='Model name to visualize')
parser$add_argument('--selected_treatment_comparisons', type='character', help='Selected treatment comparisons to visualize')

# parse arguments from command line
args <- parser$parse_args()

# define cell type
celltype <- args$cell_type
model_name <- args$model_name
selected_treatment_comparisons <- args$selected_treatment_comparisons

celltype <- "PBMC"
model_name <- "LPS_100.000_DMSO_0.025_vs_DMSO_0.100_DMSO_0.025"
selected_treatment_comparisons <- "DMSO_0.100_DMSO_0.025 vs LPS_100.000_DMSO_0.025, DMSO_0.100_DMSO_0.025 vs Thapsigargin_10.000_DMSO_0.025, LPS_100.000_DMSO_0.025 vs Thapsigargin_10.000_DMSO_0.025"


output_file_plot <- file.path(
    "..","..","figures","Binary_Classification",model_name,celltype,"pr_curves_testing_plot.png"
)
output_file_legend <- file.path(
    "..","..","figures","Binary_Classification",model_name,celltype,"pr_curves_testing_legend.png"
)

results_dir <- file.path(
    "..","..","results","Binary_Classification",model_name,celltype
)
results_file <- file.path(
    results_dir,"testing_metrics.csv"
)

# Read in the results file
df <- read.csv(results_file)
head(df,3)

unique(df$treatments_tested)

selected_treatment_comparisons <- unlist(strsplit(selected_treatment_comparisons,", "))
selected_treatment_comparisons

# filter the df by selected treatment comparisons in the treatment column
tmp_df <- df[which(df$treatment %in% selected_treatment_comparisons),]
unique(tmp_df$treatments_tested)

head(tmp_df)
tmp_df$treatments_tested <- gsub(" vs ", "\n", tmp_df$treatments_tested)
head(tmp_df)

pr_curve_gg <- (
    ggplot(tmp_df, aes(x = Recall, y = Precision))
    + geom_line(aes(color = treatments_tested, linetype = shuffled_data))
    + theme_bw()
    + xlab("Recall")
    + ylab("Precision")

    + scale_linetype_manual(
        name = "Shuffled\ntraining\ndata",
        labels = shuffled_labels,
        values = shuffled_linetypes
    )

    + guides(
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 2),
    )
    + coord_fixed()
    + figure_theme
    # Decrease spacing in legend
    + theme(
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.2, "cm"),
        legend.key.size = unit(2.5, "lines"),
        legend.key.width = unit(1, "lines")
    )
    + ggtitle(paste0("Precision-Recall Curve for ","\n", model_name, " model"))
)
# detach the legend
legend <- get_legend(pr_curve_gg)
pr_curve_gg <- pr_curve_gg + theme(legend.position = "none")



ggsave(output_file_plot, pr_curve_gg, height = 5.5, width = 8.5, dpi = 500)
ggsave(output_file_legend, legend, height = 5.5, width = 8.5, dpi = 500)
