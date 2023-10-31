library(ggplot2)
library(argparse)

# set up argument parser
parser <- ArgumentParser(description = "Plot redundancy of a set of sequences")
parser$add_argument("-c", "--cell_type", help = "Cell_type", required = TRUE)

# parse arguments
args <- parser$parse_args()

cell_type <- args$cell_type

# define output fraph paths
skree_plot_path <- file.path(paste0("../figures/", cell_type, "_skree_plot.png"))
redundancy_index_plot_path <- file.path(paste0("../figures/", cell_type, "_redundancy_index_plot.png"))

redundancy_file_path <- file.path(paste0(
    "../",
    "results/",
    cell_type,
    "_redundancy_analysis.csv"))

redundancy_df <- read.csv(redundancy_file_path)
head(redundancy_df)


# ggplot visualization for redundancy analysis line plot
line_plot <- (
    ggplot(redundancy_df, aes(x=K, y=r2, color=Shuffle))
    + geom_line()
    + theme_bw()
    + ylim(-1,1)
    + labs(x="Components (K)", y="R2")
    + ggtitle("Scree Plot for \nCanonical Correlation Analysis \nComponents (K)")
    + theme(plot.title = element_text(hjust = 0.5, size=20),
            legend.title = element_blank(),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=14),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14)
    )

)
ggsave(skree_plot_path, line_plot, width=10, height=8, dpi=600)


# read in the redundancy index results

# set path
redundancy_index_file_path <- file.path(paste0(
    "../",
    "results/",
    cell_type,
    "_redundancy_index.csv"))

# read in the redundancy index results
redundancy_index_df <- read.csv(redundancy_index_file_path)
head(redundancy_index_df)


# find the min and max redundancy index for each K
minimum_value = min(min(redundancy_index_df$u_k_RI),min(redundancy_index_df$v_k_RI))
maximum_value = max(max(redundancy_index_df$u_k_RI),max(redundancy_index_df$v_k_RI))

RI_plot <- (
    ggplot(redundancy_index_df, aes(x=u_k_RI, y=v_k_RI, color=Shuffle))
    + geom_point()
    + theme_bw()
    + xlim(minimum_value, maximum_value)
    + ylim(minimum_value, maximum_value)
    + xlab("Morphology Data Redundancy Index")
    + ylab("nELISA Data Redundancy Index")
    + geom_abline(intercept = 0, slope = 1)
    + ggtitle("Redundancy Index Plot of \nMorphology and nELISA Data")
    + theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5)

    )

)
ggsave(redundancy_index_plot_path, RI_plot, width=10, height=8, dpi=600)
