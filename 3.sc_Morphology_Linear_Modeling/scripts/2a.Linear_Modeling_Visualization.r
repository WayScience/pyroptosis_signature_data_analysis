suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparse))

# define command line arguments
parser <- ArgumentParser(description = "Visualize linear modeling results")
# add arguments
parser$add_argument('--celltype', type='character', help='Cell type to visualize')

# parse arguments from command line
args <- parser$parse_args()

# define cell type
celltype <- args$celltype


# lm file path
lm_file <- file.path(paste0("./results/", celltype, "/lm_one_beta.tsv"))

# figure paths
# scatter plot of beta values
lm_cp_fig <- file.path(paste0("./figures/", celltype, "/lm_one_beta_scatter.pdf"))
# scatter plot of absolute beta values
lm_cp_fig_abs <- file.path(paste0("./figures/", celltype, "/lm_one_beta_scatter_abs.pdf"))
# scatter plot of beta values with facet and feature type coloring
lm_facet_fig <- file.path(paste0("./figures/", celltype, "/lm_one_beta_facet_beta.pdf"))
# plot of beta values per cellular compartment and feature type
lm_coef_fig <- file.path(paste0("./figures/", celltype, "/lm_one_beta_coef_per_compartment.pdf"))

# if path does not exist, create it
if (!dir.exists(file.path(paste0("./figures/", celltype)))) {
    dir.create(file.path(paste0("./figures/", celltype)))
}

# read in linear modeling results
lm_df <- readr::read_tsv(lm_file, col_types = readr::cols(.default = "d", feature ="c", dosage_treatments_list = "c"))
head(lm_df, 2)

lm_df <- lm_df %>%
    dplyr::arrange(desc(abs(oneb_Metadata_Treatment_Dose_Inhibitor_Dose))) %>%
    tidyr::separate(
        feature,
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
    dplyr::arrange(desc(abs(oneb_Metadata_Treatment_Dose_Inhibitor_Dose)))

# Clean channel for visualization
lm_df$channel_learned <- dplyr::recode(lm_df$channel_cleaned,
        "CorrDNA" = "nuclei",
        "CorrMito" = "Mito",
        "CorrER" = "ER",
        "CorrGasdermin" = "gasdermin",
        "CorrPM" = "PM",
        .default = "other",
        .missing="other"
    )

head(lm_df, 2)

print(dim(lm_df))
head(lm_df, 2)
unique(lm_df$channel_learned)
lm_df$abs_Metadata_number_of_singlecells <- abs(lm_df$Metadata_number_of_singlecells)
lm_df$abs_oneb_Metadata_Treatment_Dose_Inhibitor_Dose <- abs(lm_df$oneb_Metadata_Treatment_Dose_Inhibitor_Dose)

loop_list <- unique(lm_df$dosage_treatments_list)
# drop 'DMSO_0.100_DMSO_0.025-DMSO_0.100_DMSO_0.025' from loop_list to avoid error in plotting
loop_list <- loop_list[!grepl('DMSO_0.100_DMSO_0.025-DMSO_0.100_DMSO_0.025', loop_list)]

pdf(file=lm_cp_fig)
for (i in 1:length(loop_list)){

    df <- lm_df[lm_df$dosage_treatments_list == loop_list[i],]
    # define the treatment without the control group (DMSO_0.025)
    treatment = strsplit(loop_list[i], "-D")[[1]][1]

    lm_fig_gg <- (
        ggplot(df, aes(x = Metadata_number_of_singlecells,, y = oneb_Metadata_Treatment_Dose_Inhibitor_Dose))

        + geom_point(aes(size = r2_score, color = channel_learned,), alpha = 0.7)

        + scale_size_continuous(range = c(2, 8), limits = c(0, 1))


        + geom_vline(xintercept = 0, linetype = "dashed", color = "red")
        + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
        + geom_density2d(color="black", show.legend = FALSE)
        + theme_bw()
        + guides(
            color = guide_legend(title = "Channel\n(if applicable)", order = 1),
            size = guide_legend(title = "R2 score")
        )
        # make legend dots bigger
        + ylab("Treatment contribution (LM beta coefficient)")
        + xlab("Cell count contribution (LM beta coefficient)")
        + ggtitle(paste0("How CellProfiler features contribute\nto ",treatment, " treatments and cell density"))
    )
    plot(lm_fig_gg)
}
dev.off()


pdf(file=lm_cp_fig_abs)
for (i in 1:length(loop_list)){

    df <- lm_df[lm_df$dosage_treatments_list == loop_list[i],]
    # define the treatment without the control group (DMSO_0.025)
    treatment = strsplit(loop_list[i], "-D")[[1]][1]
    lm_fig_gg <- (
        ggplot(df, aes(x = abs_oneb_Metadata_Treatment_Dose_Inhibitor_Dose, y = r2_score))
        + geom_point(aes(size = abs_Metadata_number_of_singlecells, color = channel_learned), alpha = 0.7)
        + scale_size_continuous(range = c(2, 8), limits = c(0, 1))


        + theme_bw()
        + guides(
            color = guide_legend(title = "Channel\n(if applicable)", order = 1),
            size = guide_legend(title = "Cell count contributution")
        )
        + ylab("R2 score of LM feature")
        + xlab("Treatment and Dose contribution (LM beta coefficient)")
        + ggtitle(paste0("How CellProfiler features contribute\nto ",treatment, " treatments and cell density"))
    )

    plot(lm_fig_gg)
}
dev.off()

pdf(file=lm_facet_fig)
for (i in 1:length(loop_list)){
    df <- lm_df[lm_df$dosage_treatments_list == loop_list[i],]
    # define the treatment without the control group (DMSO_0.025)
    treatment = strsplit(loop_list[i], "-D")[[1]][1]
    lm_facet_fig_gg <- (
        ggplot(df, aes(x = Metadata_number_of_singlecells, y = oneb_Metadata_Treatment_Dose_Inhibitor_Dose))
        + geom_point(aes(size = r2_score, color = feature_group), alpha = 0.7)
        + facet_wrap("~channel_learned")
        + geom_vline(xintercept = 0, linetype = "dashed", color = "red")
        + geom_hline(yintercept = 0, linetype = "dashed", color = "red")
        + scale_size_continuous(range = c(2, 8), limits = c(0, 1))

        + theme_bw()
        + guides(
            color = guide_legend(title = "Feature group\n(if applicable)", order = 1),
            size = guide_legend(title = "R2 score")
        )
        # set x tick labels to scientific notation
        + scale_x_continuous(labels = scales::scientific_format())
        + ylab("Treatment dose contribution (LM beta coefficient)")
        + xlab("Cell count contribution (LM beta coefficient)")
        + ggtitle("How CellProfiler features (by group) contribute to treatment and cell density")
        + scale_color_brewer(palette="Dark2")
        + theme(
            axis.text = element_text(size = 7),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            strip.text = element_text(size = 8),
            strip.background = element_rect(
                colour = "black",
                fill = "#fdfff4"
            )
        )
        + ggtitle(paste0("LM Coefficients for: ",treatment))
    )
    plot(lm_facet_fig_gg)
}
dev.off()


head(lm_df,2)
# drop rows in the feature_group column that are Location
lm_df <- lm_df[lm_df$feature_group != "Location",]
# drop rows in the feature_group column that are Neighbors or AreaShape
lm_df <- lm_df[lm_df$feature_group != "Neighbors",]
lm_df <- lm_df[lm_df$feature_group != "AreaShape",]

unique(lm_df$channel_learned)

pdf(file=lm_coef_fig)
for (i in 1:length(loop_list)){
    df <- lm_df[lm_df$dosage_treatments_list == loop_list[i],]
    df <- df %>%
    dplyr::group_by(feature_group, channel_learned, compartment) %>%
    dplyr::slice_max(order_by = oneb_Metadata_Treatment_Dose_Inhibitor_Dose, n = 1)
    # define the treatment without the control group (DMSO_0.025)
    treatment = strsplit(loop_list[i], "-D")[[1]][1]
    coef_gg <- (
        ggplot(df, aes(x = channel_learned, y = feature_group))
        + geom_point(aes(fill = abs(oneb_Metadata_Treatment_Dose_Inhibitor_Dose)), pch = 22, size = 5)
        + facet_wrap("~compartment", ncol = 1)
        + theme_bw()
        + scale_fill_continuous(
            name="Top Abs. val\ntreatment\nlinear model\ncoefficient",
            low = "darkblue",
            high = "yellow",
            limits = c(
                min(abs(df$oneb_Metadata_Treatment_Dose_Inhibitor_Dose)),
                max(abs(df$oneb_Metadata_Treatment_Dose_Inhibitor_Dose))),
        )
        + xlab("Channel")
        + ylab("Feature")
        + theme(
            axis.text = element_text(size = 7),
            axis.text.x = element_text(angle = 90, size = 7),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 9),
            legend.title = element_text(size = 10),
            strip.text = element_text(size = 8),
            strip.background = element_rect(
                colour = "black",
                fill = "#fdfff4"
            )
        )
        + ggtitle(paste0("LM Coefficients for: ",treatment))
    )
    plot(coef_gg)
}
dev.off()
