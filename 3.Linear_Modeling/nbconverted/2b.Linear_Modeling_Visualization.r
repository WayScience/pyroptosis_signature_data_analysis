suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

lm_file <- file.path("./results/lm_two_beta.tsv")

lm_cp_fig <- file.path("./figures/two_beta_results.pdf" )

lm_df <- readr::read_tsv(lm_file, col_types = readr::cols(.default = "d", feature ="c", inducer1_inhibitor_inhibitor_dose__inducer1_dose = "c"))
head(lm_df)

unique(lm_df$inducer1_inhibitor_inhibitor_dose__inducer1_dose)


# Arrange by absolute value coefficient
# Split out components of feature name for visualization
lm_df <- lm_df %>%
    dplyr::arrange(desc(abs(twob_Metadata_Treatment_Inhibitor_Dose))) %>%
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
    dplyr::mutate(channel_cleaned = channel)



unique(lm_df$channel)

unique(lm_df$channel_cleaned)

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

print(dim(lm_df))
head(lm_df)
unique(lm_df$channel_learned)
lm_df$abs_Metadata_number_of_singlecells <- abs(lm_df$Metadata_number_of_singlecells)
lm_df$abs_twob_Metadata_Treatment_Inhibitor_Dose <- abs(lm_df$twob_Metadata_Treatment_Inhibitor_Dose)
lm_df$abs_Treatment_Dose <- abs(lm_df$Treatment_Dose)

head(lm_df)

loop_list <- unique(lm_df$inducer1_inhibitor_inhibitor_dose__inducer1_dose)
x_list <- c('abs_twob_Metadata_Treatment_Inhibitor_Dose','abs_Treatment_Dose')

pdf(file=lm_cp_fig )
for (i in 1:length(loop_list)){
    df <- lm_df[lm_df$inducer1_inhibitor_inhibitor_dose__inducer1_dose == loop_list[i],]
    for (j in 1:length(x_list)){
        lm_fig_gg <- (
            ggplot(df, aes(x = .data[[x_list[j]]], y = r2_score))
            + geom_point(aes(size = abs_Metadata_number_of_singlecells, color = channel_learned), alpha = 0.7)
            + theme_bw()
            + guides(
                color = guide_legend(title = "Channel\n(if applicable)", order = 1),
                size = guide_legend(title = "Cell count contributution")
            )
            + ylab("R2 score of LM feature")
            + xlab(paste0(x_list[j]," contribution (LM beta coefficient)"))
            + ggtitle(paste0("How CellProfiler features contribute\nto ",loop_list[i], "\ntreatments and cell density"))
        )
    plot(lm_fig_gg)
    }
}
dev.off()


