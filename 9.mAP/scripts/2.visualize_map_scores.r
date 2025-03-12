suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))
suppressPackageStartupMessages(suppressWarnings(library(arrow)))
source("../../figures/utils/figure_themes.r")

# set path to the data morphology

morphology_path <- file.path("..","data","processed","mAP_scores","morphology","activity_map.parquet")
shuffled_morphology_path <- file.path("..","data","processed","mAP_scores","morphology","activity_map_shuffled.parquet")
# set path to the secretome data

secretome_path <- file.path("..","data","processed","mAP_scores","secretome","activity_map.parquet")
shuffled_secretome_path <- file.path("..","data","processed","mAP_scores","secretome","activity_map_shuffled.parquet")

df_morphology <- arrow::read_parquet(morphology_path) %>%
    dplyr::mutate(shuffled = "Non-shuffled") %>%
    dplyr::mutate(data_type = "Morphology") %>%
    # rename the mean_average_precision column to specifcy morphology
    # drop unnecessary columns
    dplyr::select(-c("Metadata_reference_index", "indices", "p_value", "corrected_p_value", "below_p", "below_corrected_p","-log10(p-value)"))

df_shuffled_morphology <- arrow::read_parquet(shuffled_morphology_path) %>%
    dplyr::mutate(shuffled = "Shuffled") %>%
    dplyr::mutate(data_type = "Morphology") %>%
    # rename the mean_average_precision column to specifcy morphology
    dplyr::select(-c("Metadata_reference_index", "indices", "p_value", "corrected_p_value", "below_p", "below_corrected_p","-log10(p-value)"))


df_secretome <- arrow::read_parquet(secretome_path) %>%
    dplyr::mutate(shuffled = "Non-shuffled") %>%
    dplyr::mutate(data_type = "Secretome") %>%
    # rename the mean_average_precision column to specifcy secretome
    dplyr::select(-c("Metadata_reference_index", "indices", "p_value", "corrected_p_value", "below_p", "below_corrected_p","-log10(p-value)"))


df_shuffled_secretome <- arrow::read_parquet(shuffled_secretome_path) %>%
    dplyr::mutate(shuffled = "Shuffled") %>%
    dplyr::mutate(data_type = "Secretome") %>%
    # rename the mean_average_precision column to specifcy secretome
    dplyr::select(-c("Metadata_reference_index", "indices", "p_value", "corrected_p_value", "below_p", "below_corrected_p","-log10(p-value)"))

df <- dplyr::bind_rows(df_morphology, df_shuffled_morphology, df_secretome, df_shuffled_secretome)
head(df)

# split out the morphology and secretome data
morphology_data <- df %>% dplyr::filter(data_type == "Morphology")
secretome_data <- df %>% dplyr::filter(data_type == "Secretome")
# rename the mean_average_precision column to specifcy morphology
morphology_data <- morphology_data %>% dplyr::rename(mAP_moprhology = mean_average_precision)
secretome_data <- secretome_data %>% dplyr::rename(mAP_secretome = mean_average_precision)
# drop the data_type column
morphology_data <- morphology_data %>% dplyr::select(-data_type)
secretome_data <- secretome_data %>% dplyr::select(-data_type)
# merge the data together to plot
df <- merge(morphology_data, secretome_data,by = c("Metadata_Treatment", "Metadata_labels", "shuffled"))

levels_list <- c(
    'Media',
    'DMSO_0.100_%_DMSO_0.025_%',
    'DMSO_0.100_%_DMSO_1.000_%',
    'DMSO_0.100_%_Z-VAD-FMK_30.000_uM',
    'DMSO_0.100_%_Z-VAD-FMK_100.000_uM',

    'Disulfiram_0.100_uM_DMSO_0.025_%',
    'Disulfiram_1.000_uM_DMSO_0.025_%',
    'Disulfiram_2.500_uM_DMSO_0.025_%',

    'Flagellin_0.100_ug_per_ml_DMSO_0.025_%',
    'Flagellin_1.000_ug_per_ml_DMSO_0.025_%',
    'Flagellin_1.000_ug_per_ml_Disulfiram_1.000_uM',

    'LPS_0.010_ug_per_ml_DMSO_0.025_%',
    'LPS_0.100_ug_per_ml_DMSO_0.025_%',
    'LPS_1.000_ug_per_ml_DMSO_0.025_%',

    'LPS_Nigericin_1.000_ug_per_ml_1.000_uM_DMSO_0.025_%',
    'LPS_Nigericin_1.000_ug_per_ml_3.000_uM_DMSO_0.025_%',
    'LPS_Nigericin_1.000_ug_per_ml_10.000_uM_DMSO_0.025_%',
    'LPS_Nigericin_1.000_ug_per_ml_10.000_uM_Disulfiram_1.000_uM',
    'LPS_Nigericin_1.000_ug_per_ml_10.000_uM_Z-VAD-FMK_100.000_uM',

    'LPS_10.000_ug_per_ml_DMSO_0.025_%',
    'LPS_10.000_ug_per_ml_Disulfiram_0.100_uM',
    'LPS_10.000_ug_per_ml_Disulfiram_1.000_uM',
    'LPS_10.000_ug_per_ml_Disulfiram_2.500_uM',
    'LPS_10.000_ug_per_ml_Z-VAD-FMK_100.000_uM',

    'LPS_100.000_ug_per_ml_DMSO_0.025_%',
    'LPS_Nigericin_100.000_ug_per_ml_1.000_uM_DMSO_0.025_%',
    'LPS_Nigericin_100.000_ug_per_ml_3.000_uM_DMSO_0.025_%',
    'LPS_Nigericin_100.000_ug_per_ml_10.000_uM_DMSO_0.025_%',

    'H2O2_100.000_nM_DMSO_0.025_%',
    'H2O2_100.000_uM_DMSO_0.025_%',
    'H2O2_100.000_uM_Disulfiram_1.000_uM',
    'H2O2_100.000_uM_Z-VAD-FMK_100.000_uM',
    'Thapsigargin_1.000_uM_DMSO_0.025_%',
    'Thapsigargin_10.000_uM_DMSO_0.025_%',

    'Topotecan_5.000_nM_DMSO_0.025_%',
    'Topotecan_10.000_nM_DMSO_0.025_%',
    'Topotecan_20.000_nM_DMSO_0.025_%'
)

df$Metadata_labels <- factor(df$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))
df$Metadata_Treatment <- factor(df$Metadata_Treatment, levels =levels_list)

# scatter plot
scatter_compare_treatment <- (
    ggplot(df, aes(x=mAP_moprhology, y=mAP_secretome, col = Metadata_labels, shape=shuffled))
    + geom_point(size=3, alpha=0.5)
    + labs(x="Morphology mAP score", y="Secretome mAP score")
    + theme_bw()
    + ggtitle("Comparison of mAP scores")
    + ylim(0,1)
    + xlim(0,1)
    # Change the legend title
    # change the legend shape
    + scale_shape_manual(
        name="Shuffle type",
        labels=c(
            "Non-shuffled",
            "Shuffled features",
            "Shuffled phenotypes"
        ),
        values=c(19, 17, 15)
    )
    + scale_color_manual(
        name="Class",
        labels=c(
            "Control",
            "Apoptosis",
            "Pyroptosis"
        ),
        values=c(
            brewer.pal(3, "Dark2")[2],
            brewer.pal(3, "Dark2")[1],
            brewer.pal(3, "Dark2")[3]
    )
)
    + figure_theme
    # add y = x line
    + geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "black")
)
scatter_compare_treatment

df <- df %>%
    mutate(Metadata_Treatment = case_when(
        Metadata_Treatment =='DMSO_0.100_%_DMSO_0.025_%' ~ "DMSO 0.1% - DMSO 0.025%",
        Metadata_Treatment =='DMSO_0.100_%_DMSO_1.000_%' ~ "DMSO 0.1% - DMSO 1.0%",
        Metadata_Treatment =='DMSO_0.100_%_Z-VAD-FMK_100.000_uM' ~ "DMSO 0.1% - Z-VAD-FMK 100.0 uM",
        Metadata_Treatment =='DMSO_0.100_%_Z-VAD-FMK_30.000_uM' ~ "DMSO 0.1% - Z-VAD-FMK 30.0 uM",
        Metadata_Treatment =='Flagellin_1.000_ug_per_ml_DMSO_0.025_%' ~ "Flagellin 1.0 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='Flagellin_1.000_ug_per_ml_Disulfiram_1.000_uM' ~ "Flagellin 1.0 ug/ml - Disulfiram 1.0 uM",
        Metadata_Treatment =='LPS_0.010_ug_per_ml_DMSO_0.025_%' ~ "LPS 0.01 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='LPS_0.100_ug_per_ml_DMSO_0.025_%' ~ "LPS 0.1 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='Flagellin_0.100_ug_per_ml_DMSO_0.0_%' ~ "Flagellin 0.1 ug/ml - DMSO 0.0%",
        Metadata_Treatment =='Flagellin_0.100_ug_per_ml_DMSO_0.025_%' ~ "Flagellin 0.1 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='Disulfiram_0.100_uM_DMSO_0.025_%' ~ "Disulfiram 0.1 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_1.000_ug_per_ml_1.000_uM_DMSO_0.025_%' ~ "LPS 1.0 ug/ml + Nigericin 1.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_1.000_ug_per_ml_10.000_uM_DMSO_0.025_%' ~ "LPS 1.0 ug/ml + Nigericin 10.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_1.000_ug_per_ml_10.000_uM_Disulfiram_1.000_uM' ~ "LPS 1.0 ug/ml + Nigericin 10.0 uM - Disulfiram 1.0 uM",
        Metadata_Treatment =='LPS_Nigericin_1.000_ug_per_ml_10.000_uM_Z-VAD-FMK_100.000_uM' ~ "LPS 1.0 ug/ml + Nigericin 10.0 uM - Z-VAD-FMK 100.0 uM",
        Metadata_Treatment =='LPS_Nigericin_1.000_ug_per_ml_3.000_uM_DMSO_0.025_%' ~ "LPS 1.0 ug/ml + Nigericin 3.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_1.000_ug_per_ml_DMSO_0.025_%' ~ "LPS 1.0 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='Flagellin_1.000_ug_per_ml_DMSO_0.0_%' ~ "Flagellin 1.0 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='Disulfiram_1.000_uM_DMSO_0.025_%' ~ "Disulfiram 1.0 uM - DMSO 0.025%",
        Metadata_Treatment =='Thapsigargin_1.000_uM_DMSO_0.025_%' ~ "Thapsigargin 1.0 uM - DMSO 0.025%",
        Metadata_Treatment =='Topotecan_10.000_nM_DMSO_0.025_%' ~ "Topotecan 10.0 nM - DMSO 0.025%",
        Metadata_Treatment =='LPS_10.000_ug_per_ml_DMSO_0.025_%' ~ "LPS 10.0 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='LPS_10.000_ug_per_ml_Disulfiram_0.100_uM' ~ "LPS 10.0 ug/ml - Disulfiram 0.1 uM",
        Metadata_Treatment =='LPS_10.000_ug_per_ml_Disulfiram_1.000_uM' ~ "LPS 10.0 ug/ml - Disulfiram 1.0 uM",
        Metadata_Treatment =='LPS_10.000_ug_per_ml_Disulfiram_2.500_uM' ~ "LPS 10.0 ug/ml - Disulfiram 2.5 uM",
        Metadata_Treatment =='LPS_10.000_ug_per_ml_Z-VAD-FMK_100.000_uM' ~ "LPS 10.0 ug/ml - Z-VAD-FMK 100.0 uM",
        Metadata_Treatment =='Thapsigargin_10.000_uM_DMSO_0.025_%' ~ "Thapsigargin 10.0 uM - DMSO 0.025%",
        Metadata_Treatment =='H2O2_100.000_nM_DMSO_0.025_%' ~ "H2O2 100.0 nM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_100.000_ug_per_ml_1.000_uM_DMSO_0.025_%' ~ "LPS 100.0 ug/ml + Nigericin 1.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_100.000_ug_per_ml_10.000_uM_DMSO_0.025_%' ~ "LPS 100.0 ug/ml + Nigericin 10.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_Nigericin_100.000_ug_per_ml_3.000_uM_DMSO_0.025_%' ~ "LPS 100.0 ug/ml + Nigericin 3.0 uM - DMSO 0.025%",
        Metadata_Treatment =='LPS_100.000_ug_per_ml_DMSO_0.025_%' ~ "LPS 100.0 ug/ml - DMSO 0.025%",
        Metadata_Treatment =='H2O2_100.000_uM_DMSO_0.025_%' ~ "H2O2 100.0 uM - DMSO 0.025%",
        Metadata_Treatment =='H2O2_100.000_uM_Disulfiram_1.000_uM' ~ "H2O2 100.0 uM - Disulfiram 1.0 uM",
        Metadata_Treatment =='H2O2_100.000_uM_Z-VAD-FMK_100.000_uM' ~ "H2O2 100.0 uM - Z-VAD-FMK 100.0 uM",
        Metadata_Treatment =='Disulfiram_2.500_uM_DMSO_0.025_%' ~ "Disulfiram 2.5 uM - DMSO 0.025%",
        Metadata_Treatment =='Topotecan_20.000_nM_DMSO_0.025_%' ~ "Topotecan 20.0 nM - DMSO 0.025%",
        Metadata_Treatment =='Topotecan_5.000_nM_DMSO_0.025_%' ~ "Topotecan 5.0 nM - DMSO 0.025%",
        Metadata_Treatment =='media_ctr_0.0_0_Media_ctr_0.0_0' ~ "Media ctr 0.0 0",
        Metadata_Treatment =='media_ctr_0.0_0_Media_0.0_0' ~ "Media ctr 0.0 0"
    ))
    # replace Media ctr 0.0 0 with Media
df$Metadata_Treatment <- gsub("Media ctr 0.0 0", "Media", df$Metadata_Treatment)

# split the Metadata_Treatment into two columns by the " - " delimiter
df <- df %>%
    separate(Metadata_Treatment, c("inducer", "inhibitor"), sep = " - ", remove = FALSE)

# replace the inhibitor NA with Media
df$inhibitor <- ifelse(is.na(df$inhibitor), "Media", df$inhibitor)

# make the group_treatment column a factor
df$inducer <- factor(
    df$inducer,
    levels = c(
        'Media',
        'DMSO 0.1%',

        'Flagellin 0.1 ug/ml',
        'Flagellin 1.0 ug/ml',

        'LPS 0.01 ug/ml',
        'LPS 0.1 ug/ml',
        'LPS 1.0 ug/ml',
        'LPS 10.0 ug/ml',
        'LPS 100.0 ug/ml',

        'LPS 1.0 ug/ml + Nigericin 1.0 uM',
        'LPS 1.0 ug/ml + Nigericin 3.0 uM',
        'LPS 1.0 ug/ml + Nigericin 10.0 uM',

        'LPS 100.0 ug/ml + Nigericin 1.0 uM',
        'LPS 100.0 ug/ml + Nigericin 3.0 uM',
        'LPS 100.0 ug/ml + Nigericin 10.0 uM',

        'H2O2 100.0 nM',
        'H2O2 100.0 uM',

        'Disulfiram 0.1 uM',
        'Disulfiram 1.0 uM',
        'Disulfiram 2.5 uM',

        'Thapsigargin 1.0 uM',
        'Thapsigargin 10.0 uM',

        'Topotecan 5.0 nM',
        'Topotecan 10.0 nM',
        'Topotecan 20.0 nM'
    )
)

# make the group_treatment column a factor
df$inhibitor <- factor(
    df$inhibitor,
    levels = c(
        'Media',
        'DMSO 0.025%',
        'DMSO 1.0%',

        'Disulfiram 0.1 uM',
        'Disulfiram 1.0 uM',
        'Disulfiram 2.5 uM',

        'Z-VAD-FMK 30.0 uM',
        'Z-VAD-FMK 100.0 uM'
    )
)
head(df)

# save the df to as parquet for plotting later on
arrow::write_parquet(df, file.path("..","data","processed","mAP_scores","morphology_secretome_comparison.parquet"))

width <- 15
height <- 15
options(repr.plot.width=width, repr.plot.height=height)
# scatter plot with fill being the treatment dose
scatter_by_treatment <- (
    ggplot(df, aes(x=mAP_moprhology, y=mAP_secretome, col = inducer, shape=inhibitor))
    + geom_point(size=3, alpha=1)
    + labs(x="Morphology mAP score", y="Secretome mAP score")
    + theme_bw()
    + ggtitle("Comparison of mAP scores")
    + ylim(0,1)
    + xlim(0,1)
    + figure_theme
    # Change the legend title
    # change the legend shape
    + scale_color_manual(
        name = "Inducer",
        labels = c(
            'Media',
            'DMSO 0.1%',

            'Flagellin 0.1 ug/ml',
            'Flagellin 1.0 ug/ml',

            'LPS 0.01 ug/ml',
            'LPS 0.1 ug/ml',
            'LPS 1.0 ug/ml',
            'LPS 10.0 ug/ml',
            'LPS 100.0 ug/ml',

            'LPS 1.0 ug/ml + Nigericin 1.0 uM',
            'LPS 1.0 ug/ml + Nigericin 3.0 uM',
            'LPS 1.0 ug/ml + Nigericin 10.0 uM',

            'LPS 100.0 ug/ml + Nigericin 1.0 uM',
            'LPS 100.0 ug/ml + Nigericin 3.0 uM',
            'LPS 100.0 ug/ml + Nigericin 10.0 uM',

            'H2O2 100.0 nM',
            'H2O2 100.0 uM',

            'Disulfiram 0.1 uM',
            'Disulfiram 1.0 uM',
            'Disulfiram 2.5 uM',

            'Thapsigargin 1.0 uM',
            'Thapsigargin 10.0 uM',

            'Topotecan 5.0 nM',
            'Topotecan 10.0 nM',
            'Topotecan 20.0 nM'
        ),
        values = colors)
    + scale_shape_manual(
        name = "Inhibitor",
        labels = c(
            'Media',
            'DMSO 0.025%',
            'DMSO 1.0%',

            'Disulfiram 0.1 uM',
            'Disulfiram 1.0 uM',
            'Disulfiram 2.5 uM',

            'Z-VAD-FMK 30.0 uM',
            'Z-VAD-FMK 100.0 uM'

        ),
        values = shapes
    )
    # make the legend 1 column
    + guides(color = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))
    + ggplot2::coord_fixed()
    + facet_grid(shuffled~.)
    # add y = x line
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
)
scatter_by_treatment
