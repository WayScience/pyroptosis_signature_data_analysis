suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(patchwork)))
suppressPackageStartupMessages(suppressWarnings(library(cowplot)))
suppressPackageStartupMessages(suppressWarnings(library(RcppTOML)))
suppressPackageStartupMessages(suppressWarnings(library(pheatmap)))
suppressPackageStartupMessages(suppressWarnings(library(lattice)))
suppressPackageStartupMessages(suppressWarnings(library("RColorBrewer")))
suppressPackageStartupMessages(suppressWarnings(library(gplots)))
suppressPackageStartupMessages(suppressWarnings(library(ComplexHeatmap)))
suppressPackageStartupMessages(suppressWarnings(library(ggplotify)))
# import the R them
figure_themes <- source(file.path("../figure_themes.r"))

# Load in the treatment list
toml_path <- file.path("..","..","1.Exploratory_Data_Analysis/utils/params.toml")
p <- parseTOML(toml_path)
list_of_treatments <- c(p$list_of_treatments$treatments)
list_of_treatments

# Figure 2A

# Load data
data_path_cytokine_values <- file.path("..","..","2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/PBMC_all_cytokine_values_per_treatment_per_well.csv")
cytokine_values <- read.csv(data_path_cytokine_values, header = TRUE, sep = ",")
# filter out the treatments that are not in the list
cytokine_values <- cytokine_values[cytokine_values$oneb_Metadata_Treatment_Dose_Inhibitor_Dose %in% list_of_treatments,]

list_of_treatments <- c(
    'DMSO_0.100_DMSO_0.025',
    'Thapsigargin_1.000_DMSO_0.025',
    'Thapsigargin_10.000_DMSO_0.025',
    'H2O2_100.000_DMSO_0.025',

    'LPS_Nigericin_1.000_1.0_DMSO_0.025',
    'LPS_Nigericin_1.000_3.0_DMSO_0.025',
    'LPS_Nigericin_1.000_10.0_DMSO_0.025',

    'Flagellin_0.100_DMSO_0.025',
    'Flagellin_1.000_DMSO_0.025',
    'LPS_0.010_DMSO_0.025',
    'LPS_1.000_DMSO_0.025',
    'LPS_100.000_DMSO_0.025')


# order the treatments
# cytokine_values$oneb_Metadata_Treatment_Dose_Inhibitor_Dose <- factor(cytokine_values$oneb_Metadata_Treatment_Dose_Inhibitor_Dose, levels = list_of_treatments)

# set plot size
options(repr.plot.width=15, repr.plot.height=5)
# Plot
cytokine_scatter_plot <- (
    ggplot(
        data = cytokine_values,
        aes(
            x = IL.1.beta..NSU.,
            y = TNF.alpha..NSU.,
            color = oneb_Metadata_Treatment_Dose_Inhibitor_Dose
        )
    )
    + geom_point()
    + xlab("IL-1 beta Abundance")
    + ylab("TNF alpha Abundance")
    + # rename legend title
    labs(color = "Treatment")
)

cytokine_scatter_plot

# set plot size
options(repr.plot.width=15, repr.plot.height=5)
# Plot
cytokine_scatter_plot1 <- (
    ggplot(
        data = cytokine_values,
        aes(
            x = IL.1.beta..NSU.,
            y = CCL24..NSU.,
            color = oneb_Metadata_Treatment_Dose_Inhibitor_Dose
        )
    )
    + geom_point()
    + xlab("IL-1 beta Abundance")
    + ylab("TNF alpha Abundance")
)

cytokine_scatter_plot1

# import melted dataframes
# Figure 2A

# Load data
data_path_cytokine_values_melted <- file.path("..","..","2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/PBMC_all_cytokine_values_per_treatment_per_well_melted.csv")
cytokine_values_melted <- read.csv(data_path_cytokine_values_melted, header = TRUE, sep = ",")
# filter out the treatments that are not in the list
cytokine_values_melted <- cytokine_values_melted[cytokine_values_melted$oneb_Metadata_Treatment_Dose_Inhibitor_Dose %in% list_of_treatments,]

head(cytokine_values_melted)

# select a few cytokines to visualize
cytokine_values_melted <- cytokine_values_melted[cytokine_values_melted$Cytokine %in% c(
    "IL-1 beta [NSU]",
    "TNF alpha [NSU]",
    "CCL24 [NSU]",
    "IL-18 [NSU]",
    "IL-6 [NSU]",
    "Osteopontin (OPN) [NSU]",
    "CCL13 [NSU]",
    "IL-2 [NSU]"
    ),]
head(cytokine_values_melted)

# set the order of the cytokines
cytokine_values_melted$Cytokine <- factor(cytokine_values_melted$Cytokine, levels = c(
    "IL-1 beta [NSU]",
    "TNF alpha [NSU]",

    "IL-18 [NSU]",
    "IL-6 [NSU]",


    "IL-2 [NSU]",
    "Osteopontin (OPN) [NSU]",
    "CCL13 [NSU]",
    "CCL24 [NSU]"
    ))

head(cytokine_values_melted)

# aggregate the data to get both the mean and the standard deviation
cytokine_values_melted_agg <- cytokine_values_melted %>%
    group_by(Cytokine, oneb_Metadata_Treatment_Dose_Inhibitor_Dose) %>%
    summarise(
        mean = mean(Cytokine_Value),
        sd = sd(Cytokine_Value)
    )
head(cytokine_values_melted_agg)

# set plot size
options(repr.plot.width=10, repr.plot.height=5)
# make the bar plots for the cytokine values
cytokine_bar_plot <- (
    ggplot(data=cytokine_values_melted_agg,
           aes(x=oneb_Metadata_Treatment_Dose_Inhibitor_Dose,
                y=mean))
    + geom_bar(stat="identity", position="dodge", color='black', fill='#980CAB')
    + theme_bw()
    + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))
    + facet_wrap(.~Cytokine, nrow =1)
    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    + ylab("Cytokine Abundacne")
    + xlab("Treatment")

)

cytokine_bar_plot


# Load data
data_path_cytokine_values_melted <- file.path("..","..","2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/PBMC_all_cytokine_values_per_treatment_per_well_melted.csv")
cytokine_values_melted <- read.csv(data_path_cytokine_values_melted, header = TRUE, sep = ",")
# filter out the treatments that are not in the list
cytokine_values_melted <- cytokine_values_melted[cytokine_values_melted$oneb_Metadata_Treatment_Dose_Inhibitor_Dose %in% list_of_treatments,]
# heatmap of the cytokine values
# get the aggregated values for the cytokine values across the treatments
cytokine_values_melted_agg <- cytokine_values_melted %>%
    group_by(oneb_Metadata_Treatment_Dose_Inhibitor_Dose, Cytokine) %>%
    summarise(Cytokine_Value = mean(Cytokine_Value))



# drop two columns by name
cytokine_values_agg <- subset(cytokine_values, select = -c(fourb_Metadata_Treatment_Dose_Inhibitor_Dose, Metadata_position_x))
# aggregate the cytokine values across the treatments and cytokine
cytokine_values_agg <- cytokine_values_agg %>%
    group_by(oneb_Metadata_Treatment_Dose_Inhibitor_Dose) %>%
    summarise_all(mean)
# cytokine_values_agg <- as.matrix(cytokine_values_agg)

head(cytokine_values_melted_agg)

# create a matrix of the cytokine values for the clustering and heatmap
# unmelt the data
cytokine_values_agg <- reshape2::dcast(cytokine_values_melted_agg, oneb_Metadata_Treatment_Dose_Inhibitor_Dose ~ Cytokine, value.var = "Cytokine_Value")# make oneb_Metadata_Treatment_Dose_Inhibitor_Dose the rownames
rownames(cytokine_values_agg) <- cytokine_values_agg$oneb_Metadata_Treatment_Dose_Inhibitor_Dose
# drop the column
cytokine_values_agg <- subset(cytokine_values_agg, select = -c(oneb_Metadata_Treatment_Dose_Inhibitor_Dose))
# convert to matrix
# cytokine_values_agg <- as.matrix(cytokine_values_agg)
#

cytokine_values_agg

row_dend <- as.dendrogram(hclust(dist(cytokine_values_agg)))
col_dend <- as.dendrogram(hclust(dist(t(cytokine_values_agg))))

# remove '[NSU]' from the column names
colnames(cytokine_values_agg) <- gsub("\\[NSU\\]", "", colnames(cytokine_values_agg))


cytokine_values_agg <- as.matrix(cytokine_values_agg)

# set plot size
options(repr.plot.width=25, repr.plot.height=5)
heatmap_plot_all <- (
  Heatmap(
  (cytokine_values_agg),
  col = brewer.pal(9, "GnBu"),
  cluster_rows = TRUE,    # Cluster rows
  cluster_columns = TRUE, # Cluster columns
  show_row_names = TRUE,  # Show row names
  show_column_names = TRUE, # Show column names
  column_names_gp = gpar(fontsize = 10), # Column name label formatting
  row_names_gp = gpar(fontsize = 12),    # Row name label formatting
  heatmap_legend_param = list(title = "Legend", at = c(0, 1)),
  # make the tiles rectangular
  rect_gp = gpar(col = NA),
  )
)
heatmap_plot_all

# import the tukey test results
tukey_results_path <- file.path("..","..","2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/tukey_test_results.csv")
# read in the data
tukey_results <- read.csv(tukey_results_path, header = TRUE, sep = ",")

# reload in the cytokine values
data_path_cytokine_values_melted <- file.path("..","..","2.Nomic_nELISA_Analysis/0.Exploratory_Analysis/PBMC/results/PBMC_all_cytokine_values_per_treatment_per_well_melted.csv")
cytokine_values_melted <- read.csv(data_path_cytokine_values_melted, header = TRUE, sep = ",")
# filter out the treatments that are not in the list
cytokine_values_melted <- cytokine_values_melted[cytokine_values_melted$oneb_Metadata_Treatment_Dose_Inhibitor_Dose %in% list_of_treatments,]

# aggregate the cytokine values across the treatments and cytokine
cytokine_values_melted_agg <- cytokine_values_melted %>%
    group_by(oneb_Metadata_Treatment_Dose_Inhibitor_Dose, Cytokine) %>%
    summarise(Cytokine_Value = mean(Cytokine_Value))

cytokine_values_melted_agg_filtered <- cytokine_values_melted_agg[cytokine_values_melted_agg$Cytokine %in% unique(tukey_results$cytokine),]
head(cytokine_values_melted_agg_filtered)

# un melt the data cytokine_values_melted_agg_filtered
cytokine_values_melted_agg_filtered <- reshape2::dcast(cytokine_values_melted_agg_filtered, oneb_Metadata_Treatment_Dose_Inhibitor_Dose ~ Cytokine, value.var = "Cytokine_Value")

# make oneb_Metadata_Treatment_Dose_Inhibitor_Dose the index
row.names(cytokine_values_melted_agg_filtered) <- cytokine_values_melted_agg_filtered$oneb_Metadata_Treatment_Dose_Inhibitor_Dose
# drop the column
cytokine_values_melted_agg_filtered <- subset(cytokine_values_melted_agg_filtered, select = -c(oneb_Metadata_Treatment_Dose_Inhibitor_Dose))

cytokine_values_melted_agg_filtered <- as.matrix(cytokine_values_melted_agg_filtered)

# remove '[NSU]' from the column names
colnames(cytokine_values_melted_agg_filtered) <- gsub("\\[NSU\\]", "", colnames(cytokine_values_melted_agg_filtered))

col_func <- colorRampPalette(brewer.pal(9, "Purples"))
# set plot size
options(repr.plot.width=25, repr.plot.height=5)
heatmap_anova_cytokines <- (
    Heatmap(
        (cytokine_values_melted_agg_filtered),
        col = brewer.pal(9, "GnBu"),
        cluster_rows = TRUE,    # Cluster rows
        cluster_columns = TRUE, # Cluster columns
        show_row_names = TRUE,  # Show row names
        show_column_names = TRUE, # Show column names
        column_names_gp = gpar(fontsize = 12), # Column name label formatting
        row_names_gp = gpar(fontsize = 12),    # Row name label formatting
        heatmap_legend_param = list(title = "Abundance", at = c(0, 1)),
        # make the tiles rectangular
        rect_gp = gpar(col = NA)
    )
)

heatmap_anova_cytokines

# read in the umap results
umap_results_path <- file.path("..","..","2.Nomic_nELISA_Analysis/1.umap/PBMC/results/nELISA_plate_430420_umap_PBMC.csv")

umap_results_selected_treatments_path <- file.path("..","..","2.Nomic_nELISA_Analysis/1.umap/PBMC/results/nELISA_plate_430420_umap_PBMC_selected_treatments.csv")
# read in the data
umap_results <- read.csv(umap_results_path, header = TRUE, sep = ",")
umap_results_selected_treatments <- read.csv(umap_results_selected_treatments_path, header = TRUE, sep = ",")


umap_results_selected_treatments$oneb_Metadata_Treatment_Dose_Inhibitor_Dose <- factor(umap_results_selected_treatments$oneb_Metadata_Treatment_Dose_Inhibitor_Dose, levels = list_of_treatments)

# set the plot size
options(repr.plot.width=15, repr.plot.height=5)
# plot the umap results
umap_plot_all_treatments <- (
    ggplot(
        data = umap_results,
        aes(
            x = umap_1,
            y = umap_2,
            color = oneb_Metadata_Treatment_Dose_Inhibitor_Dose
        )
    )
    + geom_point()
)
umap_plot_all_treatments

# plot the umap results
umap_plot_selected_treatments <- (
    ggplot(
        data = umap_results_selected_treatments,
        aes(
            x = umap_1,
            y = umap_2,
            color = oneb_Metadata_Treatment_Dose_Inhibitor_Dose
        )
    )
    + geom_point()
)
umap_plot_selected_treatments

# show the plots first as assigned names
# set plot size
options(repr.plot.width=5, repr.plot.height=5)

# cytokine_scatter_plot

# cytokine_bar_plot

# heatmap_plot_all
# heatmap_anova_cytokines

# umap_plot_all_treatments
# umap_plot_selected_treatments

new_heatmap <- as.ggplot(heatmap_anova_cytokines)
cytokine_scatter_plot_legend <- get_legend(cytokine_scatter_plot)
umap_plot_selected_treatments_legend <- get_legend(umap_plot_selected_treatments)


# remove legends
umap_plot_selected_treatments <- umap_plot_selected_treatments + theme(legend.position = "none")
cytokine_scatter_plot <- cytokine_scatter_plot + theme(legend.position = "none")
cytokine_scatter_plot1 <- cytokine_scatter_plot1 + theme(legend.position = "none")



design <- "AB#CCD#
            AB#CCD#
            EEEEEE#
            FFFFFF#
            FFFFFF#"

# set plot size
options(repr.plot.width=25, repr.plot.height=20)
fig2 <- (
    cytokine_scatter_plot
    + cytokine_scatter_plot1
    + umap_plot_selected_treatments
    + cytokine_scatter_plot_legend
    + cytokine_bar_plot
    + new_heatmap
    + plot_layout(design = design)
)
fig2
ggsave(
    filename = file.path("figure2.png"),
    plot = fig2,
    width = 25,
    height = 20,
    units = "in",
    dpi = 600
)






