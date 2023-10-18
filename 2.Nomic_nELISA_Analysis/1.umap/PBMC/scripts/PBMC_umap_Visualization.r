library(ggplot2)
library(gridExtra)
library(cowplot)
library(styler)

cell_type <- "PBMC"

# import paths

nomic_data_path <- file.path("./results/nELISA_plate_430420_umap_PBMC.csv")
nomic_data_selected_treatments_path <- file.path("./results/nELISA_plate_430420_umap_PBMC_selected_treatments.csv")

figure_output_path <- file.path("./figures/")

# read data
nomic_data <- read.csv(nomic_data_path)
nomic_data_selected_treatments <- read.csv(nomic_data_selected_treatments_path)

head(nomic_data,2)
head(nomic_data_selected_treatments,2)


## plot the data
# set the plot size
# options(repr.plot.width=10, repr.plot.height=10)
nomic_data_plot <- (
    ggplot(nomic_data, aes(x = umap_1, y = umap_2, color = oneb_Treatment_Dose_Inhibitor_Dose))
    + geom_point(size = 5)
    # detach the legend and plot it separately
    + theme_bw()
    + theme(
        legend.position = "right",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 8),
        strip.background = element_rect(
            colour = "black",
            fill = "#fdfff4"
            )
            )
    + guides(color = guide_legend(ncol = 2,title.position = 'top'))
    + labs(title = "PBMC UMAP by treatment",color="Treatment", x = "UMAP1", y = "UMAP2")
    # plot the legend separately


)
# extract the legend
legend <- get_legend(nomic_data_plot)

# remove the legend from the plot # will be saved separately
nomic_data_plot <- (
    nomic_data_plot
    + theme(legend.position = "none")
)
# show the plot + legend
legend <- plot_grid(legend)
nomic_data_plot
legend

# save the plot
# set the path
plot_path <- file.path(paste0(figure_output_path,"umap_by_all_treatment.png"))
plot_legend_path <- file.path(paste0(figure_output_path,"umap_by_all_treatment_legend.png"))
# save the plot
ggsave(plot_path, nomic_data_plot, width = 8.5, height = 5.5, dpi = 500)
ggsave(plot_legend_path, legend,dpi = 500, scale = 2)

## plot the data
# set the plot size
# options(repr.plot.width=10, repr.plot.height=10)
nomic_data_selected_treatments_plot <- (
    ggplot(nomic_data_selected_treatments, aes(x = umap_1, y = umap_2, color = oneb_Treatment_Dose_Inhibitor_Dose))
    + geom_point(size = 5)
    # detach the legend and plot it separately # save separately too
    + theme_bw()
    + theme(
        legend.position = "right",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 8),
        strip.background = element_rect(
            colour = "black",
            fill = "#fdfff4"
            )
            )
    + guides(color = guide_legend(ncol = 2,title.position = 'top'))
    + labs(title = "PBMC UMAP by selected treatments",color="Treatment", x = "UMAP1", y = "UMAP2")
    # plot the legend separately


)
# extract the legend
nomic_data_selected_treatments_plot_legend <- get_legend(nomic_data_selected_treatments_plot)

# remove the legend from the plot
nomic_data_selected_treatments_plot <- (
    nomic_data_selected_treatments_plot
    + theme(legend.position = "none")
)
# show the plot + legend
nomic_data_selected_treatments_plot_legend <- plot_grid(nomic_data_selected_treatments_plot_legend)
nomic_data_selected_treatments_plot
nomic_data_selected_treatments_plot_legend

# save the plot
# set the path
plot_path <- file.path(paste0(figure_output_path,"umap_by_selected_treatment.png"))
plot_legend_path <- file.path(paste0(figure_output_path,"umap_by_selected_treatment_legend.png"))
# save the plot
ggsave(plot_path, nomic_data_selected_treatments_plot, width = 8.5, height = 5.5, dpi = 500)
ggsave(plot_legend_path, nomic_data_selected_treatments_plot_legend, dpi = 500, scale = 2)
