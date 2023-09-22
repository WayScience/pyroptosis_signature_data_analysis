

suppressPackageStartupMessages(library(ggplot2))

platte <- (
    '#D81B60'
    '#1E88E5'
    '#FFC107'
    '#004D40'
    '#519C09'
    '#980CAB')

figure_theme <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        strip.text = element_text(size = 8),
    )
)

figure_theme_wide <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 16),
    )
)
