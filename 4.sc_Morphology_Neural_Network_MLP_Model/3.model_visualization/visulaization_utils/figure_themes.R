# Set consistent figure themes and colors
suppressPackageStartupMessages(library(ggplot2))

data_split_colors <- c(
    "train" = "#1b9e77",
    "validation" = "#d95f02",
    "test" = "#7570b3",
    "holdout" = "#044389"
)

data_split_labels <- c(
    "train" = "Training",
    "validation" = "Validation",
    "test" = "Testing",
    "holdout" = "Holdout"
)

shuffled_linetypes <- c(
    "False" = "solid",
    "True" = "dashed"
)

shuffled_labels <- c(
    "False" = "False",
    "True" = "True"
)

figure_theme <- (
    theme_bw()
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
)
