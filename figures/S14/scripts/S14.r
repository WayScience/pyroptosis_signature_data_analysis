suppressPackageStartupMessages(suppressWarnings(library(ggplot2))) # plotting
suppressPackageStartupMessages(suppressWarnings(library(dplyr))) # data manipulation
suppressPackageStartupMessages(suppressWarnings(library(argparser))) # command line arguments
suppressPackageStartupMessages(suppressWarnings(library(patchwork))) # plot patchwork
suppressPackageStartupMessages(suppressWarnings(library(reshape2))) # data manipulation
suppressPackageStartupMessages(suppressWarnings(library(ggridges))) # ridgeline plots
suppressPackageStartupMessages(suppressWarnings(library(RColorBrewer))) # color palettes
suppressPackageStartupMessages(suppressWarnings(library(cowplot))) # ggplot2 drawing
suppressPackageStartupMessages(suppressWarnings(library(ggplotify))) # ggplot2 drawing

source("../../utils/figure_themes.r")


cell_type <- "PBMC"
model_name <- "MultiClass_MLP"

# set file path for importing the data
training_metrics_file <- file.path(paste0(
    "../../../4.sc_Morphology_Neural_Network_MLP_Model/results/Multi_Class/",model_name,"/",cell_type,"/training_metrics.parquet"
))

# set output file path for graphs
f1_plot_path <- file.path(paste0(
    "../figures/f1_score.png"
))

# make the output directory if it doesn't exist
dir.create(file.path(paste0(
    "../figures/"
)), showWarnings = FALSE, recursive = TRUE)



# read in the data
training_metrics <- arrow::read_parquet(training_metrics_file)


support <- training_metrics[training_metrics$metric == "support",]
# get apoptosis, healthy, and pyroptosis support rows in one df
support <- support[support$label %in% c("apoptosis", "healthy", "pyroptosis"),]


# get the rows that contain the F1 scores
f1_scores <- training_metrics[training_metrics$metric == "f1-score",]
# remove the rows that contain the macro and weighted averages
f1_scores <- f1_scores[!grepl("macro avg", f1_scores$label),]
f1_scores <- f1_scores[!grepl("weighted avg", f1_scores$label),]
# muatate the label column for multiple cases
f1_scores$label <- gsub("healthy", "Control", f1_scores$label)
f1_scores$label <- gsub("apoptosis", "Apoptosis", f1_scores$label)
f1_scores$label <- gsub("pyroptosis", "Pyroptosis", f1_scores$label)
# mutate the data type column
f1_scores$group <- gsub("train", "Training", f1_scores$group)
f1_scores$group <- gsub("test", "Testing", f1_scores$group)
f1_scores$group <- gsub("validation", "Validation", f1_scores$group)
f1_scores$group <- gsub("treatment_holdout", "Treatment holdout", f1_scores$group)
# f1_scores$group <- gsub("holdout", "Holdout", f1_scores$group)
# factorize the group column
f1_scores$group <- factor(f1_scores$group, levels = c(
    "Training", "Validation", "Testing","Treatment holdout", "Holdout"
))
# mutate the shuffled_data column
f1_scores$shuffle <- gsub("TRUE", "Shuffled", f1_scores$shuffle)
f1_scores$shuffle <- gsub("FALSE", "Not Shuffled", f1_scores$shuffle)
# cbind the support column to the f1_scores df
f1_scores <- cbind(f1_scores, support$value)
# rename the support column
colnames(f1_scores)[colnames(f1_scores) == "support$value"] <- "support"
# dived the support by 10,000 to get the number of cells
f1_scores$support <- f1_scores$support / 10000
# round the support column to 2 decimal places
f1_scores$support <- round(f1_scores$support, 2)


# make the label a factor so that the order is preserved
f1_scores$label <- factor(
    f1_scores$label, levels = c(
        "Control", "Apoptosis", "Pyroptosis"
        )
    )


head(f1_scores, 1)

# remove the treatment holdout rows
f1_scores <- f1_scores[grepl("Treatment holdout", f1_scores$group),]
f1_scores <- f1_scores[grepl("Pyroptosis", f1_scores$label),]
head(f1_scores, 1)

# set plot size
width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
# bar plot of the F1 scores
f1_score_plot <- (
    ggplot(f1_scores, aes(x = shuffle, y = value, fill = group))
    + geom_bar(stat = "identity", position = "dodge")

    + ylim(0, 1)
    + ylab("F1 score")
    + xlab("Data split")
    # change the legend title
    + labs(fill = "Predicted class")
    # change the colours
    + scale_fill_manual(values = c(
        "Treatment holdout" = brewer.pal(3, "Dark2")[3]
    ))
    # remove legend

    + figure_theme_wide
    + theme(legend.position = "none")

)
ggsave(f1_plot_path, f1_score_plot, width = width, height = height, dpi = 600)
f1_score_plot


# load in the probabilities
treatment_holdout_probabilities_path <- file.path(
    paste0(
        "../../../4.sc_Morphology_Neural_Network_MLP_Model/results/Multi_Class/",model_name,"/",cell_type,"/probabilities.parquet"
    )
)
# read in the data from the parquet file
probabilities <- arrow::read_parquet(
    treatment_holdout_probabilities_path
)
head(probabilities,2)

unique(probabilities$data_split)
unique(probabilities$shuffle)

# replace label_true value 1 with Control
probabilities$label_true <- gsub("1", "Control", probabilities$label_true)
# replace label_true value 2 with pyroptosis
probabilities$label_true <- gsub("2", "Pyroptosis", probabilities$label_true)
# replace label_true value 0 with apoptosis
probabilities$label_true <- gsub("0", "Apoptosis", probabilities$label_true)

# replace label_pred value 1 with Control
probabilities$label_pred <- gsub("1", "Control", probabilities$label_pred)
# replace label_pred value 2 with pyroptosis
probabilities$label_pred <- gsub("2", "Pyroptosis", probabilities$label_pred)
# replace label_pred value 0 with apoptosis
probabilities$label_pred <- gsub("0", "Apoptosis", probabilities$label_pred)

# replace shuffled value TRUE with Shuffled
probabilities$shuffle <- gsub("TRUE", "Shuffled", probabilities$shuffle)
# replace shuffled value FALSE with Not Shuffled
probabilities$shuffle <- gsub("FALSE", "Not Shuffled", probabilities$shuffle)

# replace data_split value treatment_holdout with Treatment Holdout|
probabilities$data_split <- gsub("treatment_holdout", "Treatment holdout", probabilities$data_split)


head(probabilities, 2)
unique(probabilities$shuffle)
unique(probabilities$data_split)


# change the label columns to be factors
probabilities$label_true <- factor(probabilities$label_true , levels = c(
    "Control", "Apoptosis", "Pyroptosis"
))
probabilities$label_pred <- factor(probabilities$label_pred , levels = c(
    "Pyroptosis", "Apoptosis", "Control"
))

# change the shuffled_data column to be a factor
probabilities$shuffle <- factor(probabilities$shuffle, levels = c(
    "Not Shuffled", "Shuffled"
))

# remove treatment holdout rows
probabilities <- probabilities[grepl("Treatment holdout", probabilities$data_split),]

height <- 5
width <- 15
options(repr.plot.width = width, repr.plot.height = height)

ridge_plot_pyroptosis <- (
    ggplot(probabilities, aes(x = pyroptosis_prob, y = label_pred, fill = label_true, group = label_pred))
    + geom_density_ridges(
        aes(fill = label_pred), alpha = 0.7, scale = 3, rel_min_height = 0.01
    )
    + facet_grid(.~data_split)
    + scale_fill_manual(values = c(
        "Control" = brewer.pal(3, "Dark2")[2],
        "Apoptosis" = brewer.pal(3, "Dark2")[1],
        "Pyroptosis" = brewer.pal(3, "Dark2")[3]
    ))
    + geom_vline(xintercept = 1, linetype = "dashed", color = "black")
    + facet_grid(shuffle ~ ., scales = "free_y")
    + scale_x_continuous(breaks = seq(0, 1, 0.5))
    + labs(title = "Pyroptosis prediction probability", y = "Predicted class",fill = "True class")
    + labs()
    + theme_bw()
    + figure_theme

    # remove legend
    # make title larger
    + theme(plot.title = element_text(size = 20, hjust = 0.5))
    + theme(legend.position = "bottom", legend.direction = "horizontal")
    # remove x axis label
    + theme(axis.title.x = element_blank())
    # add vertical line at 1
)
ridge_plot_pyroptosis


# get the legend
ridge_plot_pyroptosis <- (
    # remove the legend
    ridge_plot_pyroptosis + theme(legend.position = "none")
)

# ridge_plot_pyroptosis
# patch the plots together via the patchwork package
layout <- c(
    area(t=1, b=2, l=1, r=1), # A
    area(t=3, b=4, l=1, r=1), # B
    area(t=5, b=6, l=1, r=1)  # C
)
# set plot size
width <- 17
height <- 17
options(repr.plot.width=width, repr.plot.height=height, units = "cm", dpi = 600)
fig5_probabilities <- (
    ridge_plot_pyroptosis
    + plot_layout(design = layout)
)
fig5_probabilities

sc_preds_path <- file.path(
    paste0(
        "../../../8.cytopick_analysis/results/PBMC/single_cell_predictions.parquet"
    )
)
# read in the data from the parquet file
sc_preds <- arrow::read_parquet(
    sc_preds_path
)
head(sc_preds,2)


# define df subsets for each class, data split
pyroptosis_correct_treatment_holdout <- sc_preds[
    sc_preds$labels == "pyroptosis" &
    sc_preds$correct == TRUE &
    sc_preds$shuffle == FALSE &
    sc_preds$data_split == "treatment_holdout",]


# make a list of the data frames for each class
correct_class_dfs <- list(
    pyroptosis_correct_treatment_holdout
)



width <- 2
height <- 2
options(repr.plot.width = width, repr.plot.height = height)
# define function to return the image object
get_image <- function(df, i){
    # Load the PNG file
    img <- png::readPNG(df$image_crop_path[i])
    # Convert the image to a raster object
    g <- grid::rasterGrob(img, interpolate=TRUE)

    # Create a ggplot
    p <- ggplot() +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_void()  # Remove axes and labels

    # Print the plot
    return(p)
}


pyroptosis_correct_treatment_holdout_image_1 <- get_image(pyroptosis_correct_treatment_holdout, 1)
# add the title
pyroptosis_correct_treatment_holdout_image_1 <- (
    pyroptosis_correct_treatment_holdout_image_1
    + theme(plot.title = element_text(size = 18, hjust = 0.5))
    + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
)
pyroptosis_correct_treatment_holdout_image_2 <- get_image(pyroptosis_correct_treatment_holdout, 2)
# add the title
pyroptosis_correct_treatment_holdout_image_2 <- (
    pyroptosis_correct_treatment_holdout_image_2
    + theme(plot.title = element_text(size = 18, hjust = 0.5))
    + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
)

pyroptosis_correct_treatment_holdout_image_3 <- get_image(pyroptosis_correct_treatment_holdout, 3)
# add the title
pyroptosis_correct_treatment_holdout_image_3 <- (
    pyroptosis_correct_treatment_holdout_image_3
    + theme(plot.title = element_text(size = 18, hjust = 0.5))
    + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
)

pyroptosis_correct_treatment_holdout_image_4 <- get_image(pyroptosis_correct_treatment_holdout, 4)
# add the title
pyroptosis_correct_treatment_holdout_image_4 <- (
    pyroptosis_correct_treatment_holdout_image_4
    + theme(plot.title = element_text(size = 18, hjust = 0.5))
    + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
)

pyroptosis_correct_treatment_holdout_image_5 <- get_image(pyroptosis_correct_treatment_holdout, 5)
# add the title
pyroptosis_correct_treatment_holdout_image_5 <- (
    pyroptosis_correct_treatment_holdout_image_5
    + theme(plot.title = element_text(size = 18, hjust = 0.5))
    + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
)


pyroptosis_correct_treatment_holdout_image_1
pyroptosis_correct_treatment_holdout_image_2
pyroptosis_correct_treatment_holdout_image_3
pyroptosis_correct_treatment_holdout_image_4
pyroptosis_correct_treatment_holdout_image_5

width <- 15
height <- 5
options(repr.plot.width = width, repr.plot.height = height)

layout <- c(
    area(t=1, b=2, l=1, r=2), # A
    area(t=1, b=2, l=3, r=4), # B
    area(t=1, b=2, l=5, r=6), # C
    area(t=1, b=2, l=7, r=8), # D
    area(t=1, b=2, l=9, r=10)  # E

)

pyroptosis_correct_treatment_holdout_image <- (

        # wrap_elements(full = ridge_plot_pyroptosis)
        pyroptosis_correct_treatment_holdout_image_1
        + pyroptosis_correct_treatment_holdout_image_2
        + pyroptosis_correct_treatment_holdout_image_3
        + pyroptosis_correct_treatment_holdout_image_4
        + pyroptosis_correct_treatment_holdout_image_5

    + plot_layout(design = layout)
)

pyroptosis_correct_treatment_holdout_image

layout <- c(
    area(t=1, b=2, l=1, r=2), # A
    area(t=1, b=2, l=3, r=4), # B
    area(t=3, b=4, l=1, r=4)  # C

)
# set plot size
width <- 17
height <- 10
options(repr.plot.width=width, repr.plot.height=height, units = "cm", dpi = 600)
figs10 <- (
    f1_score_plot
    + wrap_elements(full = ridge_plot_pyroptosis)
    + wrap_elements(full = pyroptosis_correct_treatment_holdout_image)
    + plot_layout(design = layout, widths = c(10, 10))
    # make bottom plot not align
    + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
)
figs10

# save the plot

ggsave("../figures/S14.png", figs10, width = width, height = height, dpi = 600)

