suppressPackageStartupMessages(suppressWarnings(library(ggplot2))) # plotting
suppressPackageStartupMessages(suppressWarnings(library(dplyr))) # data manipulation
suppressPackageStartupMessages(suppressWarnings(library(argparser))) # command line arguments
source("../../../figures/utils/figure_themes.r")


cell_type <- "SHSY5Y"
model_name <- "MLP_subset"


# set file path for importing the data
training_metrics_file <- file.path(paste0(
    "../../results/Multi_Class/",model_name,"/",cell_type,"/training_metrics.parquet"
))
confusion_matrix_file <- file.path(paste0(
    "../../results/Multi_Class/",model_name,"/",cell_type,"/confusion_matrices.parquet"
))

# set output file path for graphs
f1_plot_path <- file.path(paste0(
    "../../figures/Multi_Class/",model_name,"/",cell_type,"/f1_score.png"
))

confusion_matrix_plot_path <- file.path(paste0(
    "../../figures/Multi_Class/",model_name,"/",cell_type,"/confusion_matrix.png"
))
# set the path to the results
pr_curves_path <- file.path(paste0(
        "../../results/Multi_Class/",model_name,"/",cell_type,"/PR_curves.parquet"
))


# read in the data
training_metrics <- arrow::read_parquet(training_metrics_file)
confusion_matrix <- arrow::read_parquet(confusion_matrix_file)
PR_curves <- arrow::read_parquet(pr_curves_path)


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
f1_scores$group <- gsub("treatment_holdout", "Treatment Holdout", f1_scores$group)
f1_scores$group <- gsub("holdout", "Holdout", f1_scores$group)
# factorize the group column
f1_scores$group <- factor(f1_scores$group, levels = c(
    "Training", "Validation", "Testing","Treatment Holdout", "Holdout"
))
# mutate the shuffled_data column
f1_scores$shuffled_data <- gsub("TRUE", "Shuffled", f1_scores$shuffled_data)
f1_scores$shuffled_data <- gsub("FALSE", "Not Shuffled", f1_scores$shuffled_data)
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


# set plot size
width <- 10
height <- 5
options(repr.plot.width = width, repr.plot.height = height)
# bar plot of the F1 scores
f1_score_plot <- (
    ggplot(f1_scores, aes(x = shuffled_data, y = value, fill = group))
    + geom_bar(stat = "identity", position = "dodge")

    + ylim(0, 1)
    + facet_wrap(~label)
    + ylab("F1 Score")
    + xlab("Data Split")
    # change the legend title
    + labs(fill = "Predicted Class")
    # change the colours
    + scale_fill_manual(values = c(
        "Training" = "#88F2F2",
        "Validation" = "#056CF2",
        "Testing" = "#A6382E",
        "Treatment Holdout" = "#04B404",
        "Holdout" = "#F2A900"
    ))
    + figure_theme_wide

)
ggsave(f1_plot_path, f1_score_plot, width = width, height = height, dpi = 600)Healthy
f1_score_plot


confusion_matrix



# round the Recall vlaues to 2 decimal places
confusion_matrix$Recall <- round(confusion_matrix$Recall, 2)
# mutate the label column for multiple cases
confusion_matrix$True_Label <- gsub("healthy", "Control", confusion_matrix$True_Label)
confusion_matrix$True_Label <- gsub("apoptosis", "Apoptosis", confusion_matrix$True_Label)Healthy
confusion_matrix$True_Label <- gsub("pyroptosis", "Pyroptosis", confusion_matrix$True_Label)
confusion_matrix$Predicted_Label <- gsub("healthy", "Control", confusion_matrix$Predicted_Label)
confusion_matrix$Predicted_Label <- gsub("apoptosis", "Apoptosis", confusion_matrix$Predicted_Label)
confusion_matrix$Predicted_Label <- gsub("pyroptosis", "Pyroptosis", confusion_matrix$Predicted_Label)

# make the True Label and Predicted Label columns factors
confusion_matrix$True_Label <- factor(
    confusion_matrix$True_Label, levels = c(
        "Control", "Apoptosis", "Pyroptosis"
        )
    )
confusion_matrix$Predicted_Label <- factor(
    confusion_matrix$Predicted_Label, levels = c(
       "Pyroptosis", "Apoptosis", "Control"
        )
    )


# mutate the shuffled_data column
confusion_matrix$shuffled_data <- gsub("TRUE", "Shuffled", confusion_matrix$shuffled_data)
confusion_matrix$shuffled_data <- gsub("FALSE", "Not Shuffled", confusion_matrix$shuffled_data)
# mutate the data type column
confusion_matrix$data_split <- gsub("testing", "Testing", confusion_matrix$data_split)
confusion_matrix$data_split <- gsub("treatment_holdout", "Treatment Holdout", confusion_matrix$data_split)
confusion_matrix$data_split <- gsub("holdout", "Hold Out Wells", confusion_matrix$data_split)
# make the data split column a factor
confusion_matrix$data_split <- factor(confusion_matrix$data_split, levels = c(
    "Testing", "Hold Out Wells","Treatment Holdout"
))


# split the confusion matrix into testing/holdout and treatment holdout
confusion_matrix_testing <- confusion_matrix[confusion_matrix$data_split == "Testing",]
confusion_matrix_holdout <- confusion_matrix[confusion_matrix$data_split == "Hold Out Wells",]
confusion_matrix_treatment_holdout <- confusion_matrix[confusion_matrix$data_split == "Treatment Holdout",]

# combine the testing and holdout dataframes
confusion_matrix_testing_holdout <- rbind(confusion_matrix_testing, confusion_matrix_holdout)


# plot dimensions
width <- 14
height <- 11
options(repr.plot.width = width, repr.plot.height = height)
# plot a confusion matrix
confusion_matrix_plot <- (
    ggplot(confusion_matrix_testing_holdout, aes(x = True_Label, y = Predicted_Label))
    + facet_grid(data_split~shuffled_data)
    + geom_point(aes(color = Recall), size = 50, shape = 15)
    + geom_text(aes(label = Recall))
    + scale_color_gradient("Recall", low = "white", high = "dark red",limits = c(0, 1))
    + theme_bw()
    + ylab("Predicted Class")
    + xlab("True Class")
    + figure_theme


)
ggsave(confusion_matrix_plot_path, confusion_matrix_plot, width = width, height = height, dpi = 600)
confusion_matrix_plot


# read in the data
treatment_holdout_single_cell_predictions_path <- file.path(
    paste0(
        "../../results/Multi_Class/",model_name,"/",cell_type,"/treatment_holdout_single_cell_predictions.parquet"
    )
)
# read in the data from the parquet file
treatment_holdout_single_cell_predictions <- arrow::read_parquet(
    treatment_holdout_single_cell_predictions_path
)
head(treatment_holdout_single_cell_predictions)

unique(treatment_holdout_single_cell_predictions$shuffle)
# split the data into shuffled and not shuffled
treatment_holdout_single_cell_predictions_shuffled <- treatment_holdout_single_cell_predictions[
    treatment_holdout_single_cell_predictions$shuffle == "TRUE",]
treatment_holdout_single_cell_predictions_not_shuffled <- treatment_holdout_single_cell_predictions[
    treatment_holdout_single_cell_predictions$shuffle == "FALSE",]

# calculate the recall of class 2 (pyroptosis)
recall_shuffle <- treatment_holdout_single_cell_predictions_shuffled %>%
    filter(true_label == 2 & predicted_label == 2) %>%
    nrow() / treatment_holdout_single_cell_predictions_shuffled %>%
    filter(true_label == 2) %>%
    nrow()
# calculate the precision of class 2 (pyroptosis)
precision_shuffle <- treatment_holdout_single_cell_predictions_shuffled %>%
    filter(true_label == 2 & predicted_label == 2) %>%
    nrow() / treatment_holdout_single_cell_predictions_shuffled %>%
    filter(predicted_label == 2) %>%
    nrow()
# calculate the F1 score which is the harmonic mean of precision and recall
f1_score_shuffle <- 2 * (precision_shuffle * recall_shuffle) / (precision_shuffle + recall_shuffle)


# calculate the recall of class 2 (pyroptosis)
recall <- treatment_holdout_single_cell_predictions_not_shuffled %>%
    filter(true_label == 2 & predicted_label == 2) %>%
    nrow() / treatment_holdout_single_cell_predictions_not_shuffled %>%
    filter(true_label == 2) %>%
    nrow()
# calculate the precision of class 2 (pyroptosis)
precision <- treatment_holdout_single_cell_predictions_not_shuffled %>%
    filter(true_label == 2 & predicted_label == 2) %>%
    nrow() / treatment_holdout_single_cell_predictions_not_shuffled %>%
    filter(predicted_label == 2) %>%
    nrow()
# calculate the F1 score which is the harmonic mean of precision and recall
f1_score <- 2 * (precision * recall) / (precision + recall)

print(paste("F1 score for shuffled data:", f1_score_shuffle))
print(paste("F1 score for not shuffled data:", f1_score))


# replace strings in pr_curves
PR_curves$label <- gsub("apoptosis", "Apoptosis", PR_curves$label)
PR_curves$label <- gsub("healthy", "Control", PR_curves$label)
PR_curves$label <- gsub("pyroptosis", "Pyroptosis", PR_curves$label)

PR_curves$data_split <- gsub("train", "Training", PR_curves$data_split)
PR_curves$data_split <- gsub("testing", "Testing", PR_curves$data_split)
PR_curves$data_split <- gsub("validation", "Validation", PR_curves$data_split)
PR_curves$data_split <- gsub("treatment_holdout", "Treatment Holdout", PR_curves$data_split)
PR_curves$data_split <- gsub("holdout", "Holdout", PR_curves$data_split)

# factorize the data_split column
PR_curves$data_split <- factor(PR_curves$data_split, levels = c(
    "Training", "Validation", "Testing","Treatment Holdout", "Holdout"
))

unique(PR_curves$label)
unique(PR_curves$data_split)

# replace strings in pr_curves shuffle
PR_curves$shuffle <- gsub("TRUE", "Shuffled", PR_curves$shuffle)
PR_curves$shuffle <- gsub("FALSE", "Not Shuffled", PR_curves$shuffle)

# factorize the shuffled_data column
PR_curves$shuffle <- factor(PR_curves$shuffle, levels = c(
    "Not Shuffled", "Shuffled"
))



# make a line plot that has the shuffled and not shuffled lines
# with shuffled lines dashed and not shuffled lines solid
# color by label
pr_plot <- (
    ggplot(PR_curves, aes(x = recall, y = precision, color = label, linetype = label))
    + geom_line(aes(linetype = shuffle))
    # + scale_linetype_manual(values = c("solid", "dashed"))
    + facet_wrap(~data_split)
    + theme_bw()
)
pr_plot

