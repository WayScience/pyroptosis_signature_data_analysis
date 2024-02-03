library(ggplot2)

# set path to the data
class_df_path <- file.path("..","data","processed","aggregate_mAPs","mAP_scores_class.csv")

# read in the data
class_df <- read.csv(class_df_path)

head(class_df)

# declare the shuffled column as a factor
# replace the values in the shuffled column
class_df$shuffled <- gsub("features_shuffled", "Shuffled features", class_df$shuffled)
class_df$shuffled <- gsub("phenotype_shuffled", "Shuffled phenotypes", class_df$shuffled)
class_df$shuffled <- gsub("non-shuffled", "Non-shuffled", class_df$shuffled)
class_df$shuffled <- factor(class_df$shuffled, levels = c( "Non-shuffled", "Shuffled features", "Shuffled phenotypes"))
class_df$Metadata_labels <- factor(class_df$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))

# plot the data
barplot <- (
    ggplot(class_df, aes(x=Metadata_labels, y=mean_average_precision, fill=shuffled))
    + geom_bar(stat="identity", position="dodge")
    + labs(x="Class", y="mAP score")
    # legend title
    + scale_fill_discrete(name="Shuffle Type")
    + theme_bw()
    + ylim(0,1)
)
barplot


