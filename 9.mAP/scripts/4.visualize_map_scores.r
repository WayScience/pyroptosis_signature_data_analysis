library(ggplot2)
library(RColorBrewer)
source("../../figures/utils/figure_themes.r")

width <- 8
height <- 6
options(repr.plot.width=width, repr.plot.height=height)

# set path to the data morphology
class_df_morphology_path <- file.path("..","data","processed","aggregate_mAPs","morphology","mAP_scores_class.csv")
reg_df_morphology_path <- file.path("..","data","processed","mAP_scores","morphology","mAP_scores_regular_class.csv")
shuffled_df_morphology_path <- file.path("..","data","processed","mAP_scores","morphology","mAP_scores_shuffled_class.csv")
shuffled_feature_space_df_morphology_path <- file.path("..","data","processed","mAP_scores","morphology","mAP_scores_shuffled_feature_space_class.csv")
# read in the data
class_df_morphology <- read.csv(class_df_morphology_path)
reg_df_morphology <- read.csv(reg_df_morphology_path)
shuffled_df_morphology <- read.csv(shuffled_df_morphology_path)
shuffled_feature_space_df_morphology <- read.csv(shuffled_feature_space_df_morphology_path)

head(class_df_morphology)


# set path to the data secretome
class_df_secretome_path <- file.path("..","data","processed","aggregate_mAPs","secretome","mAP_scores_class.csv")
reg_df_secretome_path <- file.path("..","data","processed","mAP_scores","secretome","mAP_scores_regular_class.csv")
shuffled_df_secretome_path <- file.path("..","data","processed","mAP_scores","secretome","mAP_scores_shuffled_class.csv")
shuffled_feature_space_df_secretome_path <- file.path("..","data","processed","mAP_scores","secretome","mAP_scores_shuffled_feature_space_class.csv")
# read in the data
class_df_secretome <- read.csv(class_df_secretome_path)
reg_df_secretome <- read.csv(reg_df_secretome_path)
shuffled_df_secretome <- read.csv(shuffled_df_secretome_path)
shuffled_feature_space_df_secretome <- read.csv(shuffled_feature_space_df_secretome_path)

head(class_df_secretome)

# declare the shuffled column as a factor
# replace the values in the shuffled column
# declare the shuffled column as a factor
# replace the values in the shuffled column
class_df_morphology$shuffled <- gsub("features_shuffled", "Shuffled features", class_df_morphology$shuffled)
class_df_morphology$shuffled <- gsub("phenotype_shuffled", "Shuffled phenotypes", class_df_morphology$shuffled)
class_df_morphology$shuffled <- gsub("non-shuffled", "Non-shuffled", class_df_morphology$shuffled)
class_df_morphology$shuffled <- factor(class_df_morphology$shuffled, levels = c( "Non-shuffled", "Shuffled features", "Shuffled phenotypes"))
class_df_morphology$Metadata_labels <- factor(class_df_morphology$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))

class_df_secretome$shuffled <- gsub("features_shuffled", "Shuffled features", class_df_secretome$shuffled)
class_df_secretome$shuffled <- gsub("phenotype_shuffled", "Shuffled phenotypes", class_df_secretome$shuffled)
class_df_secretome$shuffled <- gsub("non-shuffled", "Non-shuffled", class_df_secretome$shuffled)
class_df_secretome$shuffled <- factor(class_df_secretome$shuffled, levels = c( "Non-shuffled", "Shuffled features", "Shuffled phenotypes"))
class_df_secretome$Metadata_labels <- factor(class_df_secretome$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))



# plot the data
barplot_morphology <- (
    ggplot(class_df_morphology, aes(x=Metadata_labels, y=mean_average_precision, fill=shuffled))
    + geom_bar(stat="identity", position="dodge")
    + labs(x="Class", y="mAP score")
    # legend title
    + scale_fill_discrete(name="Shuffle Type")
    + theme_bw()
    + ylim(0,1)
    + ggtitle("Morphology")
    + figure_theme
)

barplot_secretome <- (
    ggplot(class_df_secretome, aes(x=Metadata_labels, y=mean_average_precision, fill=shuffled))
    + geom_bar(stat="identity", position="dodge")
    + labs(x="Class", y="mAP score")
    # legend title
    + scale_fill_discrete(name="Shuffle Type")
    + theme_bw()
    + ylim(0,1)
    + ggtitle("Secretome")
    + figure_theme

)


barplot_morphology
barplot_secretome



# combine the dataframes
all_df_morphology <- rbind(reg_df_morphology, shuffled_df_morphology, shuffled_feature_space_df_morphology)
all_df_morphology$shuffled <- gsub("features_shuffled", "Shuffled features", all_df_morphology$shuffled)
all_df_morphology$shuffled <- gsub("phenotype_shuffled", "Shuffled phenotypes", all_df_morphology$shuffled)
all_df_morphology$shuffled <- gsub("non-shuffled", "Non-shuffled", all_df_morphology$shuffled)
all_df_morphology$shuffled <- factor(all_df_morphology$shuffled, levels = c( "Non-shuffled", "Shuffled features", "Shuffled phenotypes"))
all_df_morphology$Metadata_labels <- factor(all_df_morphology$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))
head(all_df_morphology)

all_df_secretome <- rbind(reg_df_secretome, shuffled_df_secretome, shuffled_feature_space_df_secretome)
all_df_secretome$shuffled <- gsub("features_shuffled", "Shuffled features", all_df_secretome$shuffled)
all_df_secretome$shuffled <- gsub("phenotype_shuffled", "Shuffled phenotypes", all_df_secretome$shuffled)
all_df_secretome$shuffled <- gsub("non-shuffled", "Non-shuffled", all_df_secretome$shuffled)
all_df_secretome$shuffled <- factor(all_df_secretome$shuffled, levels = c( "Non-shuffled", "Shuffled features", "Shuffled phenotypes"))
all_df_secretome$Metadata_labels <- factor(all_df_secretome$Metadata_labels, levels = c("Control", "Apoptosis", "Pyroptosis"))
head(all_df_secretome)

boxplot_morphology <- (
    ggplot(all_df_morphology, aes(x=Metadata_labels, y=average_precision, fill=shuffled))
    + geom_boxplot()
    + labs(x="Class", y="mAP score")
    # legend title
    + scale_fill_discrete(name="Shuffle Type")
    + theme_bw()
    + ylim(0,1)
    + ggtitle("Morphology")
    + figure_theme


)
boxplot_morphology

boxplot_secretome <- (
    ggplot(all_df_secretome, aes(x=Metadata_labels, y=average_precision, fill=shuffled))
    + geom_boxplot()
    + labs(x="Class", y="mAP score")
    # legend title
    + scale_fill_discrete(name="Shuffle Type")
    + theme_bw()
    + ylim(0,1)
    + ggtitle("Secretome")
    + figure_theme


)
boxplot_secretome

# cobine the dfs
head(all_df_morphology)
# get the average precision, shuffled, and Metadata_labels columns by name
subset_morphology <- all_df_morphology[,c("average_precision", "shuffled", "Metadata_labels")]
# rename the average_precision column to moprhology_ap
colnames(subset_morphology)[colnames(subset_morphology)=="average_precision"] <- "morphology_ap"

# get the average precision, shuffled, and Metadata_labels columns by name
subset_secretome <- all_df_secretome[,c("average_precision", "shuffled", "Metadata_labels")]
# rename the average_precision column to secretome_ap
colnames(subset_secretome)[colnames(subset_secretome)=="average_precision"] <- "secretome_ap"

# merge the dataframes
merged_df <- merge(subset_morphology, subset_secretome, by=c("shuffled", "Metadata_labels"))
head(merged_df)



# aggregate the data by shuffled and Metadata_labels
merged_agg <- aggregate(. ~ shuffled + Metadata_labels, data=merged_df, FUN=mean)
# combine the shuffled and Metadata_labels columns
merged_agg$group <- paste(merged_agg$shuffled, merged_agg$Metadata_labels, sep="_")
# change the text in the group column
merged_agg$group <- gsub("Non-shuffled_Control", "Non-shuffled\nControl", merged_agg$group)
merged_agg$group <- gsub("Shuffled features_Control", "Shuffled features\nControl", merged_agg$group)
merged_agg$group <- gsub("Shuffled phenotypes_Control", "Shuffled phenotypes\nControl", merged_agg$group)
merged_agg$group <- gsub("Non-shuffled_Apoptosis", "Non-shuffled\nApoptosis", merged_agg$group)
merged_agg$group <- gsub("Shuffled features_Apoptosis", "Shuffled features\nApoptosis", merged_agg$group)
merged_agg$group <- gsub("Shuffled phenotypes_Apoptosis", "Shuffled phenotypes\nApoptosis", merged_agg$group)
merged_agg$group <- gsub("Non-shuffled_Pyroptosis", "Non-shuffled\nPyroptosis", merged_agg$group)
merged_agg$group <- gsub("Shuffled features_Pyroptosis", "Shuffled features\nPyroptosis", merged_agg$group)
merged_agg$group <- gsub("Shuffled phenotypes_Pyroptosis", "Shuffled phenotypes\nPyroptosis", merged_agg$group)
# make the group column a factor
merged_agg$group <- factor(
    merged_agg$group,
    levels = c(
        "Non-shuffled\nControl",
        "Shuffled features\nControl",
        "Shuffled phenotypes\nControl",

        "Non-shuffled\nApoptosis",
        "Shuffled features\nApoptosis",
        "Shuffled phenotypes\nApoptosis",

        "Non-shuffled\nPyroptosis",
        "Shuffled features\nPyroptosis",
        "Shuffled phenotypes\nPyroptosis"))

merged_agg

# plot the data
scatter_compare <- (
    ggplot(merged_agg, aes(x=morphology_ap, y=secretome_ap, col = group, shape=group))
    + geom_point(size=3, alpha=1)
    + labs(x="Morphology mAP score", y="Secretome mAP score")
    + theme_bw()
    + ggtitle("Comparison of mAP scores")
    + ylim(0,1)
    + xlim(0,1)
    # Change the legend title
    # change the legend shape
    + scale_shape_manual(
        name="Class",
        labels=c(
            "Non-shuffled\nControl",
            "Non-shuffled\nApoptosis",
            "Non-shuffled\nPyroptosis",

            "Shuffled features\nControl",
            "Shuffled features\nApoptosis",
            "Shuffled features\nPyroptosis",

            "Shuffled phenotypes\nControl",
            "Shuffled phenotypes\nApoptosis",
            "Shuffled phenotypes\nPyroptosis"
        ),
        values=c(19, 17, 15, 19, 17, 15, 19, 17, 15)
    )
    + scale_color_manual(
        name="Class",
        labels=c(
            "Non-shuffled\nControl",
            "Non-shuffled\nApoptosis",
            "Non-shuffled\nPyroptosis",

            "Shuffled features\nControl",
            "Shuffled features\nApoptosis",
            "Shuffled features\nPyroptosis",

            "Shuffled phenotypes\nControl",
            "Shuffled phenotypes\nApoptosis",
            "Shuffled phenotypes\nPyroptosis"
        ),
        values=c(
            brewer.pal(3, "Dark2")[1],
            brewer.pal(3, "Dark2")[1],
            brewer.pal(3, "Dark2")[1],
            brewer.pal(3, "Dark2")[2],
            brewer.pal(3, "Dark2")[2],
            brewer.pal(3, "Dark2")[2],
            brewer.pal(3, "Dark2")[3],
            brewer.pal(3, "Dark2")[3],
            brewer.pal(3, "Dark2")[3]
    )
)
    + figure_theme

)
scatter_compare
