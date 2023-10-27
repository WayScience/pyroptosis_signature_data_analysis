suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))

model_name <- "DMSO_0.025_vs_LPS_100"
celltype <- "SHSY5Y"

# read in the data
output_file <- file.path(
    "..","..","figures","Regression",model_name,celltype,"regression_plots.pdf"
)
results_dir <- file.path(
    "..","..","results","Regression",model_name,celltype
)
results_file <- file.path(
    results_dir,"regression_results_training.csv"
)

# Read in the results file
df <- read.csv(results_file)
head(df,3)

print(unique(df$shuffled_data))
print(unique(df$model_type))
print(length(unique(df$cytokine)))
print((unique(df$oneb_Metadata_Treatment_Dose_Inhibitor_Dose)))

# # select rwos that are training data only
# df <- df %>% filter(model_type == "validation")

pdf(file=output_file)
for (i in 1:length(unique(df$cytokine))){
    tmpdf <- df %>% filter(cytokine == unique(df$cytokine)[i])
    # linear trend + confidence interval
    p3 <- (
    ggplot(tmpdf, aes(
        x=Actual,
        y=Average.Predicted,
        shape=shuffled_data,
        color=oneb_Metadata_Treatment_Dose_Inhibitor_Dose,
        linetype=shuffled_data
        ))
    + geom_point()
    # + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)
    + geom_smooth(method=lm , se=TRUE, formula = y ~ x)
    + labs(x="Actual", y="Predicted", title="Regression Plot")
    + theme(plot.title = element_text(hjust = 0.5))
    + theme(legend.title=element_blank())
    + ggtitle(paste0("Regression for ",unique(df$cytokine)[i])
        )
    )
    plot(p3)
}
dev.off()

# linear trend + confidence interval
p3 <- (
    ggplot(tmpdf, aes(
        x=Actual,
        y=Average.Predicted,
        color=oneb_Metadata_Treatment_Dose_Inhibitor_Dose,
        linetype=shuffled_data
        ))
    + geom_point()
    # + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE)
    + geom_smooth(method=lm , se=TRUE, formula = y ~ x)
    + labs(x="Actual", y="Predicted", title="Regression Plot")
    + theme(plot.title = element_text(hjust = 0.5))
    + theme(legend.title=element_blank())
    + ggtitle(paste0("Regression for ",unique(df$cytokine)[i])
        )
    )
p3
