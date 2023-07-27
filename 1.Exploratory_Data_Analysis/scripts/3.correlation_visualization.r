suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(corrplot)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(argparse)))

# set up parser
parser <- ArgumentParser()

# set up arguments
parser$add_argument("--cell_type", default="none",
    help="define the cell type"
    )

parser$add_argument("--level", default="none",
    help="defines the data level of aggregation"
    )

parser$add_argument("--group", default="none",
    help="defines the group to be used for correlation analysis"
    )

# parse arguments
args <- parser$parse_args()

# define vars from parsed args
cell_type = args$cell_type
level = args$level
group = args$group

project_path = getwd()
filename = paste0(group,".csv")
corr_file = file.path(project_path, "results","correlation",cell_type,level,filename)
print(corr_file)
# read in csv file
df <- read.csv(corr_file)
row.names(df) <- df$oneb_Metadata_Treatment_Dose_Inhibitor_Dose
# drop first column
df <- df[, -1]
df <- as.matrix(df)

# set path of the plot
mainDir = project_path
figure_dir = file.path(mainDir, "Figures", "corrplot", cell_type, level)
print(figure_dir)
# set path of the plot
if (file.exists(figure_dir, recursive = TRUE)){
    setwd(figure_dir)
} else {
    dir.create(figure_dir, recursive = TRUE)
    setwd(figure_dir)
}
print(figure_dir)

# plot corr plot for all wells
# set the file name
file = paste0(group,".png")
print(file)
png(height=1800, width=1800, file=(file.path(figure_dir,file)), type = "cairo")
corrplot(df,
        type = 'lower',
        order = 'hclust',
        tl.col = 'black',
        cl.ratio = 0.2,
        cl.cex=2,
        tl.cex=2,
        tl.srt = 45,
        mar=c(0,0,4,0))
dev.off()
