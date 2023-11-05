suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))


cell_type = "PBMC"


# path set
input_file_path <- file.path(paste0("../results/","regression/",cell_type))
# read in the data
output_path <- file.path(paste0("../figures/","regression/",cell_type,"/"))
# create output directory if it doesn't exist
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)


## function to process the data for visualization
process_subset_data <- function(data_path){
    # read in the data
    data <- read.csv(data_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    # get the basename of the files

    data <- data %>%
        dplyr::arrange(desc(abs(coefficients))) %>%
        tidyr::separate(
            feature_names,
            into = c(
                "compartment",
                "feature_group",
                "measurement",
                "channel",
                "parameter1",
                "parameter2"
            ),
            sep = "_",
            remove = FALSE
        ) %>%
        dplyr::mutate(channel_cleaned = channel) %>%
        dplyr::arrange(desc(abs(coefficients)))

    # Clean channel for visualization
    data$channel_learned <- dplyr::recode(data$channel,
            "CorrDNA" = "nuclei",
            "CorrMito" = "Mito",
            "CorrER" = "ER",
            "CorrGasdermin" = "gasdermin",
            "CorrPM" = "PM",
            .default = "other",
            .missing="other"
    )
    data <- data %>%
        dplyr::group_by(feature_group, channel_learned, compartment) %>%
        dplyr::slice_max(order_by = coefficients, n = 1)
    return(data)
}


# get all files in a directory
files <- list.files(path = input_file_path, pattern = "*.csv", full.names = TRUE)
coef_gg_file <- file.path(paste0(output_path,"/","top_abs_val_coefficients_enet.pdf"))
pdf(file=coef_gg_file, width=7, height=4)
for (i in files){
    filename <- basename(i)
    # split the string at the first _
    filename <- strsplit(filename, "_", fixed = TRUE)[[1]]
    cytokine <- filename[1]
    shuffle <- filename[2]
    # preprocess the data
    data <- process_subset_data(i)
    # plot the data
    coef_gg <- (
        ggplot(data, aes(x = channel_learned, y = feature_group))
        + geom_point(aes(fill = abs(coefficients)), pch = 22, size = 6)
        + facet_wrap("~compartment", ncol = 3)
        + theme_bw()
        + scale_fill_continuous(
            name="Top Abs. val\ntreatment\nlinear model\ncoefficient",
            low = "darkblue",
            high = "yellow",
        )
        + xlab("Channel")
        + ylab("Feature")
        + theme(
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            title = element_text(size = 14),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 12),
        )
        # rotate x axis labels
        + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        + ggtitle(paste0("Top Abs. val treatment ElasticNet coefficients for \n",cytokine,"\n",shuffle," model"))
        + theme(plot.title = element_text(hjust = 0.5))
    )
    plot(coef_gg)
}
dev.off()

