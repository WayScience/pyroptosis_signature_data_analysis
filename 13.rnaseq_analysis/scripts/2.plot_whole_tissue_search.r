list_of_packages <- c("ggplot2", "dplyr", "DESeq2","biomaRt","tidyr")
for(package in list_of_packages){
suppressPackageStartupMessages(suppressMessages(suppressWarnings(library(package,character.only=TRUE))))
}

gsdm_genes <- c(
    "GSDMD",
    "GSDME"
)

all_gtex_tissue_data <- file.path("../data/genes_of_interest_gtex_tissue_data.tsv")
gtex_df <- read.table(all_gtex_tissue_data, header = TRUE, sep = "\t")
# trim the decimals in Ensembl
gtex_df$Ensembl <- gsub("\\..*", "", gtex_df$Ensembl)
head(gtex_df)

gsdm_df <- gtex_df[gtex_df$geneID %in% gsdm_genes,]
gsdm_df$geneID <- factor(gsdm_df$geneID, levels = gsdm_genes)
head(gsdm_df)

# get the mean fpkm of each gene in each tissue across all samples
gsdm_df <- gsdm_df %>%
  group_by(geneID, tissue) %>%
  summarise(mean_fpkm = mean(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup()
head(gtex_df)

width <- 30
height <- 30
options(repr.plot.width = width, repr.plot.height = height)
gsdm_plot <- (
    ggplot(gsdm_df, aes(x = geneID, y = mean_fpkm, fill = geneID))
    + geom_bar(stat = "identity", position = "dodge")
    + facet_wrap(~tissue, nrow= 10, scales = "free_y")
    + theme(legend.position = "none")
    + theme(strip.text = element_text(size = 16))
    + theme(axis.text.x = element_text(size = 16)))
# save the plot
ggsave("../figures/gtex_gsdm_expression_by_tissue.png", gsdm_plot, width = width, height = height, units = "in")
gsdm_plot

tlr_and_casp_genes <- c(
 "TLR4",
 "TLR5",
 "CASP1",
 "CASP3",
 "CASP7",
 "CASP11"
)
casp_df <- gtex_df[gtex_df$geneID %in% tlr_and_casp_genes,]
head(casp_df)

# get the mean fpkm of each gene in each tissue across all samples
casp_df <- casp_df %>%
  group_by(geneID, tissue) %>%
  summarise(mean_fpkm = mean(fpkm), sd_fpkm = sd(fpkm)) %>%
  ungroup()
head(casp_df)

casp_plot <- (
    ggplot(casp_df, aes(x = geneID, y = mean_fpkm, fill = geneID))
    + geom_bar(stat = "identity", position = "dodge")
    + facet_wrap(~tissue, nrow= 10, scales = "free_y")
    # rotate the x-axis labels
    + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    # hide the legend
    + theme(legend.position = "none")
    + theme(strip.text = element_text(size = 16))
    + theme(axis.text.x = element_text(size = 16)))
casp_plot
# save the plot
ggsave("../figures/gtex_casp_tlr_expression_by_tissue.png", casp_plot, width = width, height = height, units = "in")
