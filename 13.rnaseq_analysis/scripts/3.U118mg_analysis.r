list_of_packages <- c("ggplot2", "dplyr", "DESeq2","biomaRt","tidyr", "tidyverse")
for(package in list_of_packages){
suppressPackageStartupMessages(suppressMessages(suppressWarnings(library(package,character.only=TRUE))))
}

figures_dir<-file.path("../figures")
if(!dir.exists(figures_dir)){
 dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
}

dataset1 <- read.table("../data/U118mg/GSE152291_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")
dataset2 <- read.table("../data/U118mg/GSE48865_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")

# keep only GSM4610662 and GSM4610663 columns from dataset1
dataset1 <- dataset1 %>% dplyr::select(GeneID,GSM4610662, GSM4610663)
dim(dataset1)
dim(dataset2)
# combine the two datasets
merged_dataset <- merge(dataset1, dataset2, by="GeneID", all=TRUE)
dim(merged_dataset)
head(merged_dataset)

genes <- c(
    "GSDMA",
    "GSDMB",
    "GSDMC",
    "GSDMD",
    "GSDME",
    "TLR4",
    "TLR5",
    "CASP1",
    "CASP2",
    "CASP3",
    "CASP4",
    "CASP5",
    "CASP6",
    "CASP7",
    "CASP8",
    "CASP9",
    "CASP10",
    "CASP11"
)


mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

# run and if if an error occurs, try again
while (TRUE) {
    tryCatch({
        entrezgene_id_list <- getBM(
            mart = mart,
            attributes = c('hgnc_symbol','entrezgene_id', "transcript_length"),
            filter = 'hgnc_symbol',
            values = genes,
            uniqueRows = TRUE)
        break
    }, error = function(e) {
        print(e)
    })
}
head(entrezgene_id_list)


# get the average transcript length for each gene group by hgnc_symbol and entrezgene_id
entrezgene_id_list <- entrezgene_id_list %>%
    group_by(hgnc_symbol, entrezgene_id) %>%
    summarise(transcript_length = mean(transcript_length)) %>%
    ungroup()
    # get only genes that are present in the entrezgene_id_column of the entrezgene_id_list
merged_dataset <- merged_dataset[merged_dataset$GeneID %in% entrezgene_id_list$entrezgene_id,]

# convert entrezgene_id_list from a tibble to a data.frame
entrezgene_id_list <- as.data.frame(entrezgene_id_list)
merged_dataset <- as.data.frame(merged_dataset)
merged_dataset$entrezgene_id <- as.character(merged_dataset$GeneID)

head(merged_dataset,1)
head(entrezgene_id_list,1)

merged_dataset$entrezgene_id
entrezgene_id_list$entrezgene_id

# merge the two datasets on the entrezgene_id column
merged_dataset <- merge(merged_dataset, entrezgene_id_list, by="entrezgene_id")

# get the counts
# define the fpkm for each column that begins with GSM
for (i in 1:ncol(merged_dataset)) {
if (grepl("GSM", colnames(merged_dataset)[i])) {
 merged_dataset[,colnames(merged_dataset)[i]] <- merged_dataset[,colnames(merged_dataset)[i]] / merged_dataset$transcript_length
}
}


# remove GeneID column
merged_dataset <- merged_dataset[,!grepl("GeneID", colnames(merged_dataset))]
# convert to tibble

merged_dataset <- as_tibble(merged_dataset)
head(merged_dataset)


# melt the dataset to long format
merged_dataset_long <- merged_dataset %>%
    dplyr::select(entrezgene_id, hgnc_symbol, transcript_length, starts_with("GSM")) %>%
    pivot_longer(cols = starts_with("GSM"), names_to = "sample", values_to = "fpkm")
# add a column to indicate u118mg bool for GSM4610662 and GSM4610663
merged_dataset_long$u118mg <- ifelse(grepl("GSM4610662|GSM4610663", merged_dataset_long$sample), TRUE, FALSE)
head(merged_dataset_long)

# get the mean and standard deviation for each gene for each sample
merged_dataset_long <- merged_dataset_long %>%
    group_by(hgnc_symbol, entrezgene_id,sample,u118mg)
head(merged_dataset_long)

# save the dataset to a file parquet file
arrow::write_parquet(merged_dataset_long,"../data/U118mg/merged_dataset_long.parquet")
