list_of_packages <- c("dplyr","tidyr","biomaRt")
for(package in list_of_packages){
suppressPackageStartupMessages(suppressMessages(suppressWarnings(library(package,character.only=TRUE))))
}

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

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host = "www.ensembl.org", path = "/biomart/martservice")

data_path <- file.path("../data/gtex_data")
data_files <- list.files(data_path, pattern = "gct", full.names = TRUE)

rna_df <- data.frame()
for (file in data_files){
    tmp_df <- read.table(file, header = TRUE, skip = 2, sep = "\t", stringsAsFactors = FALSE)
    # trim the decimal from the Name column

    tmp_df <- tmp_df[tmp_df$Description %in% genes,]

    # rename the Description column to geneID
    colnames(tmp_df)[2] <- "geneID"
    colnames(tmp_df)[1] <- "Ensembl"

    # drop the Name column
    tmp_df <- tmp_df[, !names(tmp_df) %in% "Name"]
    tmp_df <- tmp_df[, !names(tmp_df) %in% "Description"]
    # convert the GeneID to gene symbol
    gene_id <- unique(tmp_df$geneID)
    lookup <- getBM(
        mart = mart,
        attributes = c('transcript_length','hgnc_symbol'),
        filter = 'hgnc_symbol',
        values = gene_id,
        uniqueRows = TRUE
        )
    tmp_df <- merge(tmp_df, lookup, by.x="geneID", by.y="hgnc_symbol", all=TRUE)
    # make the df tidy long format
    tmp_df <- tmp_df %>%
        pivot_longer(cols = -c(Ensembl, geneID,transcript_length), names_to = "Sample", values_to = "expression") %>%
        mutate(tissue = gsub("../data/gtex_data/gene_reads_v10_", "", file)) %>%
        mutate(tissue = gsub(".gct", "", tissue))
    rna_df <- rbind(rna_df, tmp_df)
}
head(rna_df)

rna_df$fpkm <- rna_df$expression / rna_df$transcript_length

# save the compiled gtex data to a file

all_gtex_tissue_data <- file.path("../data/genes_of_interest_gtex_tissue_data.tsv")


write.table(rna_df, file = all_gtex_tissue_data, sep = "\t", row.names = FALSE)
