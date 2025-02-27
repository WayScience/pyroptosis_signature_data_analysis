list_of_packages <- c("ggplot2", "dplyr", "DESeq2","biomaRt","tidyr")
for(package in list_of_packages){
suppressPackageStartupMessages(suppressMessages(suppressWarnings(library(package,character.only=TRUE))))
}

figures_dir<-file.path("../figures")
if(!dir.exists(figures_dir)){
 dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
}

dataset1 <- read.table("../data/GSE74886_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")
dataset2 <- read.table("../data/GSE191270_raw_counts_GRCh38.p13_NCBI.tsv", header=TRUE, sep="\t")
head(dataset1)
head(dataset2)

# merge the two datasets
merged_dataset <- merge(dataset1, dataset2, by="GeneID", all=TRUE)
print(dim(merged_dataset))
head(merged_dataset)

# convert the GeneID to gene symbol
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
gene_id <- merged_dataset$GeneID
lookup <- getBM(
mart = mart,
attributes = c('entrezgene_id', 'ensembl_gene_id',
 'gene_biotype','hgnc_symbol','transcript_length'),
filter = 'entrezgene_id',
values = gene_id,
uniqueRows = TRUE)
# merge the two datasets
merged_dataset <- merge(merged_dataset, lookup, by.x="GeneID", by.y="entrezgene_id", all=TRUE)
# write the merged dataset to a file
head(merged_dataset)


dataset1 <- data.frame(merged_dataset$GeneID, merged_dataset$hgnc_symbol, merged_dataset$transcript_length, merged_dataset$GSM1937035, merged_dataset$GSM1937037)
dataset2 <- data.frame(merged_dataset$GeneID, merged_dataset$hgnc_symbol, merged_dataset$transcript_length, merged_dataset$GSM5742788, merged_dataset$GSM5742789, merged_dataset$GSM5742790, merged_dataset$GSM5742791, merged_dataset$GSM5742792, merged_dataset$GSM5742793, merged_dataset$GSM5742794, merged_dataset$GSM5742795, merged_dataset$GSM5742796)
# rename the columns to remove merged_dataset. prefix
colnames(dataset1) <- c("ensembl_gene_id",		"hgnc_symbol",	"transcript_length", "GSM1937035", "GSM1937037")
colnames(dataset2) <- c("ensembl_gene_id",		"hgnc_symbol",	"transcript_length", "GSM5742788", "GSM5742789", "GSM5742790", "GSM5742791", "GSM5742792", "GSM5742793", "GSM5742794", "GSM5742795", "GSM5742796")
head(dataset1)
head(dataset2)

# get the counts
# define the fpkm for each column that begins with GSM
for (i in 1:ncol(dataset1)) {
if (grepl("GSM", colnames(dataset1)[i])) {
 dataset1[,colnames(dataset1)[i]] <- dataset1[,colnames(dataset1)[i]] / dataset1$transcript_length
}
}
for (i in 1:ncol(dataset2)) {
if (grepl("GSM", colnames(dataset2)[i])) {
 dataset2[,colnames(dataset2)[i]] <- dataset2[,colnames(dataset2)[i]] / dataset2$transcript_length
}
}

# plot the FPKM for GSDMD genes
genes <- c(
 "GSDMA",
 "GSDMB",
 "GSDMC",
 "GSDMD",
 "GSDME"
)
# get only genes of interest for dataset1
dataset1 <- dataset1[dataset1$hgnc_symbol %in% genes,]
# get only genes of interest for dataset2
dataset2 <- dataset2[dataset2$hgnc_symbol %in% genes,]
# drop esemble gene id, transcript length and gene lengths
dataset1$ensembl_gene_id <- NULL
dataset1$transcript_length <- NULL
dataset1$gene_lengths <- NULL
dataset2$ensembl_gene_id <- NULL
dataset2$transcript_length <- NULL
dataset2$gene_lengths <- NULL

# make each dataset tidy long format
dataset1 <- gather(dataset1, key="sample", value="fpkm", -hgnc_symbol)
dataset2 <- gather(dataset2, key="sample", value="fpkm", -hgnc_symbol)
# concat the two datasets
dataset <- rbind(dataset1, dataset2)

# get the mean and standard deviation for each gene for each sample
dataset <- dataset %>%
group_by(hgnc_symbol, sample) %>%
summarise(mean_fpkm = mean(fpkm), sd_fpkm = sd(fpkm))
head(dataset)

plot <- (
 ggplot(dataset, aes(x=hgnc_symbol, y=mean_fpkm, fill=hgnc_symbol))
 + geom_bar(stat="identity")
 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
 + labs(y="Mean FPKM")
 # remove the x-axis label
 + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
 + geom_errorbar(aes(ymin=mean_fpkm-sd_fpkm, ymax=mean_fpkm+sd_fpkm), width=.2, position=position_dodge(.9))
 # rename legend title
+ labs(fill="Gene")
 + theme_bw()
 + facet_wrap(~sample)
 # remove x axis ticks and labels
    + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

)
ggsave(file.path(figures_dir, "GSDMD_gene_expression_SHSY5Y_dataset.png"), plot, width=12, height=8)
plot

