#________Load required packages________________
library(TCseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(cowplot)
library(ggplot2)

#_________lncRNA_tRNA count data______________________________________________
# lncRNA_tRNA list
count_lncRNA_tRNA_genes <- read.table("overlapping_tRNA_lncRNA.txt",  header = TRUE, sep = "\t", row.names = "gene_id")

# CHECK to see
index <- which(rownames(count_lncRNA_tRNA_genes) == "ENSG00000178977.4")
count_lncRNA_tRNA_genes[index, ]

#_______________ sample information ________________
countsTable <- read.table("overlapping_tRNA_lncRNA.txt",  header = TRUE, sep = "\t", row.names = "gene_id")
countsTable <- countsTable[, c(2:21)]

sample_list <- colnames(countsTable)
my_data_frame <- data.frame(x = sample_list)
colnames(my_data_frame) <- "sampleid"

# Extract timepoint from sample ID
my_data_frame$timepoint <- gsub(".*_(\\d+d).*", "\\1", my_data_frame$sampleid)

# Rename sample IDs
sample_name <- paste0("s", 1:20)
my_data_frame$sampleid <- sample_name
colnames(countsTable) <- sample_name

# Assign group values
group_values <- seq(1, 20)
my_data_frame$group <- group_values
experiment <- my_data_frame

#_____________genomic intervals_______________________
RPKM_data_analysis <- read.table("MSC.GENCODEv41.RPKM.genPos.consensus.txt", header = TRUE, sep = "\t")

# Match indices
indices <- match(rownames(count_lncRNA_tRNA_genes), RPKM_data_analysis$Geneid)
RPKM_lncRNA_tRNA_data_analysis <- RPKM_data_analysis[indices, ]

# Check
identical(RPKM_lncRNA_tRNA_data_analysis$Geneid, rownames(count_lncRNA_tRNA_genes))

genomicIntervals <- RPKM_lncRNA_tRNA_data_analysis[, 1:4]
colnames(genomicIntervals) <- c("id", "chr", "start", "end")

# Clean up non-numeric data in countsTable
row_names <- rownames(countsTable)
countsTable <- apply(countsTable, 2, function(x) as.numeric(as.character(x)))
non_numeric_rows <- which(!complete.cases(countsTable))
countsTable[non_numeric_rows, ] <- NA
rownames(countsTable) <- row_names

# Time course analysis
tca <- TCA(design = experiment, genomicFeature = genomicIntervals, counts = countsTable)
counts(tca) <- countsTable

# Load SummarizedExperiment package and create a SummarizedExperiment object
suppressWarnings(library(SummarizedExperiment))
se <- SummarizedExperiment(assays = list(counts = countsTable), colData = experiment)
tca <- TCAFromSummarizedExperiment(se = se, genomicFeature = genomicIntervals)

# Differential analysis
tca <- DBanalysis(tca)
tca <- DBanalysis(tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)

# Get differential analysis results
DBres <- DBresult(tca, group1 = "0d", group2 = c("1d", "3d", "7d", "14d"))
head(DBres$"1dvs0d")

# Filter significant differential events
DBres.sig <- DBresult(tca, group1 = "0d", group2 = c("1d", "3d", "7d", "14d"), top.sig = TRUE)
head(DBres.sig$"14dvs0d")

# Further filtering for log2-fold > 2 and p-value < 0.01
DBres.sig <- DBresult(tca, group1 = "0d", group2 = c("1d", "3d", "7d", "14d"), top.sig = TRUE, pvalue = "paj", pvalue.threshold = 0.01, abs.fold = 2, direction = "both")

# Save differential gene results
d_14_vs_0d <- as.data.frame(DBres.sig[["14dvs0d"]])
write.csv(d_14_vs_0d, file = "14dvs0d_sig_0.01.csv")

# Similar for other time points
d_1_vs_0d <- as.data.frame(DBres.sig[["1dvs0d"]])
write.csv(d_1_vs_0d, file = "1dvs0d_sig_0.01.csv")
d_3_vs_0d <- as.data.frame(DBres.sig[["3dvs0d"]])
write.csv(d_3_vs_0d, file = "3dvs0d_sig_0.01.csv")
d_7_vs_0d <- as.data.frame(DBres.sig[["7dvs0d"]])
write.csv(d_7_vs_0d, file = "7dvs0d_sig_0.01.csv")

# Time course table values (logFC)
tca <- timecourseTable(tca, value = "FC", norm.method = "rpkm", filter = TRUE)
tca <- timecourseTable(tca, value = "FC", control.group = "0d", norm.method = "rpkm", filter = TRUE)
tca <- timecourseTable(tca, value = "FC", control.group = "0d", norm.method = "rpkm", filter = TRUE, pvalue = "fdr", pvalue.threshold = 0.01, abs.fold = 2, direction = "both")
m <- tcTable(tca)
head(m)
which(rownames(m) == "ENSG00000178977.4")

# Time course table values (normalized read counts)
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE, pvalue = "fdr", pvalue.threshold = 0.01, abs.fold = 2, direction = "both")
t <- tcTable(tca)
head(t)
which(rownames(t) == "ENSG00000178977.4")

# Save differential gene expression values
timepoint_expression_value_fdr_0.01 <- as.data.frame(t)
write.csv(timepoint_expression_value_fdr_0.01, file = "timepoint_expression_value_fdr_0.01.csv", row.names = TRUE)

# Check Sameen's lncRNA "ENSG00000178977.4"
row_index <- which(rownames(t) == "ENSG00000178977.4")
print(t[row_index, ])

# Clustering analysis
set.seed(123)
tca <- timeclust(tca, algo = "cm", k = 5, standardize = TRUE)

# Visualization
p <- timeclustplot(tca, value = "z-score(RPKM)", cols = 3)
combined_plot <- cowplot::plot_grid(plotlist = p, ncol = 3)
ggsave(filename = "clustering.pdf", plot = combined_plot, device = "pdf", width = 14, height = 18)

# Save all clusters data
for (i in 1:5) {
  write.table(p[[i]]$data, file = paste0("cluster", i, ".csv"), sep = ",")
}

# Plot clusters with annotations
for (i in 1:5) {
  genes_cluster <- table(p[[i]]$data$group)
  number_of_genes <- nrow(p[[i]]$data)
  unique_genes <- length(unique(p[[i]]$data$group))
  total_genes <- paste0("N = ", unique_genes)
  
  p_with_label <- p[[i]] +
    annotate("text", x = Inf, y = Inf, label = total_genes, vjust = 1, hjust = 1, size = 5, color = "black")
  
  print(p_with_label)
  ggsave(filename = paste0("cluster", i, ".pdf"), plot = p_with_label, device = "pdf")
}
