#________Load requiered packages________________
library(TCseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(cowplot)
library(ggplot2)

#_________lncRNA_tRNA count data______________________________________________
#lncRNA_tRNA list

count_lncRNA_tRNA_genes<- read.table("overlapping_tRNA_lncRNA.txt",  header=T, sep="\t", row.names="gene_id" )

#CHECK to see
index<- which(rownames(count_lncRNA_tRNA_genes) == "ENSG00000178977.4")
count_lncRNA_tRNA_genes[index,]

#_______________ sample information ________________

countsTable <- read.table("overlapping_tRNA_lncRNA.txt",  header=T, sep="\t", row.names="gene_id")
countsTable <- countsTable[,c(2:21)]
#countsTable <- lapply(countsTable, function(x) as.numeric(x))
sample_list <- colnames(countsTable)
my_data_frame <- data.frame(x= sample_list )
colnames(my_data_frame) <- "sampleid"
my_data_frame$timepoint <- gsub(".*_(\\d+d).*","\\1",my_data_frame$sampleid)
sample_name <- paste0("s",1:20)
my_data_frame$sampleid <- sample_name
colnames(countsTable) <- sample_name
group_values <- seq(1,20)
my_data_frame$group <- group_values
experiment <- my_data_frame

#_____________genomic intervals_______________________

RPKM_data_analysis <- read.table("MSC.GENCODEv41.RPKM.genPos.consensus.txt", header=T, sep="\t")
indices<- match(rownames(count_lncRNA_tRNA_genes), RPKM_data_analysis$Geneid)
RPKM_lncRNA_tRNA_data_analysis<- RPKM_data_analysis[indices,]

#check
identical(RPKM_lncRNA_tRNA_data_analysis$Geneid ,rownames(count_lncRNA_tRNA_genes))


genomicIntervals<- RPKM_lncRNA_tRNA_data_analysis[,1:4]


#genomicIntervals <- genomicIntervals[,c(2,3,4,1)]
colnames(genomicIntervals) <- c("id","chr","start","end")

# countsTable <- apply(countsTable,2 ,function(x) as.numeric(as.character(x)))
# non_numeric_col <- sapply(countsTable,function(x) any(is.na(x)))
# non_numeric_rows <- which(!complete.cases(countsTable))
# non_numeric_rows

row_names <- rownames(countsTable)

# Step 2: Convert data to numeric
countsTable <- apply(countsTable, 2, function(x) as.numeric(as.character(x)))

# Step 3: Identify non-numeric rows
non_numeric_rows <- which(!complete.cases(countsTable))

# Step 4: Handle non-numeric rows (e.g., replace with NA or perform some other operation)
# For example, you can replace non-numeric values with NA:
countsTable[non_numeric_rows, ] <- NA

# Step 5 (optional): Reassign modified data back to the original countsTable
# Only if you want to update the original countsTable
# If you don't need to modify the original data, you can skip this step.
rownames(countsTable) <- row_names



tca <- TCA(design = experiment, genomicFeature = genomicIntervals,   counts = countsTable)

counts(tca) <- countsTable

suppressWarnings(library(SummarizedExperiment))
se <- SummarizedExperiment(assays=list(counts = countsTable), colData = experiment)
tca <- TCAFromSummarizedExperiment(se = se, genomicFeature = genomicIntervals)

# Differential Analysis

tca <- DBanalysis(tca)
tca <- DBanalysis(tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
DBres <- DBresult(tca, group1 = "0d", group2 = c("1d","3d","7d","14d"))
str(DBres, strict.width = "cut")
head(DBres$"1dvs0d")

#Signiffcant differential events (log2-fold > 2 or log2-fold < -2, adjusted p-value < 0.05)
DBres.sig <- DBresult(tca, group1 = "0d", group2 = c("1d","3d","7d","14d"), top.sig = TRUE)
str(DBres.sig, strict.width = "cut")
head(DBres$"14dvs0d")

#Signiffcant differential events (log2-fold > 2 or log2-fold < -2, adjusted p-value < 0.01)
DBres.sig <- DBresult(tca, group1 = "0d", group2 = c("1d","3d","7d","14d"), top.sig = TRUE, pvalue ="paj",pvalue.threshold = 0.01, abs.fold = 2, direction = "both")
str(DBres.sig, strict.width = "cut")
head(DBres$"14dvs0d")



d_14_vs_0d <- as.data.frame(DBres.sig[["14dvs0d"]])
write.csv(d_14_vs_0d, file = "14dvs0d_sig_0.01.csv")

d_1_vs_0d <- as.data.frame(DBres.sig[["1dvs0d"]])
write.csv(d_1_vs_0d, file = "1dvs0d_sig_0.01.csv")

d_3_vs_0d <- as.data.frame(DBres.sig[["3dvs0d"]])
write.csv(d_3_vs_0d, file = "3dvs0d_sig_0.01.csv")

d_7_vs_0d <- as.data.frame(DBres.sig[["7dvs0d"]])
write.csv(d_7_vs_0d, file = "7dvs0d_sig_0.01.csv")

# values are logFC
tca <- timecourseTable(tca, value = "FC", norm.method = "rpkm", filter = TRUE)
tca <- timecourseTable(tca, value = "FC",control.group= "0d", norm.method = "rpkm", filter = TRUE)
tca <- timecourseTable(tca, value = "FC",control.group= "0d", norm.method = "rpkm", filter = TRUE, pvalue = "fdr", pvalue.threshold = 0.01, abs.fold = 2,direction = "both")
m <- tcTable(tca)
head( m)
which(rownames(m)== "ENSG00000178977.4")
# values are normalized read counts
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE)
tca <- timecourseTable(tca, value = "expression", norm.method = "rpkm", filter = TRUE, pvalue = "fdr", pvalue.threshold = 0.01, abs.fold = 2,direction = "both")
t <- tcTable(tca)
head(t)
which(rownames(t)== "ENSG00000178977.4")
#save Differential genes
timepoint_expression_value_fdr_0.01<- as.data.frame(t)
write.csv(timepoint_expression_value_fdr_0.01, file = "timepoint_expression_value_fdr_0.01.csv",row.names=TRUE)
#check Sameen's lncRNA "ENSG00000178977.4"
row_index <- which(rownames(t) == "ENSG00000178977.4")
print(t[row_index,])

tca
set.seed(123)
tca <- timeclust(tca, algo = "cm", k = 5, standardize = TRUE)



#visualization
p <- timeclustplot(tca, value = "z-score(RPKM)", cols = 3)

combined_plot <- cowplot::plot_grid(plotlist =p,ncol =3)
ggsave(filename = "clustering.pdf",plot = combined_plot, device = "pdf",width =14,height= 18)

#save all clusters data
write.table(p[[1]]$data,file = "cluster1.csv", sep =",")
write.table(p[[2]]$data,file = "cluster2.csv", sep =",")
write.table(p[[3]]$data,file = "cluster3.csv", sep =",")
write.table(p[[4]]$data,file = "cluster4.csv", sep =",")
write.table(p[[5]]$data,file = "cluster5.csv", sep =",")



#plot cluster 1:
print(p[[1]])

genes_cluster1 <- table(p[[1]]$data$group)
number_of_genes_cluster1 <- nrow(p[[1]]$data)
number_of_genes_cluster1
unique_number_of_genes_cluster1 <- unique(p[[1]]$data$group)
unique_number_of_genes_cluster1
length(unique_number_of_genes_cluster1)
total_gene_cluster1 <- paste0("N =",length(unique_number_of_genes_cluster1))

#check
cluster1<-read.csv("cluster1.csv", header=T)
which(cluster1$group == "ENSG00000178977.4")


# p1<-(p[[1]] +geom_text(aes(label =total_gene_cluster1),
#                          x=Inf, y=-Inf,size=5,hjust =1,vjust =-49,color = "black"))
p1_with_label <- p[[1]] +
  annotate("text", x = Inf, y = Inf, label = total_gene_cluster1, vjust = 1, hjust = 1, size = 5, color = "black")
print(p1_with_label)
ggsave(filename = "cluster1.pdf",plot = p1_with_label, device = "pdf")
#plot cluster 2:
print(p[[2]])

genes_cluster2 <- table(p[[2]]$data$group)
number_of_genes_cluster2 <- nrow(p[[2]]$data)
number_of_genes_cluster2
unique_number_of_genes_cluster2 <- unique(p[[2]]$data$group)
unique_number_of_genes_cluster2
length(unique_number_of_genes_cluster2)
total_gene_cluster2 <- paste0("N =",length(unique_number_of_genes_cluster2))

#check
cluster2<-read.csv("cluster2.csv", header=T)
which(cluster2$group == "ENSG00000178977.4")


p2_with_label <- p[[2]] +
  annotate("text", x = Inf, y = Inf, label = total_gene_cluster2, vjust = 1, hjust = 1, size = 5, color = "black")
print(p2_with_label)
ggsave(filename = "cluster2.pdf",plot =p2_with_label , device = "pdf")

#plot cluster 3:
print(p[[3]])

genes_cluster3 <- table(p[[3]]$data$group)
number_of_genes_cluster3 <- nrow(p[[3]]$data)
number_of_genes_cluster3
unique_number_of_genes_cluster3 <- unique(p[[3]]$data$group)
unique_number_of_genes_cluster3
length(unique_number_of_genes_cluster3)
total_gene_cluster3 <- paste0("N =",length(unique_number_of_genes_cluster3))

#check
cluster3<-read.csv("cluster3.csv", header=T)
which(cluster3$group == "ENSG00000178977.4")



p3_with_label <- p[[3]] +
  annotate("text", x = Inf, y = Inf, label = total_gene_cluster3, vjust = 1, hjust = 1, size = 5, color = "black")
print(p3_with_label)
ggsave(filename = "cluster3.pdf",plot = p3_with_label, device = "pdf")

#plot cluster 4:
print(p[[4]])

genes_cluster4 <- table(p[[4]]$data$group)
number_of_genes_cluster4 <- nrow(p[[4]]$data)
number_of_genes_cluster4
unique_number_of_genes_cluster4 <- unique(p[[4]]$data$group)
unique_number_of_genes_cluster4
length(unique_number_of_genes_cluster4)
total_gene_cluster4 <- paste0("N = ",length(unique_number_of_genes_cluster4))

#check
cluster4<-read.csv("cluster4.csv", header=T)
which(cluster4$group == "ENSG00000178977.4")


p4_with_label <- p[[4]] +
  annotate("text", x = Inf, y = Inf, label = total_gene_cluster4, vjust = 1, hjust = 1, size = 5, color = "black")

print(p4_with_label)
ggsave(filename = "cluster4.pdf",plot = p4_with_label, device = "pdf")

#plot cluster 5:
print(p[[5]])

genes_cluster5 <- table(p[[5]]$data$group)
number_of_genes_cluster5 <- nrow(p[[5]]$data)
number_of_genes_cluster5
unique_number_of_genes_cluster5 <- unique(p[[5]]$data$group)
unique_number_of_genes_cluster5
length(unique_number_of_genes_cluster5)
total_gene_cluster5 <- paste0("N = ",length(unique_number_of_genes_cluster5))

#check
cluster5<-read.csv("cluster5.csv", header=T)
which(cluster5$group == "ENSG00000178977.4")



p5_with_label <- p[[5]] +
  annotate("text", x = Inf, y = Inf, label = total_gene_cluster5, vjust = 1, hjust = 1, size = 5, color = "black")

print(p5_with_label)
ggsave(filename = "cluster5.pdf",plot = p5_with_label, device = "pdf")


