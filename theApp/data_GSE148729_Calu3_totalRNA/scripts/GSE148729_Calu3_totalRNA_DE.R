###############################################
# This script will perform DE analysis on the samples by comparing gene expression in adjacent time points in each group
# mock infection, cov1 infection, and cov2 infection
###############################################
library(DESeq2)
library(ggplot2)
library(stringr)
library(factoextra)
setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE148729_Calu3_totalRNA")

# Load data and convert it to a numeric matrix
data <- read.table("GSE148729_Calu3_totalRNA_readcounts.txt", header = T, sep = "\t")
row.names(data) <- data[,1]
data <- data[,-1]
data <- as.matrix(data)
data <- data[which(rowSums(data) > 0),]
mode(data) <- "integer"

# Extract metadata from the column names
listMeta <- strsplit(colnames(data), "\\.")
meta <- data.frame(Sample = colnames(data),
                   Cell = sapply(listMeta, function(x) return(x[1])),
                   Infection = sapply(listMeta, function(x) return(x[2])),
                   Time = sapply(listMeta, function(x) return(x[3])),
                   Replicate = sapply(listMeta, function(x) return(x[4])))
row.names(meta) <- meta[,1]
meta <- meta[,-1]
#write.table(meta, "meta.txt", sep = "\t", row.names = T, col.names = T, quote = F)

##### EXPLORATORY ANALYSIS ON DATA #####
# PCA before normalization
preNormalizePCA <- prcomp(t(data))
autoplot(preNormalizePCA, data = meta, colour = "Replicate", size = "Time", shape = "Infection")
# PCA after normalization
dataVST <- varianceStabilizingTransformation(data)
postNormalizePCA <- prcomp(t(dataVST))
autoplot(postNormalizePCA, data = meta, colour = "Replicate", size = "Time", shape = "Infection")

setwd("./DE")

##### COMPARISONS #####
##  DE S1 4H VS 12H ##
# Select the approppriate data
data1 <- data[,1:4]
meta1 <- meta[1:4,]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$log2FoldChange)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res1 <- res
res1$log2FoldChange <- (-1) * res1$log2FoldChange
ggplot(res1, aes(x = log2FoldChange, y = -log(padj))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res1, "s1-4h-s1-12h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S1 12H VS 24H ##
# Select the approppriate data
data1 <- data[,3:6]
meta1 <- meta[3:6,]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$log2FoldChange)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res2 <- res
ggplot(res2, aes(x = log2FoldChange, y = -log(padj))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res2, "s1-12h-s1-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)


##  DE S2 4H VS 12H ##
# Select the approppriate data
data1 <- data[,7:10]
meta1 <- meta[7:10,]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$log2FoldChange)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res3 <- res
res3$log2FoldChange <- (-1) * res3$log2FoldChange
ggplot(res3, aes(x = log2FoldChange, y = -log(padj))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res3, "s2-4h-s2-12h.txt", row.names = T, col.names = T, sep = "\t", quote = F)


##  DE S2 12H VS 24H ##
# Select the approppriate data
data1 <- data[,9:12]
meta1 <- meta[9:12,]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$log2FoldChange)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res4 <- res
ggplot(res4, aes(x = log2FoldChange, y = -log(padj))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res4, "s2-12h-s2-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)


##  DE S1 4H VS MOCK 4H ##
# Select the approppriate data
data1 <- data[,c(1,2,13,14)]
meta1 <- meta[c(1,2,13,14),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Infection)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res12 <- res
ggplot(res12, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res12, "mock-4h-s1-4h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S2 4H VS MOCK 4H ##
# Select the approppriate data
data1 <- data[,c(7,8,13,14)]
meta1 <- meta[c(7,8,13,14),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Infection)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res13 <- res
ggplot(res13, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res13, "mock-4h-s2-4h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S1 24H VS MOCK 24H ##
# Select the approppriate data
data1 <- data[,c(5,6,15,16)]
meta1 <- meta[c(5,6,15,16),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Infection)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res14 <- res
ggplot(res14, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res14, "mock-24h-s1-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S2 24H VS MOCK 24H ##
# Select the approppriate data
data1 <- data[,c(11,12,15,16)]
meta1 <- meta[c(11,12,15,16),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Infection)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res15 <- res
ggplot(res15, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res15, "mock-24h-s2-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S1 4H VS 24H ##
# Select the approppriate data
data1 <- data[,c(1,2,5,6)]
meta1 <- meta[c(1,2,5,6),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res16 <- res
res16$log2FoldChange <- res16$log2FoldChange * (-1)
ggplot(res16, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res16, "s1-4h-s1-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE S2 4H VS 24H ##
# Select the approppriate data
data1 <- data[,c(7,8,11,12)]
meta1 <- meta[c(7,8,11,12),]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$pvalue)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res17 <- res
res17$log2FoldChange <- res17$log2FoldChange * (-1)
ggplot(res17, aes(x = log2FoldChange, y = -log(pvalue))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res17, "s2-4h-s2-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)

##  DE MOCK 4H VS 24H ##
# Select the approppriate data
data1 <- data[,13:16]
meta1 <- meta[13:16,]
# DE analysis
obj <- DESeqDataSetFromMatrix(countData = data1, colData = meta1, design = ~ Time)
obj <- DESeq(obj)
res <- results(obj)
res <- res[which(!is.na(res$log2FoldChange)),]
res <- as.data.frame(res)
res <- res[order(res$log2FoldChange, decreasing = T),]
res$Colour <- sapply(res$padj, function(x) if (is.na(x) | x > 0.05) return("Not") else return("Significant"))
res18 <- res
res18$log2FoldChange <- res18$log2FoldChange * (-1)
ggplot(res18, aes(x = log2FoldChange, y = -log(padj))) + geom_point(aes(color = Colour)) + xlim(-14, 8) + ylim(0, 50)
write.table(res18, "mock-4h-mock-24h.txt", row.names = T, col.names = T, sep = "\t", quote = F)
