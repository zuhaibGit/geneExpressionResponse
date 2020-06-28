library(pheatmap)
library(RColorBrewer)
setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE56677_Calu3/")

# Upload gene mappings so that we can map the Ensembl IDs we found and see how they are in the MERS data
maps <- read.table("../geneMappings/mappings_agilentHumanG3.txt", header = T, sep = "\t")
# Upload long_ data files
longMers <- read.table("./long_GSE56677_Calu3_mersCov.txt", header = T, sep = "\t")
longMock <- read.table("./long_GSE56677_Calu3_mock.txt", header = T, sep = "\t")

#################### Finding genes up-regulated in Mers AND Cov1 and Cov2, but down-regulated in the mock infections ######################
upInCov1_2 <- read.table("../data_GSE148729_Calu3_totalRNA/results/upInCov1_2.txt")
upInCov1_2 <- upInCov1_2[,1]
upInCov1_2 <- sapply(strsplit(upInCov1_2, "\\."), function(x) return(x[1]))

subsetMers <- longMers[which(longMers$Gene %in% maps$ID[which(maps$ENSEMBL_GENE_ID %in% upInCov1_2)]),]
subsetMock <- longMock[which(longMock$Gene %in% maps$ID[which(maps$ENSEMBL_GENE_ID %in% upInCov1_2)]),]

# Get lines for each gene in mock infection group
lr_mock <- lapply(unique(subsetMock[,1]), function(g) {
  timePoints <- subsetMock[grep(g, subsetMock$Gene),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mock) <- unique(subsetMock[,1])
lr_mock <- do.call(rbind, lr_mock)
colnames(lr_mock) <- c("mock_slope", "mock_yInt")

# Get lines for each gene in mers infection group
lr_mers <- lapply(unique(subsetMers[,1]), function(g) {
  timePoints <- subsetMers[grep(g, subsetMers$Gene),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mers) <- unique(subsetMers[,1])
lr_mers <- do.call(rbind, lr_mers)
colnames(lr_mers) <- c("mers_slope", "mers_yInt")

df <- cbind(lr_mers, lr_mock)

upInAllCov <- maps$ENSEMBL_GENE_ID[which(maps$ID %in% names(which(apply(df, 1, function(x) {
  if (x[1] >= 0 && x[3] < 0) {
    return(T)
  } else {
    return(F)
  }
}))))]

write.table(upInAllCov, "upInAllCov.txt", row.names = F, col.names = F, quote = F)





############### Finding genes that are down-regulated in mers AND Cov1 and Cov2, and up-regulated in mock infections
downInCov1_2 <- read.table("../data_GSE148729_Calu3_totalRNA/results/downInCov1_2.txt")
downInCov1_2 <- downInCov1_2[,1]
downInCov1_2 <- sapply(strsplit(downInCov1_2, "\\."), function(x) return(x[1]))

subsetMers <- longMers[which(longMers$Gene %in% maps$ID[which(maps$ENSEMBL_GENE_ID %in% upInCov2)]),]
subsetMock <- longMock[which(longMock$Gene %in% maps$ID[which(maps$ENSEMBL_GENE_ID %in% upInCov2)]),]



# Get lines for each gene in mock infection group
lr_mock <- lapply(unique(subsetMock[,1]), function(g) {
  timePoints <- subsetMock[grep(g, subsetMock$Gene),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mock) <- unique(subsetMock[,1])
lr_mock <- do.call(rbind, lr_mock)
colnames(lr_mock) <- c("mock_slope", "mock_yInt")

# Get lines for each gene in mers infection group
lr_mers <- lapply(unique(subsetMers[,1]), function(g) {
  timePoints <- subsetMers[grep(g, subsetMers$Gene),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mers) <- unique(subsetMers[,1])
lr_mers <- do.call(rbind, lr_mers)
colnames(lr_mers) <- c("mers_slope", "mers_yInt")

df <- cbind(lr_mers, lr_mock)

downInAll <- maps$ENSEMBL_GENE_ID[which(maps$ID %in% names(which(apply(df, 1, function(x) {
  if (x[1] <= 0 && x[3] > 0) {
    return(T)
  } else {
    return(F)
  }
}))))]
