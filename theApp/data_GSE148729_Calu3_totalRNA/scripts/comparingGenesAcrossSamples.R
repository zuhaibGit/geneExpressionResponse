###########################################
# This script takes the data in long_ files and finds the least squares line for each gene in each sample
# The slopes of these lines can be sued to compares a gene across the three samples
###########################################

setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE148729_Calu3_totalRNA")

mock <- read.table("long_GSE148729_Calu3_mockInfection.txt", sep = "\t", header = T)
cov1 <- read.table("long_GSE148729_Calu3_sarsCov1.txt", sep = "\t", header = T)
cov2 <- read.table("long_GSE148729_Calu3_sarsCov2.txt", sep = "\t", header = T)

genes <- unique(mock$Gene)

# Takes in the DE between time points of some gene, and returns x,y coordinates for the line
# as well as the color of the points based on whether it was significantly expressed.
# Note: Time points must be sorted
makeLine <- function(timePoints) {
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  sig <- c("Not", timePoints$Colour)
  return(data.frame(x = x , y = y, sig = sig))
}

numRow_cov1 <- nrow(cov1) / 2
piecewise_cov1 <- lapply(1:numRow_cov1, function(x) {
  return(makeLine(cov1[c(x, x + numRow_cov1),]))
})
names(piecewise_cov1) <- cov1$Gene[1:numRow_cov1]

numRow_cov2 <- nrow(cov2) / 2
piecewise_cov2 <- lapply(1:numRow_cov2, function(x) {
  return(makeLine(cov2[c(x, x + numRow_cov2),]))
})
names(piecewise_cov2) <- cov2$Gene[1:numRow_cov2]

numRow_mock <- nrow(mock)
piecewise_mock <- lapply(1:numRow_mock, function(x) {
  return(makeLine(mock[x,]))
})
names(piecewise_mock) <- mock$Gene[1:numRow_mock]


randGene <- sample(genes, 1)
