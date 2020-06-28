
setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE56677_Calu3/")

# Load long_ files
mers <- read.table("./long_GSE56677_Calu3_mersCov.txt", header = T, sep = "\t")
mock <- read.table("./long_GSE56677_Calu3_mock.txt", header = T, sep = "\t")

# The long_files show each gene once for each time point. I want to determine how many rows are between each of the time points for a gene
vec <- as.integer(row.names(mers[grep("A_33_P3245660", mers[,1]),]))
num <- c(vec[-1], vec[5]) - vec
num <- num[1]

# Calculate the equaiton of the least squares line for genes in mers data
lr_mers <- lapply(1:num, function(g) {
  timePoints <- mers[c(g, g + num, g + 2*num, g + 3*num, g + 4*num),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mers) <- mers[1:num,1]
lr_mers <- do.call(rbind, lr_mers)
colnames(lr_mers) <- c("mers_slope", "mers_yInt")

# Calculate the equaiton of the least squares line for genes in mock infection data
lr_mock <- lapply(1:num, function(g) {
  timePoints <- mock[c(g, g + num, g + 2*num, g + 3*num, g + 4*num),]
  minTimePoint <- timePoints[1,2]
  x <- c(minTimePoint, timePoints$T2)
  y <- c(0, cumsum(timePoints$log2FoldChange))
  theModel <- lm(y~x)
  return(c(slope = theModel$coefficients[[2]], yInt = theModel$coefficients[[1]]))
})
names(lr_mock) <- mock[1:num,1]
lr_mock <- do.call(rbind, lr_mock)
colnames(lr_mock) <- c("mock2_slope", "mock2_yInt")

# Joing them together
retDF <- cbind(lr_mers, lr_mock)

# Write to file
write.table(retDF, "leastSquaresLinesPerGenePerInfections.txt", row.names = T, col.names = T, quote = F, sep = "\t")

