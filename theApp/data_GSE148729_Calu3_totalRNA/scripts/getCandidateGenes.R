setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE148729_Calu3_totalRNA")

mock <- read.table("long_GSE148729_Calu3_mockInfection.txt", sep = "\t", header = T)
cov1 <- read.table("long_GSE148729_Calu3_sarsCov1.txt", sep = "\t", header = T)
cov2 <- read.table("long_GSE148729_Calu3_sarsCov2.txt", sep = "\t", header = T)

lines <- read.table("leastSquaresLinesPerGenePerInfections.txt", header = T, sep = )
lines <- lines[,c(1,3,5)]

genes <- unique(mock$Gene)


##### DETERMINING WHICH GENES HAVE STATISTICALLY SIGNIFICANT LOG FOLD CHANGES #####
# If any gene has a statistically significant log fold change over the course of the experiment, then it's counted as significant

# Mock
sig_mock <- mock$Colour

# Cov1 genes
numRow_cov1 <- nrow(cov1) / 2
sig_cov1 <- sapply(1:numRow_cov1, function(x) {
  if ("Significant" %in% cov1[c(x, x + numRow_cov1),5]) {
    return("Significant") 
  } else {
    return("Not")
  }
})

# Cov2 genes
numRow_cov2 <- nrow(cov2) / 2
sig_cov2 <- sapply(1:numRow_cov2, function(x) {
  if ("Significant" %in% cov2[c(x, x + numRow_cov2),5]) {
    return("Significant") 
  } else {
    return("Not")
  }
})


###### DETERMINING WHICH SAMPLE HAS A DIFFERENT TRAJECTORY PER GENE #####
# Using the least squares lines for each gene and each sample, any sample in which a gene has a 
oddSample <- t(apply(lines, 1, function(x) {
  slopes <- sort(x)
  high.mid <- atan(slopes[3]) - atan(slopes[2])
  mid.low <- atan(slopes[2]) - atan(slopes[1])
  if (high.mid >= 5 * mid.low) {
    return(c(oddOne = names(slopes)[3], Ratio = high.mid / mid.low))
  } else if (mid.low >= 5 * high.mid) {
    return(c(oddOne = names(slopes[1]), Ratio = mid.low / high.mid))
  } else {
    return(c(oddOne = "None", Ratio = NA))
  }
}))
oddSample <- as.data.frame(oddSample)
oddSample[,2] <- as.numeric(oddSample[,2])
names(oddSample) <- c("oddSample", "Ratio")

lines <- cbind(lines, oddSample)

lines <- cbind(lines, sig_mock, sig_cov1, sig_cov2)


oddCov2 <- lines[which(lines$oddSample == "cov2_slope"),]
sigOddCov2 <- oddCov2[which(sapply(1:nrow(oddCov2), function(x) {
  if ("Significant" %in% oddCov2[x,6:8]) {
    return(T)
  } else {
    return(F)
  }
})),]
sigOddCov2 <- sigOddCov2[order(sigOddCov2$Ratio, decreasing = T),]

write.table(sigOddCov2, "candidateGenes_coreCov2.txt", row.names = T, col.names = T, quote = F, sep = "\t")

oddcov1 <- lines[which(lines$oddSample == "cov1_slope"),]
sigOddcov1 <- oddcov1[which(sapply(1:nrow(oddcov1), function(x) {
  if ("Significant" %in% oddcov1[x,6:8]) {
    return(T)
  } else {
    return(F)
  }
})),]
sigOddcov1 <- sigOddcov1[order(sigOddcov1$Ratio, decreasing = T),]

write.table(sigOddcov1, "candidateGenes_coreCov1.txt", row.names = T, col.names = T, quote = F, sep = "\t")

oddmock <- lines[which(lines$oddSample == "mock_slope"),]
sigOddmock <- oddmock[which(sapply(1:nrow(oddmock), function(x) {
  if ("Significant" %in% oddmock[x,6:8]) {
    return(T)
  } else {
    return(F)
  }
})),]
sigOddmock <- sigOddmock[order(sigOddmock$Ratio, decreasing = T),]

write.table(sigOddmock, "candidateGenes_panCov1-2.txt", row.names = T, col.names = T, quote = F, sep = "\t")
