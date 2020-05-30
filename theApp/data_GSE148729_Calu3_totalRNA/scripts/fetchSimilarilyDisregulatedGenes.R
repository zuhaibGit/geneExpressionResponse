#####################################################
# This script will take the candidate genes and fetch those which have similar trajectories in the Cov1- and Cov2-infected
# cells
#####################################################

setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/data_GSE148729_Calu3_totalRNA")

candCov1 <- read.table("candidateGenes_coreCov1.txt", header = T, sep = "\t")
candCov2 <- read.table("candidateGenes_coreCov2.txt", header = T, sep = "\t")
candPan <- read.table("candidateGenes_panCov1-2.txt", header = T, sep = "\t")

theData <- read.table("GSE148729_Calu3_totalRNA_normalizedCounts.txt", header = T, sep = "\t")
row.names(theData) <- theData$geneID
theData <- theData[,-1]

candCov1$maxCount <- apply(theData[row.names(candCov1),], 1, max)
candCov2$maxCount <- apply(theData[row.names(candCov2),], 1, max)
candPan$maxCount <- apply(theData[row.names(candPan),], 1, max)

# Obtain genes that are regulated differently in the pan Cov infection then in the mock
sapply(row.names(candPan[which(apply(candPan, 1, function(x) {
  one <- as.numeric(x[1])
  two <- as.numeric(x[2])
  three <- as.numeric(x[3])
  if (one <= 0 && (two >= 0 && three >= 0)) {
    if (x[7] == "Significant" && x[8] == "Significant") {
      return(T)
    } else {
      return(F)
    }
  } else {
    return(F)
  }
})),]), print)



# Obtain genes that are differently regulated in Cov2 than in Cov1 and mock
candCov2[which(apply(candCov2, 1, function(x) {
  one <- as.numeric(x[1])
  two <- as.numeric(x[2])
  three <- as.numeric(x[3])
  if (three <= 0 && (two >= 0 && one >= 0)) {
    #return(T)
    if (x[8] == "Significant") {
      return(T)
    } else {
      return(F)
    }
  } else {
    return(F)
  }
})),]


# Obtain genes that are differently regulated in Cov1 than in Cov2 and mock
candCov1[which(apply(candCov1, 1, function(x) {
  one <- as.numeric(x[1])
  two <- as.numeric(x[2])
  three <- as.numeric(x[3])
  if (two >= 0 && (three <= 0 && one <= 0)) {
    #return(T)
    if (x[7] == "Significant") {
      return(T)
    } else {
      return(F)
    }
  } else {
    return(F)
  }
})),]
