library(genefilter)
setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp")

# Get the long_ files
longCov1 <- read.table("./data_GSE148729_Calu3_totalRNA/long_GSE148729_Calu3_sarsCov1.txt", header = T)
longCov2 <- read.table("./data_GSE148729_Calu3_totalRNA/long_GSE148729_Calu3_sarsCov2.txt", header = T)
longMock1 <- read.table("./data_GSE148729_Calu3_totalRNA/long_GSE148729_Calu3_mockInfection.txt", header = T)
longMers <- read.table("./data_GSE56677_Calu3/long_GSE56677_Calu3_mersCov.txt", header = T)
longMock2 <- read.table("./data_GSE56677_Calu3/long_GSE56677_Calu3_mock.txt", header = T)

lines1 <- read.table("./data_GSE148729_Calu3_totalRNA/leastSquaresLinesPerGenePerInfections.txt")
row.names(lines1) <- sapply(strsplit(row.names(lines1), "\\."), function(x) return(x[1]))
lines2 <- read.table("./data_GSE56677_Calu3/leastSquaresLinesPerGenePerInfections.txt")
maps <- read.table("./geneMappings/mappings_agilentHumanG3.txt", header = T, sep = "\t", quote = "")

# Map the microarray IDs to Ensembl IDs
newRowNames <- sapply(row.names(lines2), function(x) {
  return(maps$ENSEMBL_GENE_ID[which(maps$ID == x)])
})
mappings <- newRowNames[which(!is.na(newRowNames))]
mappings <- data.frame(Gene = newLines2)
newLines2 <- merge(lines2, mappings, by = 0, all.y = T)
newLines2 <- newLines2[,-1]

# Combine the ensembl and (mapped) microarray data by the gene ID column
df <- merge(lines1, newLines2, by.y = "Gene", by.x = 0, all.y = T)
# We don't need the y-intercepts, just the slopes.
df <- df[,c(1,2,4,6,8,10)]
# Remove the NAs
whichNAs <- which(is.na(df[,2]))
df <- df[-whichNAs,]


# Determine min and max slopes to determine if we should use slopes or angles
aDF <- df[,-1]
maxSlope <- max(unlist(aDF))
minSlope <- min(unlist(aDF))
# We want to determine if we should convert these slopes to angles-to-the-x-axis. As the angle theta is related to the slope m, by 
# theta = arctan(m). Because of the small angle approximation, the arctan function in the (minSlope, maxSlope) range is very similar to the linear function, except
# near the end of the range, where large slopes have closer angle values. So the angles will be harder to differentiate near the large-end of the (minSloope, maxSlope0
# range. So we'll opt to use slopes instead of angles)

#####################################################################
# First, identify the genes where the Cov cases have similar expression and the mock infections have similar expression
# We'll do t-tests on each gene to see which genes have different expression patterns in Cov than in mock infections
data <- df
groups <- factor(c("A", "B", "B", "B", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[4]], x[[5]]))) > mean(as.numeric(c(x[[2]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
  })),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[4]], x[[5]]))) <= mean(as.numeric(c(x[[2]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "G_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "G_down.txt", row.names = F, col.names = F, quote = F)




# Now we'll look at which genes have different expression in CoV1 and mers compared to CoV2 and mock infections
data <- df
groups <- factor(c("A", "B", "A", "B", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[5]]))) > mean(as.numeric(c(x[[2]], x[[4]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[5]]))) <= mean(as.numeric(c(x[[2]], x[[4]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "E_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "E_down.txt", row.names = F, col.names = F, quote = F)




# Now we'll look compare the group of Cov2 and mers to Cov1 and mock infections
data <- df
groups <- factor(c("A", "A", "B", "B", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[4]], x[[5]]))) > mean(as.numeric(c(x[[2]], x[[3]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[4]], x[[5]]))) <= mean(as.numeric(c(x[[2]], x[[3]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "F_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "F_down.txt", row.names = F, col.names = F, quote = F)





# Now, Cov1 and Cov2 compared to mers and mock infections
data <- df
groups <- factor(c("A", "B", "B", "A", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[4]]))) > mean(as.numeric(c(x[[2]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (mean(as.numeric(c(x[[3]], x[[4]]))) <= mean(as.numeric(c(x[[2]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "D_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "D_down.txt", row.names = F, col.names = F, quote = F)








# Comparing Cov1 to everything else
data <- df
groups <- factor(c("A", "B", "A", "A", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[3]]) > mean(as.numeric(c(x[[2]], x[[4]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[3]]) <= mean(as.numeric(c(x[[2]], x[[4]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "A_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "A_down.txt", row.names = F, col.names = F, quote = F)








# Comparing Cov2 to everything else
data <- df
groups <- factor(c("A", "A", "B", "A", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[4]]) > mean(as.numeric(c(x[[2]], x[[3]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[4]]) <= mean(as.numeric(c(x[[2]], x[[3]], x[[5]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "B_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "B_down.txt", row.names = F, col.names = F, quote = F)









# Comparing Cov2 to everything else
data <- df
groups <- factor(c("A", "A", "A", "B", "A"))
res <- rowttests(as.matrix(aDF), factor(groups))
data$p.value <- res$p.value
data <- data[order(data$p.value),]

# Filter out the genes where the p-value of <= 0.05
dfSig <- data[which(data$p.value <= 0.05),]

# To filter out the important significant genes, filter by genes that have a "Significant" change somewhere in their timelines.
dfSig$sig_cov1 <- apply(dfSig, 1, function(x) {
  a <- longCov1[grep(x[1], longCov1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_cov2 <- apply(dfSig, 1, function(x) {
  a <- longCov2[grep(x[1], longCov2$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSig$sig_mock <- apply(dfSig, 1, function(x) {
  a <- longMock1[grep(x[1], longMock1$Gene),]
  if ("Significant" %in% a$Colour) {
    return(T)
  } else {
    return(F)
  }
})
dfSigFiltered <- dfSig[which(apply(dfSig, 1, function(x) if (x[8]==T || x[9]==T || x[10]==T) return(T) else return(F))),]

# Separate the filtered genes by which ones are up-regulated in the CoVs compared to the mock
dfSigFilteredUp <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[5]]) > mean(as.numeric(c(x[[2]], x[[3]], x[[4]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

# Separate the filtered genes by which ones are down-regulated in the CoVs compared to the mock
dfSigFilteredDown <- dfSigFiltered[which(apply(dfSigFiltered, 1, function(x) { 
  #print(mean(x[[3]], x[[4]], x[[5]]))
  if (as.numeric(x[[5]]) <= mean(as.numeric(c(x[[2]], x[[3]], x[[4]], x[[6]])))) {
    return(T)
  } else {
    return(F)
  }
})),]

write.table(dfSigFilteredUp[,1], "C_up.txt", row.names = F, col.names = F, quote = F)
write.table(dfSigFilteredDown[,1], "C_down.txt", row.names = F, col.names = F, quote = F)
