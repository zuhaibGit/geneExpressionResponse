##########
# THIS SCRIPT WILL TAKE THE RAW COUNTS AND DO A UANTILE NORMALIZATION, TO CORRECT FOR DIFFERING LIBRARY SIZES
# NOTE: THIS ISN'T THE SAME NORMALIZATION DONE DURING THE DE ANALYSIS. THIS IS JUST TO ALLOW US TO BE ABLE TO
# ROUGHPLY COMPARE GENE COUNTS
##########

setwd("/home/zuhaib/Desktop/covid19Research/hackSeqRNA/Pan-Coronavirus-Gene-Regulatory-Networks/theApp/data_GSE148729_Calu3_totalRNA")

data <- read.table("GSE148729_Calu3_totalRNA_readcounts.txt", header = T, sep = "\t")
row.names(data) <- data$gene_id
data <- data[,-1]

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

normData <- as.data.frame(quantile_normalisation(data))

geneID <- row.names(normData)
normData <- cbind(geneID, normData)

write.table(format(normData, digits=1, scientific=F), "GSE148729_Calu3_totalRNA_normalizedCounts.txt", row.names = F, col.names = T, sep = "\t", quote = F)
