################################################
# THIS SCRIPT TAKES THE AGILENT GENE MAPPINGS, AND ADDS THE ENSEMBL GENE IDs IN ADDITION TO THE
# ENSEMBL TRANSCRIPT IDs.
################################################

library(biomaRt)

setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/geneMappings")
data <- read.table("agilentHumanG3.txt", header = T, sep = "\t", quote = "")


# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# query biomart
results <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "entrezgene_id", "uniprot_gn_symbol"),
                 filters = "ensembl_transcript_id", values = data$ENSEMBL_ID,
                 mart = mart)
results

res <- as.list(results[,2])
names(res) <- results[,1]
data$ENSEMBL_GENE_ID <- sapply(data$ENSEMBL_ID, function(x) {
  return(res[[x]])
})
data$ENSEMBL_GENE_ID <- as.character(data$ENSEMBL_GENE_ID)

res2 <- as.list(results[,3])
names(res2) <- results[,1]
data$ENTREZ_ID <- sapply(data$ENSEMBL_ID, function(x) {
  return(res2[[x]])
})
data$ENTREZ_ID <- as.character(data$ENTREZ_ID)

res3 <- as.list(results[,4])
names(res3) <- results[,1] 
data$UNIPROT_ID <- sapply(data$ENSEMBL_ID, function(x) {
  return(res3[[x]])
})
data$UNIPROT_ID <- as.character(data$UNIPROT_ID)

retDF <- apply(data, 2, function(x) {
  x[which(x == "")] <- NA
  x[which(is.na(x))] <- NA
  x[which(x == "NULL")] <- NA
  return(x)
})

write.table(retDF, "agilentHumanG3_mappings.txt", col.names = T, row.names = F, sep = "\t", quote = F)


attrs <- listAttributes(mart)
