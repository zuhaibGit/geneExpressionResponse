library(biomaRt)

setwd("/home/zuhaib/Desktop/covid19Research/geneResponseApp/geneExpressionResponse/theApp/geneMappings")
data <- read.table("../data_GSE148729_Calu3_totalRNA/normalized_GSE148729.txt", header = T, sep = "\t")
genes <- data[,1]
genes <- sapply(strsplit(data$geneID, "\\."), function(x) return(x[[1]]))

# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

filters <- listFilters(mart)
attrs <- listAttributes(mart)
res <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), filters = "ensembl_gene_id", values = sapply(strsplit(data$geneID, "\\."), function(x) return(x[[1]])), mart = mart)

# These genes weren't found in biomaRt
nfGenes <- genes[which(!(genes %in% res[,1]))]
res2 <-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), filters = "ensembl_gene_id", values = genes[which(!(genes %in% res[,1]))], mart = mart) 

# Add those genes in
adf <- as.data.frame(do.call(rbind, lapply(nfGenes, function(x) return(c(ensembl_gene_id = x, hgnc_symbol = "", entrezgene_id = "")))))
res <- rbind(res, adf)

write.table(res, "ensembl_mappings.txt", row.names = F, col.names = T, sep = "\t", quote = F)
