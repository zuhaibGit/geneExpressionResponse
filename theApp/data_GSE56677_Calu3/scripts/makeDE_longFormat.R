######################################
# This script will take in the DE files anc convert them into long format
# See README.md for details on the format
######################################
library(ggplot2)
setwd("/home/zuhaib/Desktop/covid19Research/microarrayData/GSE56677/DE")

fls <- list.files()[grep("m", list.files())]
dfs <- lapply(fls, function(x) return(read.table(x, header = T, sep = "\t")[,c(1,2,6,8)]))
names(dfs) <- fls

## DEALING WITH DUPLICATE PROBE ENTRIES
# For duplicate entries, we're averaging the LFC and adjusted p-values.
# Issues with this will be mitigated in the app, as the data will be fully
# displayed if the user selects a duplicated probeset.
dfs <- lapply(dfs, function(adf) {
  aggLFC <- aggregate(adf["logFC"], by = adf["ID"], mean)
  aggP <- aggregate(adf["adj.P.Val"], by = adf["ID"], mean)
  aggColour <- sapply(aggP[,2], function(x) {
    if (x < 0.05) {
      return("Significant")
    } else {
      return("Not")
    }
  })
  agg <- merge(aggLFC, aggP, by = "ID")
  agg$Colour <- aggColour
  return(agg)
})

probes <- dfs[[1]]$ID

setwd("../")

##### BUILDING THE MERS MATRIX #####
### 0H VS 3H DATA ###
mers.0h.mers.3h <- data.frame(ID = probes, T1 = 0, T2 = 3)
mers.0h.mers.3h <- merge(mers.0h.mers.3h, dfs$`mers-0h-mers-3h.txt`, all = T)
### 3H VS 7H DATA ###
mers.3h.mers.7h <- data.frame(ID = probes, T1 = 3, T2 = 7)
mers.3h.mers.7h <- merge(mers.3h.mers.7h, dfs$`mers-3h-mers-7h.txt`, all = T)
### 7H VS 12H DATA ###
mers.7h.mers.12h <- data.frame(ID = probes, T1 = 7, T2 = 12)
mers.7h.mers.12h <- merge(mers.7h.mers.12h, dfs$`mers-7h-mers-12h.txt`, all = T)
### 12H VS 18H DATA ###
mers.12h.mers.18h <- data.frame(ID = probes, T1 = 12, T2 = 18)
mers.12h.mers.18h <- merge(mers.12h.mers.18h, dfs$`mers-12h-mers-18h.txt`, all = T)
### 18H VS 24H DATA ###
mers.18h.mers.24h <- data.frame(ID = probes, T1 = 18, T2 = 24)
mers.18h.mers.24h <- merge(mers.18h.mers.24h, dfs$`mers-18h-mers-24h.txt`, all = T)

mers <- do.call(rbind, list(mers.0h.mers.3h, mers.3h.mers.7h, mers.7h.mers.12h, mers.12h.mers.18h, mers.18h.mers.24h))
write.table(mers, "long_GSE56677_Calu3_mersCov.txt", row.names = F, col.names = T, quote = F, sep = "\t")

##### BUILDING THE MOCK MATRIX #####
### 0H VS 3H DATA ###
mock.0h.mock.3h <- data.frame(ID = probes, T1 = 0, T2 = 3)
mock.0h.mock.3h <- merge(mock.0h.mock.3h, dfs$`mock-0h-mock-3h.txt`, all = T)
### 3H VS 7H DATA ###
mock.3h.mock.7h <- data.frame(ID = probes, T1 = 3, T2 = 7)
mock.3h.mock.7h <- merge(mock.3h.mock.7h, dfs$`mock-3h-mock-7h.txt`, all = T)
### 7H VS 12H DATA ###
mock.7h.mock.12h <- data.frame(ID = probes, T1 = 7, T2 = 12)
mock.7h.mock.12h <- merge(mock.7h.mock.12h, dfs$`mock-7h-mock-12h.txt`, all = T)
### 12H VS 18H DATA ###
mock.12h.mock.18h <- data.frame(ID = probes, T1 = 12, T2 = 18)
mock.12h.mock.18h <- merge(mock.12h.mock.18h, dfs$`mock-12h-mock-18h.txt`, all = T)

mock <- do.call(rbind, list(mock.0h.mock.3h, mock.3h.mock.7h, mock.7h.mock.12h, mock.12h.mock.18h))
write.table(mock, "long_GSE56677_Calu3_mock.txt", row.names = F, col.names = T, quote = F, sep = "\t")
