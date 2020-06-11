############################################################
# THIS SCRIPT WILL DO THE DE ANALYSIS ON THE PROBESTS. IT WILL ALSO QUANTILE NORMALIZE THE DATA
############################################################

library(limma)
library(ggplot2)
library(reshape)
setwd("/home/zuhaib/Desktop/covid19Research/microarrayData/GSE56677/arrayFiles")
filesToRead <- list.files()[grep("GSM", list.files())]

data <- read.maimages(filesToRead, source="agilent", green.only=TRUE)
exp <- data$E
row.names(exp) <- data$genes$ProbeName
setwd("../")
meta <- read.table("meta.txt", header = T, sep = "\t")
colnames(exp) <- sapply(strsplit(colnames(exp), "_"), function(x) return(x[1]))
colnames(exp) <- sapply(colnames(exp), function(x) {
  return(meta$Sample_title[which(meta$Sample_geo_accession == x)])
})

########### QUALITY CONTROL BEFORE CLEANING ###########
# logExp <- as.data.frame(log2(exp))
# names(logExp) <- sapply(strsplit(names(logExp), "_"), function(x) return(x[1]))
# longLogExp <- melt(logExp)
# ggplot(longLogExp, aes(value)) + geom_density(aes(color=variable))
# rm(logExp)
# rm(longLogExp)

########## BACKGROUND CORRECTION WITH NORMEXP ############
exp <- backgroundCorrect(exp, method = "normexp")
# logExp <- as.data.frame(log2(exp))
# names(logExp) <- sapply(strsplit(names(logExp), "_"), function(x) return(x[1]))
# longLogExp <- melt(logExp)
# ggplot(longLogExp, aes(value)) + geom_density(aes(color=variable))
# rm(logExp)
# rm(longLogExp)

############ NORMALIZATION #############
exp <- normalizeQuantiles(exp)
# logExp <- as.data.frame(log2(exp))
# names(logExp) <- sapply(strsplit(names(logExp), "_"), function(x) return(x[1]))
# longLogExp <- melt(logExp)
# ggplot(longLogExp, aes(value)) + geom_density(aes(color=variable))
# rm(logExp)
# rm(longLogExp)

#write.table(exp, "normalizedExpression.txt", row.names = T, col.names = T, quote = F, sep = "\t")

setwd("./DE")

# LoCov 0h vs 3h
meta1 <- meta[intersect(union(which(meta$Time == "0h"), which(meta$Time == "3h")), which(meta$Infection == "LoCoV")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time3h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a1 <- a
#write.table(a1, "mers-0h-mers-3h.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# LoCov 3h vs 7h
meta1 <- meta[intersect(union(which(meta$Time == "3h"), which(meta$Time == "7h")), which(meta$Infection == "LoCoV")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time7h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a2 <- a
#write.table(a2, "mers-3h-mers-7h.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# LoCov 7h vs 12h
meta1 <- meta[intersect(union(which(meta$Time == "7h"), which(meta$Time == "12h")), which(meta$Infection == "LoCoV")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time7h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a3 <- a
a3$logFC <- (-1) * a3$logFC
#write.table(a3, "mers-7h-mers-12h.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# LoCov 12h vs 18h
meta1 <- meta[intersect(union(which(meta$Time == "12h"), which(meta$Time == "18h")), which(meta$Infection == "LoCoV")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time18h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a4 <- a
#write.table(a4, "mers-12h-mers-18h.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# LoCov 18h vs 24h
meta1 <- meta[intersect(union(which(meta$Time == "18h"), which(meta$Time == "24h")), which(meta$Infection == "LoCoV")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time24h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a5 <- a
#write.table(a5, "mers-18h-mers-24h.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Mock 0h vs 3h
meta1 <- meta[intersect(union(which(meta$Time == "0h"), which(meta$Time == "3h")), which(meta$Infection == "Mock")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time3h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a6 <- a
#write.table(a6, "mock-0h-mock-3h.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# Mock 3h vs 7h
meta1 <- meta[intersect(union(which(meta$Time == "3h"), which(meta$Time == "7h")), which(meta$Infection == "Mock")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time7h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a7 <- a
#write.table(a7, "mock-3h-mock-7h.txt", sep = "\t", row.names = F, col.names = T, quote = F)

# Mock 7h vs 12h
meta1 <- meta[intersect(union(which(meta$Time == "7h"), which(meta$Time == "12h")), which(meta$Infection == "Mock")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time7h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a8 <- a
a8$logFC <- (-1) * a8$logFC
#write.table(a8, "mock-7h-mock-12h.txt", sep = "\t", row.names = F, col.names = T, quote = F)


# Mock 12h vs 18h
meta1 <- meta[intersect(union(which(meta$Time == "12h"), which(meta$Time == "18h")), which(meta$Infection == "Mock")),]
data1 <- log2(exp[,meta1$Sample_title])
desMat <- model.matrix(~Time, meta1)
fits <- lmFit(data1, design = desMat)
fitContrasts <- contrasts.fit(fits, makeContrasts("Time18h", levels = desMat))
ebys <- eBayes(fitContrasts)
a <- topTable(ebys, number = 70000)
a$Colour <- sapply(a$adj.P.Val, function(x) {if (x <= 0.05) return("Significant") else return("Not")})
a9 <- a
#write.table(a9, "mock-12h-mock-18h.txt", sep = "\t", row.names = F, col.names = T, quote = F)


