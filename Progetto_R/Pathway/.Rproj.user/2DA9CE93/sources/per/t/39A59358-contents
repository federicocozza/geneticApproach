library(Biobase)
library(RcppCNPy)
load("/Users/federicocozza/Desktop/Federico Cozza/Tagliaferri/geneticApproach/all_exprs_data.RData")
lc4 <- exprs(expr.lc4)
lc4 <- lc4[1:9446,]

genes <- rownames(lc4)

for(i in 1:length(genes)){
  genes[i] = strsplit(genes[i],split="_")[[1]][1]
}
rownames(lc4) <- genes

labelsLC4 <- pData(expr.lc4)
labelsLC4 <- labelsLC4[c(-1,-2,-3,-4,-5,-7)]
labelsLC4[which(labelsLC4 == "SCC1"),"class"] = 1
labelsLC4[which(labelsLC4 == "SCC2"),"class"] = 2
labelsLC4[which(labelsLC4 == "AC1"),"class"] = 3
labelsLC4[which(labelsLC4 == "AC2"),"class"] = 4

labelsLC4 <- as.numeric(labelsLC4$class)

source("enrichmentReactomeAndKegg.R")
enrichRes <- enrichmentReactomeAndKeggByList(rownames(lc4))
namePathway <- c(enrichRes$KEGGPathway@result$Description,enrichRes$ReactomePathway@result$Description)
mapping <- data.frame(ID =  c(paste("KEGG_",enrichRes$KEGGPathway@result$ID,sep=""),paste("REACTOME_",enrichRes$ReactomePathway@result$ID,sep="")),namePathway= c(enrichRes$KEGGPathway@result$Description,enrichRes$ReactomePathway@result$Description))
rownames(mapping) <- mapping$ID

source("getPathway.R")
minimumSizeOfPathway = 15
allPathways <- getPathways(enrichRes,RNA = lc4,minimumSizeOfPathway)

pathwaysNumber <- length(allPathways)
subjectsNumber <- length(labelsLC4)
featurePathways <- matrix(0, nrow = pathwaysNumber, ncol = subjectsNumber)

for (i in 1:pathwaysNumber) {
  pathway = allPathways[[i]]
  mediansPathway <- colMedians(pathway)
  featurePathways[i, ] <- mediansPathway
}

featurePathways <- t(featurePathways)
npySave("dataset.npy", featurePathways)
npySave("labels.npy", labelsLC4)