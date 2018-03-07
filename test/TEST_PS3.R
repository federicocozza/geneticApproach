library(Biobase)

load("~/TesiLuca/Data/all_exprs_data.RData")
PS3 <- exprs(expr.ps3)
idx <- which(rownames(PS3) == "AFFX-BioB-3_at")
PS3 <- PS3[1:(idx-1),]

rownames(PS3) <- gsub("_at","",rownames(PS3))
colnames(PS3) <- gsub(" ","",colnames(PS3))


labelsPS3 <- pData(expr.ps3)
labelsPS3 <- labelsPS3[c(-1,-2,-3,-4,-5)]
labelsPS3[which(labelsPS3 == "NormalSkin"),"class"] = 1
labelsPS3[which(labelsPS3 == "UninvolvedSkin"),"class"] = 2
labelsPS3[which(labelsPS3 == "InvolvedSkin"),"class"] = 3
rownames(labelsPS3) <- gsub(" ","",rownames(labelsPS3))

labelMapping <-  data.frame(OriginalLabel = c("NormalSkin","UninvolvedSkin","InvolvedSkin"),numberLabels = 1:length(unique(labelsPS3)))

colnames(labelsPS3) = "Classes"

library(SNFtool)
PS3 <- t(standardNormalization(t(PS3)))

source("R/enrichmentReactomeAndKegg.R")
enrichRes <- enrichmentReactomeAndKeggByList(rownames(PS3))
namePathway <- c(enrichRes$KEGGPathway@result$Description,enrichRes$ReactomePathway@result$Description)
mapping <- data.frame(ID =  c(paste("KEGG_",enrichRes$KEGGPathway@result$ID,sep=""),paste("REACTOME_",enrichRes$ReactomePathway@result$ID,sep="")),namePathway= c(enrichRes$KEGGPathway@result$Description,enrichRes$ReactomePathway@result$Description))
rownames(mapping) <- mapping$ID

source("R/getPathway.R")
minimumSizeOfPathway = 15
allPathways <- getPathways(enrichRes,RNA = PS3,minimumSizeOfPathway)

source("R/SVM.R")
pythonPath <- "/home/neuronelab/miniconda2/bin/python"
SVMResults <- SVM(pythonPath,allPathways,labelsPS3,core=8,outer_fold = 3,inner_fold = 2)
mp <- SVMResults$accuracies$accuraciesMatrix
mp$np = "NULL"
rownames(mp) = as.character(mp$pathway)

for(i in 1:nrow(mp)){
  mp[as.character(mp[i,"pathway"]),"np"] = as.character(mapping[strsplit(as.character(mp[i,"pathway"]),split="-")[[1]][1],"namePathway"])
}
library(ggplot2)
source("R/reductionPathwaysSVM.R")
SVMWeights <- SVMResults$svmWeights
dataset <- PS3
reductionPathwaysRes <- reductionPathwaysSVM(pythonPath,SVMWeights,dataset,labelsPS3,percentage = 0.8,core = 8, outer_fold = 3 ,inner_fold = 2)

source("R/permutationTest.R")

accuraciesMatrix <- reductionPathwaysRes$accuracies$accuraciesMatrix
accuraciesMatrixRed <- accuraciesMatrix[accuraciesMatrix$size >= minimumSizeOfPathway,]

pValueCutOff <- 0.05
numberPermutation = 1000

permResults <- permutationTest(pythonPath,numberPermutation,accuraciesMatrix = accuraciesMatrixRed,SVMWeighs = reductionPathwaysRes$svmWeights,dataSet = dataset,labels = labelsPS3,pValueCutOff = pValueCutOff,core= 16,outer_fold=3,inner_fold=2)
permResults$permutationResults$pValue <- p.adjust(permResults$permutationResults$pValue, method = "fdr")

permResults$significantlyPathways <- permResults$permutationResults[permResults$permutationResults$pValue < pValueCutOff,]

dir.create("PS3")
dir.create("PS3/AccuracyVsSizePathwaysPlots")
ggsave(permResults$accuracyVsPathwaysSizePlot,filename = paste("PS3/AccuracyVsSizePathwaysPlots/AccuracyVSPathwaysSizeAll.png",sep=""), width = 14, height = 12)

mp$sizeReduced = 1
mp$accuracyReduced = 1
mp$pValue = 1
mp$genes = "genes"
rownames(accuraciesMatrix) = accuraciesMatrix$pathway
rownames(permResults$permutationResults) = permResults$permutationResults$namePathway
SVMWeightsRed <- reductionPathwaysRes$svmWeights
source("R/transformInSymbols.R")
for(i in 1:nrow(mp)){
  mp[as.character(mp[i,"pathway"]),"sizeReduced"] = accuraciesMatrix[as.character(mp[i,"pathway"]),"size"]
  mp[as.character(mp[i,"pathway"]),"accuracyReduced"] = accuraciesMatrix[as.character(mp[i,"pathway"]),"accuracy"]
  mp[as.character(mp[i,"pathway"]),"pValue"] = permResults$permutationResults[as.character(mp[i,"pathway"]),"pValue"]
  genes <- ""
  for(j in 1:length(SVMWeightsRed[as.character(mp[i,"pathway"])][[1]])){
    genes <- paste(genes,transformInSymbols(names(SVMWeightsRed[as.character(mp[i,"pathway"])][[1]])[j])[[1]],"(",round(SVMWeightsRed[as.character(mp[i,"pathway"])][[1]][j]*100,2),"%),",sep="")
  }
  
  mp[as.character(mp[i,"pathway"]),"genes"] <- genes
}

ordCol = c("np","pathway","class","size","accuracy","sizeReduced","accuracyReduced","pValue","genes")
mp <- mp[,ordCol]

## ORDER BY P VALUE
mp <- mp[order(mp$pValue),]

WriteXLS::WriteXLS(mp,"PS3/PS3-AllPathwaysNOFS.xls")

SignPathways <- mp[as.character(permResults$significantlyPathways$namePathway),]
SignPathways <- SignPathways[order(SignPathways$pValue),]
WriteXLS::WriteXLS(SignPathways,"PS3/PS3-OnlySignificantlyPathwaysNOFS.xls")

PS3_NOFS_HistPlot <- ggplot(SignPathways, aes(accuracyReduced)) + geom_histogram(aes(fill=class),col="black") + ggtitle("Pathways Accuracy") +
  theme_grey() + stat_bin(geom="text",colour="black",size=3.5,aes(label = ifelse(..count..>0,..count..,""),fill=class),position = position_stack(vjust=0.5)) +scale_x_continuous(breaks = seq(0,1,0.05),limits=c(0,1.1)) + xlab("Accuracy")
ggsave(PS3_NOFS_HistPlot,filename = "PS3/HistogramAccuracySignificantPathways.png")

dir.create("PS3/RDATA")
save.image(file="PS3/RDATA/PS3-AFTER-PERMUTATION-TEST.RData")

####### ACCURATEZZA VS SIZE GRAFICO

pValueCutOff <- 0.05
significantClass <- unique(as.character(permResults$permutationResults$classes))
for(cl in significantClass){
  toPlot <- permResults$permutationResults
  toPlot <- toPlot[toPlot$classes == cl,]
  significance <- rep("significant",nrow(toPlot))
  significance[toPlot$pValue > pValueCutOff] = "not significant"
  toPlot <- cbind(toPlot,significance)
  
  percentage <- 1 - toPlot[,"pValue"]
  toPlot <- cbind(toPlot,percentage)
  
  correlationData = c()
  for(i in unique(toPlot$classes)){
    correlationData[[i]] <- cor(toPlot[toPlot$classes == i,]$accuracy,toPlot[toPlot$classes == i,]$size)
  }
  maxSize <- max(permResults$permutationResults$size) + 5
  xPos <- min(permResults$permutationResults$accuracy)+0.075
  
  positions <- seq(maxSize,maxSize+4*(length(unique(toPlot$classes))),4)
  
  correlationString <- c()
  for(i in 1:length(correlationData)){
    correlationString[[i]] <- paste(names(correlationData)[i]," r = ",round(as.numeric(correlationData[i]),2),sep="")
  }
  correlationString <- rev(correlationString)
  
  correlationString <- c(correlationString,"Correlation:")
  p <- ggplot(toPlot,aes(x=accuracy,y=size,color = classes)) + geom_point(aes(shape=significance,size=percentage)) +
    xlab("Accuracy") + ylab("Size") + xlim(min(permResults$permutationResults$accuracy),max(permResults$permutationResults$accuracy)) +
    ggtitle(paste("Accuracy vs Pathways Size class ", cl,sep="")) +
    annotate("text",x=rep(xPos,length(unique(toPlot$classes))+1),y=positions,label = correlationString) +
    geom_smooth(method = lm,se=T,fullrange=F) + scale_colour_discrete(drop=TRUE,limits = levels(toPlot$classes))
  
  ggsave(p,filename = paste("PS3/AccuracyVsSizePathwaysPlots/AccuracyVSPathwaysSize-",cl,".png",sep=""), width = 14, height = 12)
}

####CONTROLLO CONSISTENZA GENI
# for(i in 1:length(allRedPath)){
#   typePathway <- strsplit(names(allRedPath[i]),split="_")[[1]][1]
#   if(typePathway == "KEGG"){
#     namePathwayZ <- strsplit(strsplit(names(allRedPath[i]),split="_")[[1]][2],split="-")[[1]][1]
#     genesByPath <- strsplit(enrichRes$KEGGPathway@result[namePathwayZ,"geneID"],split="/")[[1]]
#     geneByMatrix <- names(allRedPath[[i]])
#     
#     if(sum(geneByMatrix %in% genesByPath) != length(geneByMatrix)){
#       print(paste(i))
#     }
#   }
# }


###BOX PLOT VALORI ESPRESSIONE###
namesSignPathways <- as.character(permResults$significantlyPathways$namePathway)
allRedPath <- reductionPathwaysRes$svmWeights[as.character(permResults$significantlyPathways$namePathway)]

SVMWeightsOnlySign <- reductionPathwaysRes$svmWeights[as.character(permResults$significantlyPathways$namePathway)]
dir.create("PS3/BoxPlot-PS3-NOFS")

library(ggplot2)
for(i in 1:length(allRedPath)){
  classComb1 <- strsplit(strsplit(names(SVMWeightsOnlySign[names(allRedPath[i])]),split="-")[[1]][2],split="vs")[[1]]
  genes <- names(SVMWeightsOnlySign[[names(allRedPath[i])]])
  prova <- as.data.frame(t(PS3[genes,]))
  prova <- stack(prova)
  prova$class <- labelsPS3$Classes
  prova <- prova[c(which(prova$class == classComb1[1]),which(prova$class == classComb1[2])),]
  prova$class <- as.factor(prova$class)
  prova$ind <- factor(prova$ind,levels=rev(genes))
  
  
  pp <- ggplot(prova) + geom_boxplot(aes(x = ind, y = values,fill = class)) + ggtitle(names(allRedPath)[i])+ theme(text = element_text(size=20)) + theme_bw() + coord_flip()
  ggsave(paste("PS3/BoxPlot-PS3-NOFS/",names(allRedPath)[i],".png",sep=""), plot = pp,width = 800, height = 1400,limitsize = FALSE,units = "mm")
  print(paste(i,"di",length(allRedPath)))
}

source("R/classifyByGenesList.R")

importantGenes <- list()
namePathways <- names(SVMWeightsOnlySign)
combinations <- as.vector(unique(permResults$significantlyPathways$classes))
for(j in combinations){
  SVMWeightsOnlySignClass <- SVMWeightsOnlySign[grepl(j,namePathways,fixed=T)]
  for(i in 1:length(SVMWeightsOnlySignClass)){
    importantGenes[[j]] <- union(importantGenes[[j]],names(SVMWeightsOnlySignClass[[i]])[1])
  }
}

classifyByGenesResult <- list()
for(i in combinations){
  classifyByGenesResult[[i]] <- classifyByGenesList(pythonPath,i,importantGenes[[i]],PS3,labelsPS3,core = 8, outer_fold = 3 ,inner_fold = 2)
}

source("R/transformInSymbols.R")
ImportantGenesDF <- data.frame(class=names(classifyByGenesResult))
ImportantGenesDF$accuracy = ""
ImportantGenesDF$genes = ""
for(i in 1:length(classifyByGenesResult)){
  ImportantGenesDF[i,"accuracy"] <- classifyByGenesResult[[i]]$accuracy
  
  genes <- ""
  for(j in 1:length(classifyByGenesResult[[i]]$svmWeights[[1]])){
    genes <- paste(genes,transformInSymbols(names(classifyByGenesResult[[i]]$svmWeights[[1]])[j])[[1]],"(",round(classifyByGenesResult[[i]]$svmWeights[[1]][j]*100,2),"%),",sep="")
  }
  ImportantGenesDF[i,"genes"] <- genes
}
rownames(ImportantGenesDF) <- ImportantGenesDF$class

write.table(ImportantGenesDF,"PS3/ImportantGenes.xls",row.names = F)

pathwaysBest <- data.frame()
for(i in combinations){
  pathComb <- SignPathways[which(SignPathways$class == i),]
  pathComb <- pathComb[which(pathComb$accuracyReduced>ImportantGenesDF[which(ImportantGenesDF$class==i),"accuracy"]),]
  pathwaysBest <- rbind(pathwaysBest,pathComb)
}

###creare boxplot dei geni più importanti
library(ggplot2)
dir.create("PS3/ImportantGenes-BOXPLOT")
for(i in combinations){
  classComb1 <- strsplit(i,split="vs")[[1]]
  
  prova <- as.data.frame(t(PS3[importantGenes[[i]],]))
  if(length(importantGenes[[i]])==1){
    # prova <- as.data.frame(prova)
    # prova <- t(prova)
    # colnames(prova) <- importantGenes[[i]]
    rownames(prova) <- importantGenes[[i]]
    prova <- as.data.frame(t(prova))
    prova$class <- labelsPS3$Classes
    prova <- prova[c(which(prova$class == classComb1[1]),which(prova$class == classComb1[2])),]
    prova$class <- as.factor(prova$class)
    
    pp <- ggplot(prova) + geom_boxplot(aes(x = class, y = prova[,1],fill = class),data = prova) + ylab(importantGenes[[i]]) + ggtitle(paste("importantGenes",i))
    ggsave(paste("PS3/ImportantGenes-BOXPLOT/ImportantGenes-",i,".png",sep=""), plot = pp,width = 800, height = 1400,limitsize = FALSE,units = "mm")
  }else{
    prova <- stack(prova)
    prova$class <- labelsPS3$Classes
    prova <- prova[c(which(prova$class == classComb1[1]),which(prova$class == classComb1[2])),]
    prova$class <- as.factor(prova$class)
    pp <- ggplot(prova) + geom_boxplot(aes(x = ind, y = values,fill = class),data = prova) +
      ggtitle(paste("importantGenes",i)) + theme_bw() + coord_flip()
    ggsave(paste("PS3/ImportantGenes-BOXPLOT/ImportantGenes-",i,".png",sep=""), plot = pp,width = 800, height = 1400,limitsize = FALSE,units = "mm")
  }
}

hmaps <- list()
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

allImportantGenes <- c()
for(i in 1:length(importantGenes)){
  allImportantGenes <- union(allImportantGenes,importantGenes[[i]])
}
PS3Norm <- range01(PS3[allImportantGenes,])

library(gplots)

dir.create("PS3/ImportantGenes-HeatMap")
for(i in 1:length(importantGenes)){
  a <- names(importantGenes)[i]
  title <- paste(labelMapping[strsplit(a,split="vs")[[1]][1],"OriginalLabel"],"vs",labelMapping[strsplit(a,split="vs")[[1]][2],"OriginalLabel"])
  
  if(length(importantGenes[[i]]) == 1){
    # X1 <- as.matrix(t(PS3Norm[importantGenes[[i]],which(labelsPS3$Classes== strsplit(a,split="vs")[[1]][1])]))
    # X2 <- as.matrix(t(PS3Norm[importantGenes[[i]],which(labelsPS3$Classes==strsplit(a,split="vs")[[1]][2])]))
    # heatData <- cbind(X1,X2)
    # rownames(heatData) <- importantGenes[[i]]
  }else{
    heatData <- cbind(PS3Norm[importantGenes[[i]],which(labelsPS3$Classes== strsplit(a,split="vs")[[1]][1])],PS3Norm[importantGenes[[i]],which(labelsPS3$Classes==strsplit(a,split="vs")[[1]][2])])
    png(file = paste("PS3/ImportantGenes-HeatMap/",names(importantGenes)[i],".png",sep=""),width = 800,height = 600)
    heatmap.2(heatData,Colv = FALSE,dendrogram = "row",trace = "none",colsep=dim(PS3Norm[importantGenes[[i]],which(labelsPS3$Classes==strsplit(a,split="vs")[[1]][1])])[2],sepcolor ="black",main=title,xlab="Patients",ylab= "Genes",labCol = F)
    dev.off()
  }
}

#### calcolo delle statistiche
PathStatData <- data.frame()
dir.create("PS3/Statistiche")
dir.create("PS3/HeatMapPathway")
for(i in 1:length(allRedPath)){
  ### statistiche (media varianza)
  classComb1 <- strsplit(strsplit(names(SVMWeightsOnlySign[names(allRedPath[i])]),split="-")[[1]][2],split="vs")[[1]]
  genes <- names(SVMWeightsOnlySign[[names(allRedPath[i])]])
  pathway <- as.data.frame(PS3[genes,])
  idxCB1 <- which(labelsPS3$Classes == classComb1[1])
  idxCB2 <- which(labelsPS3$Classes == classComb1[2])
  pathwayOrd <- pathway[,c(idxCB1,idxCB2)]
  
  pathwayComb1 <- pathwayOrd[,1:length(idxCB1)]
  pathwayComb2 <- pathwayOrd[,(length(idxCB1)+1):(length(idxCB1)+length(idxCB2))]
  
  mean1 <- rowMeans(pathwayComb1)
  mean2 <- rowMeans(pathwayComb2)
  var1 <- apply(pathwayComb1, 1, var)
  var2 <- apply(pathwayComb2, 1, var)
  
  median1 <- apply(pathwayComb1, 1, median)
  median2 <- apply(pathwayComb2, 1, median)
  
  absMean <- abs(mean1-mean2)
  weight <-   1:length(absMean)
  weight <- rev(weight)
  absMeanWeighted <- absMean*weight
  
  dataStat <- data.frame(genes = names(mean1),mean1,mean2,absMean,var1,var2,median1,median2)
  write.table(dataStat, file= paste("PS3/Statistiche/",names(allRedPath[i]),".csv",sep=""),row.names = F)
  
  PathStatData <- rbind(PathStatData,data.frame(namePathway = names(allRedPath[i]),valueNorm = sum(absMean)/length(absMean),valueWeighted = sum(absMeanWeighted)/length(absMean)))
  
  ### creazione heatmap
  jpeg(file = paste("PS3/HeatMapPathway/",names(allRedPath[i]),".jpeg",sep=""),width=1400,height=1200)
  heatmap.2(as.matrix(pathwayOrd),Colv = FALSE,dendrogram = "row",trace = "none",colsep=length(idxCB1),sepcolor ="black",main=names(allRedPath[i]),xlab="Patients",ylab= "Genes",labCol = F)
  dev.off()
}
write.table(PathStatData,file="PS3/StatisticheRiassuntive.csv",row.names = F)


### creare heatmap classi x pathway
namePathwayAll <- c()
for(i in 1:length(allRedPath)){
  namePathwayAll <- union(namePathwayAll,strsplit(names(allRedPath[i]),split="-")[[1]][1])
}
namePathwayAll <- gsub("KEGG","K",namePathwayAll)
namePathwayAll <- gsub("REACTOME","R",namePathwayAll)

allPathMatrix <- matrix(nrow=length(namePathwayAll),ncol=4,0)
rownames(allPathMatrix) <- namePathwayAll
colnames(allPathMatrix) <- 1:4
allRedPath2 <- allRedPath
names(allRedPath2) <- gsub("KEGG","K",names(allRedPath2))
names(allRedPath2) <- gsub("REACTOME","R",names(allRedPath2))

for(i in 1:length(allRedPath2)){
  genes <- names(allRedPath2[[i]])
  pp <- PS3[genes,]
  pp1 <- pp[,which(labelsPS3 == "1")]
  pp2 <- pp[,which(labelsPS3 == "2")]
  pp3 <- pp[,which(labelsPS3 == "3")]
  pp4 <- pp[,which(labelsPS3 == "4")]
  
  mc1 <- mean(pp1)
  mc2 <- mean(pp2)
  mc3 <- mean(pp3)
  mc4 <- mean(pp4)
  
  allPathMatrix[strsplit(names(allRedPath2[i]),split="-")[[1]][1],1] = mc1
  allPathMatrix[strsplit(names(allRedPath2[i]),split="-")[[1]][1],2] = mc2
  allPathMatrix[strsplit(names(allRedPath2[i]),split="-")[[1]][1],3] = mc3
  allPathMatrix[strsplit(names(allRedPath2[i]),split="-")[[1]][1],4] = mc4
}

dev.off()
jpeg(file = paste("PS3/HeatMap-PathwaysRelation.jpeg",sep=""),width=1500,height=1200)
heatmap.2(allPathMatrix,Colv = FALSE,dendrogram = "row",trace = "none",main="Pathways relation",xlab="class",ylab= "pathway",labCol = F,margins=c(2,6))
dev.off()
### Grafo

##### creare grafo
# require(igraph)

##GRAFO TOTALE
# GraphData <- matrix(0,length(names(allRedPath)),length(names(allRedPath)))
# rownames(GraphData) <- names(allRedPath)
# colnames(GraphData) <- names(allRedPath)
#
# for(i in 1:length(allRedPath)){
#   for(j in 1:length(allRedPath)){
#     GraphData[i,j] <- length(intersect(rownames(allRedPath[[i]]),rownames(allRedPath[[j]])))/length(union(rownames(allRedPath[[i]]),rownames(allRedPath[[j]])))
#   }
# }
# diag(GraphData) = 0
#
# GraphDataRed <- GraphData
# GraphDataRed[GraphDataRed < quantile(GraphDataRed,0.9)[[1]]] = 0
# sum(GraphDataRed != 0)

# GraphDataRed <- GraphData
# GraphDataRed[GraphDataRed < quantile(GraphDataRed,0.85)[[1]]] = 0
# sum(GraphDataRed != 0)
#
# GraphDataRed <- GraphData
# GraphDataRed[GraphDataRed < quantile(GraphDataRed,0.95)[[1]]] = 0
# sum(GraphDataRed != 0)

# hist(GraphData,main = "Histogram",breaks = 50,xlab = "Weight",col="grey")
# abline(v=quantile(GraphDataRed,0.9)[[1]],col="red")
#
# g <- graph.adjacency(GraphDataRed,weighted = T,mode = "undirected")
#
# bet <- betweenness(g,directed = F,weights = rep(1,length(E(g))))
# deg <- degree(g)
#
# # wtc <- walktrap.community(g,weights=E(g)$weight, steps = 4)
#
# acc <- SignPathways[names(bet),"accuracyReduced"]
# names(acc) <- as.character(SignPathways[names(bet),"pathway"])

# ## create 3-partite graph unweighted betweennes accuracy and degree
# library(MASS)
# dev.off()
# par(xpd = T, mar = par()$mar + c(0,0,0,7))
# DataForPlot <- data.frame(betweenness = bet, accuracy= acc, degree = deg)
# ii <- cut(as.numeric(acc), breaks = seq(min(as.numeric(acc)), max(as.numeric(acc)), len = 50), include.lowest = TRUE)
# colors <- colorRampPalette(c("green", "red"))(49)[ii]
# parcoord(DataForPlot,col = colors,var.label = T)
# title("Relations unweighted betweennes, accuracy and degree")
# legend(x=3.03,y=1,c("low accuracy","high accuracy"),col=c("green","red"),pch=16)


# ### create 3-partite graph weighted betweennes accuracy and strenght
# weightbet <- betweenness(g,directed = F,weights = (1-(E(g)$weight)))
# strengthI <- strength(g,weight=(1-(E(g)$weights)))
# dev.off()
# par(xpd = T, mar = par()$mar + c(0,0,0,7))
# DataForPlot2 <- data.frame(weightedBetweenness = weightbet, accuracy= acc, strength = strengthI)
# ii <- cut(as.numeric(acc), breaks = seq(min(as.numeric(acc)), max(as.numeric(acc)), len = 50), include.lowest = TRUE)
# colors <- colorRampPalette(c("green", "red"))(49)[ii]
# parcoord(DataForPlot2,col = colors,var.label = T)
# title("Relations weighted betweennes, accuracy and strength")
# legend(x=3.03,y=1,c("low accuracy","high accuracy"),col=c("green","red"),pch=16)

# ## hive plot fail
# axis = c(rep("a",length(acc)),rep("b",length(acc)),rep("c",length(acc)))
# value = c(acc[1:length(acc)],deg[1:length(deg)],bet[1:length(bet)])
# label = c(names(acc),names(deg),names(bet))
# id = paste(axis,label,sep="-")
# rCol <- c(rep(255,length(acc)),rep(0,length(acc)),rep(0,length(acc)))
# gCol <- c(rep(0,length(acc)),rep(255,length(acc)),rep(0,length(acc)))
# bCol <- c(rep(0,length(acc)),rep(0,length(acc)),rep(255,length(acc)))
# HiveData <- data.frame(id=id,axis=axis,label=label,value=value,height=1,r=rCol,g=gCol,b=bCol,a=200)
# aEdge <- paste("a",names(acc),sep="-")
# bEdge <- paste("b",names(acc),sep="-")
# cEdge <- paste("c",names(acc),sep="-")
#
# firstCol <- c(aEdge,aEdge,bEdge)
# secCol <- c(bEdge,cEdge,cEdge)
#
# EdgeHiveData <- data.frame(startnodelabel=firstCol,endnodelabel = secCol)
# names(EdgeHiveData) = c("start node label", "end node label")
#
# write.table(HiveData,file="nodes.txt",sep="\t",quote=F,row.names = F,col.names = T)
# write.table(EdgeHiveData,file="edge.txt",sep="\t",quote=F,row.names = F,col.names = T)


#########################
### Pathway combinations
########################

source("R/pathwaysCombination.R")
significantlyPathways <- permResults$significantlyPathways
probabilities <- reductionPathwaysRes$probabilities
pathCombProbabilities <- pathwaysCombination(pythonPath,probabilities,significantlyPathways,labelsClass = labelsPS3,core=32,outer_fold=3,inner_fold=2,maxNumCombination = 2)

pathCombPlot <- ggplot(pathCombProbabilities$combination, aes(accuracy)) + geom_histogram(aes(fill=class),col="black") + ggtitle("Combination Pathways Accuracies") +
  theme_grey() + stat_bin(geom="text",colour="black",size=3.5,aes(label = ifelse(..count..>0,..count..,""),fill=class),position = position_stack(vjust=0.5)) +scale_x_continuous(breaks = seq(0,1,0.05),limits=c(0,1.1))
ggsave("PS3/PathwayCombinationHistogram.png", plot = pathCombPlot,width = 800, height = 600,limitsize = FALSE,units = "mm")

pathCombProbResume <- data.frame()

cutOff1vs2 <- ImportantGenesDF[which(ImportantGenesDF$class == "1vs2"),"accuracy"]
pathComb1vs2 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "1vs2"),]
pathCombs1vs2Best <- pathComb1vs2[which(pathComb1vs2$accuracy > cutOff1vs2),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class = unique(as.character(pathCombs1vs2Best$class))[1],nComb= nrow(pathCombs1vs2Best)))

cutOff1vs3 <- ImportantGenesDF[which(ImportantGenesDF$class == "1vs3"),"accuracy"]
pathComb1vs3 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "1vs3"),]
pathCombs1vs3Best <- pathComb1vs3[which(pathComb1vs3$accuracy > cutOff1vs3),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class= unique(as.character(pathCombs1vs3Best$class))[1],nComb = nrow(pathCombs1vs3Best)))

cutOff1vs4 <- ImportantGenesDF[which(ImportantGenesDF$class == "1vs4"),"accuracy"]
pathComb1vs4 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "1vs4"),]
pathCombs1vs4Best <- pathComb1vs4[which(pathComb1vs4$accuracy > cutOff1vs4),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class = unique(as.character(pathCombs1vs4Best$class))[1],nComb = nrow(pathCombs1vs4Best)))

cutOff2vs4 <- ImportantGenesDF[which(ImportantGenesDF$class == "2vs4"),"accuracy"]
pathComb2vs4 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "2vs4"),]
pathCombs2vs4Best <- pathComb2vs4[which(pathComb2vs4$accuracy > cutOff2vs4),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class = unique(as.character(pathCombs2vs4Best$class))[1],nComb = nrow(pathCombs2vs4Best)))

cutOff3vs2 <- ImportantGenesDF[which(ImportantGenesDF$class == "3vs2"),"accuracy"]
pathComb3vs2 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "3vs2"),]
pathCombs3vs2Best <- pathComb3vs2[which(pathComb3vs2$accuracy > cutOff3vs2),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class = unique(as.character(pathCombs3vs2Best$class))[1],nComb = nrow(pathCombs3vs2Best)))

cutOff3vs4 <- ImportantGenesDF[which(ImportantGenesDF$class == "3vs4"),"accuracy"]
pathComb3vs4 <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == "3vs4"),]
pathCombs3vs4Best <- pathComb3vs4[which(pathComb3vs4$accuracy > cutOff3vs4),]
pathCombProbResume <- rbind(pathCombProbResume,data.frame(class = unique(as.character(pathCombs3vs4Best$class))[1],nComb = nrow(pathCombs3vs4Best)))

allCombinationBest <- rbind(pathCombs1vs2Best,pathCombs1vs3Best,pathCombs1vs4Best,pathCombs2vs4Best,pathCombs3vs2Best,pathCombs3vs4Best)
WriteXLS::WriteXLS(allCombinationBest,ExcelFileName = "PS3/allCombinationBest.xls",row.names = F)

save.image(file="PS3/RDATA/PS3-NO-FS-AfterPathwayCombination.RData")



################################
######## META FEATURES  ########
################################

pSign <- as.character(permResults$significantlyPathways$namePathway)
rp <- reductionPathwaysRes$svmWeights[pSign]
metaFeatures <- matrix(0,nrow = length(rp),ncol=nrow(labelsPS3))
rownames(metaFeatures) <- names(rp)
colnames(metaFeatures) <- rownames(labelsPS3)

for(i in 1:length(rp)){
  metaPathway <- matrix(0,length(rp[[i]]),nrow(labelsPS3))
  for(j in 1:length(rp[[i]])){
    ggg <- names(rp[[i]])[j]
    www <- (rp[[i]])[j]
    metaPathway[j,] <- www*PS3[ggg,]
  }
  metaFeatures[i,] <- colMeans(metaPathway)
}

combinations <- as.character(unique(permResults$significantlyPathways$classes))
metaFeaturesResults <- matrix(0,nrow = length(combinations),ncol=2)
rownames(metaFeaturesResults) <- combinations
colnames(metaFeaturesResults) <- c("Class","Accuracy")
for(comb in combinations){
  path <- paste("p",comb,sep="-")
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathwayMetaFeatures",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  cols <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  metaClasses <- metaFeatures[which(permResults$significantlyPathways$classes == comb),cols]
  
  if(class(metaClasses) == "numeric"){
    metaClasses <- as.data.frame(t(metaClasses))
    rownames(metaClasses) <- rownames(metaFeatures)[which(permResults$significantlyPathways$classes == comb)]
  }
  
  dir.create(paste(temporanyDirectory,"/",path,sep=""))
  write.table(t(metaClasses),paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)
  
  labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
  labels2classes <- labelsPS3[cols,"Classes"]
  labels2classes <- as.data.frame(labels2classes)
  rownames(labels2classes) <- rownames(labelsPS3)[cols]
  colnames(labels2classes) <- "x"
  
  write.table(labels2classes,file=labelPath ,quote=F)
  
  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core = 8,outer_fold = 3,inner_fold = 2,"C"))
  
  accuracies <- getPathwaysAccuracies(temporanyDirectory)
  
  svmWeights.pathwayMetaFeatures <- getSVMweights(temporanyDirectory)
  
  metaFeaturesResults[comb,"Class"] <- comb
  metaFeaturesResults[comb,"Accuracy"] <- accuracies$accuraciesMatrix$accuracy[1]
  
  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)
}

##################################
######## MEDIAN FEATURES  ########
##################################

pSign <- as.character(permResults$significantlyPathways$namePathway)
rpm <- reductionPathwaysRes$svmWeights[pSign]
medianFeatures <- list()

combinations <- as.character(unique(permResults$significantlyPathways$classes))
for(comb in combinations){
  ncols <- sum(c(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1],labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  medianFeatures[[comb]] <- matrix(0,nrow=sum(permResults$significantlyPathways$classes == comb),ncol = ncols)
  j=1
  for(i in 1:length(rpm)){
    if(grepl(comb,names(rpm)[i])){
      cols <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
      dataToFeatures <- PS3[names(rpm[[i]]),cols]
      medianFeatures[[comb]][j,] <- as.numeric(apply(dataToFeatures,2,FUN = median))
      j=j+1
      colnames(medianFeatures[[comb]]) <- colnames(dataToFeatures)
      rownames(medianFeatures[[comb]]) <- names(rpm)[which(grepl(comb,names(rpm)))]
    }
  }
}

medianFeaturesResults <- matrix(0,nrow = length(combinations),ncol=2)
rownames(medianFeaturesResults) <- combinations
colnames(medianFeaturesResults) <- c("Class","Accuracy")
for(comb in combinations){
  path <- paste("p",comb,sep="-")
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathwayMedianFeatures",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  metaClasses <- medianFeatures[[comb]]
  
  dir.create(paste(temporanyDirectory,"/",path,sep=""))
  write.table(t(metaClasses),paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)
  
  cols <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  
  labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
  labels2classes <- labelsPS3[cols,"Classes"]
  labels2classes <- as.data.frame(labels2classes)
  rownames(labels2classes) <- rownames(labelsPS3)[cols]
  colnames(labels2classes) <- "x"
  
  write.table(labels2classes,file=labelPath ,quote=F)
  
  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core = 8,outer_fold = 3,inner_fold = 2,"C"))
  
  accuracies <- getPathwaysAccuracies(temporanyDirectory)
  
  svmWeights.pathwayMedianFeatures <- getSVMweights(temporanyDirectory)
  
  medianFeaturesResults[comb,"Class"] <- comb
  medianFeaturesResults[comb,"Accuracy"] <- accuracies$accuraciesMatrix$accuracy[1]
  
  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)
}


#######################################
######## CORRELATION FEATURES  ########
#######################################
combinations <- as.character(unique(permResults$significantlyPathways$classes))
ncols <- list()
maxCorGenes <- list()
for(comb in combinations){
  ncols[[comb]] <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  maxCorGenes[[comb]] <- list()
  for(i in 1:length(rpm)){
    if(grepl(comb,names(rpm)[i])){
      matrixToAdd <- PS3[names(rpm[[i]]),ncols[[1]]]
      corGenes <- cor(t(matrixToAdd))
      maxCorGenes[[comb]] <- union(maxCorGenes[[comb]],names(which.max(rowMeans(corGenes))))
    }
  }
}


correlationFeaturesResults <- matrix(0,nrow = length(combinations),ncol=2)
rownames(correlationFeaturesResults) <- combinations
colnames(correlationFeaturesResults) <- c("Class","Accuracy")
for(comb in combinations){
  path <- paste("p",comb,sep="-")
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/correlationFeatures",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  corClasses <- PS3[unlist(maxCorGenes[[comb]]),ncols[[comb]]]
  
  if(length(unlist(maxCorGenes[[comb]])) == 1){
    corClasses <- as.data.frame(t(corClasses))
    rownames(corClasses) <- unlist(maxCorGenes[[comb]])
  }
  
  dir.create(paste(temporanyDirectory,"/",path,sep=""))
  write.table(t(corClasses),paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)
  
  cols <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  
  labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
  labels2classes <- labelsPS3[cols,"Classes"]
  labels2classes <- as.data.frame(labels2classes)
  rownames(labels2classes) <- rownames(labelsPS3)[cols]
  colnames(labels2classes) <- "x"
  
  write.table(labels2classes,file=labelPath ,quote=F)
  
  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core = 8,outer_fold = 3,inner_fold = 2,"C"))
  
  accuracies <- getPathwaysAccuracies(temporanyDirectory)
  
  svmWeights.pathwayCorrelationFeatures <- getSVMweights(temporanyDirectory)
  
  correlationFeaturesResults[comb,"Class"] <- comb
  correlationFeaturesResults[comb,"Accuracy"] <- accuracies$accuraciesMatrix$accuracy[1]
  
  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)
}

corDataWithGenes <- correlationFeaturesResults
corDataWithGenes <- as.data.frame(corDataWithGenes)

for(i in rownames(corDataWithGenes)){
  corDataWithGenes[i,"genesSymbol"] <- paste(transformInSymbols(unlist(maxCorGenes[[i]])),collapse=",")
  corDataWithGenes[i,"genesID"] <- paste(unlist(maxCorGenes[[i]]),collapse=",")
}

write.table(corDataWithGenes,file="PS3/CorrelatioFeaturesData.csv",row.names = F)

##########################################################
############## INFERRING PATHWAY ACTIVITY ################
##########################################################

library(SNFtool)


pSign <- as.character(permResults$significantlyPathways$namePathway)
rpm <- reductionPathwaysRes$svmWeights[pSign]

IPFeatures <- list()

for(i in 1:length(rpm)){
  classType <- strsplit(names(rpm[i]),split="-")[[1]][2]
  patientsClass1 <- which(labelsPS3$Classes == strsplit(classType,split="vs")[[1]][1])
  patientsClass2 <- which(labelsPS3$Classes == strsplit(classType,split="vs")[[1]][2])
  
  patientsIdx <- c(patientsClass1,patientsClass2)
  
  matrixPathway <- PS3[names(rpm[[i]]),patientsIdx]
  matrixPathway <- standardNormalization(matrixPathway)
  ttestPathway <- c()
  for(j in 1:nrow(matrixPathway)){
    gene <- matrixPathway[j,]
    
    geneClass1 <- gene[1:length(patientsClass1)]
    geneClass2 <- gene[(length(patientsClass1)+1):(length(patientsClass1)+length(patientsClass2))]
    
    ttres <- t.test(geneClass1,geneClass2)
    ttestPathway <- c(ttestPathway,ttres$statistic)
  }
  names(ttestPathway) <- rownames(matrixPathway)
  
  
  mStatistics <- mean(ttestPathway)
  if(mStatistics >= 0){
    idxOrdered <- order(-ttestPathway)
  }else{ ##mstatistics < 0
    idxOrdered <- order(ttestPathway)
  }
  genesOrdered <- rownames(matrixPathway)[idxOrdered]
  
  matrixPathwayOrdered <- matrixPathway[genesOrdered,]
  sg <- list()
  aps <- list()
  for(k in 1:nrow(matrixPathwayOrdered)){
    if(k==1){
      ap <- colSums(t(as.matrix(matrixPathwayOrdered[1:k,])))/sqrt(k)
    }else{
      ap <- colSums(matrixPathwayOrdered[1:k,])/sqrt(k)
    }
    apClass1 <- ap[1:length(patientsClass1)]
    apClass2 <- ap[(length(patientsClass1)+1):(length(patientsClass1)+length(patientsClass2))]
    sg[[k]] <- t.test(apClass1,apClass2)$statistic
    aps[[k]] <- ap
  }
  
  ####selezione di K
  features <- ""
  for(index in 1:length(sg)){
    if(index == length(sg)){
      features <- aps[[index]]
    }else{
      if(sg[[index+1]] < sg[[index]]){
        features <- aps[[index]]
        break
      }
    }
  }
  IPFeatures[[classType]] <- cbind(IPFeatures[[classType]],features)
  colnames(IPFeatures[[classType]])[ncol(IPFeatures[[classType]])] <- names(rpm[i])
}

combinations <- as.character(unique(permResults$significantlyPathways$classes))

IPFeaturesResults <- matrix(0,nrow = length(combinations),ncol=2)
rownames(IPFeaturesResults) <- combinations
colnames(IPFeaturesResults) <- c("Class","Accuracy")
for(comb in combinations){
  path <- paste("p",comb,sep="-")
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathwayIPFeatures",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  metaClasses <- IPFeatures[[comb]]
  
  dir.create(paste(temporanyDirectory,"/",path,sep=""))
  write.table(metaClasses,paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)
  
  cols <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]),which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
  
  labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
  labels2classes <- labelsPS3[cols,"Classes"]
  labels2classes <- as.data.frame(labels2classes)
  rownames(labels2classes) <- rownames(labelsPS3)[cols]
  colnames(labels2classes) <- "x"
  
  write.table(labels2classes,file=labelPath ,quote=F)
  
  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core = 8,outer_fold = 3,inner_fold = 2,"C"))
  
  accuracies <- getPathwaysAccuracies(temporanyDirectory)
  
  svmWeights.pathwayIPFeatures <- getSVMweights(temporanyDirectory)
  
  IPFeaturesResults[comb,"Class"] <- comb
  IPFeaturesResults[comb,"Accuracy"] <- accuracies$accuraciesMatrix$accuracy[1]
  
  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)
}


###############################################################
############## FINE INFERRING PATHWAY ACTIVITY ################
###############################################################



######################################################################
################# ALL GENES SVM ######################################
######################################################################
source("R/classifyByGenesList.R")
combinations <- as.vector(unique(permResults$permutationResults$classes))

allGenesResults <- list()
for(i in combinations){
  allGenesResults[[i]] <- classifyByGenesList(pythonPath,i,rownames(PS3),PS3,labelsPS3,core = 8, outer_fold = 3 ,inner_fold = 2)
}

######################################################################
################# FINE ALL GENES SVM #################################
######################################################################


######################################################################
################ All results #########################################
######################################################################

allRes <- matrix(0,nrow = length(unique(permResults$permutationResults$classes)) ,ncol=4)
rownames(allRes) <- unique(permResults$permutationResults$classes)

for(i in 1:length(correlationFeaturesResults[,"Accuracy"])){
  allRes[names(correlationFeaturesResults[,"Accuracy"])[i],1] <- correlationFeaturesResults[names(correlationFeaturesResults[,"Accuracy"])[i],"Accuracy"]
}

for(i in 1:length(medianFeaturesResults[,"Accuracy"])){
  allRes[names(medianFeaturesResults[,"Accuracy"])[i],2] <- medianFeaturesResults[names(medianFeaturesResults[,"Accuracy"])[i],"Accuracy"]
}

for(i in 1:length(metaFeaturesResults[,"Accuracy"])){
  allRes[names(metaFeaturesResults[,"Accuracy"])[i],3] <- metaFeaturesResults[names(metaFeaturesResults[,"Accuracy"])[i],"Accuracy"]
}

for(i in 1:nrow(ImportantGenesDF)){
  allRes[rownames(ImportantGenesDF)[i],4] <- ImportantGenesDF[rownames(ImportantGenesDF)[i],"accuracy"]
}

allRes <- cbind(allRes,0,0,0,0)
for(comb in combinations){
  permResClasses <- permResults$significantlyPathways[which(permResults$significantlyPathways$classes==comb),]
  maximum <- which(permResClasses$accuracy == max(permResClasses$accuracy))
  maxPathways <- permResClasses[maximum,]
  maxPathwaysString <- ""
  for(i in 1:nrow(maxPathways)){
    maxPathwaysString <- paste(maxPathwaysString,maxPathways[i,"namePathway"]," (",maxPathways[i,"accuracy"],")",sep="")
    if(i != nrow(maxPathways)){
      maxPathwaysString <- paste(maxPathwaysString,"-",sep=" ")
    }
  }
  allRes[comb,5] <- maxPathwaysString
}

for(comb in combinations){
  combProbabilitiesForClasses <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == comb),]
  if(nrow(combProbabilitiesForClasses) == 0){
    allRes[comb,6] <- 0
  }else{
    maximum <- which(combProbabilitiesForClasses$accuracy == max(combProbabilitiesForClasses$accuracy))
    maxCombinations <- combProbabilitiesForClasses[maximum,]
    maxCombinationsString <- ""
    for(i in 1:nrow(maxCombinations)){
      maxCombinationsString <- paste(maxCombinationsString,maxCombinations[i,"combinationPathways"]," (",maxCombinations[i,"accuracy"],")",sep="")
      if(i != nrow(maxPathways)){
        maxCombinationsString <- paste(maxCombinationsString,"-",sep=" ")
      }
    }
    allRes[comb,6] <- maxCombinationsString
  }
}

for(comb in combinations){
  allRes[comb,7] <- allGenesResults[[comb]]$accuracy
}

for(comb in rownames(IPFeaturesResults)){
  allRes[comb,8] <- IPFeaturesResults[comb,"Accuracy"]
}

allRes <- as.data.frame(allRes)
allRes$classes <- rownames(allRes)

colnames(allRes) <- c("Correlation","Median","MetaFeatures","ImportantGenes","BestPathways","BestCombinationsPathways","AllGenesSVM","InferringPathway","Classes")

allRes <- allRes[,c("Correlation","Median","MetaFeatures","InferringPathway","AllGenesSVM","ImportantGenes","BestPathways","BestCombinationsPathways","Classes")]
library(WriteXLS)
WriteXLS(as.data.frame(allRes),ExcelFileName ="PS3/FinalClassification.xls")

allResNumber <- allRes
allResNumber <- as.matrix(allResNumber)
for(i in 1:nrow(allResNumber)){
  allResNumber[i,7] <- as.numeric(strsplit(strsplit(as.character(allResNumber[i,7]),split="[(]")[[1]][2],split="[)]")[[1]][[1]])
  if(allResNumber[i,8] != 0){
    allResNumber[i,8] <- strsplit(strsplit(as.character(allResNumber[i,8]),split="[(]")[[1]][2],split="[)]")[[1]][[1]]
  }
}

allResNumber <- as.data.frame(allResNumber)

new.df<-melt(allResNumber,id.vars="Classes")
new.df$value <- as.numeric(new.df$value)
finalClassificationPlot <- ggplot(new.df,aes(Classes,value,fill=variable))+geom_histogram(stat="identity",position="dodge") + ggtitle("Accuracy") + xlab("Class combination") + ylab("Accuracy")
ggsave(finalClassificationPlot ,filename = paste("PS3/FinalClassificationHistogram.png",sep=""), width = 14, height = 12)

save.image(file="PS3/RDATA/PS3-NO-FS-DopoConfrontoClassificazione.RData")

######################################################################
################ fine all results ####################################
######################################################################


######################################################################
################ AGGIUNTA CTD INFORMATION ############################
######################################################################
library(readr)
CTD_diseases_pathways <- read_csv("../Data/CTD_diseases_pathways.csv",col_names = FALSE, skip = 29)
names(CTD_diseases_pathways) <- c("disease","mesh","namePathway","IDPathway","gene")
signP <- as.character(permResults$significantlyPathways$namePathway)

### rimuovo le classi dai pathway
for(i in 1:length(signP)){
  signP[[i]] <- strsplit(signP[[i]],split="-")[[1]][1]
}

### uniformo gli ID
for(i in 1:length(CTD_diseases_pathways$IDPathway)){
  if(startsWith(CTD_diseases_pathways$IDPathway[i],"REACT")){
    CTD_diseases_pathways$IDPathway[i] <- paste("REACTOME_",substring(CTD_diseases_pathways$IDPathway[i],13,nchar(CTD_diseases_pathways$IDPathway[i])),sep="")
  }
  
  if(startsWith(CTD_diseases_pathways$IDPathway[i],"KEGG")){
    CTD_diseases_pathways$IDPathway[i] <- paste("KEGG_",substring(CTD_diseases_pathways$IDPathway[i],6,nchar(CTD_diseases_pathways$IDPathway[i])),sep="")
  }
}

SignPathways$knowGenes = ""
SignPathways$weightsPathway = ""
for(i in 1:nrow(SignPathways)){
  toConcat <- ""
  weightPath <- 0
  namePathway <- as.character(SignPathways[i,"pathway"])
  weightsGenes <- names(reductionPathwaysRes$svmWeights[[namePathway]])
  genesPathway <- transformInSymbols(weightsGenes)
  
  namePathway <- strsplit(namePathway,split="-")[[1]][1]
  CTDPath <- CTD_diseases_pathways[which(CTD_diseases_pathways$IDPathway == namePathway),]
  if(nrow(CTDPath) == 0){
    toConcat <- "Pathway not available"
    weightPath <- 0
  }else{
    CTDPath <- CTDPath[which(grepl("lung",CTDPath$disease,ignore.case = T)),]
    if(nrow(CTDPath) == 0){
      toConcat <- "No genes connected to disease"
      weightPath <- 0
    }else{
      genesByCTD <- unique(CTDPath$gene)
      
      positions <- which(genesPathway %in% genesByCTD)
      if(length(positions) == 0){
        toConcat <- "No genes connected to disease"
        weightPath <- 0
      }else{
        weightPath <- 0
        for(j in 1:length(positions)){
          toConcat <- paste(toConcat,genesPathway[positions[j]],"(",positions[j],"/",SignPathways[i,"sizeReduced"],") ",sep="")
          weightPath <- weightPath + (1-(positions[j]/SignPathways[i,"sizeReduced"]))
        }
      }
    }
  }
  SignPathways[i,"weightsPathway"] <- as.character(weightPath)
  SignPathways[i,"knowGenes"] <- toConcat
}


WriteXLS::WriteXLS(SignPathways,"PS3/PS3-SignificantlyPathwaysWithGenesCTD.xls")

######################################################################
################ FINE CTD INFORMATION ################################
######################################################################


########################################
############## RFE #####################
########################################

valutationCouples <- function(bestCouples,pathCombClass,SignPathways){
  if(nrow(bestCouples) > 1){
    pathwayCombinationValue <- list()
    for(i in 1:nrow(bestCouples)){
      validationPathways <- list()
      pathways <- as.character(bestCouples[i,"combinationPathways"])
      pp <- strsplit(pathways,split="/")[[1]]
      pathwayCombinationValue[[pathways]] <- sum(as.numeric(SignPathways[pp,"weightsPathway"]))
      
    }
    idxBestCouples <- which.max(pathwayCombinationValue)[[1]]
    
  }else{
    idxBestCouples <- 1 
  }
  return(bestCouples[idxBestCouples,])
}

allCombinations <- pathCombProbabilities$combination
combinations <- as.character(unique(pathCombProbabilities$combination$class))
RFEFinalResults <- list()
for(comb in combinations){
  pathCombClass <- pathCombProbabilities$combination[which(pathCombProbabilities$combination$class == comb),]
  bestAccuracy <- pathCombClass[which.max(pathCombClass$accuracy),"accuracy"]
  
  bestCouples <- pathCombClass[which(pathCombClass$accuracy == bestAccuracy),]
  
  ### scelta della migliore combinazione
  bestComb <- valutationCouples(bestCouples,pathCombClass,SignPathways)
  bestComb <- cbind(bestComb,"size"=2)
  pathComb <- SignPathways[SignPathways$class == comb,]
  
  pathToDelete <- c(strsplit(as.character(bestComb[1,1]),split="/")[[1]])
  idxToDelete <- which(rownames(pathComb) %in% pathToDelete)
  
  pathCombMin <- pathComb[-idxToDelete,]
  matrixFeatures <- matrix(0,nrow = nrow(reductionPathwaysRes$probabilities[pathToDelete][[1]]),ncol=2)
  rownames(matrixFeatures) <- rownames(reductionPathwaysRes$probabilities[pathToDelete][[1]])
  matrixFeatures[,1] <- reductionPathwaysRes$probabilities[pathToDelete][[1]]$probabilities
  matrixFeatures[,2] <- reductionPathwaysRes$probabilities[pathToDelete][[2]]$probabilities
  colnames(matrixFeatures) <- pathToDelete
  
  
  fResults <- bestComb
  
  ## Ciclo per numero di classificatione da effettuare
  for(i in 1:(nrow(pathComb)-2)){
    results <- c()
    
    temporanyDirectory <- getwd()
    temporanyDirectory <- paste(temporanyDirectory,"/RFE",sep="")
    temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
    dir.create(temporanyDirectory)
    
    ## Ciclo per vedere qual è la miglior combinazione (aggiunta dei pathway alla matrice)
    for(j in 1:(nrow(pathCombMin))){
      path <- paste(j,comb,sep="-")
      
      pathToAdd <- as.character(pathCombMin[j,"pathway"])
      newMatrixFeatures <- cbind(matrixFeatures,reductionPathwaysRes$probabilities[pathToAdd][[1]]$probabilities)
      colnames(newMatrixFeatures)[ncol(newMatrixFeatures)] <- pathToAdd
      
      dir.create(paste(temporanyDirectory,"/",path,sep=""))
      write.table(newMatrixFeatures,file = paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)
      
      
      labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
      idx <- c(which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][1]), which(labelsPS3$Classes == strsplit(comb,split="vs")[[1]][2]))
      labels2Classes <- labelsPS3[idx,"Classes"]
      labels2Classes <- as.data.frame(labels2Classes)
      rownames(labels2Classes) <- rownames(labelsPS3)[idx]
      colnames(labels2Classes) <- "x"
      write.table(labels2Classes,file = labelPath ,quote=F)
    }
    
    system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core = 8,outer_fold = 3,inner_fold = 2,"C"))
    
    accuracies <- getPathwaysAccuracies(temporanyDirectory)
    accuracies <- as.matrix(accuracies$accuraciesMatrix)
    svmWeights.RFEFeatures <- getSVMweights(temporanyDirectory)
    for(k in 1:length( svmWeights.RFEFeatures)){
      accuracies[k,"pathway"] <- as.character(paste(names(svmWeights.RFEFeatures[[k]]),collapse = "/"))
    }
    
    results <- as.data.frame(accuracies)
    
    unlink(temporanyDirectory,recursive=T)
    unlink(temporanyDirectory_res,recursive = T)
    
    colnames(results) <- c("combinationPathways","class","accuracy","size")
    bestAccuracy <-  as.character(results[which.max(results$accuracy),"accuracy"])
    
    bestCouples <- results[which(results$accuracy == bestAccuracy),]
    
    ### scelta della migliore combinazione
    bestComb <- valutationCouples(bestCouples,pathCombClass,SignPathways)
    pathwaysToRemoveAll <- strsplit(x = as.character(bestComb$combinationPathways),split="/")[[1]]
    
    idxToRem <- which(rownames(pathComb) %in% pathwaysToRemoveAll)
    pathCombMin <- pathComb[-idxToRem,]
    
    fResults <- rbind(fResults,bestComb)
    
    pToAdd <- strsplit(as.character(fResults[i+1,"combinationPathways"]),split="/")[[1]][!strsplit(as.character(fResults[i+1,"combinationPathways"]),split="/")[[1]] %in% strsplit(as.character(fResults[i,"combinationPathways"]),split="/")[[1]]]
    
    matrixFeatures <- cbind(matrixFeatures,reductionPathwaysRes$probabilities[pToAdd][[1]]$probabilities)
    colnames(matrixFeatures)[ncol(matrixFeatures)] <- pToAdd
  }
  RFEFinalResults[[comb]] <- fResults
}


dir.create("PS3/RFE")
dir.create("PS3/RFE/plots")
for(i in 1:length(RFEFinalResults)){
  toPlot <- RFEFinalResults[[i]]
  
  toPlot$accuracy <- as.numeric(toPlot$accuracy) 
  toPlot$size <- as.numeric(toPlot$size)
  RFEPlot <- ggplot(data=toPlot, aes(x=size, y=accuracy, group=1)) +
    geom_line(linetype = "dashed",color ="grey")+
    geom_point()+ xlab("Size") + ylab("Accuracy")+ ylim(0.5, 1)+ 
    ggtitle(paste("RFE",names(RFEFinalResults)[i])) 
  ggsave(RFEPlot,filename = paste("PS3/RFE/plots/",names(RFEFinalResults)[i],".png",sep=""), width = 14, height = 12)
  
  write.table(RFEFinalResults[[i]],paste("PS3/RFE/",names(RFEFinalResults)[i],".xls"),row.names = F)
}

allDataForPlot <- data.frame()
for(i in 1:length(RFEFinalResults)){
  toAddNow <- RFEFinalResults[[i]][,c("accuracy","class","size")]
  rownames(toAddNow) <- NULL
  allDataForPlot <- rbind(allDataForPlot,toAddNow)
}
allDataForPlot$accuracy <- as.numeric(allDataForPlot$accuracy)
allDataForPlot$class <- as.factor(allDataForPlot$class)
allDataForPlot$size <- as.numeric(allDataForPlot$size)

RFEPlotAll <- ggplot(data=allDataForPlot, aes(x=size, y=accuracy, colour=class)) +
  geom_line(linetype = "dashed",size=1)+
  geom_point()+ xlab("Size") + ylab("Accuracy")+ ylim(0.5, 1)+ 
  ggtitle("RFE All combinations") 

ggsave(RFEPlotAll,filename = paste("PS3/RFE/plots/AllCombinations.png",sep=""), width = 14, height = 12)

save.image("PS3/RDATA/PS3-AFTER-RFE.RData")

########################################
############## FINE RFE ################
########################################