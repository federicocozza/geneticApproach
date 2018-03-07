
permutationTest <- function(pythonPath,numberPermutation,accuraciesMatrix,SVMWeighs,dataSet,labels,pValueCutOff,core= 8,outer_fold=3,inner_fold=2){
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/permutationTest",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  createPermutationDataset(accuraciesMatrix,dataset,labels,temporanyDirectory,numberPermutation)

  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core,outer_fold,inner_fold,"P"))

  permRes <- getPermutationResults(accuraciesMatrix,temporanyDirectory_res,numberPermutation)

  classes <- c()
  for(i in permRes$namePathway){
    classes <- c(classes,strsplit(as.character(i),split="-")[[1]][2])
  }
  permRes <- cbind(permRes,classes)

  significantlyPathways <- permRes[permRes$pValue < pValueCutOff,]

  #significantlyPathwaysHeatmap <- significantlyPathwaysPlot(significantlyPathways)

  accuracyVsPathwaysSizePlot <- accuracyVsPathwaysSizePlot(significantlyPathways,permRes,pValueCutOff)

  #SVMWeightsRelationPlot <- SVMWeightsRelationPlot(SVMWeights,permRes,significantlyPathways)

  #heatmapIntersection <- heatmapIntersection(SVMWeights,permRes,significantlyPathways)

  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)

  return(list(permutationResults = permRes,significantlyPathways = significantlyPathways,accuracyVsPathwaysSizePlot = accuracyVsPathwaysSizePlot))
}


createPermutationDataset <- function(dimensionsPathways,dataset,labels,path,numberPermutation){
  for(i in 1:dim(dimensionsPathways)[[1]]){
    d <- dimensionsPathways[i,]
    newPath <- paste(path,"/",d$size,"-",d$class,sep="")
    if(! dir.exists(newPath)){
      dir.create(newPath)
      for(j in 1:numberPermutation){
        permutationPath <- paste(newPath,"/permPath",j,".txt",sep="")
        labelPath <- paste(newPath,"/labels",sep="")
        idx <- c(which(labels == strsplit(x = as.character(d$class),split = "vs")[[1]][1]), which(labels == strsplit(x = as.character(d$class),split = "vs")[[1]][2]))
        labels2classes <- labels[idx,"Classes"]
        labels2classes <- as.data.frame(labels2classes)
        rownames(labels2classes) <- rownames(labels)[idx]
        colnames(labels2classes) <- "x"
        patients <- rownames(labels2classes)
        write.table(labels2classes,file=labelPath ,quote=F)

        a <- sample(rownames(dataset),d$size,replace = FALSE)

        write.table(t(dataset[a,patients]),file=permutationPath)
      }

    }
  }
}

getPermutationResults <- function(dimensionsPathways,path,nperm){
  res <- data.frame()
  for(i in 1:dim(dimensionsPathways)[[1]]){
    d <- dimensionsPathways[i,]
    finalPath <-  paste(path,"/",paste(d$size,"-",d$class,sep=""),"/test_avg_accuracy.txt",sep="")
    accuracies <- read.csv(finalPath,header = T)
    np <- length(which(d$accuracy > accuracies$avg_acc))
    res <- rbind(res,data.frame(d$pathway,1-(np/nperm),d$accuracy,d$size))
  }
  colnames(res) <- c("namePathway","pValue","accuracy","size")
  return(res)
}


significantlyPathwaysPlot <- function(significantlyPathways){
  toPlot <- significantlyPathways

  toPlot[,"namePathway"] <- as.character(as.vector(toPlot[,"namePathway"]))
  for(i in 1:nrow(toPlot)){
    toPlot[i,"namePathway"] <- strsplit(as.character(toPlot[i,"namePathway"]),split="-")[[1]][[1]]
  }

  p <- ggplot(toPlot, aes(x = classes, y = namePathway, fill = pValue)) + ggtitle("Significantly Pathways")+ ylab("Pathways") + xlab("Class combination")+ geom_tile() +scale_fill_gradientn(colors = c("white","red"))
  return(p)
}


accuracyVsPathwaysSizePlot <- function(significantlyPathways,permutationResults,pValueCutOff){
  toPlot <- permutationResults


  significance <- rep("significant",nrow(toPlot))
  significance[toPlot$pValue > pValueCutOff] = "not significant"
  toPlot <- cbind(toPlot,significance)

  percentage <- 1 - toPlot[,"pValue"]
  toPlot <- cbind(toPlot,percentage)

  correlationData = c()
  for(i in unique(toPlot$classes)){
    correlationData[[i]] <- cor(toPlot[toPlot$classes == i,]$accuracy,toPlot[toPlot$classes == i,]$size)
  }
  maxSize <- max(toPlot$size) + 5
  xPos <- min(toPlot$accuracy)+0.075

  positions <- seq(maxSize,maxSize+4*(length(unique(toPlot$classes))),4)

  correlationString <- c()
  for(i in 1:length(correlationData)){
    correlationString[[i]] <- paste(names(correlationData)[i]," r = ",round(as.numeric(correlationData[i]),2),sep="")
  }
  correlationString <- rev(correlationString)

  correlationString <- c(correlationString,"Correlation:")
  p <- ggplot(toPlot,aes(x=accuracy,y=size,color = classes)) + geom_point(aes(shape=significance,size=percentage)) +
    geom_smooth(method = lm,se=T,fullrange=F,aes(fill=classes)) +
    xlab("Accuracy") + ylab("Size") +
    ggtitle("Accuracy vs Pathways Size") +
    annotate("text",x=rep(xPos,length(unique(toPlot$classes))+1),y=positions,label = correlationString)

  return(p)
}


SVMWeightsRelationPlot <- function(SVMWeights,permResults,significantlyPathways){
  require(reshape)
  classes <- unique(permResults$classes)

  pathways <- c()
  for(i in 1:nrow(significantlyPathways)){
    pathways[[i]] <- strsplit(as.character(significantlyPathways[i,"namePathway"]),split="-")[[1]][[1]]
  }

  allPathways <- c()
  for(p in pathways){
    allPathways <- c(allPathways,paste(p,classes,sep="-"))
  }

  weights <- c()
  for(i in 1:length(SVMWeights[allPathways])){
    weights <- c(weights,SVMWeights[allPathways][[i]] * length(SVMWeights[allPathways][[i]]))
  }

  results <- list()
  for(p in pathways){
    pathwaysComb <- paste(p,classes,sep="-")
    genes <- c()
    for(pc in pathwaysComb){
      genes <- union(genes,names(SVMWeights[[pc]]))
    }

    d <- matrix(nrow=length(genes),ncol=length(classes),NA)
    colnames(d) <- pathwaysComb
    rownames(d) <- genes

    for(pc in pathwaysComb){
      d[names(SVMWeights[[pc]]),pc] = SVMWeights[[pc]] * length(SVMWeights[[pc]])
    }

    newColNames <- c()
    for(i in 1:length(colnames(d))){
      newColNames[i] <- strsplit(colnames(d)[i],split="-")[[1]][[2]]
    }
    colnames(d) <- newColNames
    toPlot <- melt(d)
    colnames(toPlot) <- c("genes","class","geneImportance")
    toPlot$genes <- as.factor(toPlot$genes)

    plot <- ggplot(toPlot, aes(x = class, y = genes, fill = geneImportance)) + ggtitle(p) + geom_tile() +
      scale_fill_gradientn(limits=c(min(weights),max(weights)),colors = c("white","red"))
    results[[p]] <- plot
  }
  return(results)
}

heatmapIntersection <- function(SVMWeights,permResults,significantlyPathways){
  require(reshape)
  classes <- unique(permResults$classes)

  allSignificantlyPathways <- SVMWeights[significantlyPathways[,"namePathway"]]

  m <- matrix(nrow = length(allSignificantlyPathways),ncol = length(allSignificantlyPathways), NA)
  colnames(m) <- names(allSignificantlyPathways)
  rownames(m) <- names(allSignificantlyPathways)

  for(i in 1:length(allSignificantlyPathways)){
    for(j in 1:length(allSignificantlyPathways))
      m[i,j] <- length(intersect(names(allSignificantlyPathways[[i]]),names(allSignificantlyPathways[[j]])))/length(union(names(allSignificantlyPathways[[i]]),names(allSignificantlyPathways[[j]])))
  }

  toPlot<- melt(m)
  colnames(toPlot) <- c("PathwayX","PathwayY","Intersection")
  p <- ggplot(toPlot, aes(x = PathwayX, y = PathwayY, fill = Intersection)) + ggtitle("Intersection between pathways") + geom_tile() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +scale_fill_gradientn(colors = c("white","red"))

  return(p)
}
