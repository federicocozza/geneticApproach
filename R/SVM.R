require(ggplot2)
source("R/getResultsPathways.R")
SVM <- function(pythonPath, allPathways, labelsSNF,core, outer_fold ,inner_fold){
  allPairPathways <- list()
  for(i in 1:length(allPathways)){
    allPairPathways <- c(allPairPathways,splitDataset(names(allPathways)[i],allPathways[[names(allPathways)[i]]],labelsSNF$Classes))
  }
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathway",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res/",sep="")

  dir.create(temporanyDirectory)
  for(i in 1:length(allPairPathways)){
    newPath <- paste(temporanyDirectory,"/",names(allPairPathways)[i],sep="")
    dir.create(newPath)
    write.table(t(allPairPathways[[names(allPairPathways)[i]]]$X),file = paste(newPath,"/",names(allPairPathways)[i],".txt",sep=""),quote=F)
    write.table(allPairPathways[[names(allPairPathways)[i]]]$y, file = paste(newPath,"/labels",sep="") ,quote=F)
  }

  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core,outer_fold,inner_fold,"C"))


  accuracies <- getPathwaysAccuracies(temporanyDirectory)

  svmWeights <- getSVMweights(temporanyDirectory)

  probabilities <- getSVMProbabilities(temporanyDirectory)

  HAplot <- heatmapAccuracyPlot(accuracies$accuraciesMatrix)


  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)

  return(list(accuracies = accuracies,svmWeights = svmWeights,probabilities = probabilities,heatmapAccuracy = HAplot))
}





splitDataset = function(namePathway,pathway,labels){
  listSplit = list()
  numberOfClasses <- length(unique(labels))
  if(numberOfClasses < 2){
    stop("Can't classify without 2 classes")
  }
  if(numberOfClasses >= 2){
    for(i in 1:(numberOfClasses-1)){
      for(j in (i+1):(numberOfClasses)){
        idx = c(which(labels == unique(labels)[i]), which(labels == unique(labels)[j]))
        pathway2classes <- pathway[,idx]
        labels2classes <- labels[idx]
        listSplit[[paste(namePathway,"-",unique(labels)[i],"vs",unique(labels)[j],sep="")]] = list(X = pathway2classes,y = labels2classes)
      }
    }
  }
  return(listSplit)
}



