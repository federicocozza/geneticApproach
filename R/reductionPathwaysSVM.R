
source("R/getResultsPathways.R")

reductionPathwaysSVM <- function(pythonPath,SVMWeights,dataset,labelsSNF,percentage = 0.8,core, outer_fold ,inner_fold){
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathwayReduced",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res/",sep="")
  dir.create(temporanyDirectory)
  for(i in 1:length(SVMWeights)){
    finalRank <- SVMWeights[[i]]
    ng <-  sum(cumsum(finalRank) < percentage)
    ggenes <-  names(finalRank)[1:ng]
    ggenes <- gsub("[.]", "-", ggenes)

    path <- names(SVMWeights)[i]

    labelPath <- paste(temporanyDirectory,"/",path,"/labels",sep="")
    idx <- c(which(labelsSNF == strsplit(x = as.character(strsplit(path,split="-")[[1]][2]),split = "vs")[[1]][1]), which(labelsSNF == strsplit(x = as.character(strsplit(path,split="-")[[1]][2]),split = "vs")[[1]][2]))
    labels2classes <- labelsSNF[idx,"Classes"]
    labels2classes <- as.data.frame(labels2classes)
    rownames(labels2classes) <- rownames(labelsSNF)[idx]
    colnames(labels2classes) <- "x"
    patients <- rownames(labels2classes)

    small_mat <-  dataset[ggenes,patients]
    dir.create(paste(temporanyDirectory,"/",path,sep=""))
    write.table(t(small_mat),paste(temporanyDirectory,"/",path,"/",path,".txt",sep=""),quote=F)

    write.table(labels2classes,file=labelPath ,quote=F)
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

