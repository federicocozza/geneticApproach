classifyByGenesList <- function(pythonPath,classType,genes,dataset,classLabels,core = 8, outer_fold = 3 ,inner_fold = 2){
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/ImportantGenes",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res/",sep="")
  
  dir.create(temporanyDirectory)

  newPath <- paste(temporanyDirectory,"/ImportantGenes-",classType,sep="")

  idx <- c(which(classLabels$Classes == strsplit(classType,split="vs")[[1]][1]),which(classLabels$Classes == strsplit(classType,split="vs")[[1]][2]))
  labels2classes <- classLabels[idx,"Classes"]
  labels2classes <- as.data.frame(labels2classes)
  rownames(labels2classes) <- rownames(classLabels)[idx]
  colnames(labels2classes) <- "x"
  patients <- rownames(labels2classes)
  
  dir.create(newPath)

  dataToWrite <- t(dataset[genes,patients])
  if(length(genes) == 1){
    dataToWrite <- as.data.frame(dataToWrite)
    rownames(dataToWrite) <- genes
    dataToWrite <- t(dataToWrite)
  }
  write.table(dataToWrite,file = paste(newPath,"/","ImportantGenes-",classType,".txt",sep=""),quote=F)
  write.table(labels2classes, file = paste(newPath,"/labels",sep="") ,quote=F)

  system(paste(pythonPath,'Python/stacking_level1_multi.py',temporanyDirectory,core,outer_fold,inner_fold,"C"))
  
  accuracies <- getPathwaysAccuraciesWithoutPlots(temporanyDirectory)
  svmWeights <- getSVMweights(temporanyDirectory)
  probabilities <- getSVMProbabilities(temporanyDirectory)
  
  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)

  return(list(accuracy = accuracies$accuraciesMatrix$accuracy,svmWeights = svmWeights,probabilities = probabilities))
}


getPathwaysAccuraciesWithoutPlots <- function(path){
  path_res <- paste(path,"_res/",sep="")
  pathways <- list.files(path_res)
  p1 <- paste(paste(path,pathways,sep="/"),"/",pathways,".txt",sep="")
  res <-  data.frame()
  for(i in 1:length(pathways)){
    x <- read.table(p1[i])
    p <-  pathways[i]
    a <- read.table(paste(path_res,p,"/","test_avg_accuracy.txt",sep=""),sep=",",header = T)
    psplit <- strsplit(p,split="-")
    res <- rbind(res,data.frame(p,psplit[[1]][2],a[["avg_acc"]],ncol(x)))
  }
  
  colnames(res) <- c("pathway","class","accuracy","size")
  
  return(list(accuraciesMatrix = res))
}

  