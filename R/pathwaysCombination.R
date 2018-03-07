
pathwaysCombination <- function(pythonPath,probabilities,significantlyPathways,labelsClass,core=8,outer_fold=3,inner_fold=2,maxNumCombination = 2){
  temporanyDirectory <- getwd()
  temporanyDirectory <- paste(temporanyDirectory,"/pathwaysCombination",sep="")
  temporanyDirectory_res <- paste(temporanyDirectory,"_res",sep="")
  dir.create(temporanyDirectory)
  createStackingLevel2Combination(significantlyPathways,temporanyDirectory,probabilities,labelsClass,maxNumCombination)
  if(identical(list.files(temporanyDirectory), character(0))){
    print("No combination possible...")
    unlink(temporanyDirectory,recursive=T)
    unlink(temporanyDirectory_res,recursive = T)
    return()
  }
  system(paste(pythonPath,'Python/stacking_level2_multi.py',temporanyDirectory,core,outer_fold,inner_fold))

  r <- getStackingLevel2Results(temporanyDirectory,significantlyPathways)

  unlink(temporanyDirectory,recursive=T)
  unlink(temporanyDirectory_res,recursive = T)
  return(r)
}


createStackingLevel2Combination <- function(significantlyPathways,path,probabilities,labelsClass,maxNumCombination){

  classi <- unique(significantlyPathways[,"classes"])

  for(j in classi){
    if(length(which(significantlyPathways[,"classes"]==j)) <= 1){
      next
    }
    curClasses <- paste(path,"/Stacking_",j,sep="")
    dir.create(curClasses)
    ptws <- significantlyPathways[which(significantlyPathways[,"classes"]==j),]
    rownames(ptws) <- ptws[,"namePathway"]
    maxComb <- min(nrow(ptws),maxNumCombination)
    for(k in 2:maxComb){
      dir.create(paste(path,"/Stacking_",j,"/size",k,sep=""))
      classCombin <- combn(ptws[,"namePathway"],k)
      for(z in 1:ncol(classCombin)){
        dir.create(paste(path,"/Stacking_",j,"/size",k,"/combination-",z,sep=""))
        newFeatures <- c()
        cnames <- c()
        for(i in 1:k){
          if(i == 1){
            labelPath <- paste(path,"/Stacking_",j,"/size",k,"/combination-",z,"/labels",sep="")
            idx <- c(which(labelsClass == strsplit(x = as.character(j),split = "vs")[[1]][1]), which(labelsClass == strsplit(x = as.character(j),split = "vs")[[1]][2]))
            labels2classes <- labelsClass[idx,"Classes"]
            labels2classes <- as.data.frame(labels2classes)
            rownames(labels2classes) <- rownames(labelsClass)[idx]
            colnames(labels2classes) <- "x"
            write.table(labels2classes,file=labelPath ,quote=F)

          }

          newFeatures <- cbind(newFeatures,probabilities[[as.character(classCombin[i,z])]][,"probabilities"])
          cnames <- c(cnames,as.character(ptws[as.character(classCombin[i,z]),"namePathway"]))
        }
        colnames(newFeatures) <- cnames
        rownames(newFeatures) <- rownames(probabilities[[as.character(classCombin[i,z])]])

        write.table(newFeatures,paste(path,"/Stacking_",j,"/size",k,"/combination-",z,"/features.txt",sep=""),quote=F)
      }
    }
  }
}

getStackingLevel2Results <- function(path,significantlyPathways){
  library(reshape2)
  library(ggplot2)
  path_res <- paste(path,"_res/",paste="",sep="")
  results <- data.frame(combinationPathways = character(),accuracy = double(),class =  character())
  stacking <- paste(path_res, list.files(path_res),sep="")
  for(s in stacking){
    sizes <- paste(s,"/",list.files(s),sep="")
    for(size in sizes){
      combinations <- paste(size,"/",list.files(size),sep="")
      for(comb in combinations){
        accuracy <- read.table(file=paste(comb,"/test_avg_accuracy.txt",sep=""),sep=",",header=T)$avg_acc
        pathways <- read.table(paste(gsub("_res","",comb),"/features.txt",sep=""),header=T,check.names = F)
        results <- rbind(results,data.frame(combinationPathways = paste(colnames(pathways),collapse="/"),accuracy = accuracy,class=strsplit(strsplit(paste(colnames(pathways),collapse="/"),split = "/")[[1]],split="-")[[1]][2]))
      }
    }
  }

  AccMat <-  matrix(0,nrow=nrow(significantlyPathways),ncol=nrow(results))
  row.names(AccMat) <- significantlyPathways$namePathway
  for(i in 1:nrow(results)){
    ptws <- unlist(strsplit(as.character(results[i,"combinationPathways"]),split="/"))
    AccMat[ptws,i] <- results[i,"accuracy"]
  }

  plist = list()

  for (comb in unique(results$class)){
    AccMat2 <- AccMat[,which(results[,"class"] %in% comb == T)]
    if(class(AccMat2) == "numeric"){
      next
    }
    toRem <- which(rowSums(AccMat2,na.rm=T) == 0)
    AccMat2 = AccMat2[-toRem,]
    mat.melted <- melt(AccMat2)
    colnames(mat.melted) <- c("Pathways","Combination","value")
    p <- ggplot(mat.melted, aes(x = Combination, y = Pathways, fill = value)) + ggtitle(paste("Class",comb)) + geom_tile() +
      scale_fill_gradientn(limits =c(min(as.vector(AccMat[AccMat>0])),max(as.vector(AccMat))),colors = c("white","red"))

    plist[[length(plist)+1]] = p

  }

  return(list(combination = results,plots = plist))
}

