getPathwaysAccuracies <- function(path){
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

  plotHist <- ggplot(res, aes(accuracy)) + geom_histogram(aes(fill=class),col="black") + ggtitle("Pathways Accuracy") +
    theme_grey() + stat_bin(geom="text",colour="black",size=3.5,aes(label = ifelse(..count..>0,..count..,""),fill=class),position = position_stack(vjust=0.5)) +scale_x_continuous(breaks = seq(0,1,0.05),limits=c(0,1.1))


  plotRegression <- ggplot(res,aes(x=accuracy,y=size,color = class,size=accuracy)) + geom_point() + geom_smooth(method =lm,se=T,fullrange=F,aes(fill=class),show.legend=T)+
    xlab("Accuracy") + ylab("Size") +  ggtitle("Accuracy vs Pathways Size")


  return(list(accuraciesMatrix = res, plotHist=plotHist,plotRegression = plotRegression))
}


getSVMweights <- function(path){
  path_res <- paste(path,"_res",sep="")
  path_res <- paste(path_res,"/",list.files(path_res),sep="")

  genes_weight <- list()

  pathways <- list.files(path)

  pathPathway <- paste(path,"/",list.files(path),"/",list.files(path),".txt",sep="")
  for(i in 1:length(path_res)){
    finalPathRes <- path_res[i]
    files <- list.files(finalPathRes)
    path_mat = read.table(pathPathway[i],check.names = F)
    foldf <- c()
    for(f in files){
      if(length(grep(pattern="_ranked_genes_",f)) == 1){
        foldt <-  read.table(paste(path_res[i],"/",f,sep=""),quote="\"")
        colnames(foldt) <- colnames(path_mat)
        foldt <- abs(foldt)
        foldt <- foldt/sum(foldt)
        foldf <- rbind(foldf,foldt)
      }
    }
    finalRank = colMeans(foldf)
    finalRank = finalRank[order(finalRank,decreasing = TRUE)]
    genes_weight[[pathways[i]]] <- finalRank
  }
  return(genes_weight)
}


getSVMProbabilities <- function(path){
  path_res <- paste(path,"_res/",sep="")
  pathway <- list.files(path_res)
  path_res <- paste(path_res,list.files(path_res),"/level1_features.txt",sep="")

  path_p <- paste(path,"/",list.files(path),"/",list.files(path),".txt",sep="")
  prob <- list()

  for(i in 1:length(path_res)){
    f <- path_res[i]
    probabilities <- as.data.frame(read.table(f,sep=",",header=T)[,2])
    colnames(probabilities) <- c("probabilities")
    rownames(probabilities) <- rownames(read.table(path_p[i]))

    prob[[pathway[[i]]]] <-  probabilities
  }

  return(prob)
}

heatmapAccuracyPlot <- function(accuraciesMatrix){
  toPlot <- accuraciesMatrix
  toPlot[,"pathway"] <- as.character(as.vector(toPlot[,"pathway"]))
  for(i in 1:nrow(toPlot)){
    toPlot[i,"pathway"] <- strsplit(as.character(toPlot[i,"pathway"]),split="-")[[1]][[1]]
  }
  p <- ggplot(toPlot, aes(x = class, y = pathway, fill = accuracy)) + ggtitle("Heatmap Accuracy") + geom_tile() +scale_fill_gradientn(colors = c("white","red"))
  return(p)
}
