getBestEnrichment <- function(dataForPlot,numberOfGenes,numberOfPathway,piVal, enrichResults,DFPresults){

  dataForPlot <- dataForPlot[which(dataForPlot[,"numberOfGenes"] == numberOfGenes),]
  dataForPlot <- dataForPlot[which(dataForPlot[,"numberOfPathway"] == numberOfPathway),]
  dataForPlot <- dataForPlot[which(as.character(dataForPlot[,"piVal"]) == as.character(piVal)),]

  if(nrow(dataForPlot) > 1){
    warning("there are more then 1 combinations")
  }
  results <- list()
  for(i in as.numeric(rownames(dataForPlot))){
    results[[length(results)+1]] <- list(enrichRes = enrichResults[[i]],disciminantGenes = rownames(DFPresults[[i]]$discriminantGenes))
  }

  for(i in 1:length(results)){
    row1 <-  cbind(results[[i]]$enrichRes$KEGGPathway@result$Description,results[[i]]$enrichRes$KEGGPathway@result$ID,results[[i]]$enrichRes$KEGGPathway@result$Count,results[[i]]$enrichRes$KEGGPathway@result$p.adjust,results[[1]]$enrichRes$KEGGPathway@result$GeneRatio,"KEGG")
    row2 <- cbind(results[[i]]$enrichRes$ReactomePathway@result$Description,results[[i]]$enrichRes$ReactomePathway@result$ID,results[[i]]$enrichRes$ReactomePathway@result$Count,results[[i]]$enrichRes$ReactomePathway@result$p.adjust,results[[1]]$enrichRes$ReactomePathway@result$GeneRatio,"Reactome")
    d <- rbind(row1,row2)
    colnames(d) <- c("namePathway","idPathway" ,"count","pValueAdjust","geneRatio","typePathway")
    d <- as.data.frame(d)
    d[,"count"] <- as.numeric(as.vector(d[,"count"]))
    d[,"pValueAdjust"] <- as.numeric(as.vector(d[,"pValueAdjust"]))
    d[,"geneRatio"] <- sapply(as.vector(d[,"geneRatio"]), function(x) eval(parse(text=x)))

    p <- ggplot(data = d,aes(geneRatio,namePathway,colour = pValueAdjust)) +
      geom_point(aes(size=count)) +
      geom_polygon(aes(geneRatio,namePathway,fill=typePathway)) +
      ggtitle("Significant Pathways") +
      xlab("Gene Ratio") +
      ylab("Pathways")+
      scale_fill_manual(values = c('Reactome'='red','KEGG'='black'))+
      theme(axis.text.y = element_text(color=d$typePathway))

      results[[i]]$plot <-  p
      IDs <- c(bestEnrich[[1]]$enrichRes$KEGGPathway@result$ID,bestEnrich[[1]]$enrichRes$ReactomePathway@result$ID)
  }

  return(results)
}



