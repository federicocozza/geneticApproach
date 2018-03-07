library(plotly)

DFPandEnrichment <- function(DFPresults,EnrichmentResults){
  allData <- matrix(0,ncol=3,nrow=length(DFPresults))
  colnames(allData) <- c("numberOfPathway","numberOfGenes","piVal")
  for(i in 1:length(DFPresults)){
    allData[i,"numberOfPathway"] <- EnrichmentResults[[i]]$NKEGG + EnrichmentResults[[i]]$NReactome
    # if(is.null(class(DFPresults[i]))){
    #   allData[i,"numberOfGenes"] <- 0
    # }else{
      allData[i,"numberOfGenes"] <- nrow(as.matrix(DFPresults[[i]]$discriminantGenes))
    # }
    allData[i,"piVal"] <- DFPresults[[i]]$piVal
  }

  allData <- as.data.frame(allData)
  p <- plot_ly(x=allData$piVal,y=allData$numberOfPathway,z=allData$numberOfGenes,color=allData$numberOfGenes) %>% add_markers() %>% layout(scene = list(xaxis = list(title = 'piVal'),
                        yaxis = list(title = 'Number of pathways'),
                        zaxis = list(title = 'Number of genes')))
  return(list(dataForPlot = allData,plot = p))
}
