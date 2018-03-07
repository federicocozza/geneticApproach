getPathways <- function(bestEnrich,DFPgenesEntrez,RNA,minimumSizeOfPathway = 5){
  pathwaysKegg <- bestEnrich$KEGGPathway[,c("ID","geneID")]
  pathwaysReactome <- bestEnrich$ReactomePathway[,c("ID","geneID")]

  namesKegg = row.names(pathwaysKegg)
  namesReactome = row.names(pathwaysReactome)

  allPathway <- list()

  for(i in 1:nrow(pathwaysKegg)){
    cur <- pathwaysKegg[i,]
    genes <- strsplit(cur$geneID,"/")[[1]]
    if(length(genes) < minimumSizeOfPathway){
      next()
    }

    # geneSymbols <- DFPgenesEntrez[which((DFPgenesEntrez[,2] %in% genes)),1]
    #
    # if (length(genes) != length(geneSymbols)) {
    #   warning("Error: there are some Entrez IDs without a corresponding gene symbol in the original dataset.")
    # }
    allPathway[paste("KEGG_",cur$ID,sep="")] <- list(RNA[as.vector(genes),])
  }

  for(i in 1:nrow(pathwaysReactome)){
    cur <- pathwaysReactome[i,]
    genes <- strsplit(cur$geneID,"/")[[1]]
    if(length(genes) < minimumSizeOfPathway){
      next()
    }

    # geneSymbols <- DFPgenesEntrez[which((DFPgenesEntrez[,2] %in% genes)),1]
    #
    # if (length(genes) != length(geneSymbols)) {
    #   warning("Error: there are some Entrez IDs without a corresponding gene symbol in the original dataset.")
    # }
    allPathway[paste("REACTOME_",cur$ID,sep="")] <- list(RNA[as.vector(genes),])
  }

  return(allPathway)
}
