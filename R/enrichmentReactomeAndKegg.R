require(parallel)
require(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(SummarizedExperiment)


enrichmentReactomeAndKegg <- function(DFPresults){
  print(paste("Numbrer of enrichment to analyzed:",length(DFPresults)))
  results <- list()
  for(i in 1:length(DFPresults)){
    cat(paste("Enrichment",i,"..."))
    if(nrow(DFPresults[[i]]$discriminantGenes) == 0 || class(DFPresults[[i]]$discriminantGenes) == "character"){
      results[[length(results)+1]] <- list(KEGGPathway = list(),ReactomePathway = list(),NKEGG = 0,NReactome = 0)
      next
    }

    DFPgenes <- as.factor(rownames(DFPresults[[i]]$discriminantGenes))
    # rownames(genesSymbolsToEntrezID) <- genesSymbolsToEntrezID[,"geneSymbol"]
    # DFPgenesEntrez <- unique(genesSymbolsToEntrezID[DFPgenes,"entrezID"])
    # DFPgenesEntrez <- DFPgenesEntrez[!is.na(DFPgenesEntrez)]

    KEGGResult <- enrichKEGG(DFPgenes)
    ReactomeResult <- ReactomePA::enrichPathway(gene = DFPgenes, pvalueCutoff = 0.05)

    results[[length(results)+1]] <- list(KEGGPathway = KEGGResult,ReactomePathway = ReactomeResult,NKEGG = ifelse(is.null(KEGGResult),0,nrow(KEGGResult)),NReactome = ifelse(is.null(ReactomeResult),0,nrow(ReactomeResult)))

    cat(" DONE!!!\n")
  }

  return(results)
}


enrichmentReactomeAndKeggByList <- function(genes){
  DFPgenes <- genes
    
  KEGGResult <- enrichKEGG(DFPgenes,pvalueCutoff = 0.05,pAdjustMethod = "BH")
  ReactomeResult <- ReactomePA::enrichPathway(gene = DFPgenes, pvalueCutoff = 0.05,pAdjustMethod = "BH")
    
  results <- list(KEGGPathway = KEGGResult,ReactomePathway = ReactomeResult,NKEGG = ifelse(is.null(KEGGResult),0,nrow(KEGGResult)),NReactome = ifelse(is.null(ReactomeResult),0,nrow(ReactomeResult)))
  return(results)
}
