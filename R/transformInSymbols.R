transformInSymbols <- function(genesId){
  library(org.Hs.eg.db)
  library(annotate)
  
  geneSymb <- getSYMBOL(genesId,data="org.Hs.eg")
  
  return(geneSymb)
}
