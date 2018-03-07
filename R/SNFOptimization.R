require(SNFtool)
require(survival)
require(survminer)
require(ggplot2)
require(GGally)

source("R/SNFProcedure.R")

SNFOptimization <- function(expData,patientsName,survivalData,iteration,kappa,alpha,maxNumberOfCluster){

  l <-  list()

  for(k in kappa){
    for(a in alpha){
      numClusters <- 2:maxNumberOfCluster
      snf <- SNFProcedure(expData,k, iteration,a,numClusters)
      for (nc in unique(unlist(snf$nc_tip))){

        cls <- spectralClustering(snf$W,nc)
        col <- getColorsForGroups(cls,colors=rainbow(nc))
        ClustGroup <- factor(paste("Cluster", cls, sep="."),levels=unique(paste("Cluster",cls,sep=".")))
        survivalData$ClustGroup=ClustGroup

        attach(survivalData)
        fit <- survfit(Surv(Survival, Death) ~ ClustGroup,data = survivalData)
        suv <- survminer::ggsurvplot(fit, risk.table = TRUE, risk.table.height = 0.5,
                                        xlim = c(0,5000), break.time.by = 500, pval = TRUE)
        pVal <- survminer:::surv_pvalue(fit,method = "survdiff",data=survivalData)
        #detach(survivalData)

        obj <- list()
        obj$k <- k
        obj$a <- a
        obj$nc <- nc
        obj$cls <- cls
        obj$table_cls <- table(cls)
        obj$survivalDataPlot <- suv
        obj$pVal <- pVal$pval
        l[[length(l)+1]] <- obj
      }
    }
  }

  pd <- list()
  for(i in l){
    if(! as.character(length(i$table_cls)) %in% names(pd)){
      pd[[as.character(length(i$table_cls))]] = i
    }
    if(pd[[as.character(length(i$table_cls))]]$pVal > i$pVal){
      pd[[as.character(length(i$table_cls))]] = i

    }
  }

  for(i in 1:length(pd)){
    optClust <- pd[[i]]
    labelsSNF <- as.data.frame(optClust$cls)
    rownames(labelsSNF) <- patientsName
    colnames(labelsSNF) <- "Classes"
    optClust$labelsSNF <- labelsSNF
    pd[[i]] <- optClust
  }

  res <- list("OptimalCluster" = pd,"AllCluster" = l)
  return(res)
}


