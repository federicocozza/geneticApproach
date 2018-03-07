library(igraph)
library(NetPathMiner)

createInteractionGraph <- function(SVMWeights,significantlyPathways,dataSet){
  all_genes_together <- SVMWeights[significantlyPathways$namePathway]

  for(i in 1:length(all_genes_together)){
    names(all_genes_together[[i]]) <- gsub("[.]","-",names(all_genes_together[[i]]))
  }

  plotlist <- list()
  for(comb in as.character(unique(significantlyPathways$classes))){
    all_genes <- all_genes_together[which(significantlyPathways$classes %in% comb)]
    name_pathways <- names(all_genes)
    names(all_genes) <- NULL

    unique_genes <- unique(names(unlist(all_genes)))

    genes_expression <- dataSet[unique_genes,]

    adj <- cor(t(genes_expression))
    adjPos = adj
    adjNeg = adj
    adjPos[adjPos<0]=0
    adjNeg[adjNeg>0]=0
    adjPos <- 1-adjPos
    adjNeg <- 1-abs(adjNeg)

    gp <- graph.adjacency(adjPos,weighted = T,mode = "undirected")
    gp_mst <- igraph::minimum.spanning.tree(gp)

    gn <- graph.adjacency(adjNeg,weighted = T,mode = "undirected")
    gn_mst <- igraph::minimum.spanning.tree(gn)

    gp_mst_adj = as.matrix(get.adjacency(graph = gp_mst,attr = "weight"))
    gn_mst_adj = -1 * as.matrix(get.adjacency(graph = gn_mst,attr = "weight"))

    ADJ = gp_mst_adj + gn_mst_adj

    g <- graph.adjacency(ADJ, mode = "undirected",weighted = TRUE)
    edge_col = ifelse(E(g)$weight<0,"red","green")

    groups = lapply(X = all_genes,FUN = names)

    for(i in 1:length(all_genes)){
      all_genes[[i]] = all_genes[[i]] * length(all_genes[[i]])
    }

    rcolors = rainbow(length(all_genes))
    node_colors = rep(rcolors[1],nrow(ADJ))
    names(node_colors) = colnames(ADJ)

    for(gene in rownames(ADJ)){
      gene_vector <- c()
      for (i in 1:length(all_genes)){
        gene_vector <- c(gene_vector,all_genes[[i]][gene])
      }

      if(sum(is.na(gene_vector)) != length(gene_vector)){
        node_colors[gene] <- rcolors[which.max(gene_vector)]
      }else{
        print(gene)
        node_colors[gene] <- "black"
      }
    }

    g <- setAttribute(g, "pathway", node_colors)
    V(g)$name= 1:vcount(g)
    l = layoutVertexByAttr(g, "pathway",layout=layout.kamada.kawai)

    #plotNetwork(g, vertex.color="pathway",layout=l)
    rownames(l) = rownames(ADJ)
    V(g)$name = rownames(ADJ)



    genes_weights = unlist(all_genes)

    vsize = c()
    for(index in rownames(ADJ)){
      vsize = c(vsize,max(genes_weights[index]))
    }


    par(mar=c(0,0,0,0))
    p <- plot(g, edge.width = abs(E(g)$weight) * 2, vertex.size = vsize * 10, edge.color = edge_col,mark.groups=groups,
         layout = l,vertex.label.cex = 0.6,vertex.color=node_colors,mark.col = rainbow(length(all_genes),alpha = 0.3))
    legend(x= "topleft",fill=rainbow(length(all_genes)),legend = name_pathways,title = paste("Classes", comb))
    plotlist[[length(plotlist)+1]] <- p
  }

  return(plotlist)
}
