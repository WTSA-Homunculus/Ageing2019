## generic functions to source into environment

map2go <- function(dataframe){
  # goseq expects a list object each time with names as gene IDs
  # assumes column names are ensembl_gene_id and GeneSet
  set.goseq <- list()
  gene_ids <- unique(dataframe$ensembl_gene_id)
  for(x in seq_along(gene_ids)){
    gid <- gene_ids[x]
    gset <- dataframe$GeneSet[dataframe$ensembl_gene_id == gid]
    set.goseq[[gid]] <- gset
  }
  return(set.goseq)
}

####################################
## Equality/specificity functions ##
####################################

Gini <- function(x){
  # absolute differences between each pair of points
  diff.vec <- numeric(0)
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      if(i != j){
        diff.vec <- c(diff.vec, abs(x[i] - x[j]))
      }
    }
  }
  sum.x <- sum(x)
  n <- length(x)
  gini.c <- sum(diff.vec)/((2*n) * sum.x)
  return(gini.c)
}


Tau <- function(x){
  # calculate specificity index over vector/expression profile
  max_x <- max(x)
  max_I <- 1 - (x/max_x)
  numerator <- sum(max_I)
  denominator <- length(x)
  
  tau <- numerator/denominator
  
  return(tau)
}

#############################
## Gene modules/eigengenes ##
#############################

find_eigengene <- function(dataframe, side="left"){
  # calculate the eigenvector through the set of genes
  # assume windows are columns and genes are rows
  # assume the first eigenvector as the eigengene
  spec <- svd(dataframe)
  if(side == "left"){
    eigengene <- abs(spec$v[, 1])
  } else if(side == "right"){
    eigengene <- abs(spec$u[, 1])
  }
  return(eigengene)
}


category_igb <- function(go.map, categories){
  # calculate the similarity between gene categories
  cat.igb <- list()
  for(i in 1:length(categories)){
    i.igb <- list()
    for(j in 1:length(categories)){
      num.igb <- length(intersect(go.map[[categories[i]]],
                                  go.map[[categories[j]]]))
      denom.igb <- length(union(go.map[[categories[i]]],
                                go.map[[categories[j]]]))
      igb <- num.igb/denom.igb
      i.igb[[categories[j]]] <- igb
    }
    cat.igb[[categories[i]]] <- i.igb
  }
  
  igb.df <- data.frame(lapply(cat.igb, unlist))
  return(igb.df)  
}


category_eigengene <- function(cat_groups, cat_clusters, go.map, dataframe, cluster_name,
                               pseudotime){
  # calculate the eigengene for each category, within a gene module
  
  cat2eigen <- list()
  genes.list <- list()
  for(i in seq_along(cat_groups)){
    clust <- cat_groups[i]
    cats <- names(cat_clusters)[cat_clusters == clust]
    cat.genes <- unique(unlist(go.map[cats]))
    cat.fit <- dataframe[cat.genes, ]
    cat.eigen <- find_eigengene(cat.fit)
    cat2eigen[[cat_groups[i]]] <- cat.eigen
    genes.list[[cat_groups[i]]] <- cat.genes
  }
  cat.eigen.df <- data.frame(do.call(cbind, cat2eigen))
  colnames(cat.eigen.df) <- cat_groups
  cat.eigen.df$DPT <- pseudotime
  
  cat.eigen.melt <- melt(cat.eigen.df, id.vars="DPT")
  colnames(cat.eigen.melt) <- c("DPT", "Cat", "Eigengene")
  cat.eigen.melt$Cat <- paste(cluster_name, cat.eigen.melt$Cat, sep=".")
  
  eigen.obj <- list("melt"=cat.eigen.melt,
                    "genes"=genes.list)
  return(eigen.obj)
}


functional_eigengenes <- function(clusters, dataframe, go.sig.all, pseudotime, gene_map){
  # retrieve the eigengenes for enriched functional terms within each cluster
  # keep a record of which categories are collapsed together within each cluster
  # return a list object containing:
  # melted eigengene dataframe for that cluster
  # enriched GO terms collapsed into each eigengene
  # genes that contribute to each eigengene
  
  storage_obj <- list()
  
  cluster_names <- names(clusters)
  for(c in seq_along(cluster_names)){
    clust.name <- cluster_names[c]
    clust.genes <- clusters[[clust.name]]
    # print(clust.name)
    # get the enriched terms first
    clust.sig.all <- go.sig.all[go.sig.all$Cluster == clust.name, ]
    
    # need to handle clusters with no or a single enrichment.
    # just calculate the eigengene for entire cluster?
    
    if(dim(clust.sig.all)[1] > 1){
      clust.all <- clust.sig.all[order(clust.sig.all$padjust, decreasing=FALSE), ]
      top.cats <- clust.all
      top.cats <- top.cats[order(top.cats$foldEnrich, decreasing=TRUE), ]
      
      # for each category, pull out the respective genes
      clust.go <- gene_map[clust.genes]
      #clust.go <- getgo(genes=clust.genes, genome="mm10", id="ensGene", fetch.cats="GO:BP")
      
      # remap categories to genes, from genes to categories 
      # using enriched categories in the cluster
      go2gene <- list()
      for(g in 1:length(clust.genes)){
        gene <- clust.genes[g]
        match.cats <- intersect(clust.go[[gene]], clust.all$category)
        if(length(match.cats) > 0){
          for(q in 1:length(match.cats)){
            go2gene[[match.cats[q]]] <- c(go2gene[[match.cats[q]]], gene)
          }
        }
        cat.hits <- names(go2gene)
      }
      
      # calculate the proportion of shared genes between enriched categories
      igb.df <- category_igb(go2gene, cat.hits)
      
      # cluster categories together, select approx. independent groups with IGB <= 0.4
      cats.clust <- flashClust(as.dist(igb.df), method="complete")
      collapse.cats <- cutree(cats.clust, h=0.6)
      cat.groups <- unique(collapse.cats)
      
      eigen.obj <- category_eigengene(cat.groups, collapse.cats, go2gene, dataframe, clust.name, pseudotime)
      cat.eigen.melt <- eigen.obj$melt
      eigen.genes <- eigen.obj$genes
      
      storage_obj[[clust.name]] <- list("eigen"=cat.eigen.melt,
                                        "eigen.load"=eigen.genes,
                                        "cat_clust"=collapse.cats,
                                        "cat_genes"=go2gene)
    }
    else{
      cat.eigen <- find_eigengene(dataframe[rownames(dataframe) %in% clust.genes, ])
      
      eigen.df <- data.frame(do.call(cbind, list(clust.name=cat.eigen)))
      eigen.df$DPT <- pseudotime
      eigen.melt <- melt(eigen.df, id.vars="DPT")
      colnames(eigen.melt) <- c("DPT", "Cat", "Eigengene")
      eigen.melt$Cat <- paste(c(clust.name, "1"), collapse=".")
      
      storage_obj[[clust.name]] <- list("eigen"=eigen.melt,
                                        "eigen.load"=clust.genes,
                                        "cat_clust"=paste(c(clust.name, "1"), collapse="."),
                                        "cat_genes"=list("None"=clust.genes))
    }
    
  }
  return(storage_obj)
}
