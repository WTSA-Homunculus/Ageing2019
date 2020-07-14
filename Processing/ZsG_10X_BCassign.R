#! /usr/bin/env Rscript
# This script takes in the HTO counts matrices and determines which is the dominant HTO barcode for each single cell

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(BiocParallel)
library(fitdistrplus)
source("/nfs/research1/marioni/mdmorgan/Thymus_ZsG/src/single_cell_functions.R")

ncores <- 4
mcparam <- MulticoreParam(workers=ncores)
register(mcparam)

data.dir <- "/nfs/research1/marioni/mdmorgan/Thymus_ZsG/"
hto_files <- paste0(data.dir, list.files(data.dir, pattern="HTO_counts.tsv"))
id_files <- paste0(data.dir, list.files(data.dir, pattern="_cellids.tsv"))
tag_files <- paste0(data.dir, list.files(data.dir, pattern="_features.tsv"))

#print(hto_files)
#print(id_files)
#print(tag_files)

norm_matrix <- function(X){
  as(apply(X, 2, function(P) P/sum(P)), "dgCMatrix")
}

hto_matrices <- bplapply(hto_files, readMM)
hto_ids <-  bplapply(id_files, function(X) read.table(X, header=TRUE, stringsAsFactors=FALSE)[, 1])
hto_tags <- bplapply(tag_files, function(X) read.table(X, header=TRUE, stringsAsFactors=FALSE)[, 1])

# add cell IDs to the matrix
intest_list <- list()
for(i in seq_along(1:length(hto_ids))){
  mat <- hto_matrices[[i]]
  colnames(mat) <- hto_ids[[i]]
  rownames(mat) <- hto_tags[[i]]
  print(dim(mat))
  intest_list[[i]] <- mat
}

# CPM normalise the HTO counts
norm_matrices <- bplapply(intest_list, FUN=function(X) norm_matrix(X))

# define the expected number of clusters
exp_clusters <- function(C){
  n <- length(C)
  combos <- nrow(expand.grid(c(1:n), c(1:n)))
  n + ((combos - n)/2)
}

exp.clusters <- bplapply(hto_tags, FUN=function(TC) exp_clusters(TC))
#lapply(hto_ids, function(X) print(head(X)))
#lapply(hto_tags, function(P) print(head(P)))
#print(exp.clusters)


# do a k-means cluster on the cells in feature space. Need to define a threshold for negative cells.
cell_class <- list()
#print(hto_tags)
for(i in seq_along(1:length(hto_tags))){
  nmat <- norm_matrices[[i]]
  exp.k <- length(hto_tags[[i]])
  k.means <- kmeans(x=t(as.matrix(nmat)), centers=exp.k)
  # fit a negative binomial model to the remaining unnormalised HTO counts for each HTO barcode
  quantile_class <- list()
  for(x in seq_along(1:length(hto_tags[[i]]))){
    x.tag <- hto_tags[[i]][x]
    xtag.counts <- intest_list[[i]][x, ]
    print(x.tag)
    print(dim(xtag.counts))
    # fit a background negative binomial, excluding the cells for a given cluster with the highest
    # average counts for that HTO tag
    print(x.tag)
    print(k.means$centers)
    max.center <- max(k.means$centers[, x.tag])
    max.clust <- which(k.means$centers[, x.tag] == max.center)
    hto.cells <- k.means$cluster == max.clust
    
    # exclude top 0.5% of cells as there are +ve cells that are mis-classified by the crude k-means
    not.cluster.counts <- xtag.counts[!hto.cells]
    top.05.counts <- not.cluster.counts[order(not.cluster.counts, decreasing=TRUE)]
    top.05 <- floor(length(top.05.counts) * 0.005)
    
    x.negbin <- fitdist(top.05.counts[top.05:length(top.05.counts)], "nbinom")
    x.99.q <- as.numeric(quantile(x.negbin, probs=0.99)$quantiles)
    # classify each cells that is >= the 99th quantile as "positive"
    pos.cells <- as.numeric(xtag.counts >= x.99.q)
    names(pos.cells) <- names(xtag.counts)
    quantile_class[[x.tag]] <- pos.cells
  }
  # classify each cell into the relevant HTO
  # if more than one then call it a multiplet
  hto.class <- do.call(cbind.data.frame,
                       quantile_class)
  hto.sums <- rowSums(hto.class)
  singlets <- names(hto.sums)[hto.sums == 1]
  multiplets <- names(hto.sums)[hto.sums > 1]
  dropout <- names(hto.sums)[hto.sums == 0]
  single.class <- apply(hto.class[singlets, ], 1, FUN=function(C) colnames(hto.class)[C == 1])
  multi.class <- apply(hto.class[multiplets, ], 1, FUN=function(Q) paste(colnames(hto.class)[Q == 1], collapse="."))
  
  cell_class[[i]] <- data.frame("ID"=c(singlets, multiplets, dropout), "Class"=c(rep("Singlet", length(singlets)),
                                                                                 rep("Multiplet", length(multiplets)), rep("Dropout", length(dropout))),
                                "HTO"=c(single.class[singlets], multi.class[multiplets], rep("None", length(dropout))))
}

all.cell.class <- do.call(rbind.data.frame,
                          cell_class)

write.table(all.cell.class,
            file=paste0("/nfs/research1/marioni/mdmorgan/Thymus_ZsG/", "HTO_class.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)





