#! /usr/bin/env Rscript

# This script takes in all of the 10X experiments called cells and normalizes them all together
# Hopefully this will help to mitigate against any major batch effects between the different experiments.
# Due to the size of the data I have run the cell calling and filtering on each sample separately.

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(BiocParallel)
library(biomaRt)
source("/nfs/research1/marioni/mdmorgan/Thymus_ZsG/src/single_cell_functions.R")

message("Running normalisation")
# read in each counts matrix and dimension names
ncores <- 4
mcparam <- MulticoreParam(workers=ncores)
register(mcparam)

base_path <- "/nfs/research1/marioni/mdmorgan/Thymus_ZsG/"
sample_list <- list.files(base_path, pattern="^Ageing")
sample_list <- sample_list[grepl(sample_list, pattern="HTO$")]

# read in colnames, rownames and counts matrix
count_locs <- paste0(base_path, list.files(base_path, pattern="_counts-HTO.tsv"))
cellid_locs <- paste0(base_path, list.files(base_path, pattern="_colnames-HTO.tsv"))
genes_locs <- paste0(base_path, list.files(base_path, pattern="_rownames-HTO.tsv"))

#print(count_locs)
message("Reading in sparse count matrices")
intest_matrices <- bplapply(count_locs, readMM)
intest_bcs <- bplapply(cellid_locs, function(X) read.table(X, sep="\t", header=TRUE, stringsAsFactors=FALSE)[, 1])
# extract the genes and features

message("Reading in features")
intest_features <- bplapply(genes_locs, function(G) read.table(G, sep="\t", header=TRUE, stringsAsFactors=FALSE)[, 1])

# apply gene names to matrix rows and subset the common genes
common.genes <- intersect(intersect(intersect(intersect(intersect(intest_features[[1]], intest_features[[2]]), intest_features[[3]]), intest_features[[4]]), intest_features[[5]]), intest_features[[6]])
print(length(common.genes))

message("Combining together counts matrices")
# smoosh all matrices together
intest_counts <- do.call(cbind, intest_matrices)
intest_ids <- do.call(c, intest_bcs)
print(head(intest_ids))
print(head(common.genes))
colnames(intest_counts) <- intest_ids
rownames(intest_counts) <- common.genes
print(dim(intest_counts))

#writeMM(intest_counts,
#	file=paste0(base_path, "ALL_counts.tsv"))
#write.table(intest_ids,
#	file=paste0(base_path, "ALL_counts_colnames.tsv"),
#	sep="\t", row.names=FALSE, quote=FALSE)
#write.table(common.genes,
#	file=paste0(base_path, "ALL_counts_rownames.tsv"),
#	sep="\t", row.name=FALSE, quote=FALSE)

######################
## Single-cell normalisation
######################
#set.seed(42)
# remove objects from environment
rm(list=c("intest_cells", "intest_list", "mat", "intest_bcs"))
gc()

message("Deconvolution size factor normalisation")
all.sf.norm <- size_factor_normalize(intest_counts, cell.sparse=0.98,
                                     gene.sparse=0.99, cluster.size=100)
colnames(all.sf.norm) <- intest_ids
print(dim(all.sf.norm))

write.table(all.sf.norm,
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_SFnorm-HTO.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

#rm(list=c("all.sf.norm"))
#gc()

message("Estimating deconvolution size factors")
# also calculate size factors separately for later use
all.size_factors <- get_size_factors(intest_counts, cell.sparse=0.98,
                                     gene.sparse=0.99, cluster.size=100, barcodes=intest_ids)

write.table(all.size_factors,
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_sizeFactors-HTO.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
