#! /usr/bin/env Rscript

# This script takes in all of the intestine 10X experiments to data and normalizes them all together
# Hopefully this will help to mitigate against any major batch effects between the different experiments.

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(BiocParallel)
library(biomaRt)
source("/nfs/research1/marioni/mdmorgan/Thymus_ZsG/src/single_cell_functions.R")

# rather than relying on the CellRanger to call cells I'll use emptyDrops and remove poor quality cells at a later point.
ncores <- 12
mcparam <- MulticoreParam(workers=ncores)
register(mcparam)

base_path <- "/nfs/research1/marioni/mdmorgan/Thymus_ZsG/"
sample_list <- list.files(base_path, pattern="^Ageing")
sample_list <- sample_list[grepl(sample_list, pattern="HTO")]

# locations of hdf5 molecule information across all samples
mol_locs <- paste0(base_path, sample_list, "/outs/molecule_info.h5")
out_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrices/matrix.mtx.gz")
bc_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrices/barcodes.tsv.gz")
genes_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrices/features.tsv.gz")

intest_matrices <- bplapply(out_locs, readMM)
intest_bcs <- bplapply(bc_locs, function(X) read.table(X, header=FALSE, stringsAsFactors=FALSE)[, 1])
intest_genes <- bplapply(genes_locs, function(G) read.table(G, header=FALSE, stringsAsFactors=FALSE)[, 1])

# add sample ID to barcodes
intest_barcode <- list()
for(i in seq_along(sample_list)){
  samp <- sample_list[i]
  intest_barcode[[i]] <- paste0(samp, "_", intest_bcs[[i]])
}

# call cells with emptyDrops
intest_calls <- lapply(intest_matrices, emptyDrops, niters=20000, ignore=4999, BPPARAM=mcparam, lower=100, retain=Inf)

# identify cells using a FDR 1%
sig_cells <- lapply(intest_calls, function(X) X$FDR <= 0.01 & !is.na(X$FDR))

# subset the called cells - the numbers are very similar to what CellRanger gives
intest_cells <- lapply(1:length(intest_barcode), function(i) intest_matrices[[i]][, sig_cells[[i]]])
cell_barcodes <- lapply(1:length(intest_barcode), function(x) intest_barcode[[x]][sig_cells[[x]]])

# apply gene names to matrix rows and subset the common genes
common.genes <- intersect(intersect(intersect(intersect(intest_genes[[1]],
                                                                                                                                                                                                                                                                                                                            intest_genes[[2]]),
                                                                                                                                                                                                                                                                                                                  intest_genes[[3]]),
                                                                                                                                                                                                                                                                                                        intest_genes[[4]]),
                                                                                                                                                                                                                                                                                              intest_genes[[5]]),
                                                                                                                                                                                                                                                                                    intest_genes[[6]]),
                                                                                                                                                                                                                                                                          intest_genes[[7]]),
                                                                                                                                                                                                                                                                intest_genes[[8]]),
                                                                                                                                                                                                                                                      intest_genes[[9]]),
                                                                                                                                                                                                                                            intest_genes[[10]]),
                                                                                                                                                                                                                                  intest_genes[[11]]),
                                                                                                                                                                                                                        intest_genes[[12]]),
                                                                                                                                                                                                              intest_genes[[13]]),
                                                                                                                                                                                                    intest_genes[[14]]),
                                                                                                                                                                                          intest_genes[[15]]),
                                                                                                                                                                                intest_genes[[16]]),
                                                                                                                                                                      intest_genes[[17]]),
                                                                                                                                                            intest_genes[[18]]),
                                                                                                                                                  intest_genes[[19]]),
                                                                                                                                        intest_genes[[20]]),
                                                                                                                              intest_genes[[21]]),
                                                                                                                    intest_genes[[22]]),
                                                                                                          intest_genes[[23]]),
                                                                                                intest_genes[[24]]),
                                                                                      intest_genes[[25]]),
                                                                            intest_genes[[26]]),
                                                                  intest_genes[[27]]),
                                                        intest_genes[[28]]),
                                              intest_genes[[29]]),
                                    intest_genes[[30]]),
                          intest_genes[[31]])


intest_list <- list()

for(i in seq_along(1:length(cell_barcodes))){
  mat <- intest_cells[[i]]
  rownames(mat) <- intest_genes[[i]]
  intest_list[[i]] <- mat[common.genes, ]
}

# smoosh all matrices together
intest_counts <- do.call(cbind, intest_list)
intest_ids <- do.call(c, cell_barcodes)

# remove supefluous matrices
rm(list=c("intest_matrices", "intest_cells"))
gc()

######################
## Single-cell QC
######################
## Remove very low complexity cell libraries, i.e. ones that express < 1000 genes
lib.sizes <- Matrix::colSums(intest_counts)
nonzero.genes <- Matrix::colSums(intest_counts > 0)

# remove cells with < 1000 genes
#sum(nonzero.genes > 1000)

intest_counts <- intest_counts[, nonzero.genes > 1000]
intest_ids <- intest_ids[nonzero.genes > 1000]

## remove cells with excessively high mitochondrial content
biomaRt.connection <- useMart("ensembl", "mmusculus_gene_ensembl")
gene.df <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters="ensembl_gene_id", 
                 values=rownames(intest_counts), mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

mt.counts <- intest_counts[which(common.genes %in% gene.df$ensembl_gene_id[gene.df$chromosome_name == "MT"]), ]
mt.fraction <- Matrix::colSums(mt.counts)/Matrix::colSums(intest_counts)

# fit a median-centred, MAD variance model normal to dervive p-values for outlier proportions
# this isn't the best calibrated model, so might be a little liberal
mt.p <- pnorm(mt.fraction, mean=median(mt.fraction), sd=mad(mt.fraction)*2, lower.tail=FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method="fdr") < 0.05)])

# cells with a high mitochondrial fraction > `mt.lim` are removed as outliers.
intest_counts <- intest_counts[, mt.fraction < mt.lim]
intest_ids <- intest_ids[mt.fraction < mt.lim]

# Now we have fairly decent cells I'll proceed to the normalisation using deconvolution-estimated size factors
writeMM(intest_counts,
        file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_counts.tsv")
write.table(intest_ids,
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_counts_colnames.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(rownames(intest_counts),
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_counts_rownames.tsv",
            sep="\t", row.names=FALSE, quote=FALSE)

set.seed(42)
# remove objects from environment
# rm(list=c("intest_cells", "intest_list", "mat", "intest_bcs"))
# gc()

all.sf.norm <- size_factor_normalize(intest_counts, cell.sparse=0.98,
                                     gene.sparse=0.99, cluster.size=100)
colnames(all.sf.norm) <- intest_ids

write.table(all.sf.norm,
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_SFnorm.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

# also calculate size factors separately for later use
all.size_factors <- get_size_factors(intest_counts, cell.sparse=0.98,
                                     gene.sparse=0.99, cluster.size=100, barcodes=intest_ids)

write.table(all.size_factors,
            file="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/ALL_sizeFactors.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
