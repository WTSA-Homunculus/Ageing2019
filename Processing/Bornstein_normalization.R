## process and normalized MARS-seq data from Bornstein et al Nature (2018)
## Tuft-like mTECs

library(SingleCellExperiment)
library(scran)
library(igraph)
library(scater)
library(biomaRt)
source("src/SingleCellFunctions/single_cell_functions.R")

meta.df <- read.table("GSE103967_metadata.txt.gz",
                      sep="\t", h=TRUE, stringsAsFactors=FALSE)
bad.cells <- meta.df$well[meta.df$readcount == 0]

in.path <- "Bornstein/"
in.files <- list.files(in.path, pattern="GSM")

mars.list <- list()

for(i in seq_along(in.files)){
  fle <- paste0(in.path, in.files[i])
  in.df <- read.table(fle, sep="\t",
                      h=TRUE, stringsAsFactors=FALSE)
  in.df$gene_id <- rownames(in.df)
  
  # remove the bad cells
  in.df <- in.df[, !colnames(in.df) %in% bad.cells]
  mars.list[[in.files[i]]] <- in.df
}

all.counts <- Reduce(x=mars.list,
                     f=function(x, y) merge(x, y, by='gene_id'))
rownames(all.counts) <- all.counts$gene_id
all.counts <- all.counts[, -1]

#############
#### QC #####
#############
# do some basic QC to remove poor quality cells, i.e. < 1,000 UMIs
low.umi <- log10(colSums(all.counts)) < 3.0
counts.nz <- all.counts[, !low.umi]

cell_sparsity <- apply(counts.nz == 0, 2, sum)/nrow(counts.nz)
sparse.cells <- cell_sparsity < 0.98
counts.nz <- counts.nz[, sparse.cells]

gene_sparsity <- apply(counts.nz == 0, 1, sum)/ncol(counts.nz)
sparse.genes <- gene_sparsity > 0.99
counts.nz <- counts.nz[!sparse.genes, ]

### Get the gene information
# setup the biomaRt connection to Ensembl using the correct species genome 
ensembl <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl')

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',
                                  'chromosome_name', 'start_position',
                                  'end_position', 'strand'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=rownames(all.counts))

counts.nz <- as(as.matrix(counts.nz), "dgCMatrix")
mars.sce <- SingleCellExperiment(list(counts=counts.nz))

clusters <- quickCluster(mars.sce, min.size=150,
                         method="igraph")
max.size <- floor(150/2)

mars.sce <- computeSumFactors(mars.sce,
                         max.cluster.size=max.size,
                         positive=FALSE,
                         get.spikes=FALSE,
                         assay.type='counts', clusters=clusters)

mars.sce <- normalise(mars.sce)
sf.norm <- as.data.frame(as(exprs(mars.sce), "matrix"))
sf.norm$gene_id <- rownames(counts.nz)

write.table(sf.norm,
            file="Bornstein/SF_norm.tsv",
            row.names=FALSE, quote=FALSE, sep="\t")







