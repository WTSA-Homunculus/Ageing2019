#! /usr/bin/env Rscript

library(lettercase)
library(reshape2)
library(vegan)
library(biomaRt)

# read in the normalised isoform counts
thymus.norm <- read.table("/nfs/research1/marioni/mdmorgan/Thymus/Quant.dir/ThymusQC_TranscriptNorm.tsv.gz",
                          sep="\t", header=TRUE, stringsAsFactors=FALSE)
rownames(thymus.norm) <- thymus.norm$transcript_id
thymus.norm <- thymus.norm[, -ncol(thymus.norm)]

# retrieve the gene info for all isoforms
all.transcripts <- rownames(thymus.norm)
ensembl <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id',  'external_gene_name', 'chromosome_name'),
                     filters='ensembl_transcript_id', mart=ensembl,
                     values=all.transcripts)
rownames(gene_symbol) <- gene_symbol$ensembl_transcript_id

# retieve the meta data to pull out the different clusters
thymus.meta <- read.table("/nfs/research1/marioni/mdmorgan/Thymus/Quant.dir/ThymusMarker_tSNE_PCA_meta.tsv",
                          sep="\t", h=TRUE, stringsAsFactors=FALSE)

# create a variable that represents the ages as a single continuous variable
thymus.meta$AgeInt <- as.numeric(gsub(thymus.meta$Age, pattern="wk", replacement=""))
thymus.meta$AgeFactor <- ordered(thymus.meta$AgeInt,
                                 levels=c("1", "4", "16", "32", "52"))

thymus.meta$SampleCol <- unlist(lapply(strsplit(thymus.meta$Sample, split=".", fixed=TRUE),
                                       FUN=function(X) paste(X[1:4], collapse=".")))

# only test the valid genes, genes that have more than 1 isoform
valid.genes <- unique(names(table(gene_symbol$ensembl_gene_id)[table(gene_symbol$ensembl_gene_id) > 1]))
permanova.meta <- thymus.meta[thymus.meta$SampleCol %in% colnames(thymus.norm), c("Sample", "SampleCol", "AgeFactor", "TFIDF.Cluster")] 
permanova.meta$TFIDF.Cluster <- as.factor(permanova.meta$TFIDF.Cluster)
rownames(permanova.meta) <- permanova.meta$SampleCol

# subset the gene expression data to match the meta data                                                                                                                                          
thymus.norm <- thymus.norm[, intersect(colnames(thymus.norm), permanova.meta$SampleCol)]
permanova.results <- list()

# being PERMANOVA across all clusters and ages
for(i in seq_along(valid.genes)){
  gene <- valid.genes[i]
  
  gene.isoforms <- gene_symbol[gene_symbol$ensembl_gene_id %in% c(gene), ]$ensembl_transcript_id
  gene.exprs <- thymus.norm[gene.isoforms, colnames(thymus.norm) %in% thymus.meta$SampleCol]
  gene.dist <- vegdist(as.matrix(t(gene.exprs)), method="bray")
  gene.dist[is.na(gene.dist)] <- 0
  
  gene.permanova <- adonis(gene.dist ~ TFIDF.Cluster + AgeFactor, data=permanova.meta)
  
  permanova.results[[gene]] <- gene.permanova$aov.tab
}

permanova.results <- lapply(permanova.results, as.matrix)
perm.df <- do.call(rbind.data.frame,
                   permanova.results)
perm.df$Gene <- unlist(lapply(strsplit(rownames(perm.df), split=".", fixed=TRUE),
                              FUN=function(G) paste0(G[1])))
perm.df$Param <- unlist(lapply(strsplit(rownames(perm.df), split=".", fixed=TRUE),
                               FUN=function(G) paste0(G[2])))
perm.df$`Pr(>F)`[is.na(perm.df$`Pr(>F)`)] <- 1

write.table(perm.df,
            file="/nfs/research1/marioni/mdmorgan/Thymus/Quant.dir/Transcript_PERMANOVA_results.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")
