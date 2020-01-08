### Prepping E-MTAB-8560 submission
qc.meta <- read.table("~/Dropbox/AgeingExperiment/Meta_preQC.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
sub.sdrf <- readLines("~/Dropbox/AgeingExperiment/submission9455_annotare_v1.sdrf.txt")

sdrf.list <- list()
for(x in seq_along(sub.sdrf)){
  if(x == 1){
    head.line <- unlist(strsplit(sub.sdrf[x], split="\t"))
  } else{
    line.x <- unlist(strsplit(sub.sdrf[x], split="\t"))
    sdrf.list[[paste0(x)]] <- line.x
  }
}

sdrf.df <- do.call(rbind.data.frame,
                   sdrf.list)
colnames(sdrf.df) <- head.line

# get all fastq filenames
all.fastq <- read.table("~/Dropbox/AgeingExperiment/all_fastq.md5",
                        header=FALSE, stringsAsFactors=FALSE)
all.fastq$V2 <- gsub(all.fastq$V2, pattern="\\S+/fastqs/rawdata/\\S+/", replacement="")

colnames(all.fastq) <- c("Comment[MD5]", "Array Data File")

# merge and keep all and fill down all meta-data
fastq.merge <- merge(sdrf.df, all.fastq, by=c('Array Data File', "Comment[MD5]"), all=TRUE)

fastq.merge$`Source Name` <- as.character(gsub(fastq.merge$`Array Data File`, pattern="\\.fastq\\.gz", replacement=""))
fastq.merge$`Characteristics[organism]` <- "Mus musculus"
fastq.merge$`Characteristics[strain]` <- "C57BL/6"
fastq.merge$`Unit[time unit]` <- "week"
fastq.merge$`Term Source REF` <- "EFO"
fastq.merge$`Term Accession Number` <- "UO_0000034"
fastq.merge$`Characteristics[sex]` <- "female"
fastq.merge$`Characteristics[organism part]` <- "thymus"
fastq.merge$`Characteristics[genotype]` <- "wild type genotype"
fastq.merge$`Material Type` <- "cell"
fastq.merge$`Extract Name` <- gsub(fastq.merge$`Array Data File`, pattern="\\.fastq\\.gz", replacement="")
fastq.merge$`Comment[LIBRARY_LAYOUT]` <- "PAIRED"
fastq.merge$`Comment[LIBRARY_SELECTION]` <- "PolyA"
fastq.merge$`Comment[LIBRARY_SOURCE]` <- "TRANSCRIPTOMIC SINGLE CELL"
fastq.merge$`Comment[LIBRARY_STRAND]` <- "not applicable"
fastq.merge$`Comment[LIBRARY_STRATEGY]` <- "RNA-Seq"
fastq.merge$`Comment[NOMINAL_LENGTH]` <- 500
fastq.merge$`Comment[NOMINAL_SDEV]` <- 1
fastq.merge$`Comment[ORIENTATION]` <- "5'-3'-5'-3'"
fastq.merge$`Comment[end bias]` <- "none"
fastq.merge$`Comment[input molecule]` <- "polyA RNA"
fastq.merge$`Comment[library construction]` <- "Smart-seq2"
fastq.merge$`Comment[primer]` <- "oligo-dT"
fastq.merge$`Comment[single cell isolation]` <- "FACS"
fastq.merge$`Comment[spike in]` <- "ERCC mix1"
fastq.merge$`Comment[spike in dilution]` <- "1:600000"
fastq.merge$`Assay Name` <-  gsub(fastq.merge$`Array Data File`, pattern="\\.fastq\\.gz", replacement="")
fastq.merge$`Technology Type` <- "sequencing assay"
fastq.merge$`Factor Value[single cell identifier]` <- unlist(lapply(strsplit(gsub(fastq.merge$`Array Data File`, pattern="\\.fastq\\.gz", replacement=""),
                                                                             split="_", fixed=TRUE),
                                                                    FUN=function(G) paste(G[1:3], collapse="_")))
fastq.merge$`Unit[time unit]` <- "week"
fastq.merge$`Term Source REF` <- "EFO"
fastq.merge$`Term Accession Number` <- "UO_0000034"
fastq.merge$PlateID <- unlist(lapply(strsplit(unlist(lapply(strsplit(fastq.merge$`Assay Name`, split="_", fixed=TRUE),
                                                            FUN=function(X)(paste0(X[1])))),
                                              split="-", fixed=TRUE), FUN=function(P) paste0(P[2])))
fastq.merge$Position <- unlist(lapply(strsplit(unlist(lapply(strsplit(fastq.merge$`Assay Name`, split="_", fixed=TRUE),
                                                             FUN=function(X)(paste0(X[1])))),
                                               split="-", fixed=TRUE), FUN=function(P) paste0(P[1])))
fastq.merge$Row <- gsub(fastq.merge$Position, pattern="([A-Z])([0-9]+)", replacement="\\1")
fastq.merge$Column <- as.numeric(as.character(gsub(fastq.merge$Position, pattern="([A-Z])([0-9]+)", replacement="\\2")))

# # read in post-QC meta data to confirm which ones are good QC cells
# thymus.meta <- read.table("~/Dropbox/`Characteristics[age]`ingExperiment/Data/ThymusMarker_tSNE_PCA_meta.tsv",
#                             sep="\t", h=TRUE, stringsAsFactors=FALSE)


fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "16"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "4"

# spec2
# add in mouse age information
fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "32"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "52"


# spec3
# add in mouse age information
fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "52"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "32"


# spec4
# add in mouse age information
fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "16"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "32"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "1"


# spec5
# add in mouse age information
fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "52"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "16"


# spec6
# add in mouse age information
fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(1:8)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(9:16)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Column %in% c(17:24)) & (fastq.merge$Row != "A") & 
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "16"

# add sorting day information
fastq.merge$SortDay[fastq.merge$PlateID %in% paste0("B00", c(2297, 2229, 2285, 2287, 2289, 2286))] <- 1
fastq.merge$SortDay[fastq.merge$PlateID %in% paste0("B00", c(2294, 2284, 2283, 2288, 2295, 2296))] <- 2
fastq.merge$SortDay[fastq.merge$PlateID %in% paste0("B00", c(2290, 2291, 2292, 2293, 2298, 2300))] <- 3
fastq.merge$SortDay[fastq.merge$PlateID %in% paste0("B00", c(2301, 1631, 1630, 1629, 1632))] <- 4
fastq.merge$SortDay[fastq.merge$PlateID %in% paste0("B00", c(2336, 2345, 2450, 2449, 2337, 2451, 2601, 2299))] <- 5

# exclusions
# some plates aren't full, and some of the minibulks could not be generated due to insufficient material
# add in the minibulk age and cell type designation
fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "16"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2289, 2286, 2295, 2296, 2298, 
                                                          2300, 1629, 1632, 2337, 2451)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "32"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2297, 2294, 2290, 2301, 2336, 2299)))] <- "52"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "52"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2229, 2284, 2291, 2345, 2601)))] <- "32"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "16"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "32"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2287, 2288, 2293, 1630, 2449)))] <- "1"


fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "52"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(2285, 2283, 2292, 2450)))] <- "16"


fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(9, 12, 15, 18)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "4"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(10, 13, 16, 19)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "1"

fastq.merge$`Characteristics[age]`[(fastq.merge$Row == "A") & (fastq.merge$Column %in% c(11, 14, 17, 20)) &
                 (fastq.merge$PlateID %in% paste0("B00", c(1631)))] <- "16"

fastq.merge$`Characteristics[individual]` <- fastq.merge$SortDay
fastq.merge$`Factor Value[age]` <- fastq.merge$`Characteristics[age]`
fastq.merge$`Unit[time unit]` <- "week"
fastq.merge$`Term Accession Number` <- "UO_0000034"
fastq.merge$`Term Source REF` <- "EFO"



## cell type is assigned based on plate position
fastq.merge$`Characteristics[cell type]` <- "CNTRL"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("A") & fastq.merge$Column %in% c(9:11)] <- "mTEClo.MiniBulk"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("A") & fastq.merge$Column %in% c(12:14)] <- "mTEChi.MiniBulk"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("A") & fastq.merge$Column %in% c(15:17)] <- "Dsg3+TEC.MiniBulk"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("A") & fastq.merge$Column %in% c(18:20)] <- "cTECMiniBulk"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("A") & fastq.merge$Column %in% c(21:24)] <- "EMPTY"
fastq.merge$`Characteristics[cell type]`[is.na(fastq.merge$`Characteristics[age]`)] <- "EMPTY"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("B", "C", "D") & fastq.merge$Column %in% c(1:24)] <- "mTEClo"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("E") & fastq.merge$Column %in% c(1:18)] <- "mTEClo"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("E") & fastq.merge$Column %in% c(19:24)] <- "mTEChi"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("F", "G", "H") & fastq.merge$Column %in% c(1:24)] <- "mTEChi"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("I") & fastq.merge$Column %in% c(1:12)] <- "mTEChi"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("I") & fastq.merge$Column %in% c(13:24)] <- "Dsg3+TEC"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("J", "K", "L") & fastq.merge$Column %in% c(1:24)] <- "Dsg3+TEC"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("M") & fastq.merge$Column %in% c(1:6)] <- "Dsg3+TEC"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("M") & fastq.merge$Column %in% c(7:24)] <- "cTEC"
fastq.merge$`Characteristics[cell type]`[fastq.merge$Row %in% c("N", "O", "P") & fastq.merge$Column %in% c(1:24)] <- "cTEC"
fastq.merge$`Factor Value[cell type]` <- fastq.merge$`Characteristics[cell type]`

fastq.merge$`Characteristics[developmental stage]` <- "adult"
fastq.merge$`Characteristics[developmental stage]`[fastq.merge$`Characteristics[age]` %in% c(1, 4)] <- "juvenile"

## add QC pass/fail information
fastq.merge$Sample <- gsub(fastq.merge$`Factor Value[single cell identifier]`, pattern="-", replacement=".")
fastq.merge$`Characteristics[post analysis well quality]`[fastq.merge$Sample %in% thymus.meta$Sample] <- "pass"
fastq.merge$`Characteristics[post analysis well quality]`[!fastq.merge$Sample %in% thymus.meta$Sample] <- "fail"
fastq.merge$`Characteristics[single cell well quality]` <- "Fail"
fastq.merge$`Characteristics[single cell well quality]`[fastq.merge$Sample %in% thymus.meta$Sample] <- "OK"

# add FACS marker information
fastq.merge$`Characteristics[FACS marker]` <- as.character(fastq.merge$`Characteristics[FACS marker]`)
fastq.merge$`Characteristics[FACS marker]`[fastq.merge$`Characteristics[cell type]` %in% c("cTEC")] <- "CD45-EpCAM+Ly51+CD80-"
fastq.merge$`Characteristics[FACS marker]`[fastq.merge$`Characteristics[cell type]` %in% c("mTEChi")] <- "CD45-EpCAM+Ly51-CD80hiMHCIIhi"
fastq.merge$`Characteristics[FACS marker]`[fastq.merge$`Characteristics[cell type]` %in% c("mTEClo")] <- "CD45-EpCAM+Ly51-CD80loMHCIIloDsg3-"
fastq.merge$`Characteristics[FACS marker]`[fastq.merge$`Characteristics[cell type]` %in% c("Dsg3+TEC")] <- "CD45-EpCAM+Ly51-CD80loMHCIIloDsg3+"

fastq.merge$`Characteristics[inferred cell type]` <- as.character(fastq.merge$`Characteristics[inferred cell type]`)
fastq.merge$`Characteristics[inferred cell type]` <- "Unassigned"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(1)]] <- "Post-Aire.mTEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(2)]] <- "Intertypical.TEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(3)]] <- "Mature.cTEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(5)]] <- "Tuft-like.mTEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(6)]] <- "Proliferating.TEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(7)]] <- "Mature.mTEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(8)]] <- "eTEC1"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(9)]] <- "Perintal.cTEC"
fastq.merge$`Characteristics[inferred cell type]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(10)]] <- "eTEC2"

fastq.merge$`Characteristics[single cell well quality]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(4)]] <- "Fail"
fastq.merge$`Characteristics[post analysis well quality]`[fastq.merge$Sample %in% thymus.meta$Sample[thymus.meta$TFIDF.Cluster %in% c(4)]] <- "fail"

write.table(fastq.merge,
            file="~/Dropbox/AgeingExperiment/ArrayExpress_submission.sdrf",
            sep="\t", row.names=FALSE, quote=FALSE)
