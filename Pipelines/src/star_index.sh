#! /usr/bin/bash

# create a star index from the mm10-ERCC92 genome
#BSUB-n 12
#BSUB-q research
#BSUB-M 128000
#BSUB-R "rusage[mem=128000]"
#BSUB-o /nfs/leia/research/marioni/genomes/indices/mm10_ERCC92.stdout
#BSUB-e /nfs/leia/research/marioni/genomes/indices/mm10_ERCC92.stderr

/nfs/software/marioni/STAR-2.5.3a/bin/Linux_x86_64_static/STAR --runThreadN 12 --runMode genomeGenerate  --genomeDir /nfs/leia/research/marioni/genomes/indices  --genomeFastaFiles /nfs/leia/research/marioni/genomes/mm10_ERCC92.fasta  --sjdbGTFfile /nfs/leia/research/marioni/mikemorgan/Thymus/genes/genes.gtf --sjdbOverhang 99
