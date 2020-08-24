#! /usr/bin/bash

## perform per-cluster isoform analysis across age for TECs
## don't run cluster 4 - it's an artifact
CLUSTERS=( 1 2 3 5 6 7 8 9 10 )
OUTPUT_DIR="/nfs/research1/marioni/mdmorgan/Thymus/Splice.dir"

for i in ${CLUSTERS[@]};
do
    OUTFILE=$(eval "echo '/nfs/research1/marioni/mdmorgan/Thymus/Splice.dir/Transcript_PERMANOVA_results-Cluster$i.tsv'")
    JOB="bsub -R "rusage[mem=5000]" -n 1 -T 5 -M 10000 -J PERMANOVA-Cluster$i -q research-rh7 -e $OUTPUT_DIR/Cluster$i.err -o $OUTPUT_DIR/Cluster$i.out R CMD /nfs/research1/marioni/mdmorgan/Thymus/scripts/permutation_isoforms-by_cluster.R $i $OUTFILE"
    
    #echo $JOB
    eval $JOB
done
