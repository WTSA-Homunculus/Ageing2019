#! /usr/bin/bash

#BSUB-n 1
#BSUB-q research-rh7
#BSUB-M 54000                                                                                                                                                                                      
#BSUB-R "rusage[mem=54000]"
#BSUB-e /nfs/research1/marioni/mdmorgan/Thymus/merge_counts.err
#BSUB-o /nfs/researach1/marini/mdmorgan/Thymus/merge_counts.out

INPUT_DIR="/nfs/research1/marioni/mdmorgan/Thymus/FC_QUANT/"
INPUT_REGEX="gene_counts.txt$"
#ALL_FILES=$(eval "find /nfs/research1/marioni/mdmorgan/Thymus/FC_QUANT/*gene_counts.txt | tr -s '\n' ','")
#FILE_LIST=${ALL_FILES%?}

python3 /nfs/research1/marioni/mdmorgan/Thymus/scripts/merge_counts.py --input-directory=$INPUT_DIR --file-regex=$INPUT_REGEX  > /nfs/research1/marioni/mdmorgan/Thymus/FC_QUANT/Merged_counts.txt
