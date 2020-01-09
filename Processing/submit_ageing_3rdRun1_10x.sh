#! /usr/bin/bash

## submit a 10X cell ranger job for the ageing ZsGreen-lineage tracing experiment 3rd sample run 1
LOG="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/logs/Ageing_10X_ZsG"
SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/src"

HTO_LIBFILE="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/lib_files/AgeingExperiment_3rd_Run1_HTO.tsv"  
ADT_LIBFILE="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/lib_files/AgeingExperiment_3rd_Run1_ADT.tsv"

ADT_REFLIB="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/refs/Ageing_3rd_Run1_ADT_refs.tsv"
HTO_REFLIB="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/refs/Ageing_3rd_Run1_HTO_refs.tsv"

# I need to run the feature barcoding for the CITE-seq antibodies and hashtags separately

ADT_ID="Ageing_ZsG_3rdRun1_ADT"
HTO_ID="Ageing_ZsG_3rdRun1_HTO"

HTO_JOB="bsub -R "rusage[mem=24000]" -M 36000 -T 50000 -q research-rh74 -J $HTO_ID -e $LOG.$HTO_ID.err -o $LOG.$HTO_ID.out bash $SRCDIR/submit_10X_cellranger.sh $HTO_LIBFILE $HTO_REFLIB $HTO_ID"
ADT_JOB="bsub -R "rusage[mem=24000]" -M 36000 -T 50000 -q research-rh74 -J $ADT_ID -e $LOG.$ADT_ID.err -o $LOG.$ADT_ID.out bash $SRCDIR/submit_10X_cellranger.sh $ADT_LIBFILE $ADT_REFLIB $ADT_ID"

echo $HTO_JOB
#$HTO_JOB

echo $ADT_JOB
#exec $ADT_JOB
