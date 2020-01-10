#! /usr/bin/bash

## loop over sample and run emptyDrops + cell filtering on each sample separately
SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/src"
LOGDIR="/nfs/research1/marioni/mdmorgan/Thymus_ZsG/logs"
SAMPDIR="/nfs/research1/marioni/mdmorgan/Thymus_ZsG"

for fle in `ls $SAMPDIR | grep "^Ageing*" | grep "ADT" | grep -v "tsv"`;
do
    echo $fle
    LOGERROR=$LOGDIR/$fle.err
    LOGOUT=$LOGDIR/$fle.out

    JOB="bsub -R "rusage[mem=32000]" -M 64000 -T 15000 -q research-rh74 -J $fle -e $LOGERROR -o $LOGOUT R CMD $SRCDIR/ZsG_emptyDrops_ADT.R $fle"
    echo $JOB
    eval $JOB
done
