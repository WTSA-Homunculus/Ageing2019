#! /usr/bin/bash

## Run multiple iterations of TCR simulations - 10 runs?
TCR_FASTA="/nfs/research1/marioni/mdmorgan/TCR_simulations/TCR_FASTA"
N_CELLS=1000000
OUTPUT_STUB="/nfs/research1/marioni/mdmorgan/TCR_simulations/results/TCR_simulations"
N_RUNS=(1 2 3 4 5 6 7 8 9 10)


for X in ${N_RUNS[@]};
do
    JOB="bsub -R"rusage[mem=12000]" -M 16000 -T 70 -q research-rh74 -J TCR_sim-Run$X -e $OUTPUT_STUB-Run$X.err -o $OUTPUT_STUB-Run.$X.out   python3 /nfs/research1/marioni/mdmorgan/TCR_simulations/src/TCR_generator.py --fasta-directory=$TCR_FASTA --n-cells=$N_CELLS --output-filename=$OUTPUT_STUB-Run$X.txt"
    echo $JOB
    eval $JOB
done
