#!/bin/bash

# source augur

refseq_path=/app/inputs/references/H5N1_ha.fasta
target_file=/app/inputs/H5N1/gisaid_epiflu_sequence.fasta
aligned_file=/app/alignments/aligned_H5N1_ha.fasta

nthreads=3

augur align --reference-sequence $refseq_path --sequences $target_file --output $aligned_file --nthreads $nthreads --remove-reference

