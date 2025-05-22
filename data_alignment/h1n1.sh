#!/bin/bash

# source augur

refseq_path=inputs/references/h1n1_ha.fasta
target_file=inputs/h1n1/gisaid_epiflu_sequence.fasta
aligned_file=alignments/aligned_h1n1_ha.fasta

nthreads=3

mkdir alignments

augur align --reference-sequence $refseq_path --sequences $target_file --output $aligned_file --nthreads $nthreads --remove-reference True

