#!/bin/sh

## get relevant splice events
Rscript 01-identify-and-create-splice-variant-matrix.R

## convert to bed
cat results/mixed_events.tsv | awk '{print $2"\t"$3"\t"$4"\t"$6"\t0\t"$5}' | grep "_" > results/mixed_events.bed

## get nucleotide sequence
bedtools getfasta -fi ~/Downloads/hg38.fa -bed results/mixed_events.bed > results/mixed_events.fasta

perl convert_to_peptide.pl results/mixed_events.fasta
