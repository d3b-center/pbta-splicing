#!/bin/sh

## pre-process histology files to only RNA-seq and 8 histologies (HGAT, LGAT, Med, Oligo, Cranio, ATRT, Epend, Gang)
cat input/pbta-histologies.tsv | grep RNA-Seq | awk '{if ( ($NF~/HGAT/) || ($NF~/LGAT/) || ($NF~/Oligo/) || ($NF~/Med/) || ($NF~/Gang/) || ($NF~/Epend/) || ($NF~/ATRT/) || ($NF~/Crani/) ){ print $0}}' > input/pbta-histologies.RNA-Seq.initial.tsv

## create PSI matrix keeping only one patient sample for downstream clustering
perl create_matrix_of_PSI.pl input/pbta-histologies.RNA-Seq.initial.tsv <rMATs_file_paths.txt> input/compiled_molecular_subtypes_with_clinical_pathology_feedback.tsv

## run clustering method
Rscript consensus_clustering.R

## find cluster contributors
python vtest_calc.py
