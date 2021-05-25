#!/bin/sh

echo "processing PSIs and generating tables"
perl generate_splicing_index_tab.pl ../psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/pan_cancer_rmats/filtered_samples_files.v2.txt

echo "plotting"
Rscript splicing_index.R
