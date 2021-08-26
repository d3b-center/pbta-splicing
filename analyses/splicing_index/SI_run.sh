#!/bin/sh

echo "processing PSIs and generating tables"
perl generate_splicing_index_tab_using_tumors.pl ../psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv ../../data/merge_rMATS_splicing.SE.single.tsv

echo "plotting"
Rscript splicing_index_tumors.R
