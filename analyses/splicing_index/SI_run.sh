#!/bin/sh

echo "processing PSIs and generating tables"
perl ../psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv ~/Desktop/AS-DMG/analyses/merge_rMATS/merge_rMATS_splicing.SE.single.tsv

echo "plotting"
Rscript splicing_index_tumors.R
