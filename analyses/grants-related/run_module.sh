#!/bin/sh

#SE
perl process_rmats_SE.pl ../input/v19_plus_20210311_pnoc_rna.primary.tsv
bedtools intersect -wo -a splicing_neoepitope.pri.dpsi.se.bed -b ../data/unipLocExtra.hg38.col.txt > splicing_neoepitope.primary.dpsi.se.unipLocExtra.bed
perl generate_tab.se.pl splicing_neoepitope.primary.dpsi.se.unipLocExtra.bed splicing_neoepitope.primary.dpsi.se.unipLocExtra.bed ../input/hgg_combined_tab.tpm.csv ../input/gene_counts_normals_final.csv > 
splicing_neoepitope.primary.dpsi.se.unipLocExtra.summary.txt

#MXE
#perl process_rmats_MXE.pl ~/Desktop/pbta-splicing/data/v19_plus_20210311_pnoc_rna.tsv
#bedtools intersect -wo -a splicing_neoepitope.dpsi.mxe.bed -b data/unipLocExtra.hg38.col.txt > splicing_neoepitope.dpsi.mxe.unipLocExtra.bed
#perl generate_tab.mxe.pl splicing_neoepitope.dpsi.mxe.unipLocExtra.bed splicing_neoepitope.dpsi.mxe.tsv input/hgg_combined_tab.tpm.csv input/gene_counts_normals_final.csv > splicing_neoepitope.dpsi.mxe.unipLocExtra.summary.txt
