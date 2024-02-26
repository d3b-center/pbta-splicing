#!/bin/sh

set -e
set -o pipefail

## perform diff expression on ctrl vs clk1-morph
Rscript 01-diffExpr-ctrl_vs_morph.R
echo -e "gene\tlog2FoldChange\tpadj" > results/ctrl_vs_treated.de.formatted.tsv
cat results/ctrl_vs_treated.de.tsv | awk -F "\t" '{print $8,"\t"$2"\t"$6}' | awk -F "_" '{print $2}' >> results/ctrl_vs_treated.de.formatted.tsv

## perform GSEA on diff expression results
Rscript 02-gsea-analysis.R

## plot differential splicing events between untreated vs treated
Rscript 03-plot_diff-splice-events.R

## perform ORA of mis-spliced genes after morpholino treatment
Rscript 04-ora-analysis.R

## Peform GSVA on all cell lines
Rscript --vanilla 05-conduct-gsva-analysis.R

