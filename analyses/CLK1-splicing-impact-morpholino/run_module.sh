#!/bin/sh

set -e
set -o pipefail

## perform diff expression on ctrl vs clk1-morph
Rscript --vanilla 01-diffExpr-ctrl_vs_morph.R
echo -e "gene\tlog2FoldChange\tpadj" > results/ctrl_vs_treated.de.formatted.tsv
cat results/ctrl_vs_treated.de.tsv | awk -F "\t" '{print $8,"\t"$2"\t"$6}' | awk -F "_" '{print $2}' >> results/ctrl_vs_treated.de.formatted.tsv

## plot differential splicing events between untreated vs treated
echo "DE analysis"
Rscript --vanilla 02-plot_diff-splice-events.R

## find functional sites 
echo "analyze functional sites"
bash 03-run_bedtools_intersect-morpho.sh
Rscript --vanilla 04-plot_diff-func-splice-events.R

## perform ORA of mis-spliced genes after morpholino treatment
echo "ORA analysis"
Rscript --vanilla 05-ora-analysis.R

## Peform GSVA on all cell lines
echo "GSVA analysis"
Rscript --vanilla 06-conduct-gsva-analysis.R

