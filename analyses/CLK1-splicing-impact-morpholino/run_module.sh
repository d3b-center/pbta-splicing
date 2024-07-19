#!/bin/sh

set -e
set -o pipefail

## Identify the transcripts of interest per splice variants of interest
Rscript --vanilla 00-get-splice-transcripts.R

## perform diff expression on ctrl vs clk1-morph
Rscript --vanilla 01-diffExpr-ctrl_vs_morph.R
echo -e "gene\tlog2FoldChange\tpadj" > results/ctrl_vs_treated.de.formatted.tsv
cat results/ctrl_vs_treated.de.tsv | awk -F "\t" '{print $9"\t"$2"\t"$6}' | grep -v log   >> results/ctrl_vs_treated.de.formatted.tsv

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

## Compares GSVA scores for CLK1 morpholino treated vs non-targeting morpholino and plot
echo "GSVA comparisons and plots"
Rscript -e "rmarkdown::render('07-run-gsva-comparisons.Rmd', clean = TRUE)"

## lolliplot of splicing cases
Rscript --vanilla 08-plot_total-splicing-cases.R

## intersect DE and DS
Rscript --vanilla 09-intersection-dex-des.R

## intersect with crispr scores from CCMA
Rscript --vanilla 10-crispr-screen-intersection.R
