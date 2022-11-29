#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo $script_directory
input_file="$script_directory"/""

## make plots
Rscript 01-plot_CLK1_EI_vs_ES_PSI_volcano.R
Rscript 02-plot_highExon4_vs_lowExon4_diffExpr_volcano.R

## run expression vs splicing correlation analyses
## get splicing information from input list and extract PSI for each sample
cat input/diffSplicing_cand.filterDiffExpr.txt | awk '{print $2}' | perl extract_psi_from_list.pl

## for each gene and its splicing change, make table of PSI and expression
ls results/scr/*txt | xargs -n 1 echo "perl 03-extract_psi_for_splice_variants.pl input/stranded_trans_rsem_counts_tab.tsv  " | bash

## compute correlation between PSI vs expression for each gene
perl corr_calc.pl | grep PSI | grep -v ^PSI > results/expr_vs_psi_corr_res.txt
