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
Rscript 02-plot_CLK1_high_v_low_Exon4_volcano.R

## run expression vs splicing correlation analyses ##
## get splicing information from input list and extract PSI for each sample
cat input/diffSplicing_cand.filterDiffExpr.txt | awk '{print $2}' | perl 03-extract_psi_for_splice_variants.pl

## for each gene and its splicing change, make table of PSI and expression
ls results/results_diff/*txt  | xargs -n 1 echo "perl 04-generate_expr_vs_psi_table.pl input/stranded_trans_rsem_counts_tab.tsv  " | bash

## compute correlation between PSI vs expression for each gene
perl 05-compute_corr_calc.pl  | grep PSI | grep -v ^PSI > results/expr_vs_psi_corr_res.txt
rm results/results_diff/*txt

## plot correlations
Rscript 06-plot_psi_vs_expr.R
