#!/bin/sh

## get splicing information from input list and extract PSI for each sample
cat input/diffSplicing_cand.filterDiffExpr.txt | awk '{print $2}' | perl extract_psi_from_list.pl

## for each gene and its splicing change, make table of PSI and expression
ls results/results_diff/*txt | xargs -n 1 echo "perl quant_expr_vs_psi_PCT.pl ~/Desktop/AS-DMG/data/stranded_trans_rsem_counts_tab.tsv ~/Desktop/AS-DMG/data/stranded_trans_rsem_counts_tab.tsv " | bash

## compute correlation between PSI vs expression for each gene
perl corr_calc.pl | grep PSI | grep -v ^PSI > results/expr_vs_psi_corr_res.txt
