#!/bin/sh

## histology and rmats file
hist_file="../cohort_summary/results/histologies-plot-group.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"


## process PSI and generate SBI tables
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples SE
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples A3SS
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples A5SS
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples RI

## plot values (SBI) generated from above script in CDF plot
echo "creating SBI CDF plots"
Rscript --vanilla 02-plot_splicing_burden_index.R

## plot tumors based on high vs low SBI
echo "plotting SBI high/low"
Rscript --vanilla 03-identify_and_plot_histologies_by_SBI.R

## plot distrubution of splicing types/cases
echo "plotting distributions"
Rscript --vanilla 04-plot_total-splicing-cases.R

## plot tmb vs sbi
echo "plotting SBI vs TMB"
Rscript --vanilla 05-plot-tmb-vs-sbi.R
