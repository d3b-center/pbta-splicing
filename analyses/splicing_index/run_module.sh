#!/bin/sh

## histology and rmats file
hist_file="../../data/histologies.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"
indep_samples_plus="../../data/independent-specimens.rnaseqpanel.primary-plus.tsv"

## process PSI and generate SBI tables
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples $indep_samples_plus SE
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples $indep_samples_plus A3SS
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples $indep_samples_plus A5SS
perl 01-generate_splicing-index_and_diff-events_table.pl $hist_file $rmats_file $indep_samples $indep_samples_plus RI

## plot values (SBI) generated from above script in CDF plot
echo "plotting splicing burden indices"
Rscript 02-plot_splicing_burden_index.R

## perform differential gene expression analyses on high vs low SBI tumors (HGGs)
Rscript 03-plot_diffExp_highlowSBI.R

## plot distrubution of splicing types/cases
Rscript 04-plot_total-splicing-cases.R
