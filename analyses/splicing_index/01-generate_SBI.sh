#!/bin/sh

## histology and rmats file
hist_file="../../data/histologies.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"
indep_samples_plus="../../data/independent-specimens.rnaseqpanel.primary-plus.tsv"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## process PSI and generate SBI tables
perl 03-generate_splicing-index_and_diff-events_table.SE.pl $hist_file $rmats_file $indep_samples $indep_samples_plus
perl 04-generate_splicing-index_and_diff-events_table.RI.pl $hist_file $rmats_file $indep_samples $indep_samples_plus
perl 05-generate_splicing-index_and_diff-events_table.A5SS.pl $hist_file $rmats_file $indep_samples $indep_samples_plus
perl 06-generate_splicing-index_and_diff-events_table.A3SS.pl $hist_file $rmats_file $indep_samples $indep_samples_plus

## plot values (SBI) generated from above script in CDF plot
echo "plotting splicing burden indices"
Rscript 07-plot_splicing_burden_index_SE.R
Rscript 08-plot_splicing_burden_index_RI.R
Rscript 09-plot_splicing_burden_index_A5SS.R
Rscript 10-plot_splicing_burden_index_A3SS.R
