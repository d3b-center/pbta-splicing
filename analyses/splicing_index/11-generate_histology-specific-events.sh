#!/bin/sh

## histology and rmats file
hist_file="../../data/histologies.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"
indep_samples_plus="../../data/independent-specimens.rnaseqpanel.primary-plus.tsv"

## 10% histology specific splicing based on splicing index computations
perl 12-generate_hist_spec_events_tab.pl $hist_file $rmats_file $indep_samples $indep_samples_plus 
Rscript 13-plot_histology-specific_splicing_events.R
