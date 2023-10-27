#!/bin/sh

## histology and rmats file
hist_file="input/histologies-plot-group.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples_plus="../../data/independent-specimens.rnaseqpanel.primary-plus.tsv"

## 2% histology specific splicing based on splicing index computations for SE
perl 01-generate_hist_spec_events_tab.pl $hist_file $rmats_file $indep_samples_plus
Rscript --vanilla 02-plot_histology-specific_splicing_events.R

