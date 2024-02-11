#!/bin/sh

## histology and rmats file
hist_file="../cohort_summary/results/histologies-plot-group.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"

## 2% histology specific splicing based on splicing index computations for SE
perl 01-generate_hist_spec_events_tab.pl $hist_file $rmats_file $indep_samples
Rscript --vanilla 02-plot_histology-specific_splicing_events.R

rm Rplots.pdf
