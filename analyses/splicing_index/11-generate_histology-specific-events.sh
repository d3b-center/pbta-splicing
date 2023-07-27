#!/bin/sh

## histology and rmats file
hist_file="../../data/histologies.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"
indep_samples="../../data/independent-specimens.rnaseqpanel.primary.tsv"
indep_samples_plus="../../data/independent-specimens.rnaseqpanel.primary-plus.tsv"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## 10% histology specific splicing based on splicing index computations
perl 12-generate_splicing-index_and_diff-events_table_histology-spec.SE.pl $hist_file $rmats_file
Rscript 13-plot_histology-specific_splicing_events.R
