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

## histology and rmats file
hist_file="../../data/histologies.tsv"
rmats_file="../../data/splice-events-rmats.tsv.gz"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## process PSIs and generate tables for splicing index values for each tumor
bash generate_SBI.sh

## plot values (SBI) generated from above script in CDF plot
echo "plotting splicing burden indices"
Rscript 02-plot_splicing_burden_index.R

## 5% histology specific splicing based on splicing index computations
perl 04-generate_hist_spec_events_tab.pl ../../data/histologies.tsv ../../data/splice-events-rmats.tsv.gz ../../data/independent-specimens.rnaseqpanel.primary.tsv ../../data/independent-specimens.rnaseqpanel.primary-plus.tsv
Rscript 05-plot_histology-specific_splicing_events.R

## differential gene expression
# generate volcano plot of high vs low splicing burden tumors
Rscript 05-plot_diffExp_highlowSBI.R

## plot histology make up per SBI group (high vs low)
Rscript 07-identify_and_plot_histologies_by_SBI.R

## plot TMB vs SBI
Rscript 07-plot_tmb_by_sbi.R
