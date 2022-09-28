#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

## histology input file (column orders important)
input_file = "../../data/v19_plus_20210311_pnoc_rna.tsv"

echo "input file:" $input_file
echo "process rMATS with .20 dPSI and 10 junction read counts...";

## Process rMATS files given histologies file. Keep only HGG midlines samples and storng splicing events
# ../../data/v19_plus_20210311_pnoc_rna.tsv
perl 01-extract_recurrent_splicing_events_hgg.pl $input_file

echo "bedtools intersect...";
bash 01-run_bedtools_intersect.sh

echo "make tab for ggplot ...";
bash 02-format_for_ggplot.sh

## Remove intermediatery files / off for now
rm results/splicing_events.total.*intersectUnipMod.wo.txt

## make plots
Rscript 01-plot_splicing_across_functional_sites.R 
Rscript 02-plot-flip_mixed_events.R
