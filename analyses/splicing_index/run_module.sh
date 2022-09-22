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

## hsitolgy and rmats file
hist_file="../../data/pbta-histologies.RNA-Seq.initial.tsv"
rmats_file="../../data/merge_rMATS_splicing.SE.single.tsv"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## process PSIs and generate tables for splicing index values for each tumor
echo "processing PSIs and generating tables"
perl generate_splicing_index_tab_using_tumors.pl $hist_file $rmats_file

echo "plotting SI"
Rscript 01-splicing_index_tumors.R

## 10% histology specific splicing based on splicing index computations
perl generate_hist_spec_events_tab_using_tumors.pl $hist_file $rmats_file
Rscript 02-hist_specific_splicing_tumors.R

