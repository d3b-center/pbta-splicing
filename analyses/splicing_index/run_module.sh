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
hist_file="$script_directory"/"$1" # ../psi_clustering/input/pbta-histologies.RNA-Seq.initial.tsv
rmats_file="$script_directory"/"$2" #/../data/merge_rMATS_splicing.SE.single.tsv

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

echo "processing PSIs and generating tables"
perl generate_splicing_index_tab_using_tumors.pl $hist_file $rmats_file

echo "plotting"
Rscript splicing_index_tumors.R
