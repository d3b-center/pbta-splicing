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

##default file paths for histology and rmats output
hist_file="../../data/v19_plus_20210311_pnoc_rna.tsv"
rmats_file="../../data/merge_rMATS_splicing.SE.single.tsv"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## create PSI matrix keeping only one patient sample for downstream clustering
# ../../data/v19_plus_20210311_pnoc_rna.tsv ../../data/merge_rMATS_splicing.SE.single.tsv
echo "create PSI matrix keeping only one patient sample for downstream clustering..."
perl create_matrix_of_PSI_SE_gene.pl $hist_file $rmats_file

## run clustering method, save in plots folder
echo "run clustering method, save in plots folder..."
Rscript consensus_clustering.R results/pan_cancer_splicing_SE.txt
