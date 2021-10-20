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
hist_file="$script_directory"/"$1"
rmats_file="$script_directory"/"$2"

echo "hist file:" $hist_file
echo "rmats file:" $rmats_file

## create PSI matrix keeping only one patient sample for downstream clustering
# ../../data/v19_plus_20210311_pnoc_rna.tsv ../../data/merge_rMATS_splicing.SE.single.tsv
echo "create PSI matrix keeping only one patient sample for downstream clustering..."
perl create_matrix_of_PSI_SE.pl $hist_file $rmats_file

## run clustering method, save in plots folder
echo "run clustering method, save in plots folder..."
Rscript consensus_clustering.R results/pan_cancer_splicing_SE.txt

## stacked barplots of cluster memberships, results in plots folder // not done yet
echo "make stacked barplots of cluster memberships, results in plots folder..."
Rscript stacked_barplots.R results/CC_groups.txt $hist_file

## find cluster contributors, files saved in results folder
echo "find cluster contributors, files saved in results folder..."

##delete first column name in first line only
cat results/pan_cancer_splicing_SE.txt | perl -pe 's/Splice_ID\t//g' > results/pan_cancer_splicing_SE.vtest_tmp.txt
python vtest_calc.py -i results/pan_cancer_splicing_SE.vtest_tmp.txt -c results/CC_groups.txt -t Cluster -o plots/vtest_res. -v results/vtest_calc.tsv
