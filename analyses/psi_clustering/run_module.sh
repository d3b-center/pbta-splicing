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
input_file="$script_directory"/"$1"

echo "input file:" $input_file

## pre-process histology files to only RNA-seq and 8 histologies (HGAT, LGAT, Med, Oligo, Cranio, ATRT, Epend, Gang)
echo "pre-process histology files to only RNA-seq and 8 histologies (HGAT, LGAT, Med, Oligo, Cranio, ATRT, Epend, Gang)..."
#cat input/v19_plus_20210311_pnoc.tsv | grep RNA-Seq | awk '{if ( ($NF~/HGAT/) || ($NF~/LGAT/) || ($NF~/Oligo/) || ($NF~/Med/) || ($NF~/Gang/) || ($NF~/Epend/) || ($NF~/ATRT/) || ($NF~/Crani/) ){ print $0}}' > input/pbta-histologies.RNA-Seq.filtered.tsv
#cat $input_file | grep RNA-Seq | awk '{if ( ($NF~/HGAT/) || ($NF~/LGAT/) || ($NF~/Oligo/) || ($NF~/Med/) || ($NF~/Gang/) || ($NF~/Epend/) || ($NF~/ATRT/) || ($NF~/Crani/) ){ print $0}}' > input/pbta-histologies.RNA-Seq.filtered.v2.tsv

## create PSI matrix keeping only one patient sample for downstream clustering
echo "create PSI matrix keeping only one patient sample for downstream clustering..."
perl create_matrix_of_PSI_SE.pl ../../data/v19_plus_20210311_pnoc_rna.tsv ../../data/merge_rMATS_splicing.SE.single.tsv

## run clustering method, save in plots folder
echo "run clustering method, save in plots folder..."
Rscript consensus_clustering.R results/pan_cancer_splicing_SE.txt

## stacked barplots of cluster memberships, results in plots folder // not done yet
echo "make stacked barplots of cluster memberships, results in plots folder..."
Rscript stacked_barplots.R results/CC_groups.txt ../../data/v19_plus_20210311_pnoc_rna.tsv

## find cluster contributors, files saved in results folder
echo "find cluster contributors, files saved in results folder..."

##delete first column name in first line only
cat results/pan_cancer_splicing_SE.txt | perl -pe 's/Splice_ID\t//g' > results/pan_cancer_splicing_SE.vtest_tmp.txt

python vtest_calc.py -i results/pan_cancer_splicing_SE.vtest_tmp.txt -c results/CC_groups.txt -t Cluster -o plots/vtest_res. -v results/vtest_calc.tsv
