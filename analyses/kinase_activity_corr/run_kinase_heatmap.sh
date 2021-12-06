#!/bin/bash

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Define directory
data_dir="../../data"
ref_dir="../../references"
plots_dir="plots"
input_dir="input"
pki_input_dir="input/individual_files"

# Define files used
histology_file="${data_dir}/pbta-histologies.tsv"
kinase_activity_file="${input_dir}/supplementarytable_kinase.xlsx"
clk1_psi_group_file="${input_dir}/CLK1_PSI.pan-cancer.groupings.txt"


Rscript -e "rmarkdown::render('kinase_activity_corr.Rmd', clean = TRUE)"

# Run kinase activity per gene for all genes of interest
Rscript --vanilla kinase_activity_per_gene.R \
--histology $histology_file \
--kinase $kinase_activity_file \
--group $clk1_psi_group_file \
--pki_input_dir $pki_input_dir 
