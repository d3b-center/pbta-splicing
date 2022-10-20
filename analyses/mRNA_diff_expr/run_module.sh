#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

## generate volcano plot of SFs between ctrl vs hgat
Rscript 01-volcano_plot_mRNA.R

## SRSF11 specific plots
Rscript 02-SRSF11_plots.R

## correlation plots for SRSF22 and RBM5
Rscript 03-plot_prot_vs_rna.R
