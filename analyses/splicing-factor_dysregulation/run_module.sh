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

## prepare tables of SRSF11 PSI values
bash 02-identify_preprocess_SRSF11_psi.sh

## SRSF11 specific plots
Rscript 03-SRSF11_plots.R

## correlation plots for SRSF22 and RBM5
#Rscript 04-plot_prot_vs_rna.R

## plot heatmap of select splicing factor RNA vs proteo obtained from CPTAC
Rscript 04-plot_SFs_rna_vs_prot.R
