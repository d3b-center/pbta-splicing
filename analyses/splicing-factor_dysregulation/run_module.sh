#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit


## generate volcano plot of SFs within HGG tumors
Rscript --vanilla 01-plot-diffExp_highlowSBI.R

## plot heatmap of select splicing factor RNA vs proteo obtained from CPTAC
Rscript --vanilla 02-plot_SFs_rna_vs_prot.R
