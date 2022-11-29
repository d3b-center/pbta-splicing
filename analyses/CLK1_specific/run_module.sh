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
input_file="$script_directory"/""

## make plots
Rscript 01-plot_CLK1_EI_vs_ES_PSI_volcano.R
Rscript 02-plot_CLK1_high_v_low_Exon4_volcano.R
