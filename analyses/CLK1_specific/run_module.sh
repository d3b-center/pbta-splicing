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
Rscript splicing_functional_sites.R results/splicing_events.total.pos.intersectUnip.ggplot.txt results/splicing_events.total.neg.intersectUnip.ggplot.txt
Rscript CLK1_EI_vs_ES_PSI_volcano.R   
Rsript PSI_ggstatplot_across_types.R 
Rscript highExon4_vs_lowExon4_diffExpr_volcano.R

