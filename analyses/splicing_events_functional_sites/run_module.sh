#!/bin/sh

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

#echo $script_directory
#input_file="$script_directory"/"$1"
input_file = "../../data/v19_plus_20210311_pnoc_rna.tsv"

echo "input file:" $input_file

echo "process rMATS with .20 dPSI and 10 junction read counts...";

## Process rMATS files given histologies file. Keep only HGG midlines samples and storng splicing events
# ../../data/v19_plus_20210311_pnoc_rna.tsv
perl extract_recurrent_splicing_events_hgg.pl $input_file

echo "bedtools intersect...";
./bedtools_intersect.sh

echo "make tab for ggplot ...";
./generate_ggplot_tab.sh

## Remove intermediatery files / off for now
rm results/splicing_events.total.*intersectUnipMod.wo.txt

## make plots
Rscript splicing_functional_sites.R results/splicing_events.total.pos.intersectUnip.ggplot.txt results/splicing_events.total.neg.intersectUnip.ggplot.txt
