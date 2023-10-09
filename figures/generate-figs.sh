#!/bin/bash
#
# Run all figure making scripts.

# enviroment settings
set -e
set -o pipefail

# If RUN_LOCAL is used, the time-intensive steps are skipped because they cannot
# be run on a local computer -- the idea is that setting RUN_LOCAL=1 will allow for
# local testing running/testing of all the other figures
RUN_LOCAL=${RUN_LOCAL:-0}

# Find current directory based on this script
WORKDIR=$(dirname "${BASH_SOURCE[0]}")
cd "$WORKDIR"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"
data_dir="$BASEDIR/data"

# Make output folders for all figures
mkdir -p tiffs/fig1

############################# Figure panels ##########################################

# Below, we establish figure panels in `figures/pdf/` for all manuscript figures.
# For each, we first run any `figures/scripts/` scripts that generate figure panels,
#  and then we copy any additional panels that were generated in `analyses/` modules.

##### Figure 1: Workflow and sample distribution ------------------------------------
# Copy the main figure to final directory - panels A,B

##Figs 1
cp ${analyses_dir}/splicing_index/plots/splice-types.pdf pdfs/fig1/.
cp ${analyses_dir}/splicing_index/plots/SBI-plot.SE.tiff pdfs/fig1/.
cp ${analyses_dir}/splicing_index/plots/hist_by_sbi-level_lolliplot.pdf pdfs/fig1/.
cp ${analyses_dir}/survival/plots/km_hgg_k28_OS_EFS_SIburden.pdf pdfs/fig1/.
cp ${analyses_dir}/survival/plots/km_hgg_k28_OS_EFS_SIburden.pdf pdfs/fig1/.
cp ${analyses_dir}/survival/plots/forest_hgg_OS_subtype_SI.pdf pdfs/fig1/.
cp ${analyses_dir}/splicing_index/plots/enhancedVolcano_hgg_sbi.pdf pdfs/fig1/.
cp ${analyses_dir}/histology-specific-splicing/plots/upsetR_histology-specific.es.tiff pdfs/fig1/.
cp ${analyses_dir}/histology-specific-splicing/plots/upsetR_histology-specific.ei.tiff pdfs/fig1/.
