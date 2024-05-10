#!/bin/sh
#
# Script to run all analysis modules
# Author: Your Name
# Date: 2024/05/10
#


set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

## cohort summary
echo "cohort summary"
#Rscript cohort_summary/01-generate-cohort-summary-circos-plot.R

## stranded-polyA-assessment
echo stranded-polyA-assessment
cd ${analyses_dir}/stranded-polyA-assessment
bash run_module.sh

## histology-specific splice events
echo "histology-specific splice events"
cd ${analyses_dir}/histology-specific-splicing

bash run_module.sh
