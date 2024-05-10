#!/bin/sh
#
# Script to run all analysis modules
# Author: Your Name
# Date: 2024/05/10
#

## cohort summary
echo "cohort summary"
Rscript cohort_summary/01-generate-cohort-summary-circos-plot.R

## stranded-polyA-assessment
echo stranded-polyA-assessment
bash stranded-polyA-assessment/run_module.sh

## histology-specific splice events
echo "histology-specific splice events"
bash histology-specific-splicing/run_module.sh

## splicing index
echo "splicing index"
bash splicing_index/run_module.sh

## survival
echo "survival"
bash survival/run-survival-module.sh

## splicing-factor_dysregulation
echo "splicing factor dysregulation"
bash splicing-factor_dysregulation/run_module.sh

## splicing events functional sites
echo "splicing events functional sites"
bash splicing_events_functional_sites/run_module.sh

## CLK1 splicing correlations
echo "CLK1 splicing correlations"
bash CLK1-splicing_correlations/run_module.sh

## oncoprint
echo "oncoprint"
bash oncoprint/01-oncoprint.R

## long-read-CLK1-validation
echo "long-read-CLK1-validation"
bash long-read-CLK1-validation/run_module.sh

## KNS42-cell-line
echo "KNS42-cell-line"
bash KNS42-cell-line/run_module.sh

## CLK1-splicing-impact-morpholino
echo "CLK1-splicing-impact-morpholino"
bash CLK1-splicing-impact-morpholino/run_module.sh
