#!/bin/sh

## DepMap portal cell-lines analyses
Rscript --vanilla 01-plot-and-process-depmap-data.R

## cell prolif assay
Rscript --vanilla 02-plot_cell-proliferation-assay-res.R

## PCR assay
Rscript --vanilla 03-plot-qPCR-results.R
