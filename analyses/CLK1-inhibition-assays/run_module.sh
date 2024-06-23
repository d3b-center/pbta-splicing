#!/bin/bash

## plot cell proliferation results
Rscript --vanilla 01-plot_cell-proliferation-assay.R

## plot cell viability results
Rscript --vanilla 02-plot_ctg-assay.R
