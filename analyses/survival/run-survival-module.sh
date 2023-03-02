#!/bin/bash

set -e
set -o pipefail

# run prepare-survival script
Rscript -e "rmarkdown::render('01-prepare-survival.Rmd')"

# run run-survival script
Rscript -e "rmarkdown::render('02-run-survival.Rmd')"

# run plot-survival script
Rscript 03-plot-survival.R
