#!/bin/bash

set -e
set -o pipefail

# run prepare-survival script
Rscript -e "rmarkdown::render('01-prepare-survival.Rmd')"

# run run-survival scripts
Rscript -e "rmarkdown::render('02-run-survival-SIgroup.Rmd')"
Rscript -e "rmarkdown::render('03-run-survival-SI.Rmd')"

# run plot-survival script
Rscript 04-plot-survival.R
