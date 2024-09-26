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

# run HGG survival by CLK1 status
Rscript -e "rmarkdown::render('05-survival-hgg-clk1-status.Rmd')"

# survival analyses by splicing cluster
Rscript -e "rmarkdown::render('06-survival_by_cluster.Rmd')"

# run survival by NF1 PSI status
Rscript -e "rmarkdown::render('07-run-survival-nf1-psi.Rmd')"

# plot survival by NF1 PSI status
Rscript 08-plot-survival-nf1-psi.R