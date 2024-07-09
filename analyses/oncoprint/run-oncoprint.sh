#!/bin/sh

## plot oncoprint
echo "plotting oncoprint"
Rscript --vanilla 01-oncoprint.R

Rscript --vanilla 02-oncoprint-SFs.R
