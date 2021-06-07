# Proliferative index

Module authors: Komal S. Rathi (@komalsrathi)

## Introduction

Purpose of this module is to run OCLR based Stemness profiling. 
Paper: https://www.cell.com/cell/pdf/S0092-8674(18)30358-1.pdf
Tutorial: http://tcgabiolinks.fmrp.usp.br/PanCanStem/mRNAsi.html and https://bioinformaticsfmrp.github.io/PanCanStem_Web/.

## Input

1. Expression: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds` from release hgg-dmg 20201109
2. Metadata: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds`

## Utils

1. `util/genes2hugo.R`: convert entrez ids to hugo symbols
2. `util/main.train.R`: train using in-built Progenitor Cell Biology Consortium (PCBC) dataset. The output is generated in `input/pcbc-stemsig.tsv` that has genes and associated scores.
3. `util/main.predict.R`: predict stemness using scores generated in main.train.R

## Output

Boxplots of plots with Proliferative index on the y axis and CLK1-pos vs CLK1-negative samples on the x-axis: `output/stemness_index_plot.pdf`

Table of scores: `output/stemness_scores.tsv`