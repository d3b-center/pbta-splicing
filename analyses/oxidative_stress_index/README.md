# Oxidative stress index

Module authors: Komal S. Rathi (@komalsrathi)

## Introduction

Purpose of this module is to do geneset variation analysis using R package GSVA with the following C5 BP genesets:

1. GO_RESPONSE_TO_OXIDATIVE_STRESS
2. GO_RESPONSE_TO_HYDROGEN_PEROXIDE
3. GO_CELLULAR_RESPONSE_TO_HYDROGEN_PEROXIDE

## Input

1. Expression: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds` from release hgg-dmg 20201109
2. Metadata: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds`

## Output

Boxplots of plots with Oxidative stress repair index on the y axis and CLK1-pos vs CLK1-negative samples on the x-axis: `output/oxidative_stress_index_plot.pdf`