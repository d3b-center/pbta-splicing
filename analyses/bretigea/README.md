# BRETIGEA: Brain Cell Type Specific Gene Expression Analysis

Module authors: Komal S. Rathi (@komalsrathi)

## Introduction

Purpose of this module is to do Brain cell type proportion analysis using the R package BRETIGEA. The package computes proportions for the following brain cell types:

1. Astrocytes
2. Endothelial Cells
3. Microglia
4. Neurons
5. Oligodendrocytes
6. Oligodendrocyte PC

## Input

1. Expression: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds` from release hgg-dmg 20201109
2. Metadata: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds`

## Output

Boxplots of plots with brain cell type proportions on the y axis and CLK1-pos vs CLK1-negative samples on the x-axis: `output/bretigea_plot.pdf`

