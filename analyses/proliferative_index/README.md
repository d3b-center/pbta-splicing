# Proliferative index

Module authors: Komal S. Rathi (@komalsrathi)

## Introduction

Purpose of this module is to do geneset variation analysis using R package GSVA with `Table S2. Cell cycle control gene set` obtained from:
Yuan, J., Levitin, H.M., Frattini, V. et al. Single-cell transcriptome analysis of lineage diversity in high-grade glioma. Genome Med 10, 57 (2018). https://doi.org/10.1186/s13073-018-0567-9
PMID: 30041684

## Input

1. Expression: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds` from release hgg-dmg 20201109
2. Metadata: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds`
3. Cell-cycle genelist: `input/cell-cycle-genes.txt`

## Output

Boxplots of plots with Proliferative index on the y axis and CLK1-pos vs CLK1-negative samples on the x-axis: `output/proliferative_index_plot.pdf`

