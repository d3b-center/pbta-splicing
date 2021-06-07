# DNA repair index

Module authors: Komal S. Rathi (@komalsrathi)

## Introduction

Purpose of this module is to do geneset variation analysis using R package GSVA with the following C2 CP KEGG genesets:

1. KEGG_MISMATCH_REPAIR
2. KEGG_HOMOLOGOUS_RECOMBINATION
3. KEGG_NON_HOMOLOGOUS_END_JOINING
4. KEGG_BASE_EXCISION_REPAIR

## Input

1. Expression: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-rsem-tpm-uncorrected.rds` from release hgg-dmg 20201109
2. Metadata: `../../data/pbta-hgat-dx-prog-pm-gene-expression-subset-metadata.rds`

## Output

Boxplots of plots with DNA repair index on the y axis and CLK1-pos vs CLK1-negative samples on the x-axis: `output/dna_repair_index_plot.pdf`