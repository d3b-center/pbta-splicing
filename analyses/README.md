### Get rMATS data files
<br>**Run shell script to get merged rMATS result tables**
```
./download_data.sh
```
## Analysis Modules
This directory contains various analysis modules in the pbta-splicing-hgat project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance
The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
This is in service of documenting interdependent analyses.

Note that _nearly all_ modules use the harmonized clinical data file (`histologies.tsv`) even when it is not explicitly included in the table below.

| Module |Brief Description |
|--------|------------------|
| [`splicing-burden-index`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/splicing_index) | Compute and compare splicing burden indices for each tumor sample across pediatric brain histologies. Module also performs related analyses assessing global patterns based on SBI, including total make-up of splicing types, differential gene expression between high vs low SBI tumors, survival associations, and histology-specific/shared splicing signatures.
| [`histology-specific-splicing`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/histology-specific-splicing) | Identify and visualize specific splicing signatures based on clinical-annotated histologies.
| [`survival`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/survival) | Survival associations stratified by splicing burden, histologies, and subtypes
| [`splicing-derived-clustering`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/clustering_analysis) | Consensus clustering and relevant associations of tumor splicing quantifications (PSI), including cluster-specific expression differences, gene-set enrichments and splicing burden associations.
| [`splicing-factor-dysregulation`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/splicing_index) | Differential gene expression of known splicing factors comparing controls vs midline HGG tumors
| [`splicing-events-functional-sites`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/splicing_events_functional_sites) | Identify signficiant recurrent aberrant splicing events that result in loss/gain of functional sites in midline HGGs
| [`CLK1-splicing-correlations`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/CLK1-splicing_correlations) | Determine significant correlations of CLK1 splicing with  expression, SRSF phospoproteomics, and SBI
| [`ONT_visualization`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/CLK1-splicing_correlations) | Plots and comparisons related to ONT (long-RNA-seq)
| [`KNS42_cell-line`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/KNS42_cell-line`) | Analysis assessing CLK1 impact on KNS42 cell lines
| [`CLK1-splicing-impact`](https://github.com/d3b-center/pbta-splicing/tree/main/analyses/CLK1_splicing_impact) | Identify impact of CLK1 splicing by identifying gene expression and splicing changes of control vs morpholino treated cell lines
