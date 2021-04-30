## Analysis Modules
This directory contains various analysis modules in the pbta-splicing-hgat project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance
The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
This is in service of documenting interdependent analyses.

Note that _nearly all_ modules use the harmonized clinical data file (`pbta-histologies.tsv`) even when it is not explicitly included in the table below.

| Module | Input Files | Brief Description | Output Files Consumed by Other Analyses |
|--------|-------|-------------------|--------------|
| [`pan_cancer`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/pan_cancer) | `pbta-histologies.RNA-Seq.initial.tsv` <br> `rmats_output.gz` <br> | -Consensus clustering of samples into cluster for downstream analyses (eg. survival). <br> -Generates splicing index table and plot to assess aberrant splicing (compared to n=9 healthy samples). <br> -Oncoplot generation to help visualize mutations, fusions and splicing. | N/A
| [`global patterns of aberrant splicing`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/global_patterns) | `pbta-histologies.RNA-Seq.initial.tsv` <br> `rmats_output.gz` <br> `majiq_output.gz` | Generate tables and plots to assess global patterns of aberrant splicing, including splicing types, splicing index, and splicing heterogeneity | N/A
| [`mRNA diff. expression`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/mRNA_diff_expr) | `rsem.gz`  | Generate differential expression table and volcano plot between two conditions | N/A
| [`motif_binding`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/motif_binding) | `rmats_output.gz` <br> | Identify enrichment of all known RBP motifs within mis-spliced exons and flanking introns | N/A
| [`overlap_of_spliced_events`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/overlap_of_spliced_events) | `rmats_output.gz` <br> `majiq_output.gz`  | Generate overlap of identified mis-splicing events given conditions | N/A
| [`visualize_splicing_event`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/visualize_splicing_event) | `rmats_output.gz` <br> `majiq_output.gz`  | Generate plots to help visualize splicing events | N/A
| [`splicing_events_functional_sites`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/splicing_events_functional_sites) | `rmats_output.gz` <br> `majiq_output.gz` <br> `*uniprot.bed` | Identify signficiant aberrant splicing events that result in loss/gain of functional domain/site | `dominant_events_lsvs.total.*wo.txt` <br> `dominant_events_lsvs.total.*ggplot.txt`
| [`visualizing_splicing_events_across_sites`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/visualizing_splicing_events_across_sites) | `rmats_output.gz` <br> `majiq_output.gz`  | Generate plots and compare distribution of PSI values across given set of groups (gene types, functional site types, etc) | N/A
| [`pairwise_correlation`](https://github.com/naqvia/pbta-splicing-hgat/tree/main/analyses/splicing_events_functional_sites) | `rmats_output.gz` <br> `majiq_output.gz` | Compute and plot pairwise correlations between given splicing events | NA
