## Analysis Modules
This directory contains various analysis modules in the pbta-splicing-hgat project.
See the README of an individual analysis modules for more information about that module.

### Modules at a glance
The table below is intended to help project organizers quickly get an idea of what files (and therefore types of data) are consumed by each analysis module, what the module does, and what output files it produces that can be consumed by other analysis modules.
This is in service of documenting interdependent analyses.

Note that _nearly all_ modules use the harmonized clinical data file (`pbta-histologies.tsv`) even when it is not explicitly included in the table below.

| Module | Input Files | Brief Description | Output Files Consumed by Other Analyses |
|--------|-------|-------------------|--------------|
| [`pan_cancer`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/chromosomal-instability) | `pbta-histologies.tsv` <br> `pbta-sv-manta.tsv.gz` <br> `pbta-cnv-cnvkit.seg.gz` | Evaluates chromosomal instability by calculating chromosomal breakpoint densities and by creating circular plot visuals | N/A
| [`cnv-chrom-plot`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-chrom) | `pbta-cnv-consensus-gistic.zip` <br> `analyses/copy_number_consensus_call/results/pbta-cnv-consensus.seg` | Plots genome wide visualizations relating to copy number results | N/A
| [`cnv-comparison`](https://github.com/AlexsLemonade/OpenPBTA-analysis/tree/master/analyses/cnv-comparison) | Earlier version of SEG files | *Deprecated*; compared earlier version of the CNV methods. | N/A
